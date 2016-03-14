#include <vector_types.h>
#include <cuda_runtime.h>

#include "cuda_declarations.h"
#include "alignment_score_gpu.cuh"
#include "data_management.cuh"
#include "sequences.h"
#include "main_cu.h"
#include "exceptions.h"

using Data::Sequences;
using Exceptions::Exception;

#define Y_STEPS     ALIGNMENT_SCORE_Y_STEPS
#define WIN_SIZE    winSize
#define MEM_OFFSET  memOffset
#define BLOCK_SIZE  ALIGNMENT_SCORE_BLOCK_SIZE

#define seqXNo      (blockIdx.x * blockDim.x + threadIdx.x)
#define seqYNo      (blockIdx.y * blockDim.y + threadIdx.y)
#define startX      (tex1Starts[seqXNo])
#define startY      (tex2Starts[seqYNo])

__global__ void NeedlemanWunschGlobalScoreKernel(short2* AF, int* scores, bool border = false)
{
    // SUBSTITUTION MATRIX GOES TO SHARED MEMORY
    __shared__ char shmSM[SUBSTITUTION_MATRIX_SIZE];
    short idx = threadIdx.y * blockDim.x + threadIdx.x;
    shmSM[idx] = substitutionMatrix[idx];
    idx += BLOCK_SIZE;
    shmSM[idx] = substitutionMatrix[idx];
    idx += BLOCK_SIZE;
    if(idx < 576)
        shmSM[idx] = substitutionMatrix[idx];
    
    
    __syncthreads();
    
    // |\xxx    we do not compute x
    // | \xx
    // |  \x
    // |___\

    if(border && (seqXNo > seqYNo))
        return;

    //we determine our position in the window
    //int seqXNo = blockIdx.x * blockDim.x + threadIdx.x; //there is an idea to do not hold this but
    //int seqYNo = blockIdx.y * blockDim.y + threadIdx.y; //recalculate this every time to save registers
    int blockThread = threadIdx.x + threadIdx.y * blockDim.x; //0...(BLOCK_SIZE-1)

    //int startX = tex1Starts[seqXNo];
    //int startY = tex2Starts[seqYNo];
    short2 lengthXY;
    lengthXY.x = tex1Starts[seqXNo + 1] - startX;
    lengthXY.y = tex2Starts[seqYNo + 1] - startY;

    if((lengthXY.x == 0) || (lengthXY.y == 0))//if there is nothing to do -> quit
        return;

    //startPosA == thread number within whole grid
    int startPosA = seqYNo * WIN_SIZE + seqXNo;

////for me - it will not give us expected acceleration and is expensive
//    //A -> our individual A and F matrix
//    int* A = &Atabs[startPosA];

    //initialization of the -1 row in A matrix
    // - 2 bytes for element of A matrix
    // - 2 bytes for element of F matrix
    for(short x = 0; x < lengthXY.x; x++)
    {
        //PACK(A element, F element)
        //(x + 1) because the first element should be -gapEx
        //AandF[startPosA + x * MEM_OFFSET] = PACK(-gapEx * (x + 1), SHORT_MIN + gapEx);
        short2 tmp;
        tmp.x = -gapEx * (x + 1);
        tmp.y = SHORT_MIN + gapEx;
        AF[startPosA + x * MEM_OFFSET] = tmp;
    }

    //one element of sharedA consist of:
    // - one A element
    // - one E element
    __shared__ short2 AE_shared[Y_STEPS][BLOCK_SIZE];
    //elements of Y sequence go to sharedYSeq
    __shared__ int sharedYSeq[Y_STEPS/4][BLOCK_SIZE];


    short2 AE_current;
    AE_current.x = 0;

    // |
    // |
    // |
    // V
    for (short y = 0; y < lengthXY.y; y += Y_STEPS)
    {
        short2 A_init_upleft;
        A_init_upleft.x = -gapEx * y;

        //initialialization of the -1 column in A matrix
        // - one element of A matrix
        // - one element of E matrix
        for (short i = 0; i < Y_STEPS; i++)
        {
            short2 tmp;
            tmp.x = -gapEx * (y + i + 1);
            tmp.y = SHORT_MIN + gapEx;
            AE_shared[i][blockThread] = tmp;
        }


        //we read elements of the Y sequence
        for (short i = 0; i < Y_STEPS/4; i++)
        {
            sharedYSeq[i][blockThread] = PACK_BYTES(tex1Dfetch(texSeqsY, startY + y + i*4 + 0),
                                                    tex1Dfetch(texSeqsY, startY + y + i*4 + 1),
                                                    tex1Dfetch(texSeqsY, startY + y + i*4 + 2),
                                                    tex1Dfetch(texSeqsY, startY + y + i*4 + 3));
        }


        //------>
        for (short x = 0; x < lengthXY.x; x++)
        {
            //actual up_left gets a value of recent read from the global memory
            //and actual read value is stored in first two bites of A_upleft
            A_init_upleft.y = A_init_upleft.x;

            char2 XYSeq;
            XYSeq.x = tex1Dfetch(texSeqsX, startX + x);
            
            //read from global memory
            short2 AF_up = AF[startPosA + x * MEM_OFFSET];
            AE_current.x = AF_up.x;
            

            //A_init -> up element read in previous iteration from global memory (up-left)
            A_init_upleft.x = AF_up.x;

            int F_current = AF_up.y;
            int similarity;
            short ymin = min(Y_STEPS, lengthXY.y - y); //(i < Y_STEPS) && (i + y < lengthY)
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            for(short i = 0; i < ymin; i++)
            {
                XYSeq.y = (sharedYSeq[i/4][blockThread] >> (((15-i)%4) * 8)) & 0xFF;
                
                similarity = shmSM[XYSeq.y*lettersCount + XYSeq.x];
                
                F_current = max(F_current - gapEx, AE_current.x - gapOp);
                AE_current.y = max(AE_shared[i][blockThread].y - gapEx, AE_shared[i][blockThread].x - gapOp);

                AE_current.x = max(AE_current.y, F_current);
                AE_current.x = max(AE_current.x, similarity + A_init_upleft.y);
                
                A_init_upleft.y = AE_shared[i][blockThread].x;
                AE_shared[i][blockThread] = AE_current;
            }
            //write variables to global memory for next loop
            short2 AF_tmp;
            AF_tmp.x = AE_current.x;
            AF_tmp.y = F_current;
            AF[startPosA + x * MEM_OFFSET] = AF_tmp;

        }
    }

    //here write result (A_current) to global memory
    scores[startPosA] = AE_current.x;
}

#undef Y_STEPS
#undef WIN_SIZE
#undef MEM_OFFSET
#undef BLOCK_SIZE
#undef seqXNo
#undef seqYNo
#undef startX
#undef startY

const char* NeedlemanWunschGlobalScoreGpu::getAlgorithmName()
{
    return "NeedlemanWunschGlobalScoreGpu";
}

void NeedlemanWunschGlobalScoreGpu::kernelInvocation(void* kernelParams)
{
    KernelScoreInvokingParams* params = (KernelScoreInvokingParams*) kernelParams;
    NeedlemanWunschGlobalScoreKernel<<<params->gridShape,params->blockShape>>>(params->AF, params->scores, params->border);
}
