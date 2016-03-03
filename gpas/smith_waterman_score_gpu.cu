#include <vector_types.h>
#include <cuda_runtime.h>

#include "cuda_declarations.h"
#include "alignment_score_gpu.cuh"
#include "data_management.cuh"
#include "thread_manager.h"
#include "sequences.h"
#include "main_cu.h"
#include "exceptions.h"
#include "hi_res_timer.h"

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

__global__ void SmithWatermanScoreKernel(short2* AF, int* scores, bool border = false)
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
    
    // |\xxx
    // | \xx
    // |  \x
    // |___\    we do not compute x

    if(border && (seqXNo > seqYNo))
        return;

    //we determine our position in the window
    int blockThread = threadIdx.x + threadIdx.y * blockDim.x; //0...(BLOCK_SIZE-1)

    short2 lengthXY;
    lengthXY.x = tex1Starts[seqXNo + 1] - startX;
    lengthXY.y = tex2Starts[seqYNo + 1] - startY;

    if((lengthXY.x == 0) || (lengthXY.y == 0))//if there is nothing to do -> quit
        return;

    //startPosA == thread number within whole grid
    int startPosA = seqYNo * WIN_SIZE + seqXNo;

    //initialization of the -1 row in A and F matrices
    // - 2 bytes for element of A matrix
    // - 2 bytes for element of F matrix
    for(short x = 0; x < lengthXY.x; x++)
    {
        short2 tmp;
        tmp.x = 0; //A
        tmp.y = SHORT_MIN + gapEx; //F
        AF[startPosA + x * MEM_OFFSET] = tmp;
    }

    //one element of AE_shared consist of:
    // - one A element
    // - one E element
    __shared__ short2 AE_shared[Y_STEPS][BLOCK_SIZE];
    //elements of Y sequence go to sharedYSeq
    __shared__ int sharedYSeq[Y_STEPS/4][BLOCK_SIZE];


    short2 AE_current;
    AE_current.x = 0;

    short2 ymin_score; //stores ymion and score
    ymin_score.y = 0;

    // |
    // |
    // |
    // V
    for (short y = 0; y < lengthXY.y; y += Y_STEPS)
    {
        short2 A_init_upleft;
        A_init_upleft.x = 0;

        //initialialization of the -1 column in A matrix
        // - one element of A matrix
        // - one element of E matrix
        for (short i = 0; i < Y_STEPS; i++)
        {
            short2 tmp;
            tmp.x = 0;
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



        ymin_score.x = min(Y_STEPS, lengthXY.y - y); //(i < Y_STEPS) && (i + y < lengthY)

        //------>
        for (short x = 0; x < lengthXY.x; x++)
        {
            //actual up_left gets a value of recent read from the global memory
            //and actual read value is stored in first two bites of A_upleft
            A_init_upleft.y = A_init_upleft.x;

            char2 XYSeq;
            XYSeq.x = tex1Dfetch(texSeqsX, startX + x);

            //read from global memory
            short2 SimF_current;
            SimF_current = AF[startPosA + x * MEM_OFFSET];
            AE_current.x = SimF_current.x;
            //A_init -> up element read in previous iteration from global memory (up-left)
            A_init_upleft.x = SimF_current.x;

            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            for(short i = 0; i < ymin_score.x; i++)
            {
                XYSeq.y = (sharedYSeq[i/4][blockThread] >> (((15-i)%4) * 8)) & 0xFF;

                SimF_current.x = shmSM[XYSeq.y*lettersCount + XYSeq.x];

                SimF_current.y = max(SimF_current.y - gapEx, AE_current.x - gapOp);
                AE_current.y = max(AE_shared[i][blockThread].y - gapEx, AE_shared[i][blockThread].x - gapOp);

                AE_current.x = max(AE_current.y, SimF_current.y);
                AE_current.x = max(AE_current.x, SimF_current.x + A_init_upleft.y);
                AE_current.x = max(AE_current.x, 0);

                ymin_score.y = max(ymin_score.y, AE_current.x);

                //initialize variables for next iterations
                A_init_upleft.y = AE_shared[i][blockThread].x;
                AE_shared[i][blockThread] = AE_current;
            }
            //write variables to global memory for next loop
            short2 AF_tmp;
            AF_tmp.x = AE_current.x;
            AF_tmp.y = SimF_current.y;
            AF[startPosA + x * MEM_OFFSET] = AF_tmp;


        }

    }

    //here write result (A_current) to global memory
    scores[startPosA] = ymin_score.y;
}

#undef Y_STEPS
#undef WIN_SIZE
#undef MEM_OFFSET
#undef BLOCK_SIZE
#undef seqXNo
#undef seqYNo
#undef startX
#undef startY

const char* SmithWatermanScoreGpu::getAlgorithmName()
{
    return "SmithWatermanScoreGpu";
}

void SmithWatermanScoreGpu::kernelInvocation(void* kernelParams)
{
    KernelScoreInvokingParams* params = (KernelScoreInvokingParams*) kernelParams;
    SmithWatermanScoreKernel<<<params->gridShape,params->blockShape>>>(params->AF, params->scores, params->border);
}
