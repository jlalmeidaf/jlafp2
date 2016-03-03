#include <vector_types.h>
#include <cuda_runtime.h>

#include "cuda_declarations.h"
#include "alignment_match_gpu.cuh"
#include "data_management.cuh"
#include "thread_manager.h"
#include "sequences.h"
#include "main_cu.h"
#include "exceptions.h"
#include "hi_res_timer.h"

using Data::Sequences;
using Data::MatchesManager;
using Exceptions::Exception;

#define Y_STEPS     ALIGNMENT_MATCH_Y_STEPS
#define WIN_SIZE    winSize
#define MEM_OFFSET  memOffset
#define BLOCK_SIZE  ALIGNMENT_MATCH_BLOCK_SIZE

#define seqXNo      (blockIdx.x * blockDim.x + threadIdx.x)
#define seqYNo      (blockIdx.y * blockDim.y + threadIdx.y)
#define startX      (tex1Starts[seqXNo])
#define startY      (tex2Starts[seqYNo])

/*******************************************************************************
 * "back" consist of 4bits x 8 (=32bits):                                      *
 * 4bits:                                                                      *
 * -> 0 and 1 bits:                                                            *
 *      - 00 -> stop                                                           *
 *      - 01 -> up                                                             *
 *      - 10 -> left                                                           *
 *      - 11 -> crosswise                                                      *
 * -> 2 bit:                                                                   *
 *      - 0 -> not continueUp                                                  *
 *      - 1 -> continueUp                                                      *
 * -> 3 bit:                                                                   *
 *      - 0 -> not continueLeft                                                *
 *      - 1 -> continueLeft                                                    *
 *                                                                             *
 * BACK:                                                                       *
 * back[startPosA + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET]    *
 * BACK(-1,-1) => the firs element of -1 row (and -1 column)                   *
 *******************************************************************************/
#define BACK(x,y)   back[startPosA + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET]

__global__ void NeedlemanWunschGlobalMatchKernel(short2* AF, unsigned int* back, int* scores, short rowWidth, bool border = false)
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
    
    /***************************************************************************
     * |\xxx                                                                   *
     * | \xx    we do not compute x                                            *
     * |  \x                                                                   *
     * |___\                                                                   *
     ***************************************************************************/
    if(border && (seqXNo > seqYNo))
        return;

    int blockThread = threadIdx.x + threadIdx.y * blockDim.x; //0...(BLOCK_SIZE-1)

    short2 lengthXY;
    lengthXY.x = tex1Starts[seqXNo + 1] - startX;
    lengthXY.y = tex2Starts[seqYNo + 1] - startY;

    if((lengthXY.x == 0) || (lengthXY.y == 0))//if there is nothing to do -> quit
        return;

    //startPosA == thread number within whole grid
    int startPosA = seqYNo * WIN_SIZE + seqXNo;

    //initialization of the -1 row in A matrix
    // - 2 bytes for element of A matrix
    // - 2 bytes for element of F matrix
    for(short x = 0; x < lengthXY.x; x++)
    {
        short2 tmp;
        //(x + 1) because the first element should be -gapEx
        tmp.x = -gapEx * (x + 1);
        tmp.y = SHORT_MIN + gapEx;
        AF[startPosA + x * MEM_OFFSET] = tmp;

        //fill the -1 row of "back" array
        BACK(x,-1) = 9; //0000 0000 0000 0000 0000 0000 0000 1001 == 9
    }

    //fill the -1 column of "back" array
    for(short y = 0; y < lengthXY.y; y+=Y_STEPS)
    {
        BACK(-1,y) = 1717986918; //0110 0110 0110 0110 0110 0110 0110 0110 = 1717986918
    }
    BACK(-1,-1) = 0; //stop element

    //one element of AE_shared consist of:
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
            //actual up_left gets a value of recent read value from the global memory
            //and actual read value is stored in first two bites of A_upleft
            A_init_upleft.y = A_init_upleft.x;

            char2 XYSeq;
            XYSeq.x = tex1Dfetch(texSeqsX, startX + x);

            //read from global memory
            short2 AF_up = AF[startPosA + x * MEM_OFFSET];

            //A_init -> up element read in previous iteration from global memory (up-left)
            A_init_upleft.x = AF_up.x;
            int F_up;// = AF_up.y;
            AE_current.x = AF_up.x;

            //short2 AE_left;
            int F_current = AF_up.y;
            //int F_up;
            int similarity;
            unsigned int back8 = 0;
            short ymin = min(Y_STEPS, lengthXY.y - y); //(i < Y_STEPS) && (i + y < lengthY)
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            for(short i = 0; i < ymin; i++)
            {
                //AE_left = AE_shared[i][blockThread];

                XYSeq.y = (sharedYSeq[i/4][blockThread] >> (((15-i)%4) * 8)) & 0xFF;

                similarity = max(F_current - gapEx, AE_current.x - gapOp);
                F_up = (similarity==F_current-gapEx);
                F_current = similarity;
                //similarity = substitutionMatrix[XYSeq.y*lettersCount + XYSeq.x];
                similarity = shmSM[XYSeq.y*lettersCount + XYSeq.x];
                similarity += A_init_upleft.y;

                AE_current.y = max(AE_shared[i][blockThread].y - gapEx, AE_shared[i][blockThread].x - gapOp);
                //F_current = max(F_up - gapEx, AE_current.x - gapOp);

                AE_current.x = max(AE_current.y, F_current);
                AE_current.x = max(AE_current.x, similarity);

                //"back" array
                back8 <<= 1;
                back8 |= ((AE_current.x==AE_current.y) && (AE_current.x!=F_current)) || (AE_current.x==similarity); //if go left
                back8 <<= 1;
                back8 |= (AE_current.x==F_current) || (AE_current.x==similarity); //if go up
                back8 <<= 1;
                back8 |= F_up;//(F_current == (F_up - gapEx)); //if continue up
                back8 <<= 1;
                back8 |= (AE_current.y == (AE_shared[i][blockThread].y - gapEx)); //if continue left

                //initialize variables for next iterations
                //short2 AE_tmp;
                //AE_tmp.x = AF_current.x;
                //AE_tmp.y = E_current;
                A_init_upleft.y = AE_shared[i][blockThread].x;
                AE_shared[i][blockThread] = AE_current;
                //A_init_upleft.y = AE_left.x;
                //F_up = F_current;

            }

            //we want the last row of back8 to be completed
            back8 <<= 4 * (Y_STEPS - ymin);


            short2 AF_tmp;
            AF_tmp.x = AE_current.x;
            AF_tmp.y = F_current;
            //write variables to global memory for next loop
            AF[startPosA + x * MEM_OFFSET] = AF_tmp;

            BACK(x,y) = back8;

        }
    }

    //here write result (AF_current) to global memory
    scores[startPosA] = AE_current.x;
}


/*******************************************************************************
 * If you don't want to use the texture when backtracing then                  *
 * comment the two lines below.                                                *
 *                                                                             *
 * UNFORTUNATELLY (OR NOT) WE CAN'T USE TEXTURE FOR THIS PURPOSE               *
 * BECAUSE THE MAXIMAL NUMBER OF TEXTURE ELEMENTS (NOT BYTES) IS 2^27          *
 * WHICH WE CAN EASLY EXCEED.                                                  *
 *******************************************************************************/
//#undef  BACK
//#define BACK(x,y)   tex1Dfetch(texBack, startPosA + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET)

/*******************************************************************************
 * "back" consist of 4bits x 8 (=32bits):                                      *
 * 4bits:                                                                      *
 * -> 0 and 1 bits:                                                            *
 *      - 00 -> stop                                                           *
 *      - 01 -> up                                                             *
 *      - 10 -> left                                                           *
 *      - 11 -> crosswise                                                      *
 * -> 2 bit:                                                                   *
 *      - 0 -> not continueUp                                                  *
 *      - 1 -> continueUp                                                      *
 * -> 3 bit:                                                                   *
 *      - 0 -> not continueLeft                                                *
 *      - 1 -> continueLeft                                                    *
 *                                                                             *
 *******************************************************************************/

#define STOP         0
#define UP           4
#define LEFT         8
#define CROSSWISE   12
#define DIRECTION   12
#define CONTIN_UP    2
#define CONTIN_LEFT  1
#define ELEMENT     15

__global__ void NeedlemanWunschGlobalBackKernel(unsigned int* back, short rowWidth, unsigned int* matchesX, unsigned int* matchesY, bool border = false)
{
    if(border && (seqXNo > seqYNo))
        return;

    short2 lengthXY;
    lengthXY.x = tex1Starts[seqXNo + 1] - startX;
    lengthXY.y = tex2Starts[seqYNo + 1] - startY;

    if((lengthXY.x == 0) || (lengthXY.y == 0))//if there is nothing to do -> quit
        return;

    //startPosA == thread number within whole grid
    int startPosA = seqYNo * WIN_SIZE + seqXNo;


    short2 indexXY;
    indexXY.x = lengthXY.x - 1; //lengthX (-1 because of addressing in BACK(x,y))
    indexXY.y = lengthXY.y - 1; //lengthY
    
    unsigned int back8 = BACK(indexXY.x, indexXY.y);

    short carret = 0;
    unsigned char prevDirection = CROSSWISE;// 1100 == 12 =>crosswise
    unsigned char back1; //current element of back array
    unsigned char todo;

    unsigned int tmpMatchX;
    unsigned int tmpMatchY;
 
    back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4;

    back1 = back8 & ELEMENT;
    back8 >>= 4;

    
    while(back1 & DIRECTION) //while(direction != STOP)
    {

        if( ((prevDirection & DIRECTION) == UP) && (prevDirection & CONTIN_UP) )
        {
            todo = UP;
        }
        else if( ((prevDirection & DIRECTION) == LEFT) && (prevDirection & CONTIN_LEFT) )
        {
            todo = LEFT;
        }
        else if ((back1 & DIRECTION) == UP)
        {
            todo = UP;
        }
        else if ((back1 & DIRECTION) == LEFT)
        {
            todo = LEFT;
        }
        else //if (back1 & DIRECTION == CROSSWISE)
        {
            todo = CROSSWISE;
        }

        tmpMatchY <<= 8;
        tmpMatchX <<= 8;


        if (todo == LEFT)
        {
            tmpMatchY |= (unsigned char)'-';
            tmpMatchX |= revConvConst[tex1Dfetch(texSeqsX, startX + indexXY.x)];


            indexXY.x--;
            back8 = BACK(indexXY.x, indexXY.y);
            back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4; //because of the last row of back array
        }
        else if (todo == UP)
        {
            tmpMatchX |= (unsigned char)'-';
            tmpMatchY |= revConvConst[tex1Dfetch(texSeqsY, startY + indexXY.y)];

            indexXY.y--;
            if((indexXY.y % 8) == 7)
                back8 = BACK(indexXY.x, indexXY.y);
        }
        else //if (todo == CROSSWISE)
        {
            tmpMatchX |= revConvConst[tex1Dfetch(texSeqsX, startX + indexXY.x)];
            tmpMatchY |= revConvConst[tex1Dfetch(texSeqsY, startY + indexXY.y)];

            indexXY.x--;
            indexXY.y--;

            back8 = BACK(indexXY.x, indexXY.y);
            back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4; //because of the last row of back array
        }
        

        prevDirection = todo | back1&3;
        back1 = back8 & ELEMENT;
        back8 >>= 4;


        carret++;
        if( !(carret%4) )
        {
            //save results to global memory
            matchesX[startPosA + (carret/4 - 1) * MEM_OFFSET] = tmpMatchX;
            matchesY[startPosA + (carret/4 - 1) * MEM_OFFSET] = tmpMatchY;
            tmpMatchX = 0;
            tmpMatchY = 0;
        }

    }
    

    tmpMatchX <<= 8;
    tmpMatchY <<= 8;
    tmpMatchX |= (unsigned char)0;//end of match
    tmpMatchY |= (unsigned char)0;//end of match
    tmpMatchX <<= ((4 - ((carret + 1) % 4)) % 4) * 8;
    tmpMatchY <<= ((4 - ((carret + 1) % 4)) % 4) * 8;

    
    carret+=4;
    matchesX[startPosA + (carret/4 - 1) * MEM_OFFSET] = tmpMatchX;
    matchesY[startPosA + (carret/4 - 1) * MEM_OFFSET] = tmpMatchY;
 
}

#undef STOP
#undef UP
#undef LEFT
#undef CROSSWISE
#undef DIRECTION
#undef CONTIN_UP
#undef CONTIN_LEFT
#undef ELEMENT

#undef Y_STEPS
#undef WIN_SIZE
#undef MEM_OFFSET
#undef BLOCK_SIZE
#undef seqXNo
#undef seqYNo
#undef startX
#undef startY

#undef BACK

const char* NeedlemanWunschGlobalMatchGpu::getAlgorithmName()
{
    return "NeedlemanWunschGlobalMatchGpu";
}

void NeedlemanWunschGlobalMatchGpu::kernelInvocation(void* kernelParams)
{
    KernelMatchInvokingParams* params = (KernelMatchInvokingParams*) kernelParams;
    NeedlemanWunschGlobalMatchKernel<<<params->gridShape,params->blockShape>>>(params->AF, params->back, params->scores, params->rowWidth, params->border);
}

void NeedlemanWunschGlobalMatchGpu::kernelBacktraceInvocation(void* kernelParams)
{
    KernelMatchInvokingParams* params = (KernelMatchInvokingParams*) kernelParams;
    NeedlemanWunschGlobalBackKernel<<<params->gridShape,params->blockShape>>>(params->back, params->rowWidth, params->matchesX, params->matchesY, params->border);
}
