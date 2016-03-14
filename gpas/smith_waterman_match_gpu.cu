#include "hi_res_timer.h"
#include "alignment_match_gpu.cuh"

using Data::Sequences;
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
#define BACK(x,y)   back[startPosA[blockThread] + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET]

__global__ void SmithWatermanMatchKernel(short2* AF_maxXY, unsigned int* back, int* scores, short rowWidth, bool border = false)
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
    int __shared__ startPosA[BLOCK_SIZE];
    startPosA[blockThread] = seqYNo * WIN_SIZE + seqXNo;

    //initialization of the -1 row in A matrix
    // - 2 bytes for element of A matrix
    // - 2 bytes for element of F matrix
    for(short x = 0; x < lengthXY.x; x++)
    {
        short2 tmp;
        //(x + 1) because the first element should be -gapEx
        tmp.x = 0;
        tmp.y = SHORT_MIN + gapEx;
        AF_maxXY[startPosA[blockThread] + x * MEM_OFFSET] = tmp;

        //fill the -1 row of "back" array
        BACK(x,-1) = 0; //in case of SW -> STOP when get to the end of the sequence
    }

    //fill the -1 column of "back" array
    for(short y = 0; y < lengthXY.y; y+=Y_STEPS)
    {
        BACK(-1,y) = 0; //in case of SW -> STOP when get to the end of the sequence
    }
    BACK(-1,-1) = 0; //stop element

    //one element of AE_shared consist of:
    // - one A element
    // - one E element
    __shared__ short2 AE_shared[Y_STEPS][BLOCK_SIZE];
    //elements of Y sequence go to sharedYSeq
    __shared__ int sharedYSeq[Y_STEPS/4][BLOCK_SIZE];


    short2 AF_current;
    AF_current.x = 0;

    __shared__ short2 ymin_score[BLOCK_SIZE]; //stores ymin and score
    ymin_score[blockThread].y = 0;

    __shared__ short2 maxXY[BLOCK_SIZE];
    maxXY[blockThread].x = lengthXY.x - 1;
    maxXY[blockThread].y = 0;

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

        ymin_score[blockThread].x = min(Y_STEPS, lengthXY.y - y); //(i < Y_STEPS) && (i + y < lengthY)

        //------>
        for (short x = 0; x < lengthXY.x; x++)
        {
            //actual up_left gets a value of recent read value from the global memory
            //and actual read value is stored in first two bites of A_upleft
            A_init_upleft.y = A_init_upleft.x;

            char2 XYSeq;
            XYSeq.x = tex1Dfetch(texSeqsX, startX + x);

            //read from global memory
            short2 AF_up = AF_maxXY[startPosA[blockThread] + x * MEM_OFFSET];

            //A_init -> up element read in previous iteration from global memory (up-left)
            A_init_upleft.x = AF_up.x;

            short2 AE_left;
            int E_current;
            int similarity;
            unsigned int back8 = 0;
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            //  |  /|  /|
            //  | / | / |
            //  |/  |/  V
            for(short i = 0; i < ymin_score[blockThread].x; i++)
            {
                AE_left = AE_shared[i][blockThread];

                XYSeq.y = (sharedYSeq[i/4][blockThread] >> (((15-i)%4) * 8)) & 0xFF;

                //similarity = substitutionMatrix[XYSeq.y*lettersCount + XYSeq.x];
                similarity = shmSM[XYSeq.y*lettersCount + XYSeq.x];
                similarity += A_init_upleft.y;

                E_current = max(AE_left.y - gapEx, AE_left.x - gapOp);
                AF_current.y = max(AF_up.y - gapEx, AF_up.x - gapOp);

                AF_current.x = max(E_current, AF_current.y);
                AF_current.x = max(AF_current.x, similarity);
                AF_current.x = max(AF_current.x, 0);

                //"back" array
                back8 <<= 1;
                back8 |= (((AF_current.x==E_current) && (AF_current.x!=AF_current.y)) || (AF_current.x==similarity)) && AF_current.x; //if go left
                back8 <<= 1;
                back8 |= ((AF_current.x==AF_current.y) || (AF_current.x==similarity)) && AF_current.x; //if go up
                back8 <<= 1;
                back8 |= (AF_current.y == (AF_up.y - gapEx)); //if continue up
                back8 <<= 1;
                back8 |= (E_current == (AE_left.y - gapEx)); //if continue left


                if (AF_current.x > ymin_score[blockThread].y)
                {
                    maxXY[blockThread].x = x;
                    maxXY[blockThread].y = y + i;
                }

                ymin_score[blockThread].y = max(ymin_score[blockThread].y, AF_current.x);

                //initialize variables for next iterations
                short2 AE_tmp;
                AE_tmp.x = AF_current.x;
                AE_tmp.y = E_current;
                AE_shared[i][blockThread] = AE_tmp;
                A_init_upleft.y = AE_left.x;
                AF_up = AF_current;

            }

            //we want the last row of back8 to be completed
            back8 <<= 4 * (Y_STEPS - ymin_score[blockThread].x);

            //write variables to global memory for next loop
            AF_maxXY[startPosA[blockThread] + x * MEM_OFFSET] = AF_current;
            BACK(x,y) = back8;

        }
    }

    //here write result (AF_current) to global memory
    scores[startPosA[blockThread]] = ymin_score[blockThread].y;
    AF_maxXY[startPosA[blockThread]] = maxXY[blockThread];
}

#undef BACK
#define BACK(x,y)   back[startPosA + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET]









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

__global__ void SmithWatermanBackKernel(short2* maxXY, unsigned int* back, short rowWidth, unsigned int* matchesX, unsigned int* matchesY, bool border = false)
{
    if(border && (seqXNo > seqYNo))
        return;

    short2 lengthXY;
    lengthXY.x = tex1Starts[seqXNo + 1] - startX;
    lengthXY.y = tex2Starts[seqYNo + 1] - startY;

    //THIS CONDITION IS NECESARRY IN SW!!!
    if((lengthXY.x == 0) || (lengthXY.y == 0))//if there is nothing to do -> quit
        return;
    
    

    //startPosA == thread number within whole grid
    int startPosA = seqYNo * WIN_SIZE + seqXNo;

    //short2 lengthXY;
    //lengthXY.x = tex1Starts[seqXNo + 1] - startX;
    //lengthXY.y = tex2Starts[seqYNo + 1] - startY;
    short2 myMaxXY = maxXY[startPosA];
    short2 indexXY = myMaxXY;

    //indexXY.x = lengthXY.x - 1; //lengthX (-1 because of addressing in BACK(x,y))
    //indexXY.y = lengthXY.y - 1; //lengthY


    short carret = 0;

    unsigned int back8 = BACK(indexXY.x, indexXY.y);

    back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4;

    unsigned char back1 = back8 & ELEMENT; //current element of back array
    back8 >>= 4;

    unsigned char prevDirection = CROSSWISE;// 1100 == 12 =>crosswise
    unsigned char todo;

    unsigned int tmpMatchX;
    unsigned int tmpMatchY;


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

        if ((indexXY.x == myMaxXY.x) && (indexXY.y > myMaxXY.y))
            back1 = 0;
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

const char* SmithWatermanMatchGpu::getAlgorithmName()
{
    return "SmithWatermanMatchGpu";
}

void SmithWatermanMatchGpu::kernelInvocation(void* kernelParams)
{
    KernelMatchInvokingParams* params = (KernelMatchInvokingParams*) kernelParams;
    SmithWatermanMatchKernel<<<params->gridShape,params->blockShape>>>(params->AF, params->back, params->scores, params->rowWidth, params->border);
}

void SmithWatermanMatchGpu::kernelBacktraceInvocation(void* kernelParams)
{
    KernelMatchInvokingParams* params = (KernelMatchInvokingParams*) kernelParams;//
    SmithWatermanBackKernel<<<params->gridShape,params->blockShape>>>(params->AF, params->back, params->rowWidth, params->matchesX, params->matchesY, params->border);
}
