#ifndef _ALIGNMENT_MATCH_GPU_CUH_
#define _ALIGNMENT_MATCH_GPU_CUH_

#include "alignment_score_gpu.cuh"
#include "data_management.cuh"
#include "matches_manager.h"

using namespace Data;


/*******************************************************************************
 * BLOCK_SHAPE defines the width and the height of the one block of threads.   *
 * Each thread corresponds to one pair of sequences and calculates one         *
 * one alignment. It's a good practice to have it multiple of 8. But the best  *
 * results should give 16.                                                     *
 *******************************************************************************/

#define ALIGNMENT_MATCH_BLOCK_SHAPE    16

/*******************************************************************************
 * BLOCK_SIZE defines the count of the threads runs in the multiprocessor in   *
 * the one block. We will threat it as a square and so it should have a value: *
 * BLOCK_SIZE = BLOCK_SHAPE^2                                                  *
 *******************************************************************************/

#define ALIGNMENT_MATCH_BLOCK_SIZE     (ALIGNMENT_MATCH_BLOCK_SHAPE*ALIGNMENT_MATCH_BLOCK_SHAPE)

/*******************************************************************************
 * Y_STEPS defines the number of steps within the sequence Y which will be     *
 * computed in one iteration without an access to the global memory.           *
 * Maximum Y_STEPS value is 12 (as long as size of the shared memory size is   *
 * 16KB per multiprocessor - see CUDA documentation).                          *
 *                                                                             *
 *                             !!!! WARNING !!!!                               *
 *    In contrast to the ALIGNMENT_SCORE_Y_STEPS the ALIGNMENT_MATCH_Y_STEPS   *
 *            CANNOT be modify without changes in the implementation!          *
 *******************************************************************************/

#define ALIGNMENT_MATCH_Y_STEPS         8



/*******************************************************************************
 * The RESHAPE_MEMORY_BLOCK_SIZE refere to the width and height of the shared  *
 * memory in one half-warp which should be the same value as thread count in   *
 * kernel invocation. This value shouldn't be modified.                        *
 *******************************************************************************/

#define ALIGNMENT_MATCH_BLOCK_X_SIZE   16

/*******************************************************************************
 * Function reorderMatchesMemory should be used to change order of global      *
 * memory allocated to the matches in the way that the consecutive letters of  *
 * the sequence to be in successive global memory cells.                       *
 *******************************************************************************/


#define BLOCK_X_SIZE  ALIGNMENT_MATCH_BLOCK_X_SIZE
#define MEM_OFFSET    memOffset

__global__ void reorderMatchesMemory(unsigned int* inMatches, unsigned int* outMatches, int maxMatchLength)
{

    // blocks must be launched in a grid with shape: dim3(n,1,1)
    //  ____   ____   ____   ____
    // |____| |____| |____| |____| ...
    //
    // each block transcripts 16 sequences using shared memory

    //the number of sequence that is transcripted by this thread
    int seqNoRead = blockIdx.x * BLOCK_X_SIZE + threadIdx.x;
    int seqNoWrite = blockIdx.x * BLOCK_X_SIZE + threadIdx.y;

    //BLOCK_X_SIZE + 1 -> to avoid bank conflicts
    __shared__ unsigned int shmMatches[BLOCK_X_SIZE][BLOCK_X_SIZE + 1];

    unsigned int fetch;
    unsigned int fetch2;

    //14 -> 0, 15 -> 0, 16 -> 16, 17 -> 16 ...
    int end = (maxMatchLength / BLOCK_X_SIZE) * BLOCK_X_SIZE;

    //main loop
    for (int i = 0; i < end; i += BLOCK_X_SIZE)
    {
        fetch = inMatches[seqNoRead + (i + threadIdx.y) * MEM_OFFSET];
        //changing the order of bytes in int
        fetch2 = fetch & 0xFF;
        fetch2 <<= 8;
        fetch >>= 8;
        fetch2 |= fetch & 0xFF;
        fetch2 <<= 8;
        fetch >>= 8;
        fetch2 |= fetch & 0xFF;
        fetch2 <<= 8;
        fetch >>= 8;
        fetch2 |= fetch & 0xFF;

        shmMatches[threadIdx.y][threadIdx.x] = fetch2;

        __syncthreads();

        outMatches[seqNoWrite * maxMatchLength + i + threadIdx.x] = shmMatches[threadIdx.x][threadIdx.y];

        __syncthreads();

    }

    //transcripting the end of sequecne (if maxMatchLength % BLOCK_X_SIZE != 0)
    if (end + threadIdx.y < maxMatchLength)
    {
        fetch = inMatches[seqNoRead + (end + threadIdx.y) * MEM_OFFSET];
        fetch2 = fetch & 0xFF;
        fetch2 <<= 8;
        fetch >>= 8;
        fetch2 |= fetch & 0xFF;
        fetch2 <<= 8;
        fetch >>= 8;
        fetch2 |= fetch & 0xFF;
        fetch2 <<= 8;
        fetch >>= 8;
        fetch2 |= fetch & 0xFF;

        shmMatches[threadIdx.y][threadIdx.x] = fetch2;
    }

    __syncthreads();

    if (end + threadIdx.x < maxMatchLength)
    {
        outMatches[seqNoWrite * maxMatchLength + end + threadIdx.x] = shmMatches[threadIdx.x][threadIdx.y];
    }
}

#undef BLOCK_X_SIZE
#undef MEM_OFFSET


void AlignmentMatchAlgorithm::run(ThreadManager* tm, Sequences* seqs, int gapOpen, int gapExt, int maxMultiprocessorCount, unsigned int windowSize, const char* outfile)
{
    matches1Manager = new MatchesManager(seqs->getSequenceNumber(), windowSize, ((seqs->getMaxSeqLen() - 1) / 4 + 1) * 4 * 2);
    matches2Manager = new MatchesManager(seqs->getSequenceNumber(), windowSize, ((seqs->getMaxSeqLen() - 1) / 4 + 1) * 4 * 2);
    calculateMatches = true;
    AlignmentScoreAlgorithm::run(tm, seqs, gapOpen, gapExt, maxMultiprocessorCount, windowSize, outfile);
}

AlignmentMatchAlgorithm::~AlignmentMatchAlgorithm()
{
    delete matches1Manager;
    delete matches2Manager;
    AlignmentScoreAlgorithm::~AlignmentScoreAlgorithm();
}

AlignmentScoreKernelInvokerParams* AlignmentMatchAlgorithm::getParams()
{
    AlignmentMatchKernelInvokerParams* result = new AlignmentMatchKernelInvokerParams();
    //default parameters' values (for AlignmentMatchAlgorithm class)
    result->blockShape = ALIGNMENT_MATCH_BLOCK_SHAPE;
    result->matches1Manager = matches1Manager;
    result->matches2Manager = matches2Manager;
    return result;
}

void AlignmentMatchAlgorithm::actualInvokedMethod(AlignmentScoreKernelInvokerParams* sParams, int gpuNo)
{
    //THIS FUNCTION IS CALLED BY THREAD MANAGER

    AlignmentMatchKernelInvokerParams* params = (AlignmentMatchKernelInvokerParams*)sParams;

    params->startSeqs1No = params->windowX * params->windowSize;
    params->startSeqs2No = params->windowY * params->windowSize;
    try
    {
        int maxSeqLength = params->seqs->getMaxSeqLen();



        //one element in AF matrix:
        // - 2 bytes for element of A matrix
        // - 2 bytes for element of F matrix
        short2* AF;
        if (AFs.find(gpuNo) == AFs.end()) {
            cudaMalloc(&AF, sizeof(int) * maxSeqLength * params->windowSize * params->windowSize);
            AFs[gpuNo] = AF;
        }
        AF = AFs[gpuNo];
        //sizeof(int) - one element in A matrix:
        // - 2 bytes for element of A matrix
        // - 2 bytes for element of F matrix
        int* scoresDevPtr;
        if (scoresDevPtrs.find(gpuNo) == scoresDevPtrs.end()) {
            cudaMalloc(&scoresDevPtr, sizeof(int) * params->windowSize * params->windowSize);
            scoresDevPtrs[gpuNo] = scoresDevPtr;
        }
        scoresDevPtr = scoresDevPtrs[gpuNo];
        
        //sizoeof(back) => sizeof(int) * maxSeqLength * (maxSeqLength/8) * windowSize * windowSize
        //(height+8) * (maxSeqLength+1) -> additional -1 row and -1 column in "back" array
        int height = ((maxSeqLength-1) / 8 + 1) * 8; //8->8, 9->16, 10->16 ...
        //hight of this array must be dividable by 8 (ALIGNMENT_MATCH_Y_STEPS)
        unsigned int* back;
        unsigned int  backSize = sizeof(unsigned int) * (height+8) * (maxSeqLength+1) * params->windowSize * (params->windowSize/ALIGNMENT_MATCH_Y_STEPS);
        cudaMalloc(&back, backSize);
        //texture can't be too big -> since CUDA 3.0 the line below causes errors
        //cudaBindTexture(0, texBack, back, backSize);
        
        //memory for temporary (intermediate) results (alignments/matches)
        //we need: 2x maxSeqLength*2*windowSize*windowSize
        int maxSeqLengthAlignedTo4 = ((maxSeqLength - 1) / 4 + 1) * 4;
        unsigned int* matchesSeqXDevPtr; //this array will have 4 characters packed in one int
        unsigned int* matchesSeqYDevPtr;
        cudaMalloc(&matchesSeqXDevPtr, sizeof(char) * maxSeqLengthAlignedTo4 * 2 * params->windowSize * params->windowSize);
        cudaMalloc(&matchesSeqYDevPtr, sizeof(char) * maxSeqLengthAlignedTo4 * 2 * params->windowSize * params->windowSize);

        //memory for final results (alignments/matches)
        unsigned int* outMatchesSeqXDevPtr; //this array will have 4 characters packed in one int
        unsigned int* outMatchesSeqYDevPtr;
        cudaMalloc(&outMatchesSeqXDevPtr, sizeof(char) * maxSeqLengthAlignedTo4 * 2 * params->windowSize * params->windowSize);
        cudaMalloc(&outMatchesSeqYDevPtr, sizeof(char) * maxSeqLengthAlignedTo4 * 2 * params->windowSize * params->windowSize);


        //copying input sequences to texture memory
        TexVariablesAddresses addr = copySeqsToTex(params->seqs, params->startSeqs1No, params->startSeqs2No, params->windowSize);

        

        /***********************************************************************
         * KERNEL 1                                                            *
         * score calculation and "back" matrix fill                            *
         ***********************************************************************/

        dim3 blockShape(params->blockShape,params->blockShape);
        dim3 gridShape((params->windowSize-1)/params->blockShape + 1,(params->windowSize-1)/params->blockShape +1);

        HiResTimer timer1;
        HiResTimer timer2;

        KernelMatchInvokingParams* kernelParams = new KernelMatchInvokingParams();
        kernelParams->blockShape = blockShape;
        kernelParams->gridShape = gridShape;
        kernelParams->AF = AF;
        kernelParams->back = back;
        kernelParams->scores = scoresDevPtr;
        //maxSeqLength+1 => +1 because we have to take into account the -1 column
        kernelParams->rowWidth = maxSeqLength+1;
        kernelParams->border = (params->windowX == params->windowY);

        timer1.start();
        kernelInvocation(kernelParams);
        cudaThreadSynchronize();
        timer1.stop();


        /***********************************************************************
         * KERNER 2                                                            *
         * backtracing - alignment matches generation                          *
         ***********************************************************************/

        //short2 x;
        //cudaMemcpy(&x, AF, sizeof(short2), cudaMemcpyDeviceToHost);
        //printf("xxx%d %d %d %d\n", x.x, params->seqs->getLengths()[params->windowX*params->windowSize], x.y, params->seqs->getLengths()[params->windowY*params->windowSize]);


        kernelParams->matchesX = matchesSeqXDevPtr;
        kernelParams->matchesY = matchesSeqYDevPtr;
        timer2.start();
        kernelBacktraceInvocation(kernelParams);
        cudaThreadSynchronize();
        timer2.stop();
        
        printf("Kernel[%5d] %5dms %5dms\n", params->partId, (int)timer1.getElapsedTime(), (int)timer2.getElapsedTime()); 

        //reading scores from GPU
        int* scoresHostPtr = new int[params->windowSize * params->windowSize];
        cudaMemcpy(scoresHostPtr, scoresDevPtr, sizeof(int)*params->windowSize*params->windowSize, cudaMemcpyDeviceToHost);
        

        /***********************************************************************
         * KERNER 3                                                            *
         * changing order of the results in GPU memory                         *
         ***********************************************************************/
        
        dim3 blockShape2(params->blockShape, params->blockShape);
        dim3 gridShape2( (params->windowSize*params->windowSize) / params->blockShape);
        
        /***********************************************************************
         * maxSeqLengthAlignedTo4*2/4                                          *
         *    -> *2 because alignment can be 2x as long as the longest         *
         *      sequence                                                       *
         *    -> /4 because we packed chars to int                             *
         ***********************************************************************/
        reorderMatchesMemory<<<gridShape2, blockShape2>>>(matchesSeqXDevPtr, outMatchesSeqXDevPtr, maxSeqLengthAlignedTo4*2/4);
        reorderMatchesMemory<<<gridShape2, blockShape2>>>(matchesSeqYDevPtr, outMatchesSeqYDevPtr, maxSeqLengthAlignedTo4*2/4);

        cudaMemcpy(params->matches1Manager->getWindow(params->windowX, params->windowY), outMatchesSeqXDevPtr, sizeof(char) * maxSeqLengthAlignedTo4 * 2 * params->windowSize * params->windowSize, cudaMemcpyDeviceToHost);
        cudaMemcpy(params->matches2Manager->getWindow(params->windowX, params->windowY), outMatchesSeqYDevPtr, sizeof(char) * maxSeqLengthAlignedTo4 * 2 * params->windowSize * params->windowSize, cudaMemcpyDeviceToHost);



        /***********************************************************************
         * READING SCORES                                                      *
         ***********************************************************************/

        int minX = MIN(params->windowSize, params->seqs->getSequenceNumber() - params->windowX*params->windowSize);
        int minY = MIN(params->windowSize, params->seqs->getSequenceNumber() - params->windowY*params->windowSize);
        for(int i = 0; i<minY; i++)//Y
            for(int j = 0; j<minX; j++)//X
                params->scores[params->windowY*params->windowSize + i][params->windowX*params->windowSize + j] = scoresHostPtr[i*params->windowSize + j];

        /***********************************************************************
         * FINISHING                                                           *
         ***********************************************************************/

        //dealocating memory on GPU
        //cudaFree(AF);
        cudaFree(back);
        cudaFree(matchesSeqXDevPtr);
        cudaFree(matchesSeqYDevPtr);
        cudaFree(outMatchesSeqXDevPtr);
        cudaFree(outMatchesSeqYDevPtr);
        //cudaFree(scoresDevPtr);
        cudaFree(addr.texSeqs1DevPtr);
        cudaFree(addr.texSeqs2DevPtr);

        checkError("errors in NeedlemanWunschGlobalMatchKernel");

        //dealocation memory on CPU
        delete[] scoresHostPtr;
        delete params;
        delete kernelParams;
    }
    catch (Exception* ex)
    {
        printf("%s\n",ex->getMessage());
    }

}

#endif
