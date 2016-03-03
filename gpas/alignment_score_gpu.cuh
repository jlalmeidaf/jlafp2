#ifndef _ALIGNMENT_SCORE_GPU_CUH_
#define _ALIGNMENT_SCORE_GPU_CUH_

#include "data_management.cuh"
#include "hi_res_timer.h"
#include "thread_manager.h"
#include "main_cu.h"

using std::string;
using namespace Exceptions;


/*******************************************************************************
 * BLOCK_SHAPE defines the width and the height of the one block of threads.   *
 * Each thread corresponds to one pair of sequences and calculates one         *
 * one alignment. It's a good practice to have it multiple of 8. But the best  *
 * results should give 16.                                                     *
 *******************************************************************************/

#define ALIGNMENT_SCORE_BLOCK_SHAPE    16

/*******************************************************************************
 * BLOCK_SIZE defines the count of the threads runs in the multiprocessor in   *
 * the one block. We will threat it as a square and so it should have a value: *
 * BLOCK_SIZE = BLOCK_SHAPE^2                                                  *
 *******************************************************************************/

#define ALIGNMENT_SCORE_BLOCK_SIZE    (ALIGNMENT_SCORE_BLOCK_SHAPE * ALIGNMENT_SCORE_BLOCK_SHAPE)

/*******************************************************************************
 * Y_STEPS defines the number of steps within the sequence Y which will be     *
 * computed in one iteration without an access to the global memory.           *
 * Maximum Y_STEPS value is 12 (as long as size of the shared memory size is   *
 * 16KB per multiprocessor - see CUDA documentation).                          *
 *******************************************************************************/

#define ALIGNMENT_SCORE_Y_STEPS        12

/*******************************************************************************
 * AlignmentInvokerParams::getEstimatedCompexity() implements the              *
 * functionality of estimation of the complexity of the problem contained in   *
 * the, corresponding to the AlignmentInvokerParams window.                    *
 * Good approximation of the value is a product of the max length on X         *
 * position and on Y position.                                                 *
 *                                                                             *
 * Cpx = max[L(X)]*max[L(Y)]                                                   *
 *                                                                             *
 * Where:                                                                      *
 * L(n) is the length of sequence n                                            *
 *                                                                             *
 *******************************************************************************/

long long AlignmentScoreKernelInvokerParams::getEstimatedComplexity()
{
    long long result;
    result = (long long) seqs->getWindowSum(windowSize, windowX, blockShape);
    result *= (long long) seqs->getWindowSum(windowSize, windowY, blockShape);

    //        if (windowX == windowY)
    //        {
    //            result = (result + seqs->getSquaredSum(windowSize, windowX))/2;
    //        }

    int multiprocessorCorrection;


    if (partId == 1)
    {
        multiprocessorCorrection = seqs->getSequenceNumber() % windowSize;
        multiprocessorCorrection += blockShape - 1;
        multiprocessorCorrection /= blockShape;
        multiprocessorCorrection *= multiprocessorCorrection;

        if ((multiprocessorCorrection < maxMultiprocessorCount) && (multiprocessorCorrection != 0))
        {
            result *= maxMultiprocessorCount;
            result /= multiprocessorCorrection;
        }
    }
    //        else if (seqs->getSequenceNumber() / windowSize == windowY)
    //        {
    //            multiprocessorCorrection = seqs->getSequenceNumber() % windowSize;
    //            multiprocessorCorrection += blockShape - 1;
    //            multiprocessorCorrection /= blockShape;
    //            multiprocessorCorrection *= windowSize / blockShape;
    //            if (multiprocessorCorrection < maxMultiprocessorCount)
    //            {
    //                result *= maxMultiprocessorCount;
    //                result /= multiprocessorCorrection;
    //            }
    //        }

    //        multiprocessorCorrection = seqs->getSequenceNumber();
    //        multiprocessorCorrection += windowSize - 1;
    //        multiprocessorCorrection /= windowSize;
    //        multiprocessorCorrection *= multiprocessorCorrection;
    //        multiprocessorCorrection /= maxMultiprocessorCount;

    result = MAX(result, seqs->getWindowMax(windowSize, windowX) * seqs->getWindowMax(windowSize, windowY) * maxMultiprocessorCount); //16 multiprocessors count

    return result;
}

/*******************************************************************************
 * compareAlignmentScoreKernelInvokerParams compares this datatype by the      *
 * estimated complexity.                                                       *
 *******************************************************************************/

int AlignmentScoreKernelInvokerParams::compareAlignmentScoreKernelParams(const void* first, const void* second)
{
    AlignmentScoreKernelInvokerParams* sFirst = *((AlignmentScoreKernelInvokerParams**) first);
    AlignmentScoreKernelInvokerParams* sSecond = *((AlignmentScoreKernelInvokerParams**) second);
    if (sFirst->getEstimatedComplexity() < sSecond->getEstimatedComplexity())
        return 1;
    if (sFirst->getEstimatedComplexity() == sSecond->getEstimatedComplexity())
        return 0;
    return -1;
}

/*******************************************************************************
 *
 *******************************************************************************/

void AlignmentScoreAlgorithm::run(ThreadManager* tm, Sequences* seqs, int gapOpen, int gapExt, int maxMultiprocessorCount, unsigned int windowSize, const char* outfile)
{
    //HERE TASKS ARE DEFINED


    HiResTimer timer;
    timer.start();

    int gpus = tm->getThreadsNumber();
    smParams.sm = seqs->getSubtitutionMatrix();
    smParams.gapOpen = gapOpen;
    smParams.gapExtension = gapExt;
    smParams.windowSize = windowSize;
    smParams.memoryOffset = windowSize * windowSize;

    for (int i = 0; i < gpus; i++)
    {
        tm->request(copySMToConstInThread, (void*) & smParams, tm->threadsIDs[i]);
    }

    int sequenceNumber = seqs->getSequenceNumber();
    //[windowsNumber*(windowsNumber + 1)] / 2 - the number of data parts to be processed
    int windowsNumber = (sequenceNumber - 1) / windowSize + 1;

    //here we allocate memory for results
    scores = new int*[sequenceNumber];
    scores[0] = new int[sequenceNumber * sequenceNumber];
    for (int i = 1; i < sequenceNumber; i++)
        scores[i] = &scores[0][i * sequenceNumber];
    //scores[Y][X]

    int partId = 1;
    int parts = windowsNumber * (windowsNumber + 1) / 2;
    AlignmentScoreKernelInvokerParams* params;

    AlignmentScoreKernelInvokerParams** jobs = new AlignmentScoreKernelInvokerParams*[parts];


    //long jobs first and short jobs last is better for load balancing
    for (int j = windowsNumber - 1; j >= 0; j--) //we iterate through all the windows
    {
        for (int i = j; i >= 0; i--)
        {
            params = getParams();
            params->windowX = i;
            params->windowY = j;
            params->partId = partId++;
            params->parts = parts;
            params->seqs = seqs;
            params->scores = scores;
            params->maxMultiprocessorCount = maxMultiprocessorCount;
            params->windowSize = windowSize;

            jobs[params->partId - 1] = params;
            //tm->request(NeedlemanWunschGlobalScoreKernelInvoker, (void*)params, -1);
        }
    }

    qsort(jobs, parts, sizeof (AlignmentScoreKernelInvokerParams*), AlignmentScoreKernelInvokerParams::compareAlignmentScoreKernelParams);

    for (int i = 0; i < parts; i++)
    {
        //printf("%d: %lld %d %d\n", jobs[i]->partId, jobs[i]->getEstimatedComplexity(), jobs[i]->seqs->getWindowSum(jobs[i]->windowSize, jobs[i]->windowX), jobs[i]->seqs->getWindowSum(jobs[i]->windowSize, jobs[i]->windowY));//
        InvokingParams* invokingParams = new InvokingParams();
        invokingParams->algorithm = this;
        invokingParams->params = jobs[i];
        tm->request(invoker, (void*)invokingParams, -1);
    }

    tm->wait();

    timer.stop();
    printf("%s total: %dms\n", getAlgorithmName(), (int) timer.getElapsedTime());


    //writing results to a file
    if(outfile != NULL)
    {
        FILE* file = fopen(outfile, "w");

        char* seqBuffer = new char[seqs->getMaxSeqLen() + 1];
        fprintf(file, "File produced by gpu-pairAlign\n");
        fprintf(file, "%d\n", sequenceNumber); // number of sequences
        for (int i = 0; i < sequenceNumber; i++)
        {
            // retrieving original input sequences
            for (int a = 0; a < (seqs->getLengths()[i]); a++)
            {
                seqBuffer[a] = seqs->getSubtitutionMatrix()->getRevConv()
                              [ seqs->getSequences()
                              [ seqs->getStarts()[i] + seqs->getLengths()[i] - a - 1] ]
                              - 'A' + 'a'; // to convert to small letters
            }
            seqBuffer[seqs->getLengths()[i]] = 0;

            // TODO: nazwy sekwencji
            fprintf(file, "%3d %s %d %s\n", i, seqs->getSeqName(i), seqs->getLengths()[i], seqBuffer);
        }
        delete[] seqBuffer;
        fprintf(file, "\n");


        for (int i = 0; i < sequenceNumber; i++)//Y
        {
            for (int j = 0; j <= i; j++)//X
            {
                fprintf(file, "#%d %d, score: %d\n", i, j, scores[i][j]);
                if(calculateMatches)
                {
                    fprintf(file, "%s\n",   ((AlignmentMatchAlgorithm*)this)->matches1Manager->getSequence(i, j));
                    fprintf(file, "%s\n\n", ((AlignmentMatchAlgorithm*)this)->matches2Manager->getSequence(i, j));
                }
            }
        }
        fclose(file);
    }

}

AlignmentScoreKernelInvokerParams* AlignmentScoreAlgorithm::getParams()
{
    AlignmentScoreKernelInvokerParams* result = new AlignmentScoreKernelInvokerParams();
    //default parameters' values (for AlignmentScoreAlgorithm class)
    result->blockShape = ALIGNMENT_SCORE_BLOCK_SHAPE;
    return result;
}

AlignmentScoreAlgorithm::AlignmentScoreAlgorithm()
{
    calculateMatches = false;
}

AlignmentScoreAlgorithm::~AlignmentScoreAlgorithm()
{
    delete[] scores[0];
    delete[] scores;
}

void AlignmentScoreAlgorithm::invoker(ThreadManager* tm, void* data)
{
    InvokingParams* iParams = (InvokingParams*)data;
    iParams->algorithm->actualInvokedMethod(iParams->params, tm->getThreadsInfo().gpuNo);
    delete iParams;
}

void AlignmentScoreAlgorithm::actualInvokedMethod(AlignmentScoreKernelInvokerParams* params, int gpuNo)
{
    //THIS FUNCTION IS CALLED BY THREAD MANAGER
    
    params->startSeqs1No = params->windowX * params->windowSize;
    params->startSeqs2No = params->windowY * params->windowSize;
    try
    {
        int maxSeqLength = params->seqs->getMaxSeqLen();


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
        //printf("Copying sequences part %d/%d ... ", params->partId, params->parts);



        //copying sequences to texture memory
        TexVariablesAddresses addr = copySeqsToTex(params->seqs, params->startSeqs1No, params->startSeqs2No, params->windowSize);
        //printf("OK\n");


        HiResTimer timer;
        timer.start();
        //KERNEL INVOCATION
        dim3 blockShape(params->blockShape,params->blockShape);
        dim3 gridShape((params->windowSize-1)/params->blockShape + 1,(params->windowSize-1)/params->blockShape +1);
//        NeedlemanWunschGlobalScoreKernel<<<gridShape,blockShape>>>(AF, scoresDevPtr, params->windowX == params->windowY);
        KernelScoreInvokingParams* kernelParams = new KernelScoreInvokingParams();
        kernelParams->gridShape = gridShape;
        kernelParams->blockShape = blockShape;
        kernelParams->AF = AF;
        kernelParams->scores = scoresDevPtr;
        kernelParams->border = (params->windowX == params->windowY);
        kernelInvocation((void*) kernelParams);


        cudaThreadSynchronize();
        timer.stop();
        printf("Kernel [%5d] %dms\n", params->partId, (int)timer.getElapsedTime());

        //checkes for the errors accuring in the kernel
        //(doesn't really work, but what a heck)
        string s;
        s = "errors in ";
        s += getAlgorithmName();
        checkError(s.c_str());
        //printf("%d\n", params->partId);

        //reading results from GPU
        int* scoresHostPtr = new int[params->windowSize * params->windowSize];
        cudaMemcpy(scoresHostPtr, scoresDevPtr, sizeof(int)*params->windowSize*params->windowSize, cudaMemcpyDeviceToHost);



        int minX = MIN(params->windowSize, params->seqs->getSequenceNumber() - params->windowX*params->windowSize);
        int minY = MIN(params->windowSize, params->seqs->getSequenceNumber() - params->windowY*params->windowSize);
        for(int i = 0; i<minY; i++)//Y
            for(int j = 0; j<minX; j++)//X
                params->scores[params->windowY*params->windowSize + i][params->windowX*params->windowSize + j] = scoresHostPtr[i*params->windowSize + j];


        //dealocating memory on GPU
        //cudaFree(AF);
        //cudaFree(scoresDevPtr);
        cudaFree(addr.texSeqs1DevPtr);
        cudaFree(addr.texSeqs2DevPtr);

        delete kernelParams;
    }
    catch (Exception* ex)
    {
        printf("%s\n",ex->getMessage());
    }
}


#endif
