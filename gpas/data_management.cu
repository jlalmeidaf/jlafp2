#include "data_management.cuh"
#include "sequences.h"
#include <cuda.h>

#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <math.h>
#include "thread_manager.h"
#include "main_cu.h"

using Exceptions::IndexOutOfRangeException;
using Exceptions::Exception;
using Data::Sequences;
using Data::SubstitutionMatrix;

__global__ void test(char* tab)
{
    tab[0] = substitutionMatrix[0];
    tab[1] = gapOp;
    tab[2] = gapEx;
}

__global__ void test2(char* tab)
{
    tab[0] = tex1Dfetch(texSeqsY, 1);
}

__global__ void test3(int* tab)
{
//    tab[0] = signbit(0);
//    tab[1] = signbit(1);
//    tab[2] = signbit(-1);

//    tab[0] = PACK(-23, 48);
//    tab[1] = UNPACK_FROM_LEFT(tab[0]);
//    tab[2] = UNPACK_FROM_RIGHT(tab[0]);

    short xx = 1234;
    int yy = 1234;
    tab[0] = (1==0);
    tab[1] = (1432000000==4320) || 2;
    tab[2] = (xx==yy);
}


void copySMToConst(SubstitutionMatrix* sm, char gapOpen, char gapExtension, unsigned int windowSize, unsigned int memoryOffset)
{
    if(sm->getMatrixSize() > SUBSTITUTION_MATRIX_SIZE)
        throw new IndexOutOfRangeException("Too big substitution matrix to copy to constant memory.\n");



    int lettersNumber = sm->getLettersNumber();

    cudaMemcpyToSymbol(substitutionMatrix, sm->getMatrix(), sm->getMatrixSize());
    cudaMemcpyToSymbol(lettersCount, &lettersNumber, sizeof(int));
    cudaMemcpyToSymbol(revConvConst, sm->getRevConv(), MAX_LETTERS_COUNT);
    cudaMemcpyToSymbol(gapOp, &gapOpen, sizeof(char));
    cudaMemcpyToSymbol(gapEx, &gapExtension, sizeof(char));
    cudaMemcpyToSymbol(winSize, &windowSize, sizeof(unsigned int));
    cudaMemcpyToSymbol(memOffset, &memoryOffset, sizeof(unsigned int));



    char host[] = {0,0,0,0};
    char* devPtr;
    cudaMalloc((void**)&devPtr, 4);
    test<<<1,1>>>(devPtr);
    cudaMemcpy(host, devPtr, 4, cudaMemcpyDeviceToHost);
    cudaFree(devPtr);

    if ((host[0] == sm->getMatrix()[0]) && (host[1] == gapOpen) && (host[2]==gapExtension))
    {
        //here everything is OK
    }
    else
    {
        printf("failed\n");
        printf("%c%4d%4d\n",host[0], host[1], host[2]);
        printf("%c%4d%4d\n",sm->getMatrix()[0], gapOpen, gapExtension);
    }

}

void copySMToConstInThread(ThreadManager* This, void* data)
{
    try
    {
        printf("Declaring GPU(%d) constant memory\n", This->getThreadsInfo().gpuNo);

        SubstitutionMatrixParams* params = (SubstitutionMatrixParams*) data;
        copySMToConst(params->sm, params->gapOpen, params->gapExtension, params->windowSize, params->memoryOffset);
    }
    catch (Exception* ex)
    {
        printf("%s\n",ex->getMessage());
    }

}


TexVariablesAddresses copySeqsToTex(Sequences* s, int startSeqs1No, int startSeqs2No, int windowSize)
{
    if (windowSize > MAX_ALGORITHM_WINDOW_SIZE)
        throw new IndexOutOfRangeException("Window size greater than MAX_ALGORITHM_WINDOW_SIZE.\n");

    int* lengths1 = new int[windowSize];
    int* lengths2 = new int[windowSize];
    int* starts1 = new int[windowSize + 1];
    int* starts2 = new int[windowSize + 1];
    int size1; //the size of allocated memory for texSeqs1
    int size2; //the size of allocated memory for texSeqs2

    int* starts = s->getStarts();
    int* lengths = s->getLengths();


    int maxCount = startSeqs1No + windowSize;
    maxCount = MIN(maxCount, s->getSequenceNumber());
    maxCount = MAX(maxCount-startSeqs1No, 0);
    for (int i = 0; i < maxCount; i++)
    {
        lengths1[i] = lengths[i + startSeqs1No];
        starts1[i]  = starts[i + startSeqs1No] - starts[startSeqs1No];//we substract the offset to make
                                                                      //starts1[0] == 0
    }
    for (int i = maxCount; i <= windowSize; i++)
    {
        lengths1[i] = 0;
        starts1[i] = ((maxCount>0)?(starts1[maxCount - 1] + lengths1[maxCount - 1]):0);
    }


    maxCount = startSeqs2No + windowSize;
    maxCount = MIN(maxCount, s->getSequenceNumber());
    maxCount = MAX(maxCount-startSeqs2No,0);
    for (int i = 0; i < maxCount; i++)
    {
        lengths2[i] = lengths[i + startSeqs2No];
        starts2[i]  = starts[i + startSeqs2No] - starts[startSeqs2No];//we substract the offset to make
                                                                      //starts2[0] == 0
    }
    for (int i = maxCount; i <= windowSize; i++)
    {
        lengths2[i] = 0;
        starts2[i] = ((maxCount>0)?(starts2[maxCount - 1] + lengths2[maxCount - 1]):0);
    }


    //COPYING ARRAYS OF SEQUENCES STARTS INTO CONST MEMORY
    cudaMemcpyToSymbol(tex1Starts, starts1, (windowSize + 1)*sizeof(int));
    cudaMemcpyToSymbol(tex2Starts, starts2, (windowSize + 1)*sizeof(int));


    //COPYING X SEQUENCES TO TEXTURE
    char* firstSeqHost = s->getSequences();
    firstSeqHost += s->getStarts()[startSeqs1No];

    char* firstSeqDev;
    size1 = starts1[windowSize - 1] + lengths1[windowSize - 1];//length of all sequences within the window

    cudaMalloc((void**)&firstSeqDev, sizeof(char)*size1);
    cudaMemcpy(firstSeqDev, firstSeqHost, sizeof(char)*size1, cudaMemcpyHostToDevice);
    cudaBindTexture(0, texSeqsX, firstSeqDev, size1);


    //COPYING Y SEQUENCES TO TEXTURE
    char* secondSeqHost = s->getSequences();
    secondSeqHost += s->getStarts()[startSeqs2No];

    char* secondSeqDev;
    size2 = starts2[windowSize - 1] + lengths2[windowSize - 1];

    cudaMalloc((void**)&secondSeqDev, sizeof(char)*size2);
    cudaMemcpy(secondSeqDev, secondSeqHost, sizeof(char)*size2, cudaMemcpyHostToDevice);
    cudaBindTexture(0, texSeqsY, secondSeqDev, size2);

////TESTS:
//    char host[] = {0,0,0,0};
//    char* devPtr;
//    cudaMalloc((void**)&devPtr, 4);
//    test2<<<1,1>>>(devPtr);
//    cudaMemcpy(host, devPtr, 4, cudaMemcpyDeviceToHost);
//    printf("%c\n", s->getSubtitutionMatrix()->revConvert(host[0]));

//    int host[] = {0,0,0,0};
//    int* devPtr;
//    cudaMalloc((void**)&devPtr, 4*sizeof(int));
//    test3<<<1,1>>>(devPtr);
//    cudaMemcpy(host, devPtr, 4*sizeof(int), cudaMemcpyDeviceToHost);
//    printf("%4d %4d %4d\n", host[0], host[1], host[2]);


    delete[] lengths1;
    delete[] lengths2;
    delete[] starts1;
    delete[] starts2;


    TexVariablesAddresses result;
    result.texSeqs1DevPtr = firstSeqDev;
    result.texSeqs2DevPtr = secondSeqDev;

    return result;
}
