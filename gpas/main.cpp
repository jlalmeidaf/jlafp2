#include <stdio.h>
#include <unistd.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <fstream>
#include "sequences.h"
#include "substitution_matrix.h"
#include "needleman_wunsch_global.h"
#include "needleman_wunsch_semiglobal.h"
#include "smith_waterman_local.h"
#include "main_cu.h"
#include "thread_manager.h"
#include "arguments_manager.h"

using Exceptions::Exception;
using namespace Data;
using namespace Algorithms::Sequential;
using std::map;
using namespace std;

//GLOBAL DATA
char gapOpen = 5;
char gapExtension = 2;
int maxGpu = 128;
unsigned int windowSize = 64;
char* submatrix = (char*)"data/substitution_matrices/EBLOSUM30";
char* infile = (char*)"data/test.fax";
char* outfile = NULL;
char* algName = (char*)"sws";//Smith Waterman Score
int seqX = -1;
int seqY = -1;

int maxDevMultiprocessor = 1;
int gpuNumber;
int* gpus;
int* desiredGpus = NULL;

bool deviceQuery()
{
    cudaGetDeviceCount(&gpuNumber);
    printf("Device count: %d\n", gpuNumber);

    if(gpuNumber == 0)
        return false;
    
    gpus = new int[gpuNumber];
    if(desiredGpus)
    {
        printf("Using devices(s): ");

        int i=0;
        int selectedGpuNumber = 0;
        while(desiredGpus[i] > -1)
        {
            if(desiredGpus[i] < gpuNumber)
            {
                gpus[selectedGpuNumber++] = desiredGpus[i];
                printf("%d ", desiredGpus[i]);
                if(selectedGpuNumber == gpuNumber)
                    break; // no more gpus available on the system
            }
            i++;
        }
        printf("\n");
        gpuNumber = selectedGpuNumber; // not all gpus have to be used
    }
    else
    {
        if(gpuNumber > maxGpu)
        {
            gpuNumber = maxGpu;
            printf("Using %d first device(s)\n", gpuNumber);
        }
        for(int i=0; i<gpuNumber; i++)
            gpus[i] = i;
 
    }
    if(gpuNumber == 0)
        return false;

    cudaDeviceProp devProp;
    for (int i = 0; i < gpuNumber; i++)
    {
        cudaGetDeviceProperties(&devProp, 0);
        maxDevMultiprocessor = MAX(maxDevMultiprocessor, devProp.multiProcessorCount);
    }

    return true;
}

void initCUDA(ThreadManager* This, void* data)
{
    //cudaSetDevice(This->threadsInfos[pthread_self()].gpuNo);
    int gpuNo = This->getThreadsInfo().gpuNo;
    cudaSetDevice(gpuNo);

    char* devPtr;
    cudaMalloc((void**)&devPtr, sizeof(int));
    cudaMemcpy(devPtr, &gpuNo, sizeof(int), cudaMemcpyHostToDevice);

    gpuNo = -1;
    cudaMemcpy(&gpuNo, devPtr, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(devPtr);
    
    printf("GPU: %d\n",gpuNo);
}

bool readArguments(int argc, char** argv)
{
    char* argBuf;
    ArgumentsManager arguments(argc, argv);
    
    argBuf = (char*)arguments.getParam("help", "h");
    if (argBuf != NULL)
    {
        ifstream fin("HELP");
        string temp;
        while(getline(fin, temp))
            cout << temp << endl;
        fin.close();
        return false;
    }
    argBuf = (char*)arguments.getParam("version", "v");
    if (argBuf != NULL) {
        cout << "G-PAS (GPU-based Pairwise Alignment Software) 2.0" << endl;
        return false;
    }

    argBuf = (char*)arguments.getParam("infile", "i");
    if( (argBuf != NULL) && (strcmp(argBuf,"true")) )
        infile = strdup(argBuf);

    argBuf = (char*)arguments.getParam("outfile", "o");
    if( (argBuf != NULL) && (strcmp(argBuf,"true")) )
        outfile = strdup(argBuf);

    argBuf = (char*)arguments.getParam("submatrix", "sm");
    if( (argBuf != NULL) && (strcmp(argBuf,"true")) )
        submatrix = strdup(argBuf);
    
    argBuf = (char*)arguments.getParam("winsize", "ws");
    if ((argBuf != NULL) && (atoi(argBuf)>0) && (atoi(argBuf)%16 == 0))
        windowSize = atoi(argBuf);

    argBuf = (char*)arguments.getParam("gapopen", "go");
    if ((argBuf != NULL) && (atoi(argBuf)>0))
        gapOpen = atoi(argBuf);

    argBuf = (char*)arguments.getParam("gapext", "ge");
    if ((argBuf != NULL) && (atoi(argBuf)>0))
        gapExtension = atoi(argBuf);

    argBuf = (char*)arguments.getParam("alg", "a");
    if( (argBuf != NULL) && (strcmp(argBuf,"true")) )
        algName = strdup(argBuf);

    argBuf = (char*)arguments.getParam("ngpus", "ng");
    if ((argBuf != NULL) && (atoi(argBuf) > 0))
        maxGpu = atoi(argBuf);

    argBuf = (char*)arguments.getParam("seqX", "x");
    if ((argBuf != NULL) && (atoi(argBuf) > 0))
        seqX = atoi(argBuf) - 1;

    argBuf = (char*)arguments.getParam("seqY", "y");
    if ((argBuf != NULL) && (atoi(argBuf) > 0))
        seqY = atoi(argBuf) - 1;

    argBuf = (char*)arguments.getParam("gpus", "g");
    if( (argBuf != NULL) )
    {
        desiredGpus = new int[129];
        for(int i=0; i< 129; i++)
            desiredGpus[i] = -1;
        int i, a, idx;
        bool aNumber;
        aNumber = i = a = idx = 0;
        while(argBuf[i] != 0)
        {
            int c = (int)argBuf[i] - 48;
            if(c>=0 && c<=9)
            {
                a *= 10;
                a += c;
                aNumber=true;
            }
            else if (aNumber) // probably a separator like colon
            {
                desiredGpus[idx++] = a;
                a = aNumber = 0;
                if(idx>=127)
                    return false;
            }
            i++;
        }
        if(aNumber)
            desiredGpus[idx++] = a;
        //for(int k=0; k<idx; k++)
        //    printf("%d\n", desiredGpus[k]);
    }

    return true;
    
}

int main(int argc, char** argv)
{
    
    if(!readArguments(argc, argv))
        return 1;
    
    try
    {
//        SubstitutionMatrix sm(submatrix);
//        Sequences s(infile, &sm);
//        s.load();
//        s.sortSequences();

        //Homo_sapiens.GRCh37.55.pep.all.fa
        SubstitutionMatrix sm(submatrix);
        Sequences s(infile, &sm);
        s.load();
        s.sortSequences();

//        s.writeToFile("data/5_xxx.fax",100,420);
        
//        char fileName[512];
//        for(int i=16; i<1025; i+=16)
//        {
//            sprintf(fileName, "data/test/16-1024/5_%d.fax", i);
//            s.writeToFile(fileName, i);
//        }

//        char fileName[512];
//        for(int i=0; i<6; i++)
//        {
//            sprintf(fileName, "data/test/5_4000_%d.fax", i);
//            s.writeToFile(fileName, 4000);
//        }



//        s.writeToFile("data/3_xxx.fax", 100, 420);

//        NeedlemanWunschGlobal nwg;
//        nwg.Run(&s, 234, 660, gapOpen, gapExtension);
//        nwg.PrintResults("results/2.txt");
//
//        NeedlemanWunschSemiglobal nws;
//        nws.RunAll(&s, gapOpen, gapExtension);
//        nws.Run(&s, 201, 650, gapOpen, gapExtension);
//        nws.PrintResults("results/2.txt");
//
//        SmithWatermanLocal swl;
//        swl.Run(&s, 234, 660, gapOpen, gapExtension);
//        swl.PrintResults("results/2.txt");

        if(!deviceQuery())
            return 1;

        ThreadManager* tm;
        tm = new ThreadManager(gpus,gpuNumber);
        for(int i=0; i<gpuNumber; i++)
            tm->request(initCUDA,NULL, gpus[i]);
        tm->wait();

        AlignmentScoreAlgorithm* alg = NULL;
        if(!strcmp(algName, "nwgs")) alg = new NeedlemanWunschGlobalScoreGpu();
        if(!strcmp(algName, "nwss")) alg = new NeedlemanWunschSemiglobalScoreGpu();
        if(!strcmp(algName, "sws")) alg = new SmithWatermanScoreGpu();

        if(!strcmp(algName, "nwg")) alg = new NeedlemanWunschGlobalMatchGpu();
        if(!strcmp(algName, "nws")) alg = new NeedlemanWunschSemiglobalMatchGpu();
        if(!strcmp(algName, "sw")) alg = new SmithWatermanMatchGpu();
        if(alg != NULL) alg->run(tm, &s, 5, 2, maxDevMultiprocessor, windowSize, outfile);

        if((seqX>=0) && (seqY>=0))
        {
            printf("%s\n", ((AlignmentMatchAlgorithm*)alg)->matches1Manager->getSequence(seqX, seqY));
            printf("%s\n", ((AlignmentMatchAlgorithm*)alg)->matches2Manager->getSequence(seqX, seqY));
        }
        delete alg;
      
    }
    catch (Exception* ex)
    {
        printf("%s\n",ex->getMessage());
        return 1;
    }
    return 0;
}
