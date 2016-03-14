#ifndef _MAIN_CU_H_
#define _MAIN_CU_H_

#include "sequences.h"
#include "thread_manager.h"
#include "matches_manager.h"

struct TexVariablesAddresses; //from "data_management.cuh" file

using Data::Sequences;
using Data::SubstitutionMatrix;
using Data::MatchesManager;

/*******************************************************************************
 * AlignmentScoreKernelInvokerParams is the type of the second parameter of    *
 * the functions that implements alignment algorithm calculating only scalar   *
 * score result. (e.g. NeedlemanWunschGlobalScoreGpu).                         *
 *******************************************************************************/

class AlignmentScoreKernelInvokerParams
{
public:
    Sequences* seqs;
    int startSeqs1No, startSeqs2No, windowSize;
    int blockShape;
    int windowX, windowY, partId, parts;
    int** scores;
    int maxMultiprocessorCount;
    virtual long long getEstimatedComplexity();

    static int compareAlignmentScoreKernelParams(const void* first, const void* second);
};

/*******************************************************************************
 * AlignmentMatchKernelInvokerParams is the type of the second parameter of    *
 * the functions that implements algorithm calculating only the   the match    *
 * and the score result of the alignment. (e.g.                                *
 * NeedlemanWunschGlobalMatchGpu).                                             *
 *******************************************************************************/

class AlignmentMatchKernelInvokerParams : public AlignmentScoreKernelInvokerParams
{
public:
    MatchesManager* matches1Manager;
    MatchesManager* matches2Manager;
};



struct KernelScoreInvokingParams
{
    short2* AF;
    int* scores;
    bool border;
    dim3 gridShape,blockShape;
};


struct KernelMatchInvokingParams : public KernelScoreInvokingParams
{
    unsigned int* back;
    short rowWidth;
    unsigned int* matchesX;
    unsigned int* matchesY;
};


/*******************************************************************************
 *
 *******************************************************************************/


class AlignmentScoreAlgorithm
{
public:
    virtual void run(ThreadManager* tm, Sequences* seqs, int gapOpen, int gapExt, int maxMultiprocessorCount, unsigned int windowSize, const char* outfile);
    virtual AlignmentScoreKernelInvokerParams* getParams();
    virtual const char* getAlgorithmName() = 0;
    AlignmentScoreAlgorithm();
    ~AlignmentScoreAlgorithm();
    static void invoker(ThreadManager* tm, void* data);
    virtual void actualInvokedMethod(AlignmentScoreKernelInvokerParams* params, int gpuNo);
    virtual void kernelInvocation(void* kernelParams) = 0;

    int** scores;
    bool calculateMatches;
    map<int, short2*> AFs;
    map<int, int*> scoresDevPtrs;
};

struct InvokingParams
{
    AlignmentScoreAlgorithm* algorithm;
    AlignmentScoreKernelInvokerParams* params;
};



class NeedlemanWunschGlobalScoreGpu : public AlignmentScoreAlgorithm
{
public:
    virtual const char* getAlgorithmName();
    virtual void kernelInvocation(void* kernelParams);
};

class NeedlemanWunschSemiglobalScoreGpu : public AlignmentScoreAlgorithm
{
public:
    virtual const char* getAlgorithmName();
    virtual void kernelInvocation(void* kernelParams);
};

class SmithWatermanScoreGpu : public AlignmentScoreAlgorithm
{
public:
    virtual const char* getAlgorithmName();
    virtual void kernelInvocation(void* kernelParams);
};

/*******************************************************************************
 * ALIGNMENT WITH BACKTRACKING                                                 *
 *******************************************************************************/

class AlignmentMatchAlgorithm : public AlignmentScoreAlgorithm
{
public:
    virtual void run(ThreadManager* tm, Sequences* seqs, int gapOpen, int gapExt, int maxMultiprocessorCount, unsigned int windowSize, const char* outfile);
    virtual AlignmentScoreKernelInvokerParams* getParams();
    ~AlignmentMatchAlgorithm();
    virtual void actualInvokedMethod(AlignmentScoreKernelInvokerParams* params, int gpuNo);
    virtual void kernelBacktraceInvocation(void* kernelParams) = 0;

    MatchesManager* matches1Manager;
    MatchesManager* matches2Manager;
};

class NeedlemanWunschGlobalMatchGpu : public AlignmentMatchAlgorithm
{
public:
    virtual const char* getAlgorithmName();
    virtual void kernelInvocation(void* kernelParams);
    virtual void kernelBacktraceInvocation(void* kernelParams);
};

class NeedlemanWunschSemiglobalMatchGpu : public AlignmentMatchAlgorithm
{
public:
    virtual const char* getAlgorithmName();
    virtual void kernelInvocation(void* kernelParams);
    virtual void kernelBacktraceInvocation(void* kernelParams);
};

class SmithWatermanMatchGpu : public AlignmentMatchAlgorithm
{
public:
    virtual const char* getAlgorithmName();
    virtual void kernelInvocation(void* kernelParams);
    virtual void kernelBacktraceInvocation(void* kernelParams);
};

//extern "C" void NeedlemanWunschGlobalMatchGpu(ThreadManager* tm, Sequences* seqs, int gapOpen, int gapExt, int maxMultiprocessorCount);

extern "C" TexVariablesAddresses copySeqsToTex(Sequences* s, int startSeqs1No, int startSeqs2No, int count);

#endif
