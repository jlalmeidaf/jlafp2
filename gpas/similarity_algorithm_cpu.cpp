#include "similarity_algorithm_cpu.h"

using namespace Algorithms::Sequential;
using Data::Sequences;
using Data::SubstitutionMatrix;

SimilarityAlgorithmCpu::SimilarityAlgorithmCpu()
{
    A = NULL;
    E = NULL;
    F = NULL;
    B = NULL;
    result1 = NULL;
    result2 = NULL;

    scores = NULL;
    matches1 = NULL;
    matches2 = NULL;
}

void SimilarityAlgorithmCpu::Run(Sequences* s, int seq1No, int seq2No, int gapOp, int gapEx)
{
    int* starts = s->getStarts();
    char* sequences = s->getSequences();
    int* lengths = s->getLengths();
    seq1 = &(sequences[starts[seq1No]]) - 1;
    seq2 = &(sequences[starts[seq2No]]) - 1;
    seq1Length = lengths[seq1No];
    seq2Length = lengths[seq2No];
    sm = s->getSubtitutionMatrix();

    this->gapOp = gapOp;
    this->gapEx = gapEx;

    DeallocateMemoryForSingleRun();
    AllocateMemoryForSingleRun();
    InitializeMatrices();
    FillMatrices();
    BackwardMoving();
}

void SimilarityAlgorithmCpu::AllocateMemoryForSingleRun()
{
    A = new int*[seq1Length + 1];
    E = new int*[seq1Length + 1]; //left matrix
    F = new int*[seq1Length + 1]; //up matrix
    B = new BackUpStruct*[seq1Length + 1];


    A[0] = new int[(seq1Length + 1) * (seq2Length + 1)];
    E[0] = new int[(seq1Length + 1) * (seq2Length + 1)];
    F[0] = new int[(seq1Length + 1) * (seq2Length + 1)];
    B[0] = new BackUpStruct[(seq1Length + 1) * (seq2Length + 1)];


    for (int i = 1; i < seq1Length + 1; i++)
    {
        A[i] = A[0] + (seq2Length + 1)*i;
        E[i] = E[0] + (seq2Length + 1)*i;
        F[i] = F[0] + (seq2Length + 1)*i;
        B[i] = B[0] + (seq2Length + 1)*i;
    }
    
    result1 = new char[seq1Length + seq2Length];
    result2 = new char[seq1Length + seq2Length];
}


void SimilarityAlgorithmCpu::RunAll(Sequences* s, int gapOp, int gapEx)
{
    int n = s->getSequenceNumber();
    int len = s->getMaxSeqLen()*2;

    DeallocateMemoryForAllRuns();

    matches1 = new MatchesManager(n, n, len);
    matches2 = new MatchesManager(n, n, len);
    scores = new int*[n];
    scores[0] = new int[n*n];
    for(int i=1; i<n; i++)
        scores[i] = &scores[0][i*n];

    HiResTimer timer;
    timer.start();

                           //  |
    for(int j=0; j<n; j++) // \|/
    {
        for(int i=0; i<=j; i++) // -->
        {
            char* winPtrX = matches1->getWindow(0, 0);
            char* winPtrY = matches2->getWindow(0, 0);
            winPtrX = winPtrX + len*(n*j + i);
            winPtrY = winPtrY + len*(n*j + i);

            Run(s, i, j, gapOp, gapEx);

            for(int k=0; result1[k]!=0; k++)
            {
                winPtrX[k] = result1[k];
                winPtrY[k] = result2[k];
                
            }
            scores[j][i] = score;
            //printf("%5d %5d\n", j, i);
        }
        printf("%5d\n", j);
    }

    timer.stop();
    printf("Sequential alg total time: %d\n", (int)timer.getElapsedTime());

}

void SimilarityAlgorithmCpu::PrintResults(const char* fileName)
{
    FILE* testF = fopen(fileName, "w");

    fprintf(testF, "MACIERZ A:\n\n");
    fprintf(testF, "       _");
    for (int i = 1; i <= seq2Length; i++)
        fprintf(testF, "   %c", sm->revConvert(seq2[i]));
    fprintf(testF, "\n");
    fprintf(testF, "   _");
    for (int j = 0; j <= seq2Length; j++)
    {
        fprintf(testF, "%4d", A[0][j]);
    }
    fprintf(testF, "\n");
    for (int i = 1; i <= seq1Length; i++)
    {
        fprintf(testF, "   %c", sm->revConvert(seq1[i]));
        for (int j = 0; j <= seq2Length; j++)
        {
            fprintf(testF, "%4d", A[i][j]);
        }
        fprintf(testF, "\n");
    }

//    fprintf(testF, "\n\n\n");
//
//    fprintf(testF, "MACIERZ kontynuacji left:\n\n");
//    fprintf(testF, "    ");
//    for (int i = 1; i <= seq2Length; i++)
//        fprintf(testF, "   %c", sm->revConvert(seq2[i]));
//    fprintf(testF, "\n");
//    for (int i = 1; i <= seq1Length; i++)
//    {
//        fprintf(testF, "   %c", sm->revConvert(seq1[i]));
//        for (int j = 1; j <= seq2Length; j++)
//        {
//            fprintf(testF, "   %c", (B[i][j].continueLeft)?('t'):('n'));
//        }
//        fprintf(testF, "\n");
//    }
//
//    fprintf(testF, "\n\n\n");
//
//    fprintf(testF, "MACIERZ kontynuacji up:\n\n");
//    fprintf(testF, "    ");
//    for (int i = 1; i <= seq2Length; i++)
//        fprintf(testF, "   %c", sm->revConvert(seq2[i]));
//    fprintf(testF, "\n");
//    for (int i = 1; i <= seq1Length; i++)
//    {
//        fprintf(testF, "   %c", sm->revConvert(seq1[i]));
//        for (int j = 1; j <= seq2Length; j++)
//        {
//            fprintf(testF, "   %c", (B[i][j].continueUp)?('t'):('n'));
//        }
//        fprintf(testF, "\n");
//    }

    fprintf(testF, "\n\n\n");

    fprintf(testF, "MACIERZ kierunku:\n\n");
    const char* direction = "s^<\\";
    fprintf(testF, "    ");
    for (int i = 1; i <= seq2Length; i++)
        fprintf(testF, "   %c", sm->revConvert(seq2[i]));
    fprintf(testF, "\n");
    for (int i = 1; i <= seq1Length; i++)
    {
        fprintf(testF, "   %c", sm->revConvert(seq1[i]));
        for (int j = 1; j <= seq2Length; j++)
        {
            fprintf(testF, "[%c%c%c%4d]", direction[B[i][j].backDirection], (B[i][j].continueUp)?('u'):(' '), (B[i][j].continueLeft)?('l'):(' '), A[i][j]);
        }
        fprintf(testF, "\n");
    }


    //RESULT
    fprintf(testF, "\n\nscore: %d\n", score);
    fprintf(testF, "seq1: %s\n", result1);
    fprintf(testF, "seq2: %s\n", result2);



    fclose(testF);
}

void SimilarityAlgorithmCpu::DeallocateMemoryForSingleRun()
{
    if(A != NULL)
    {
        if(A[0] != NULL)
            delete[] A[0];
        delete[] A;
    }
    if(E != NULL)
    {
        if(E[0] != NULL)
            delete[] E[0];
        delete[] E;
    }
    if(F != NULL)
    {
        if(F[0] != NULL)
            delete[] F[0];
        delete[] F;
    }
    if(B != NULL)
    {
        if(B[0] != NULL)
            delete[] B[0];
        delete[] B;
    }
    if(result1 != NULL)
        delete[] result1;
    if(result2 != NULL)
        delete[] result2;

    A = NULL;
    E = NULL;
    F = NULL;
    B = NULL;
    result1 = NULL;
    result2 = NULL;
}

void SimilarityAlgorithmCpu::DeallocateMemoryForAllRuns()
{
    if(scores != NULL) //array with pointers to the actual array with scores
    {
        if (scores[0] != NULL) //scores[0] - pointer to memory with scores
        {
            delete[] scores[0];
        }
        delete[] scores;
        scores = NULL;
    }

    if(matches1 != NULL)
        delete matches1;
    if(matches2 != NULL)
        delete matches2;
    matches1 = NULL;
    matches2 = NULL;

    

}

SimilarityAlgorithmCpu::~SimilarityAlgorithmCpu()
{
    DeallocateMemoryForSingleRun();
    DeallocateMemoryForAllRuns();
}
