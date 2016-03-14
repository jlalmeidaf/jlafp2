#ifndef _SIMILARITY_ALGORITHM_CPU_H_
#define _SIMILARITY_ALGORITHM_CPU_H_

#include "back_up_struct.h"
#include "sequences.h"
#include "substitution_matrix.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "matches_manager.h"
#include "hi_res_timer.h"

using Data::Sequences;
using Data::SubstitutionMatrix;
using Data::MatchesManager;

namespace Algorithms
{
    namespace Sequential
    {
        class SimilarityAlgorithmCpu
        {
        public:
            virtual void Run(Sequences* s, int seq1No, int seq2No, int gapOp, int gapEx);
            virtual void RunAll(Sequences* s, int gapOp, int gapEx);
            virtual void PrintResults(const char* fileName);
            SimilarityAlgorithmCpu();
            ~SimilarityAlgorithmCpu();

            MatchesManager* matches1;
            MatchesManager* matches2;
            int** scores;

        protected:
            virtual void AllocateMemoryForSingleRun();
            virtual void InitializeMatrices() = 0;
            virtual void FillMatrices() = 0;
            virtual void BackwardMoving() = 0;
            virtual void DeallocateMemoryForAllRuns();
            virtual void DeallocateMemoryForSingleRun();

            //INPUT DATA
            // seq1 and seq2 are 1-based array of sequence characters
            char* seq1;
            char* seq2;
            int seq1Length;
            int seq2Length;
            SubstitutionMatrix* sm;
            int gapOp;
            int gapEx;

            //MATRICES
            int **A;
            int **E; //left matrix
            int **F; //up matrix
            BackUpStruct **B;

            //RESULTS
            char* result1;
            char* result2;
            int score;
        };
    }
}

#endif
