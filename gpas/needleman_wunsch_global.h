#ifndef _NEEDLEMAN_WUNSCH_H_
#define	_NEEDLEMAN_WUNSCH_H_

    #include "back_up_struct.h"
    #include "sequences.h"
    #include "substitution_matrix.h"
    #include "similarity_algorithm_cpu.h"

    using namespace Data;

    namespace Algorithms
    {
        namespace Sequential
        {
            class NeedlemanWunschGlobal: public SimilarityAlgorithmCpu
            {
            protected:
                virtual void InitializeMatrices();
                virtual void FillMatrices();
                virtual void BackwardMoving();
            };
        }
    }

#endif	/* _NEEDLEMAN_WUNSCH_H */

