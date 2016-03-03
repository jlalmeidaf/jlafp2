#ifndef _NEEDLEMAN_WUNSCH_SEMIGLOBAL_H_
#define _NEEDLEMAN_WUNSCH_SEMIGLOBAL_H_

    #include "back_up_struct.h"
    #include "sequences.h"
    #include "substitution_matrix.h"
    #include "similarity_algorithm_cpu.h"
    #include "needleman_wunsch_global.h"

    namespace Algorithms
    {
        namespace Sequential
        {
            class NeedlemanWunschSemiglobal: public NeedlemanWunschGlobal
            {
            protected:
                virtual void InitializeMatrices();
                virtual void BackwardMoving();
            };
        }
    }

#endif
