#ifndef _DATA_MANAGEMENT_H_
#define _DATA_MANAGEMENT_H_

    #include "cuda_declarations.h"
    #include "substitution_matrix.h"
    #include "sequences.h"

    using Data::SubstitutionMatrix;
    using Data::Sequences;

/*******************************************************************************
 * SUBSTITUTION_MATRIX_SIZE defines maximum size of the substitution matrix    *
 * which resides in constant memory of the card and contains the rate values   *
 * at which one character in a sequence changes to other character in second   *
 * sequence. 24*24=576                                                         *
 *******************************************************************************/

    #define SUBSTITUTION_MATRIX_SIZE 576
  
/*******************************************************************************
 * MAX_ALGORITHM_WINDOW_SIZE defines maximum value of the WINDOW_SIZE and      *
 * determines the count of start position of the sequences in the continous    *
 * block of the memory in stored in a texture (See CUDA documentation).        *
 *******************************************************************************/

    #define MAX_ALGORITHM_WINDOW_SIZE 512

/*******************************************************************************
 * substitutionMatrix resides in constant memory of the graphic card and       *
 * contains the rate values at which one character in a sequence changes to    *
 * other character (e.g. BLOSUM30).                                            *
 *******************************************************************************/

    __constant__ char substitutionMatrix[SUBSTITUTION_MATRIX_SIZE];

/*******************************************************************************
 * lettersCount contains the count of the different characters in substitution *
 * matrix.                                                                     *
 *******************************************************************************/

    __constant__ int  lettersCount;
    
/*******************************************************************************
 * WINDOW_SIZE defines the size of the window within the alignment algorithms. *
 * The value indicates the width as well as the height of the space in the     *
 * sequences pair matrix that is computed in one graphic card in one           *
 * iteration.                                                                  *
 *                                                                             *
 * NOTICE:                                                                     *
 * WINDOW_SIZE must be multiple of BLOCK_SHAPE e.g. 128, 144, 160              *
 * to satisfy kernel's needs.                                                  *
 *                                                                             *
 * NOTICE2:                                                                    *
 * WINDOW_SIZE should be smaller than MAX_ALGORITHM_WINDOW_SIZE                *
 * (see data_management.cuh)                                                   *
 *                                                                             *
 * ALIGNMENT_SCORE_WINDOW_SIZE referes to the size of the window in algorithms *
 * that as a result has a scalar value of the alignment score.                 *
 *                                                                             *
 * ALIGNMENT (SCORE ONLY) best WINDOW_SIZE:                                    *
 *  - GTS250:   128                                                            *
 *  - 9600GT:    96??                                                          *
 *  - GTX280:   240                                                            *
 *  - 8600MGS:  128                                                            *
 *                                                                             *
 * ALIGNMENT WITH BACKTRACKING best WINDOW_SIZE (3_667.fax):                   *
 *  - GTS250:      64                                                          *
 *  - 9600GT:      64                                                          *
 *  - GTX280:      96                                                          *
 *  - 8600MGS:     32                                                          *
 *  - TeslaC1060: 112                                                          *
 *******************************************************************************/

    __constant__ unsigned int  winSize;

/*******************************************************************************
 * MEMORY_OFFSET defines offset between the cells of global memory resides on  *
 * the graphic card. It helps in setting and retrieving consecutive values by  *
 * consecutive threads and speed it up significantly (see also CUDA            *
 * documentation).                                                             *
 * MEMORY_OFFSET == WINDOW_SIZE * WINDOW_SIZE                                  *
 *******************************************************************************/

    __constant__ unsigned int  memOffset;

/*******************************************************************************
 * MAX_LETTERS_COUNT defines maximum count of different characters in          *
 * substitutionMatrix.                                                         *
 *******************************************************************************/

    #define MAX_LETTERS_COUNT   48

/*******************************************************************************
 *
 *******************************************************************************/

    __constant__ unsigned char revConvConst[MAX_LETTERS_COUNT];

/*******************************************************************************
 * gapOp and gapEx contain the value of penalty for respectively opening and   *
 * extension of the gap in alignment algorithm. Resides in constant memory.    *
 *******************************************************************************/

    __constant__ char gapOp;
    __constant__ char gapEx;

/*******************************************************************************
 * tex1Starts and tex2Starts contain start positions' sequences in the raw     *
 * memory of the texture.                                                      *
 *******************************************************************************/

    __constant__ int  tex1Starts[MAX_ALGORITHM_WINDOW_SIZE];
    __constant__ int  tex2Starts[MAX_ALGORITHM_WINDOW_SIZE];

/*******************************************************************************
 * texSeqs1 and texSeqs2 contain sequences in a raw memory of the texture.     *
 *                                                                             *
 * t  -------------                                                            *
 * e |             |                                                           *
 * x |  alignment  |                                                           *
 * S |             |                                                           *
 * e |             |                                                           *
 * q |   matrix    |                                                           *
 * s |             |                                                           *
 * Y  ------------                                                             *
 *  t e x S e q s X                                                            *
 *                                                                             *
 * texBack - texture with data needed to backtrace in NW & SW algorithms       *
 *******************************************************************************/

    texture<char, 1, cudaReadModeElementType> texSeqsX;
    texture<char, 1, cudaReadModeElementType> texSeqsY;
    texture<unsigned int, 1, cudaReadModeElementType> texBack;

/*******************************************************************************
 * TexVariablesAddresses is a result type of the function copySeqsToTex and    *
 * contains pointers (to free) of a global memory of the graphic card.         *
 *******************************************************************************/

    struct TexVariablesAddresses
    {
        char* texSeqs1DevPtr;
        char* texSeqs2DevPtr;
    };


/*******************************************************************************
 * smScoreParams contains a params of the copySMToConstInThread function,      *
 * which is invoked dynamically by the ThreadInvoker                           *
 *******************************************************************************/

    struct SubstitutionMatrixParams
    {
        SubstitutionMatrix* sm;
        char gapOpen, gapExtension;
        unsigned int windowSize, memoryOffset;
    } smParams;

#endif
