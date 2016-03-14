#ifndef _SUBSTITUTION_MATRIX_H_
#define	_SUBSTITUTION_MATRIX_H_ 1

    #include <stdio.h>
    #include <ctype.h>
    #include <string.h>
    #include <string>
    #include <limits.h>
    #include <stdlib.h>

    #include "exceptions.h"


    namespace Data
    {
        class SubstitutionMatrix
        {
        public:
            SubstitutionMatrix(const char* filename);
            ~SubstitutionMatrix();

            char convert(unsigned char c);//A->0
            unsigned char revConvert(char i);//0->A
            char getScore(char i_1, char i_2);//getScore(0, 2);

            unsigned int getMatrixSize();
            char* getMatrix();

            unsigned char* getRevConv();

            int getLettersNumber();

            //methods for tests only (not converted strings e.g. ACCGTG...)
            //if seq1 is uchar type, must be finished with 0
            //if seq1 is char type, must be finished with -1
            int computeScore(unsigned char* seq1, unsigned char* seq2, int gapOp, int gapEx);
            int computeScore(char* seq1, char* seq2, int gapOp, int gapEx);

        private:
            char conversion[256];   //conversion of letters to numbers (e.g. A->0, R->1, N->2)
            unsigned char revConversion[128]; //revert conversion (e.g. 0->A, 1->R, 2->N)
            char* matrix;
            int lettersNumber;      //lettersNumber^2 == number of elements in 'matrix'
        };
    }


#endif	/* _SUBSTITUTION_MATRIX_H */
