#ifndef _SEQUENCES_H_
#define _SEQUENCES_H_ 1

    #include <stdio.h>
    #include <stdlib.h>
    #include <ctype.h>
    #include <string.h>
    #include <limits.h>

    #include "substitution_matrix.h"

    #define MIN(x, y) (((x)>(y))?(y):(x))
    #define MAX(x, y) (((x)<(y))?(y):(x))
    #define MAX3(x, y, z) (((x)>(y))?(((x)>(z))?(x):(z)):(((y)>(z))?(y):(z)))
    #define MIN3(x, y, z) (((x)<(y))?(((x)<(z))?(x):(z)):(((y)<(z))?(y):(z)))

    namespace Data
    {

        class Sequences
        {
        public:
            int getSequenceNumber();
            int getEntireLength();
            int* getStarts();
            int* getLengths();
            char* getSequences();
            int getMinSeqLen();
            int getMaxSeqLen();
            float getAvgSeqLen();
            void load();
            SubstitutionMatrix* getSubtitutionMatrix();
            int getSquaredSum(int windowSize, int windowNo, int blockShape);
            int getWindowSum(int windowSize, int windowNo, int blockShape);
            int getWindowMax(int windowSize, int windowNo);
            long long getCellsCount();
            char* getSeqName(int seqNo);


            Sequences(const char* filename, SubstitutionMatrix* sm);
            void writeToFile(const char* filename, int minLength, int maxLength);
            void writeToFile(const char* filename, int howMany);
            ~Sequences();
            void sortSequences();

        private:
            int* sequenceNumber;    //number of sequences
            int* entireLength;      //number of bytes that are used for storing all sequences
            int* starts;            //indicies of sequences beginnings
            int* lengths;
            char* seqs;             //all sequences
            char* filename;
            char** seqNames;        //names of sequences from file
            int* readStarts();

            int*   minSeqLen;
            int*   maxSeqLen;
            float* avgSeqLen;
            long long cellsCount;

            SubstitutionMatrix* sm;

        };

    }
    
#endif

