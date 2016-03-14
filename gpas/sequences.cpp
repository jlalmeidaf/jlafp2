#include <math.h>
#include <vector_types.h>
#include <stdlib.h>
#include <time.h>

#include "sequences.h"
#include "exceptions.h"


#define BUFFER_SIZE 255

using namespace Data;
using namespace Exceptions;
using std::string;


Sequences::Sequences(const char* filename, SubstitutionMatrix* sm)
{
    sequenceNumber = NULL;
    entireLength = NULL;
    starts = NULL; //indicies of sequences beginnings
    seqs = NULL; //all sequences
    cellsCount = -1;
    this->filename = strdup(filename);

    printf("%d seqs\n", getSequenceNumber());
    getEntireLength();
    getStarts();
    //printf("%d seqlen\n", getEntireLength());
    //printf("%d secondstart\n", getStarts()[1]);
    printf("%d min\n", getMinSeqLen());
    printf("%d max\n", getMaxSeqLen());
    printf("%.1f avg\n", getAvgSeqLen());
    printf("%lld cells\n", getCellsCount());
    //printf("seqences:\n%s\n", getSequences());

    this->sm = sm;

    srand((unsigned)time(0));
}

Sequences::~Sequences()
{
    if (filename != NULL)
        delete[] filename;
    if (sequenceNumber != NULL)
        delete sequenceNumber;
    if (entireLength != NULL)
        delete entireLength;
    if (starts != NULL)
        delete[] starts;
    if (seqs != NULL)
        delete[] seqs;
}

int Sequences::getSequenceNumber()
{
    if (sequenceNumber == NULL)
    {
        FILE* fd;
        fd = fopen(filename, "rb");
        if (fd == NULL)
            throw new FileNotFoundException(filename);
        int result = 0;
        int actualLength;
        char buffer[BUFFER_SIZE];

        while (!feof(fd))
        {
            actualLength = (int) fread(buffer, sizeof (char), BUFFER_SIZE, fd);
            for (int i = 0; i < actualLength; i++)
            {
                if (buffer[i] == '>')
                {
                    result++;
                }
            }
        }
        fclose(fd);

        sequenceNumber = new int(result);
    }
    return *sequenceNumber;
}

int Sequences::getEntireLength()
{
    if (entireLength == NULL)
    {
        FILE* fd;
        fd = fopen(filename, "rb");
        if (fd == NULL)
            throw new FileNotFoundException(filename);

        char buffer[BUFFER_SIZE];
        int actualLength;
        int sequenceNumber = getSequenceNumber();
        int *sequencesLengths = new int[sequenceNumber];
        lengths = new int[sequenceNumber];

        for (int i = 0; i < sequenceNumber; i++)
        {
            sequencesLengths[i] = 0;
        }

        int actualSequence = -1;

        bool title = false;

        while (!feof(fd))
        {
            actualLength = fread(buffer, sizeof(char), BUFFER_SIZE, fd);
            for (int i = 0; i < actualLength; i++)
            {
                if (buffer[i] == '>')
                {
                    title = true;
                    actualSequence++;
                    continue;
                }
                if (buffer[i] == '\n')
                {
                    title = false;
                    continue;
                }
                if ((buffer[i] == '\t') || (buffer[i] == ' '))
                {
                    continue;
                }
                if (!title)
                {
                    sequencesLengths[actualSequence]++;
                }
            }
        }
        if (sequenceNumber > 0)
        {
            int min = sequencesLengths[0];
            int max = sequencesLengths[0];
            float avg;
            
            int result = sequencesLengths[0];
            avg = sequencesLengths[0];
            lengths[0] = sequencesLengths[0];
            for (int i = 1; i < sequenceNumber; i++)
            {
                lengths[i] = sequencesLengths[i];
                result += sequencesLengths[i];
                min = MIN(min, sequencesLengths[i]);
                max = MAX(max, sequencesLengths[i]);
                sequencesLengths[i] += sequencesLengths[i - 1];
            }
            avg = ((float)result) / sequenceNumber;
            minSeqLen = new int(min);
            maxSeqLen = new int(max);
            avgSeqLen = new float(avg);

            entireLength = new int(result);
            starts = sequencesLengths;
            for (int i = getSequenceNumber() - 1; i >= 1; i--)
            {
                starts[i] = starts[i - 1];
            }
            starts[0] = 0;
        }
    }
    return *entireLength;
}

int Sequences::getMinSeqLen()
{
    if (minSeqLen == NULL)
    {
        getEntireLength();
    }
    return *minSeqLen;
}

int Sequences::getMaxSeqLen()
{
    if (maxSeqLen == NULL)
    {
        getEntireLength();
    }
    return *maxSeqLen;
}

float Sequences::getAvgSeqLen()
{
    if (avgSeqLen == NULL)
    {
        getEntireLength();
    }
    return *avgSeqLen;
}

int* Sequences::getStarts()
{
    if (starts == NULL)
    {
        getEntireLength();
    }
    return starts;
}

int* Sequences::getLengths()
{
    if (lengths == NULL)
    {
        getEntireLength();
    }
    return lengths;
}

char* Sequences::getSequences()
{
    if (seqs == NULL)
    {
        int length = getEntireLength();
        seqs = new char[length + 1];
        seqs[length] = 0;

        FILE* fd;
        fd = fopen(filename, "rb");
        if (fd == NULL)
            throw new FileNotFoundException(filename);

        char buffer[BUFFER_SIZE];
        int actualLength;
        int actualSequence = -1;
        int offset = 0; //przesunięcie wewnątrz konkretnej sekwencji
        int* starts = getStarts();
        char actualChar;

        bool title = false;
        char* titleBuf = new char[1024];
        int titleCarret = 0;
        seqNames = new char*[getSequenceNumber()];

        while (!feof(fd))
        {
            actualLength = fread(buffer, sizeof(char), BUFFER_SIZE, fd);
            for (int i = 0; i < actualLength; i++)
            {
                if(title && (buffer[i] != '\n'))
                {
                    if( (buffer[i]==' ') || (buffer[i]=='\t') )
                        continue;

                    titleBuf[titleCarret] = buffer[i];
                    titleCarret++;
                    if(titleCarret==1024)
                        throw IndexOutOfRangeException("The name of sequence in input file is too long.");
                    continue;
                }
                if (buffer[i] == '>')
                {
                    title = true;
                    actualSequence++;
                    titleCarret = 0;
                    offset = 0;
                    continue;
                }
                if (buffer[i] == '\n')
                {
                    title = false;
                    titleBuf[titleCarret] = 0;
                    seqNames[actualSequence] = new char[titleCarret+1];
                    strcpy(seqNames[actualSequence], titleBuf);
                    continue;
                }
                if ((buffer[i] == '\t') || (buffer[i] == ' '))
                {
                    continue;
                }
                if (!title)
                {
                    if ((buffer[i] >= 'a') && (buffer[i] <= 'z'))
                    {
                        buffer[i] = buffer[i] - 'a' + 'A';
                    }

                    actualChar = sm->convert((unsigned char) buffer[i]);
                    //printf("%d\n", actualChar);

                    if (actualChar < 0)
                    {
                        string s;
                        s = "In sequences character '";
                        s += (unsigned char)buffer[i];
                        s += "' which is not present in substitution matrix.";
                        throw new ConversionException(s.c_str());
                    }
                    seqs[starts[actualSequence] + offset] = actualChar;
                    offset++;
                }

            }
        }

        //WE REVERT ALL SEQUENCES TO NOT TO REVERT RESULTS
        //REVERTING RESULTS IS MORE COMPLEX
        int numberOfSequences = getSequenceNumber();
        int* lengths = getLengths();
        for (int seqNo = 0; seqNo < numberOfSequences; seqNo++)
        {
            char* actualSeq = &seqs[starts[seqNo]];
            int seqLength = lengths[seqNo];
            char tmp;
            for (int i = 0; i < seqLength / 2; i++)
            {
                tmp = actualSeq[i];
                actualSeq[i] = actualSeq[seqLength - 1 - i];
                actualSeq[seqLength - 1 - i] = tmp;
            }
        }

    }
    return seqs;
}

void Sequences::load()
{
    getSequenceNumber();
    getEntireLength();
    getAvgSeqLen();
    getMaxSeqLen();
    getMinSeqLen();
    getStarts();
    getSequences();
}

SubstitutionMatrix* Sequences::getSubtitutionMatrix()
{
    return sm;
}

class SortSequence
{
public:
    int index;
    Sequences* seqs;
};

int compareSortSequences(const void* first, const void* second)
{
    SortSequence* sortSequence1 = *(SortSequence**)first;
    SortSequence* sortSequence2 = *(SortSequence**)second;
    int* lengths = sortSequence1->seqs->getLengths();
    return lengths[sortSequence2->index] - lengths[sortSequence1->index];
}

void Sequences::sortSequences()
{
    int seqNum = getSequenceNumber();
    SortSequence** sortSequences = new SortSequence*[seqNum];
    for (int i = 0; i < seqNum; i++)
    {
        sortSequences[i] = new SortSequence();
        sortSequences[i]->index = i;
        sortSequences[i]->seqs = this;
    }
    qsort(sortSequences, seqNum, sizeof(SortSequence*), compareSortSequences);
////DEBUG
//    for (int i = 0; i < seqNum; i++)
//    {
//        printf("%d\n", sortSequences[i]->index);
//    }

    int entireLength = getEntireLength();
    char* sortedSeqs = new char[entireLength];
    char* sortingCarret = sortedSeqs;
    char* globalCarret;
    int sortingLength;

    for (int i = 0; i < seqNum; i++)
    {
        int index = sortSequences[i]->index;
        sortingLength = lengths[index];
        globalCarret = seqs + starts[index];
        for (int j = 0; j < sortingLength; j++)
        {
            *sortingCarret = *globalCarret;
            sortingCarret++;
            globalCarret++;
        }
    }

    delete[] seqs;
    seqs = sortedSeqs;

    int* sortedLengths = new int[seqNum];
    char** sortedSeqNames = new char*[seqNum];
    for (int i = 0; i < seqNum; i++)
    {
        int index = sortSequences[i]->index;
        sortedLengths[i] = lengths[index];
        sortedSeqNames[i] = seqNames[index];
    }

    delete[] lengths;
    lengths = sortedLengths;

    delete[] seqNames;
    seqNames = sortedSeqNames;

    for (int i = 1; i < seqNum; i++)
    {
        starts[i] = starts[i - 1] + lengths[i - 1];
    }

////DEBUG
//    for (int i = 0; i < lengths[0]; i++)
//    {
//        printf("%c", sm->revConvert(seqs[i])); //reverted!!!!
//    }
//    printf("\n");
}

void Sequences::writeToFile(const char* filename, int minLength, int maxLength)
{
    int n = getSequenceNumber();

    int* starts = getStarts();
    int* lengths = getLengths();
    char* seqs = getSequences();

    FILE* file = fopen(filename, "w");
    for(int j = 0; j<=n; j++)
    {
        if((lengths[j] > maxLength) || (lengths[j] < minLength))
            continue;

        fprintf(file, ">sequence%d", j);
        for (int i = 0; i < lengths[j]; i++)
        {
            if(i%60==0)
                fprintf(file, "\n");
            fprintf(file, "%c", sm->revConvert(seqs[starts[j] + i]));
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void Sequences::writeToFile(const char* filename, int howMany)
{
    int n = getSequenceNumber();
    int seqIdx;
    
    int* starts = getStarts();
    int* lengths = getLengths();
    char* seqs = getSequences();

    FILE* file = fopen(filename, "w");
    for(int j = 0; j<howMany; j++)
    {
        seqIdx = rand()%n;

        fprintf(file, ">sequence%d", j);
        for (int i = 0; i < lengths[seqIdx]; i++)
        {
            if(i%60==0)
                fprintf(file, "\n");
            fprintf(file, "%c", sm->revConvert(seqs[starts[seqIdx] + i]));
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


int Sequences::getWindowSum(int windowSize, int windowNo, int blockShape)
{
    int sum = 0;
    int offset = windowSize * windowNo;
    for (int i = 0;  (i < windowSize) && (offset + i < getSequenceNumber() ); i+=blockShape)//for (int i = blockShape - 1; (i < windowSize) && (offset + i - blockShape < getSequenceNumber() ); i+=blockShape)
    {
        int startPos = offset + i;//MIN(offset + i, getSequenceNumber() - 1);
        sum += lengths[startPos];
    }
    return sum;
}


int Sequences::getSquaredSum(int windowSize, int windowNo, int blockShape)
{
    int sum = 0;
    int offset = windowSize * windowNo;
    for (int i = 0;  (i < windowSize) && (offset + i < getSequenceNumber() ); i+=blockShape)//for (int i = blockShape - 1; (i < windowSize) && (offset + i - blockShape < getSequenceNumber() ); i+=blockShape)
    {
        int startPos = offset + i;//, getSequenceNumber() - 1);
        sum += lengths[startPos] * lengths[startPos];
    }
    return sum;
}

int Sequences::getWindowMax(int windowSize, int windowNo)
{
    //int sum = 0;
    int offset = windowSize * windowNo;
    //int endPos = MIN(offset + windowSize - 1, getSequenceNumber() - 1);
    return lengths[offset];//lengths[endPos];
}

long long Sequences::getCellsCount()
{
    if(cellsCount == -1)
    {
        cellsCount = 0;
        int n = getSequenceNumber();
        int* len = getLengths();
        
        for(int i=0; i<n; i++)
            for(int j=0; j<=i; j++)
            {
                cellsCount += len[i] * len[j];
            }
	}
    return cellsCount;
}

char* Sequences::getSeqName(int seqNo)
{
    if( (seqNo >= getSequenceNumber()) || (seqNo<0) )
        throw new IndexOutOfRangeException("No such sequence.");

    return seqNames[seqNo];
}
