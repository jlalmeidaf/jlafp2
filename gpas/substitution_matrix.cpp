#include "substitution_matrix.h"


using namespace Data;
using Exceptions::FileNotFoundException;
using Exceptions::IndexOutOfRangeException;
using std::string;


SubstitutionMatrix::SubstitutionMatrix(const char* filename)
{
    for(int i=0; i<256; i++)
        conversion[i] = -1;
    for(int i=0; i<128; i++)
        revConversion[i] = 0;

    lettersNumber = 0;

    FILE* fd;
    fd = fopen(filename, "rb");
    if (fd == NULL)
        throw new FileNotFoundException(filename);

    bool comment = false;
    bool dataRead = false; //indicate whether some data has been read

    unsigned char c;
    while (!feof(fd))
    {
        c = fgetc(fd);
        if (c == '#')
        {
            comment = true;
            continue;
        }
        if(comment && c == '\n')
        {
            comment = false;
            continue;
        }
        if(comment || c == ' ')
        {
            continue;
        }
        if(dataRead && c == '\n')
        {
            break;
        }
        if(c == '\n')
        {
            continue;
        }
        
        if((c >= 'a') && (c <= 'z'))
        {
            c = c - 'a' + 'A'; // to convert small letters to big letters
        }
        conversion[c] = lettersNumber; //e.g. conversion['R'] == 1
        revConversion[lettersNumber] = c;
        lettersNumber++;
        dataRead = true;
        
    }//end of while
//    printf("%d\n", lettersNumber);
//    for(unsigned char i=0; i<255; i++)
//        printf("%c  %d\n",i, conversion[i]);
//    for(int i=0; i<128; i++)
//        printf("%d  %c\n",i, revConversion[i]);

    matrix = (char*)malloc(lettersNumber*lettersNumber*sizeof(char));

    int number;

    for(int i=0; i< lettersNumber; i++)
    {
        c = fgetc(fd);//a letter
        
        for(int j=0; j< lettersNumber; j++)
        {
            if(fscanf(fd,"%d",&number)>0)
            	matrix[i*lettersNumber + j] = (char)number;
            //printf("%3d", number);
        }
        while(!feof(fd) && fgetc(fd) != '\n');
        //printf("\n");
    }
    

    fclose(fd);
}

SubstitutionMatrix::~SubstitutionMatrix()
{
    if(matrix!=NULL)
        free(matrix);
}



char SubstitutionMatrix::convert(unsigned char c)
{
    return conversion[c];
}

unsigned char SubstitutionMatrix::revConvert(char i)
{
    if(i<0)
    {
        string s;
        s = "An error occured while doing revert conversion. ";
        s += "Function: revConvert, array: revConversion, index: ";
        s += "less than 0.";
        throw new IndexOutOfRangeException(s.c_str());
    }
    return revConversion[i];
}

char SubstitutionMatrix::getScore(char i_1, char i_2)
{
    if(i_1<0 || i_2<0)
        throw new IndexOutOfRangeException("Index out of range: SubstitutionMatrix::getScore\n");
    
    int addr = i_1*lettersNumber + i_2;
    if(addr >= lettersNumber*lettersNumber)
        throw new IndexOutOfRangeException("Index out of range: SubstitutionMatrix::getScore\n");

    return matrix[addr];
}

unsigned int SubstitutionMatrix::getMatrixSize()
{
    return lettersNumber*lettersNumber;
}

int SubstitutionMatrix::getLettersNumber()
{
    return lettersNumber;
}

char* SubstitutionMatrix::getMatrix()
{
    return matrix;
}



int SubstitutionMatrix::computeScore(unsigned char* seq1, unsigned char* seq2, int gapOp, int gapEx)
{
    int result = 0;
    bool gapSeq1Started = false;
    bool gapSeq2Started = false;
    while ((*seq1) && (*seq2))
    {
        if (*seq1=='-')
        {
            if (gapSeq1Started)
                result -= gapEx;
            else
            {
                gapSeq1Started = true;
                result -= gapOp;
            }
            gapSeq2Started = false;
        }
        else if (*seq2=='-')
        {
            if (gapSeq2Started)
                result -= gapEx;
            else
            {
                gapSeq2Started = true;
                result -= gapOp;
            }
            gapSeq1Started = false;
        }
        else
        {
            result += getScore(convert(*seq1), convert(*seq2));
            gapSeq1Started = false;
            gapSeq2Started = false;
        }
        seq1++;
        seq2++;
    }
    return result;
}

int SubstitutionMatrix::computeScore(char* seq1, char* seq2, int gapOp, int gapEx)
{
    int result = 0;
    bool gapSeq1Started = true;
    bool gapSeq2Started = true;
    while ((*seq1 != -1) && (*seq2 != -1))
    {
        if (*seq1==-2)
        {
            if (gapSeq1Started)
                result -= gapEx;
            else
            {
                gapSeq1Started = true;
                result -= gapOp;
            }
        }
        else if (*seq2==-2)
        {
            if (gapSeq2Started)
                result -= gapEx;
            else
            {
                gapSeq2Started = true;
                result -= gapOp;
            }
        }
        else
        {
            result += getScore(*seq1, *seq2);
        }
        seq1++;
        seq2++;
    }
    return result;
}

unsigned char* SubstitutionMatrix::getRevConv()
{
    return revConversion;
}
