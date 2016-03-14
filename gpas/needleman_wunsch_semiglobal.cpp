#include "needleman_wunsch_semiglobal.h"

using namespace Algorithms::Sequential;
using Data::Sequences;
using Data::SubstitutionMatrix;

void NeedlemanWunschSemiglobal::InitializeMatrices()
{
    for (int i = 0; i < seq1Length + 1; i++)
    {
        A[i][0] = 0;
        E[i][0] = INT_MIN + gapEx; //jednak musi być bardzo niska wartość zamiast i*gapEx;
        F[i][0] = INT_MIN + gapEx; //(INT_MIN - gapEx) MINUS JEST CELOWY!!!
        B[i][0].backDirection = up;
    }
    for (int j = 0; j < seq2Length + 1; j++)
    {
        A[0][j] = 0;
        E[0][j] = INT_MIN + gapEx;
        F[0][j] = INT_MIN + gapEx;
        B[0][j].backDirection = left;
    }
    B[0][0].backDirection = stop;
}

void NeedlemanWunschSemiglobal::BackwardMoving()
{
    int y = seq1Length;
    int x = seq2Length;
    score = INT_MIN;

    //SEARCHING FOR MAXIMAL SCORE IN THE LAST ROW AND THE LAST COLUMN
    for(int i=1; i<=seq2Length; i++)
    {
        if(A[seq1Length][i] > score)
        {
            score = A[seq1Length][i];
            y = seq1Length;
            x = i;
        }
    }
    for(int i=1; i<=seq1Length; i++)
    {
        if(A[i][seq2Length] > score)
        {
            score = A[i][seq2Length];
            y = i;
            x = seq2Length;
        }
    }


    
    //BACKWARD MOVING
    int carret = 0;

    //GO TO ELEMENT WITH MAXIMAL SCORE
    for(int i=seq2Length; i>x; i--)
    {
        result1[carret] = '-';
        result2[carret] = sm->revConvert(seq2[i]);
        carret++;
    }
    for(int i=seq1Length; i>y; i--)
    {
        result1[carret] = sm->revConvert(seq1[i]);
        result2[carret] = '-';
        carret++;
    }


    BackDirection prev = crosswise;
    while(B[y][x].backDirection != stop)
    {
        if (prev == up && B[y][x].continueUp) //CONTINUE GOING UP
        {                                          //GAP EXTENSION
            result1[carret] = sm->revConvert(seq1[y]);
            result2[carret] = '-';
            carret++;
            y--;
        }
        else if (prev == left && B[y][x].continueLeft) //CONTINUE GOING LEFT
        {                                         //GAP EXTENSION
            result1[carret] = '-';
            result2[carret] = sm->revConvert(seq2[x]);
            carret++;
            x--;
        }
        else
        {
            prev = B[y][x].backDirection;
            if(prev == up)
            {
                result1[carret] = sm->revConvert(seq1[y]);
                result2[carret] = '-';
                carret++;
                y--;
            }
            else if(prev == left)
            {
                result1[carret] = '-';
                result2[carret] = sm->revConvert(seq2[x]);
                carret++;
                x--;
            }
            else //prev == crosswise
            {
                result1[carret] = sm->revConvert(seq1[y]);
                result2[carret] = sm->revConvert(seq2[x]);
                carret++;
                x--;
                y--;
            }
        }
    }
    result1[carret] = 0;
    result2[carret] = 0;
}
