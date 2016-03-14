#include "needleman_wunsch_global.h"

using namespace Algorithms::Sequential;
using Data::Sequences;
using Data::SubstitutionMatrix;

void NeedlemanWunschGlobal::InitializeMatrices()
{
    for (int i = 0; i < seq1Length + 1; i++)
    {
        A[i][0] = -i*gapEx;
        E[i][0] = INT_MIN + gapEx; //jednak musi być bardzo niska wartość zamiast i*gapEx;
        F[i][0] = INT_MIN + gapEx; //(INT_MIN - gapEx) MINUS JEST CELOWY!!!
        B[i][0].backDirection = up;
    }
    for (int j = 0; j < seq2Length + 1; j++)
    {
        A[0][j] = -j*gapEx;
        E[0][j] = INT_MIN + gapEx;
        F[0][j] = INT_MIN + gapEx;
        B[0][j].backDirection = left;
    }
    B[0][0].backDirection = stop;
}


void NeedlemanWunschGlobal::FillMatrices()
{
    /*
     *   s e q 2
     * s
     * e
     * q
     * 1
     */

    //ALGORYTM
    //E - responsible for left direction
    //F - responsible for up   direction
    for (int i = 1; i <= seq1Length; i++)
    {
        for (int j = 1; j <= seq2Length; j++)
        {
            E[i][j] = MAX(E[i][j - 1] - gapEx, A[i][j - 1] - gapOp);
            B[i][j - 1].continueLeft = (E[i][j] == E[i][j - 1] - gapEx);
            F[i][j] = MAX(F[i - 1][j] - gapEx, A[i - 1][j] - gapOp);
            B[i - 1][j].continueUp = (F[i][j] == F[i - 1][j] - gapEx);

            A[i][j] = MAX3(E[i][j], F[i][j], A[i - 1][j - 1] + sm->getScore(seq1[i], seq2[j]));


            if(A[i][j] == (A[i - 1][j - 1] + sm->getScore(seq1[i], seq2[j])))
                B[i][j].backDirection = crosswise;
            else if(A[i][j] == E[i][j])
                B[i][j].backDirection = left;
            else //if(A[i][j] == F[i][j])
                B[i][j].backDirection = up;


//            if(A[i][j] == E[i][j])
//                B[i][j].backDirection = left;
//            else if(A[i][j] == F[i][j])
//                B[i][j].backDirection = up;
//            else
//                B[i][j].backDirection = crosswise;

        }
    }    
}

void NeedlemanWunschGlobal::BackwardMoving()
{
    score = A[seq1Length][seq2Length];

    //BACKWARD MOVING
    int y = seq1Length;
    int x = seq2Length;
    int carret = 0;
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
