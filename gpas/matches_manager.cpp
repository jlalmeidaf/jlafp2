#include "matches_manager.h"
#include "exceptions.h"
#include <stdio.h>

using namespace Data;
using namespace Exceptions;

MatchesManager::MatchesManager(int sequencesCount, int windowSize, int maxSequenceLength)
{
    windowsXYCount = (sequencesCount - 1)/windowSize + 1;
    windows = new char**[windowsXYCount];
    this->windowSize = windowSize;
    windows[0] = new char*[windowsXYCount*(windowsXYCount + 1)/2];
    for (int j = 1; j < windowsXYCount; j++)
        windows[j] = windows[j - 1] + j;

    unsigned long long memToAllocate = 0;
    for (int j = 0; j < windowsXYCount; j++)
        for (int i = 0; i <= j; i++)
            memToAllocate += windowSize*windowSize*maxSequenceLength;

    try
    {
        for (int j = 0; j < windowsXYCount; j++)
            for (int i = 0; i <= j; i++)
                windows[j][i] = new char[windowSize*windowSize*maxSequenceLength];
    }
    catch(...)
    {
        printf("Allocation of %.0fMB memory failed.\n", ((double)memToAllocate)/(1024.0*1024.0));
        throw new IndexOutOfRangeException("Not enough memory! Terminating.");
    }
    

    this->maxSequenceLength = maxSequenceLength;
}

//char& MatchesManager::operator ()(int x, int y, int cnum)
//{
//    if (x > y)
//        return (*this)(y, x, cnum);
//    char* currentWindow = windows[y/windowSize][x/windowSize];
//    return currentWindow[((y%windowSize)*windowSize + x%windowSize)*maxSequenceLength + cnum];
//}
//
//char* MatchesManager::operator ()(int windowX, int windowY)
//{
//
//}

char* MatchesManager::getWindow(int windowX, int windowY)
{
    return windows[windowY][windowX];
}

char* MatchesManager::getSequence(int x, int y)
{
    if (x > y)
        return getSequence(y, x);
    char* currentWindow = windows[y/windowSize][x/windowSize];
    return &currentWindow[((y%windowSize)*windowSize + x%windowSize)*maxSequenceLength];
}

MatchesManager::~MatchesManager()
{
    for (int j = 0; j < windowsXYCount; j++)
        for (int i = 0; i <= j; i++)
            delete[] windows[j][i];
    delete[] windows[0];
    delete[] windows;
}
