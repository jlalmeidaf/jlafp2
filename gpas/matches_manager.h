#ifndef _MATCHES_MANAGER_H_
#define _MATCHES_MANAGER_H_

namespace Data
{
    class MatchesManager
    {
    public:
        MatchesManager(int sequencesCount, int windowSize, int maxSequenceLength);
        ~MatchesManager();
//        char& operator() (int x, int y, int cnum);
//        char* operator() (int windowX, int windowY);
        char* getSequence(int x, int y);
        char* getWindow(int windowX, int windowY);
    private:
        char*** windows;
        int windowSize;
        int windowsXYCount;
        int maxSequenceLength;
    };
}

#endif
