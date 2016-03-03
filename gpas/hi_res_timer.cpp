#include "hi_res_timer.h"
#include <stdio.h>

void HiResTimer::start()
{
    gettimeofday(&t1, NULL);
}

void HiResTimer::stop()
{
    gettimeofday(&t2, NULL);
}

unsigned long HiResTimer::getElapsedTime()
{
    double elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

    return (unsigned long)elapsedTime;
}
