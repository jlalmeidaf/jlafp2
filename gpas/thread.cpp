#include "thread.h"

//Constructor
Thread::Thread(void*(*func)(void*), void* data)
{
    this->func = func;
    this->data = data;
}

void Thread::run()
{
    pthread_create(&thread, NULL, func, data);
}

//Wait for a thread to finish
void Thread::join()
{
    pthread_join(thread, NULL);
}

//Destroy the thread
void Thread::kill()
{
    pthread_cancel(thread);
}

pthread_t Thread::getThreadId()
{
    return thread;
}


//Wait for multiple threads
void Thread::joinAll(Thread** threads, int num)
{
    for(int i = 0; i < num; i++)
        threads[i]->join();
}

