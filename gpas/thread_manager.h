#ifndef _THREAD_MANAGER_H_
#define _THREAD_MANAGER_H_

#include <list>
#include <map>
#include <cuda_runtime_api.h>
#include <unistd.h>
#include "thread.h"
#include "exceptions.h"
#include "sem.h"
#include "mutex.h"

using std::list;
using std::map;

class ThreadManager;

typedef void (*ThreadManagerRequest)(ThreadManager*, void*);

struct Request
{
    //void(*func)(ThreadManager*, void*);
    ThreadManagerRequest func;
    void* data;
    int gpuNo; //-1 -> runs on any card
               // 1 -> runs on 1 card
};

struct ThreadsInfo
{
    int gpuNo; //the ID of the GPU that is used by specyfic thread
    Mutex* anyNewRequest;
};

class ThreadManager
{
public:
    ThreadManager(int* gpus,unsigned int threadNumber);
    ~ThreadManager();

    //arguments: the function to be called, pointer to input data, ID of GPU to execute (-1 -> any GPU)
    void request(ThreadManagerRequest func, void* data, int gpuNo);//void(*func)(ThreadManager*, void*), void* data, int gpuNo);

    //waits for all requests to end
    void wait();

    //can be called only from singleThread method or request function
    ThreadsInfo getThreadsInfo();

    int getThreadsNumber();
    int* threadsIDs; //numbers associated with threads, here: which devices are used, read only array

private:
    //METHODS
    static void* singleThread(void* data);//thanks to "data" we have access to all data in this class
    Request* findRequest();

    //DATA
    unsigned int threadsNumber;
    Thread** threads;
    list<Request*> requests; //whenever reading/writing - use mutex requestsMutex!
    Semaphore* requestsMutex;
    Semaphore* requestCountSemaphore;
    Semaphore* threadsInfosSemaphore;//threads have to wait for threadsInfos map to be filled with the data
    Semaphore* requestsDoneSemaphore;   //used in wait()
    unsigned int requestsNumber;        //used in wait()
    map<pthread_t, ThreadsInfo> threadsInfos;//maps pthread_t to a GPU
};


#endif
