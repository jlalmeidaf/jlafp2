#include "thread_manager.h"
#include "thread.h"

using Exceptions::IndexOutOfRangeException;

ThreadManager::ThreadManager(int* gpus,unsigned int threadNumber)
{
    requestsMutex = new Semaphore(1);
    requestCountSemaphore = new Semaphore(0);
    threadsInfosSemaphore = new Semaphore(0);
    requestsDoneSemaphore = new Semaphore(0);
    requestsNumber = 0;
    threadsNumber = threadNumber;
    threadsIDs = new int[threadNumber];

    threads = new Thread*[threadNumber];
    for(int i=0; i<threadNumber; i++)
    {
        threadsIDs[i] = gpus[i];

        threads[i] = new Thread(singleThread, this);
        threads[i]->run();

        ThreadsInfo ti;
        ti.gpuNo = gpus[i];
        ti.anyNewRequest = new Mutex();
        threadsInfos.insert(std::pair<pthread_t, ThreadsInfo>(threads[i]->getThreadId(), ti));
    }

    for(int i=0; i<threadNumber; i++)
    {
        threadsInfosSemaphore->signal();
    }
}

int ThreadManager::getThreadsNumber()
{
    return threadsNumber;
}

ThreadManager::~ThreadManager()
{
    for(int i=0; i<threadsNumber; i++)
        threads[i]->kill();
}

void ThreadManager::wait()
{
    while(requestsNumber > 0)
    {
        requestsDoneSemaphore->wait();
        requestsMutex->wait();
        requestsNumber--;
        requestsMutex->signal();
    }
}


void ThreadManager::request(ThreadManagerRequest func, void* data, int gpuNo = -1)//void(*func)(ThreadManager*, void*), void* data, int gpuNo = -1)
{
    requestsMutex->wait();

    Request* r = new Request();
    r->func = func;
    r->data = data;
    r->gpuNo = gpuNo;
    requests.push_back(r);


    //wake up potentially waiting threads
    map<pthread_t, ThreadsInfo>::iterator i;
    for(i=threadsInfos.begin(); i!=threadsInfos.end(); ++i)
    {
        if((gpuNo == -1) || (gpuNo == i->second.gpuNo))
        {
            i->second.anyNewRequest->unlock();
        }
    }

    requestsNumber++;
    requestCountSemaphore->signal();//increase the number of requests
    requestsMutex->signal();
}

void* ThreadManager::singleThread(void *data)
{
    ThreadManager* This = (ThreadManager*)data;
    This->threadsInfosSemaphore->wait(); //we have to wait for threadsInfos map to be filled with the data
    
    while(1)
    {
        try
        {
            Request* request = This->findRequest();
            request->func(This, request->data);

            This->requestsDoneSemaphore->signal();

        }
        catch(IndexOutOfRangeException* ex)
        {
            This->threadsInfos[pthread_self()].anyNewRequest->lock();
        }
    }

}

Request* ThreadManager::findRequest()
{
    requestCountSemaphore->wait();
    requestsMutex->wait();


    int myGpuNo = threadsInfos[pthread_self()].gpuNo;
    
    list<Request*>::iterator i;
    for(i=requests.begin(); i!=requests.end(); ++i)
    {
        if(((*i)->gpuNo == -1) || ((*i)->gpuNo == myGpuNo))
        {
            Request* result = *i; //copying pointer to result
            requests.remove(*i);
            
            requestsMutex->signal();
            return result;
        }
    }

    //WE HAVEN'T GOT ANY TASK
    requestCountSemaphore->signal();
    requestsMutex->signal();


    throw new IndexOutOfRangeException("No new tasks.");

}

ThreadsInfo ThreadManager::getThreadsInfo()
{
    return threadsInfos[pthread_self()];
}
