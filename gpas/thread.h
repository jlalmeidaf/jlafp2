#ifndef _THREADS_H_
#define _THREADS_H_

    #include <pthread.h>

    class Thread
    {
    public:
        //Create thread.
        Thread(void*(*func)(void*), void* data);

        //Wait for thread to finish.
        void join();

        //Runs the thread
        void run();

        //Destroy thread.
        void kill();

        pthread_t getThreadId();

        //Wait for multiple threads.
        static void joinAll(Thread** threads, int num);
        
    private:
        void*(*func)(void*);
        void* data;
        pthread_t thread;
    };

#endif
