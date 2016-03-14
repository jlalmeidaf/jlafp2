#ifndef _SEMAPHORE_H_
#define _SEMAPHORE_H_

    #include <semaphore.h>

    class Semaphore
    {
    public:
        Semaphore(unsigned int initValue);
        void wait();
        void signal();
        
    private:
        sem_t semaphore;
    };

#endif
