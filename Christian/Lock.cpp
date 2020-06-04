#include <assert.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <omp.h>
#include <atomic>
#include <chrono>
#include <stdlib.h>
#include <thread>
#include <boost/thread/tss.hpp>

///////////////////////////////////////////////// TAS Lock

class TAS_lock {
    private:
    std::atomic<bool> state;
    
    public:

    TAS_lock(){
        state = false;
    }
    
    void lock(){
        while (state.exchange(true))
        {}
    }

    void unlock(){
        state.exchange(false);
    }

};
        
//////////////////////////////////////////////// TTAS Lock

class TTAS_lock {
    private:
    std::atomic<bool> state;
    
    public:
    TTAS_lock(){
        state = false;
    }
    void lock(){
        while (true) {
            while (state.load()) 
            {}
            if (!state.exchange(true))
                return;
        }
    }
    void unlock(){
        state.exchange(false);
    }
};

/////////////////////////////////////////////// Ticket Lock

class Ticket_lock {
    std::atomic<int> ticket;
    int served;

    public:
    
    Ticket_lock() {
        ticket = 0;
        served = 0;
    }

    void lock() {
        int next = ticket.fetch_add(1);
        volatile int* my_served = &served;
        while (*my_served < next)
        {}
    }

    void unlock() {
        served++;
    }

};

///////////////////////////////////////////////// Array Lock
/*
class Array_lock {
    bool* flag;
    std::atomic<int> tail;
    std::atomic<int>* mySlot;
    int numthreads;

    public:

    Array_lock(int n) : flag(new bool[n]), mySlot(new std::atomic<int>[n]) {    ///////////////////////////// false sharing noch korrigieren
        for (int i = 0; i < n; ++i)
            flag[i] = false;
        flag[0] = true;
        tail = 0;
        numthreads = n;
    }

    void lock(int tid) {
        mySlot[tid] = tail.fetch_add(1)%numthreads;
        while (!flag[mySlot[tid]])
        {}
    }

    void unlock(int tid) {
        flag[mySlot[tid]] = false;
        flag[(mySlot[tid]+1)%numthreads] = true;
    }

};*/

class Array_lock {
    private:
    volatile bool* flag;
    std::atomic<int> tail;
    int numthreads;

    public:

    Array_lock(int n) : flag(new volatile bool[n]) {    
        for (int i = 0; i < n; ++i)
            flag[i] = false;
        flag[0] = true;
        tail = 0;
        numthreads = n;
    }

    void lock(int* mySlot) {
        *mySlot = tail.fetch_add(1)%numthreads;
        while (!flag[*mySlot])
        {}
    }

    void unlock(int* mySlot) {
        flag[*mySlot] = false;
        flag[(*mySlot+1)%numthreads] = true;
    }

};


class Array_lock_padded {
    private:
    volatile bool* flag;
    std::atomic<int> tail;
    int numthreads;

    public:

    Array_lock_padded(int n) : flag(new volatile bool[(n-1)*64+1]) {    
        for (int i = 0; i < n; ++i)
            flag[i*64] = false;
        flag[0] = true;
        tail = 0;
        numthreads = n;
    }

    void lock(int* mySlot) {
        *mySlot = tail.fetch_add(1)%numthreads;
        while (!flag[*mySlot*64])
        {}
    }

    void unlock(int* mySlot) {
        flag[*mySlot*64] = false;
        flag[((*mySlot*64)+64)%(numthreads*64)] = true;
    }

};

///////////////////////////////////////////////////////////// CLH Lock

class QNode {

    public:
    volatile bool locked;
    QNode* pred;
    
    QNode() {
        locked = false;
        pred = nullptr;
    }   
};
    
class CLH_lock {

    std::atomic<QNode*> tail;
    
    public:
    
    CLH_lock() {
        tail = new QNode;
    }

    void lock(QNode** pointerToNode) {
        QNode* node = new QNode;
        *pointerToNode = node;
        node->locked = true;
        node->pred = std::atomic_exchange(&tail,node);
        while (node->pred->locked) 
        {}
    }

    void unlock(QNode* node) {  
        delete node->pred;
        node->locked = false;
    }

};

///////////////////////////////////////////////////////////////////////// MCS Lock

class Node {
    
    private:
     std::atomic<bool> locked;
     std::atomic<Node*> next;

    public:

    Node() {
        next = nullptr;
        locked = false;
    }

    void setLocked(bool val) {
        this->locked = val;
    }
    void setNext(Node* val) {
        this->next = val;
    }
    bool getLocked() {
        return this->locked;
    }
    Node* getNext() {
        return this->next;
    }
};

class MCS_lock {

    public:
    std::atomic<Node*> tail;

    MCS_lock() {
        tail = nullptr;
    }

    void lock(Node* node) {
        Node* my = node;
        Node* pred = tail.exchange(my, std::memory_order_acquire);
        if (pred != nullptr) {
            my->setLocked(true);
            pred->setNext(my);
            while (my->getLocked())
            {}
        }
    }

    void unlock(Node* node) {
        Node* my = node;
        if (my->getNext() == nullptr) {
            Node* p = my;
            if (tail.compare_exchange_strong(p, nullptr, std::memory_order_release, std::memory_order_relaxed)) {
                return;
            }

            while (my->getNext() == nullptr)
            {}
        }
        my->getNext()->setLocked(false);
        my->setNext(nullptr);
    }

};

class execute {

public:
///////////////////////////////////////////////////////////// TAS

void run_TAS_lock(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    TAS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[i*8] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock.lock();
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[tid*8]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock();
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[i*8] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

///////////////////////////////////////////////////////////// TTAS

void run_TTAS_lock(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    TTAS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock.lock();
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock();
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in TTAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

///////////////////////////////////////////////////////////// Ticket

void run_Ticket_lock(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    Ticket_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock.lock();
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock();
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in Ticket Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

///////////////////////////////////////////////////////////// Array

void run_Array_lock(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    Array_lock mylock(numthreads);

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        thread_local int mySlot;
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock.lock(&mySlot);
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock(&mySlot);
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in Array Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

///////////////////////////////////////////////////////////// CLH

void run_CLH_lock(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    CLH_lock* mylock = new CLH_lock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        thread_local QNode* pointerToNode;
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock->lock(&pointerToNode);
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock->unlock(pointerToNode);
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in CLH Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

///////////////////////////////////////////////////////////// MCS

void run_MCS_lock(int numthreads, int iterations) {

    std::cout << "Hello" << std::endl;
    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    MCS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        thread_local Node my;
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock.lock(&my);
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock(&my);
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in MCS_Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

////////////////////////////////////////// Array Lock padded

void run_Array_lock_padded(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    Array_lock_padded mylock(numthreads);

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        thread_local int mySlot;
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            mylock.lock(&mySlot);
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock(&mySlot);
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in Array Lock padded with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

void run_Native_lock(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    omp_lock_t mylock;
    omp_init_lock(&mylock);

    for (int i = 0; i < numthreads; ++i)
        turns[i*8] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            omp_set_lock(&mylock);
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[tid*8]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            omp_unset_lock(&mylock);
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran in Native lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[i*8] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

void run_omp_critical(int numthreads, int iterations) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing

    for (int i = 0; i < numthreads; ++i)
        turns[i*8] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		    
        thread_local int mySlot;
        tid = omp_get_thread_num();
        while(counter < iterations)
        {
            #pragma omp critical 
            {
                try {
                    if(counter < iterations)
                    {
                        counter++;
                        turns[tid*8]++;
                        std::cout << "Thread " << tid << " is in critical section" << std::endl;
                    }
                }
                catch (int j) {
                    std::cout << "Some error occured while in CS" << std::endl;
                }
            }
        }
        
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << std::endl << "Program ran omp critical " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[i*8] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;

}

};


///////////////////////////////////////////////////////////// main starts here //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) 
{

    int mode = std::atoi(argv[1]);
    int numthreads = std::atoi(argv[2]);    // number of threads, command line input
    long int iterations = std::atoi(argv[3]);              // number of iterations in CS
		
    execute Locker;

    if (mode == 1) {
        Locker.run_TAS_lock(numthreads, iterations);
    }
    if (mode == 2) {
        Locker.run_TTAS_lock(numthreads, iterations);
    }
    if (mode == 3) {
        Locker.run_Ticket_lock(numthreads, iterations);
    }
    if (mode == 4) {
        Locker.run_Array_lock(numthreads, iterations);
    }
    if (mode == 5) {
        Locker.run_CLH_lock(numthreads, iterations);
    }
    if (mode == 6) {
        Locker.run_MCS_lock(numthreads, iterations);
    }
    if (mode == 7) {
        Locker.run_Array_lock_padded(numthreads, iterations);
    }
    if (mode == 8) {
        Locker.run_Native_lock(numthreads, iterations);
    }
    if (mode == 9) {
        Locker.run_omp_critical(numthreads, iterations);
    }
}


