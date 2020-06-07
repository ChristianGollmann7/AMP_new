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
#include <string>
#include <fstream>

#define SETUP double totruntime = 0; \
              double Timevar = 0; \
              double HowFair = 0; \
              int tid; \
              long int counter = 0; \
              double Verteilung[numthreads]; \
              double Time[numofiter]; \
              int len = 0; \
              double var = 0; \
              double turns[(numthreads-1)*8+1];

#define INIT_TURNS for (int i = 0; i < numthreads; ++i) \
                   turns[i*8] = 0;

#define POSTPROCESSING for (int i = 0; i < numthreads; ++i) \
                       { \
                           Verteilung[i] = turns[i*8]; \
                       } \
                       len = sizeof(Verteilung)/sizeof(Verteilung[0]); \
                       var = Varianz(Verteilung,len); \
                       HowFair = HowFair + var; \
                       totruntime = totruntime + runtime; 


///////////////////////////////////////////////// TAS Lock

class TAS_lock 
{
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

class TTAS_lock 
{
    private:    
    std::atomic<bool> state;
    
    public:
    TTAS_lock(){
        state = false;
    }  
    void lock(){
        while (true) {
            while (state) 
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

class Ticket_lock 
{
    private:
    std::atomic<int> ticket;
    volatile int served;

    public:  
    Ticket_lock() {
        ticket = 0;
        served = 0;
    }
    void lock() 
    {
        int next = ticket.fetch_add(1);
        while (served < next)
        {}
    }
    void unlock() 
    {
        served++;
    }
};

///////////////////////////////////////////////// Array Lock

class Array_lock 
{
    private:
    volatile bool* flag;
    std::atomic<int> tail;
    int numthreads;

    public:
    Array_lock(int n) : flag(new volatile bool[n]) 
    {    
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

class Array_lock_padded 
{
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

class QNode 
{
    public:
    std::atomic<bool> locked;
    QNode* pred;
    
    QNode() {
        locked = false;
        pred = nullptr;
    }   
};
    
class CLH_lock 
{
    private:
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

class Node 
{
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

class MCS_lock 
{
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

//////////////////////////////////////////////////////////////////// execute

class execute {

public:
    std::string file_name;

long int CS(long int counter, int iterations, double *turns,int tid)
{   
    double k = 0;
    try {
            if(counter < iterations)
            {
                counter++;
                turns[tid*8]++;
                /*for (int i = 0; i < 10000;i++)
                {
                    k = log(i);
                }*/
            }
        }
        catch (int j) {
            std::cout << "Some error occured while in CS" << std::endl;
        }

return counter;
}   

double Varianz(double *field,int len)
{
        double mean = 0;
        double var = 0;
        double tmp = 0;
        for (int i = 0; i < len; ++i)
        {
            mean = mean + field[i];
        }
        mean = mean/len;
        for (int i = 0; i < len; ++i)
        {
            tmp = tmp + (field[i] - mean) * (field[i] - mean);         
        }
        var = sqrt(tmp/len);
    return var;
}

void output(std::string Lockname, int iterations, int numthreads, int numofiter, double totruntime, double HowFair, double Timevar)
{
    file_name = "data/"+Lockname+"_"+std::to_string(numthreads)+"_"+std::to_string(iterations)+".csv";
    std::ofstream myfile;
    myfile.open (file_name);
    myfile << "NameOfLock" << ";" <<"NumberOfIterations" << ";" << "NumberOfThreads"<< ";" << "RunTime" << ";" << "HowFair" << ";" << "Timevar" <<"\n";
    myfile << Lockname << ";" << iterations << ";" << numthreads << ";" << totruntime/numofiter << ";" << HowFair/numofiter << ";" << Timevar/numofiter;	  
    myfile << std::endl;
    myfile.close();
}
///////////////////////////////////////////////////////////// TAS

void run_TAS_lock(int numthreads, int iterations, int numofiter) 
{   
    SETUP
    
    TAS_lock mylock;

    INIT_TURNS
    

    omp_set_num_threads(numthreads);
	for (int n = 0; n < numofiter; n++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock();
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock();
            }     
            if (tid==0)  
                end = std::chrono::high_resolution_clock::now();    
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;
     
        POSTPROCESSING
    }   
    Timevar = Varianz(Time,numofiter);
    output("TAS_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// TTAS

void run_TTAS_lock(int numthreads, int iterations, int numofiter) 
{  
    SETUP

    TTAS_lock mylock;

    INIT_TURNS
    
    omp_set_num_threads(numthreads); 
	for (int n = 0; n < numofiter; n++)
    {	
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock();
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock();
            }
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();
            
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;

        POSTPROCESSING
    }   
    Timevar = Varianz(Time,numofiter);
    output("TTAS_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// Ticket

void run_Ticket_lock(int numthreads, int iterations, int numofiter) 
{  
    SETUP

    Ticket_lock mylock;

    INIT_TURNS
    

    omp_set_num_threads(numthreads);
	for (int n = 0; n < numofiter; n++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock();
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock();
            }
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();
            
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;
            
        POSTPROCESSING   
    }   
    Timevar = Varianz(Time,numofiter);
    output("Ticket_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// Array

void run_Array_lock(int numthreads, int iterations, int numofiter) 
{    
    SETUP

    Array_lock mylock(numthreads);

    INIT_TURNS
    

    omp_set_num_threads(numthreads);       
	for (int n = 0; n < numofiter; n++)
    {	
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            thread_local int mySlot;
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock(&mySlot);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(&mySlot);
            }
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();
            
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;
        
        POSTPROCESSING 
    } 
    Timevar = Varianz(Time,numofiter);
    output("Array_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// Array Lock padded

void run_Array_lock_padded(int numthreads, int iterations, int numofiter) 
{    
    SETUP

    Array_lock_padded mylock(numthreads);

    INIT_TURNS
    

    omp_set_num_threads(numthreads);    
	for (int n = 0; n < numofiter; n++)
    {	
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            thread_local int mySlot;
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock(&mySlot);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(&mySlot);
            }
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();
            
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;
        
        POSTPROCESSING 
    } 
    Timevar = Varianz(Time,numofiter);
    output("Array_lock_padded",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// CLH

void run_CLH_lock(int numthreads, int iterations, int numofiter) 
{
    SETUP

    CLH_lock mylock;

    INIT_TURNS
    

    omp_set_num_threads(numthreads);
	for (int n = 0; n < numofiter; n++)
    {	
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            thread_local QNode* pointerToNode;
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock(&pointerToNode);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(pointerToNode);
            }  
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();    
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime; 
   
        POSTPROCESSING
        
    }   
    Timevar = Varianz(Time,numofiter);
    output("CLH_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// MCS

void run_MCS_lock(int numthreads, int iterations, int numofiter) 
{   
    SETUP

    MCS_lock mylock;

    INIT_TURNS
    

    omp_set_num_threads(numthreads);
	for (int n = 0; n < numofiter; n++)
    {   
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            thread_local Node my;
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                mylock.lock(&my);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(&my);
            } 
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();         
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;      
 
        POSTPROCESSING
        
    }  
    Timevar = Varianz(Time,numofiter);
    output("MCS_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

////////////////////////////////////////////////////////////// Native Lock

void run_Native_lock(int numthreads, int iterations, int numofiter) 
{   
    SETUP
    
    omp_lock_t mylock;
    omp_init_lock(&mylock);

    INIT_TURNS
    

    omp_set_num_threads(numthreads);    
	for (int n = 0; n < numofiter; n++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                omp_set_lock(&mylock);
                counter = CS(counter,iterations,turns,tid);
                omp_unset_lock(&mylock);
            }  
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();         
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;

        POSTPROCESSING
        
    }   
    Timevar = Varianz(Time,numofiter);
    output("Native_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// omp critical

void run_omp_critical(int numthreads, int iterations, int numofiter) 
{   
    SETUP
    
    INIT_TURNS
    
    omp_set_num_threads(numthreads);     
	for (int n = 0; n < numofiter; n++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter,start,end)
        {		    
            tid = omp_get_thread_num();
            #pragma omp barrier
            if (tid==0)
                start = std::chrono::high_resolution_clock::now();
            while(counter < iterations)
            {
                # pragma omp critical
                {
                    counter = CS(counter,iterations,turns,tid);
                }
            } 
            if (tid==0)
                end = std::chrono::high_resolution_clock::now();          
        }
        double runtime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Time[n] = runtime;

        POSTPROCESSING
        
    }   
    Timevar = Varianz(Time,numofiter);
    output("omp_critical",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

};


///////////////////////////////////////////////////////////// main starts here //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) 
{
    int maxmode = std::atoi(argv[1]);
    int numthreadsmax = std::atoi(argv[2]);    // number of threads, command line input
    long int iterations = std::atoi(argv[3]); // number of iterations in CS
    int numofiter = std::atoi(argv[4]); // number of iterations in CS
    int numthreads = 0;
    int mode = 0;
		
    execute Locker;
    for(int j = 1; j < maxmode+1; j++)
    {   
        mode = j;
        for (int i = 1; i < numthreadsmax+1; i++)
        {
            std::cout << "LockNumber: " << j << " ThreadNumber: " << i << std::endl;
            numthreads = i;
            if (mode == 1) {
                Locker.run_TAS_lock(numthreads, iterations, numofiter);
            }
            if (mode == 2) {
                Locker.run_TTAS_lock(numthreads, iterations, numofiter);
            }
            if (mode == 3) {
                Locker.run_Ticket_lock(numthreads, iterations, numofiter);
            }
            if (mode == 4) {
                Locker.run_Array_lock(numthreads, iterations, numofiter);
            }
            if (mode == 5) {
                Locker.run_CLH_lock(numthreads, iterations, numofiter);
            }
            if (mode == 6) {
                Locker.run_MCS_lock(numthreads, iterations, numofiter);
            }
            if (mode == 7) {
                Locker.run_Array_lock_padded(numthreads, iterations, numofiter);
            }
            if (mode == 8) {
                Locker.run_Native_lock(numthreads, iterations, numofiter);
            }
            if (mode == 9) {
                Locker.run_omp_critical(numthreads, iterations, numofiter);
            }
        }
        std::cout << std::endl;
    }
}

