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


///////////////////////////////////////////////// TAS Lock

class TAS_lock 
{
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


///////////////////////////////////////////////////////////// CLH Lock

class QNode {

    public:
    std::atomic<bool> locked;
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

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    std::string file_name;
    double totruntime = 0;

long int CS(long int counter, int iterations, double *turns,int tid)
{   
    double k = 0;
    //std::cout << "Thread " << tid << " is in critical section" << std::endl;
    try {
            if(counter < iterations)
            {
                counter++;    // function
                turns[std::max(tid*8-1,0)]++;
                //std::cout << "Thread " << tid << " is in critical section" << std::endl;
                for (int i = 0; i < 10000;i++)
                {
                    k = log(i);
                }
            }
        }
        catch (int j) {
            std::cout << "Some error occured while in CS" << std::endl;
        }

return counter;
}   

double Varianz(double *field,int len) // nur field // arrays für varianzen erstellen
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
    std::cout << std::endl;
    double runtime = 0;
    double Timevar = 0;
    double HowFair = 0;
    // setup timer variables
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    //std::cout << "numthreads" << numthreads << std::endl;
    double Verteilung[numthreads];
    double Time[numofiter];
    int len = 0;
    double var = 0;
    double turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    
    TAS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                mylock.lock();
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock();
            }           
        }
        end = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        Time[n] = runtime;

        for (int i = 0; i < numthreads; ++i)
        {
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }
        
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;

        len = sizeof(Verteilung)/sizeof(Verteilung[0]);
        var = Varianz(Verteilung,len);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
        
    }   
    Timevar = Varianz(Time,numofiter);
    output("TAS_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
    //std::cout << std::endl;
}

///////////////////////////////////////////////////////////// TTAS

void run_TTAS_lock(int numthreads, int iterations, int numofiter) 
{  
    double runtime = 0;
    double totruntime = 0;
    double HowFair = 0;
    double tmp = 0;
    double Timevar = 0;
    // setup timer variables

    //std::vector<int> Verteilung (numthreads,0);
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    double var = 0;
    int len = 0;
    double Verteilung[numthreads];
    double Time[numofiter];
    double turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    TTAS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {	
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                mylock.lock();
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock();
            }
            
        }
        end = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        Time[n] = runtime;
        
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;


        for (int i = 0; i < numthreads; ++i)
        {
            //std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }
        
        len = sizeof(Verteilung)/sizeof(Verteilung[0]);
        var = Varianz(Verteilung,len);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
        
    }
    for (int i = 0; i < numofiter; ++i)
    {
        //std::cout << "After Time n: " << i << " " << Time[i] << std::endl; 
    }
    
    Timevar = Varianz(Time,numofiter);
    output("TTAS_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// Ticket

void run_Ticket_lock(int numthreads, int iterations, int numofiter) 
{  
    double runtime = 0;
    double totruntime = 0;
    double HowFair = 0;
    double Timevar = 0;
    // setup timer variables

    //std::vector<int> Verteilung (numthreads,0);
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    double Verteilung[numthreads];
    double Time[numofiter];
    double var = 0;
    int len = 0;
    double turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    Ticket_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                mylock.lock();
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock();
            }
            
        }
            end = std::chrono::high_resolution_clock::now();
            runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
            Time[n] = runtime;
        
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;
        for (int i = 0; i < numthreads; ++i)
        {
            //std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }
        len = sizeof(Verteilung)/sizeof(Verteilung[0]);
        var = Varianz(Verteilung,len);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;     
    }   
    Timevar = Varianz(Time,numofiter);
    output("Ticket_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// Array

void run_Array_lock(int numthreads, int iterations, int numofiter) 
{    
    double runtime = 0;
    double totruntime = 0;
    double HowFair = 0;
    double Timevar = 0;
    // setup timer variables
    //std::vector<int> Verteilung (numthreads,0);
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    double Verteilung[numthreads];
    double Time[numofiter];
    double var = 0;
    int len = 0;
    double turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    Array_lock mylock(numthreads);

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {	
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            thread_local int mySlot;
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                mylock.lock(&mySlot);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(&mySlot);
            }
            
        }
            end = std::chrono::high_resolution_clock::now();
            runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
            Time[n] = runtime;
        
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;


        for (int i = 0; i < numthreads; ++i)
        {
            //std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }
        
        len = sizeof(Verteilung)/sizeof(Verteilung[0]);
        var = Varianz(Verteilung,len);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;  
    } 
    Timevar = Varianz(Time,numofiter);
    output("Array_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// CLH

void run_CLH_lock(int numthreads, int iterations, int numofiter) 
{
    double runtime = 0;
    double totruntime = 0;
    double HowFair = 0;
    double Timevar = 0;
    // setup timer variables
    //std::vector<int> Verteilung (numthreads,0);
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    double Verteilung[numthreads];
    double Time[numofiter];
    double var = 0;
    int len = 0;
    double turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    CLH_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {	
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            thread_local QNode* pointerToNode;
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                mylock.lock(&pointerToNode);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(pointerToNode);
            }      
        }
            end = std::chrono::high_resolution_clock::now();
            runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
            Time[n] = runtime;    
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;
        for (int i = 0; i < numthreads; ++i)
        {
            //std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }       
        len = sizeof(Verteilung)/sizeof(Verteilung[0]);
        var = Varianz(Verteilung,len);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
        
    }   
    Timevar = Varianz(Time,numofiter);
    output("CLH_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
}

///////////////////////////////////////////////////////////// MCS

void run_MCS_lock(int numthreads, int iterations, int numofiter) 
{   
    double runtime = 0;
    double totruntime = 0;
    double HowFair = 0;
    double tmp = 0;
    double Timevar = 0;
    // setup timer variables
    //std::vector<int> Verteilung (numthreads,0);
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    double Verteilung[numthreads];
    double Time[numofiter];
    double var = 0;
    int len = 0;
    double turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    MCS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {   
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            thread_local Node my;
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                mylock.lock(&my);
                counter = CS(counter,iterations,turns,tid);
                mylock.unlock(&my);
            }          
        }
        end = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        Time[n] = runtime;       
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;
        for (int i = 0; i < numthreads; ++i)
        {
            //std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }   
        len = sizeof(Verteilung)/sizeof(Verteilung[0]);
        var = Varianz(Verteilung,len);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
        
    }  
    Timevar = Varianz(Time,numofiter);
    output("MCS_lock",iterations,numthreads,numofiter,totruntime,HowFair,Timevar);
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
        }
        std::cout << std::endl;
    }
}

