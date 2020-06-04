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
///////////////////////////////////////////////////////////// TAS

void run_TAS_lock(int numthreads, int iterations, int numofiter) {

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime = 0;
    double totruntime = 0;
    std::string file_name;

    double HowFair = 0;
    double mean = 0;
    double var = 0;
    double tmp = 0;
    std::vector<int> Verteilung (numthreads,0);


    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    TAS_lock mylock;


    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {
        tmp = 0;
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
                        counter++;    // function
                        turns[std::max(tid*8-1,0)]++;
                        std::cout << "Thread " << tid << " is in critical section" << std::endl;
                    }
                }
                catch (int j) {
                    //std::cout << "Some error occured while in CS" << std::endl;
                }
                mylock.unlock();
            }
            
        }
        end = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;
        for (int i = 0; i < numthreads; ++i)
        {
            std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            //std::cout << i <<  std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }

        mean = iterations/numthreads;

        for (int i = 0; i < numthreads; ++i)
        {
            tmp = tmp + (Verteilung[i] - mean) * (Verteilung[i] - mean);           
        }

        var = sqrt(tmp/numthreads);
        std::cout << var << std::endl;
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
        
    }
    file_name = "data/Max"+std::to_string(iterations)+"_"+std::to_string(numthreads)+"TAS_lock.csv";
        //return results;
    std::ofstream myfile;
    myfile.open (file_name);
    myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair" << "; \n";
    myfile << iterations << ";" << numthreads << ";" << totruntime/numofiter << ";" << HowFair/numofiter << ";";	  
    myfile << std::endl;

    myfile.close();
}

///////////////////////////////////////////////////////////// TTAS

void run_TTAS_lock(int numthreads, int iterations) 
{

    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    std::vector<double> results; // results array
    double runtime;
    std::string file_name;
    double HowFair;
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    TTAS_lock mylock;
    results = std::vector<double>(3, 0);  // Creates a result array

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

       results[0] = numthreads;
    results[1] = iterations;
    results[2] = results[numthreads];
    results[3] = HowFair;

    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
    results[numthreads] = runtime;
    file_name = "data/"+std::to_string(numthreads)+"_"+std::to_string(iterations)+"_"+"TTAS_lock.csv";
    std::ofstream myfile;
    myfile.open (file_name);
    myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair" << "; \n";
    myfile << numthreads << ";" << iterations << ";" << results[numthreads] << ";" << HowFair << ";";	  
    myfile << std::endl;
    myfile.close();

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

    double HowFair;
    double mean = 0;
    double var = 0;
    double tmp = 0;
    std::vector<int> Verteilung (numthreads,0);

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
    CLH_lock mylock;

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
            mylock.lock(&pointerToNode);
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
            mylock.unlock(pointerToNode);
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

void run_MCS_lock(int numthreads, int iterations, int numofiter) {

    std::cout << "Hello" << std::endl;
    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;
    double czero = 0;
    std::string file_name;
    double HowFair;

    double mean = 0;
    double var = 0;
    double tmp = 0;
    std::vector<int> Verteilung (numthreads,0);

    std::vector <std::vector<double>> results(numofiter, std::vector<double>(4, 0));
    int tid;                                // thread ID
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    MCS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	for (int n = 0; n < numofiter; n++)
    {   
        czero = 0;
        tmp = 0;
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
                    //std::cout << "Some error occured while in CS" << std::endl;
                }
                mylock.unlock(&my);
            }
            
        }
        end = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        for (int i = 0; i < numthreads; ++i)
        {
            //std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            //std::cout << i <<  std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }

        mean = iterations/numthreads;

        for (int i = 0; i < numthreads; ++i)
        {
            tmp = tmp + (Verteilung[i] - mean) * (Verteilung[i] - mean);           
        }

        var = sqrt(tmp/numthreads);
        
        HowFair = var;

        results[n][0] = numthreads;
        results[n][1] = runtime;
        results[n][2] = var;
        //std::cout << std::endl << "Program ran in MCS_Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;
        for (int i = 0; i < numthreads; ++i)
        {   
            std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
        }
        std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
       // results[iterations][2] = 1/(1 + czero);
        //std::cout << "iterations " << iterations << "HowFair " << results[iterations][2] << "zeros " << czero << std::endl;
    }
    file_name = "data/Max"+std::to_string(iterations)+"_"+std::to_string(numthreads)+"MCS_lock.csv";
    std::ofstream myfile;
    myfile.open (file_name);
    myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair" << "; \n";
    for (int n = 0; n < numofiter; n++)
    {
        myfile << iterations << ";" << results[n][0] << ";" << results[n][1] << ";" << results[n][2] << ";";	  
        myfile << std::endl;
    }
    myfile.close();
}

};


///////////////////////////////////////////////////////////// main starts here //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) 
{

    std::string file_name;
    std::vector<double> results;
    int mode = std::atoi(argv[1]);
    int numthreads = std::atoi(argv[2]);    // number of threads, command line input
    //long int maxiterations = std::atoi(argv[3]);              // number of iterations in CS
    long int iterations = std::atoi(argv[3]);
    results = std::vector<double>(3, 0);
    int numofiter = 50;
		
    execute Locker;
    //for (int iterations = 0; iterations < maxiterations; iterations = iterations + 1)
    //{
        if (mode == 1) {
            Locker.run_TAS_lock(numthreads, iterations, numofiter);
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
            Locker.run_MCS_lock(numthreads, iterations, numofiter);
        }
    //}
    /*file_name = "data/"+std::to_string(numthreads)+"TAS_lock.csv";
    std::ofstream myfile;
    myfile.open (file_name);
    myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair \n";
    myfile << numthreads << ";" << iterations << ";" << results[numthreads] << ";" << HowFair;	  
    myfile << std::endl;
    myfile.close();
    */
}

