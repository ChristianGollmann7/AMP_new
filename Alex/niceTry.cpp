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

class LOCK 
{
    std::atomic<bool> state;
    volatile bool* flag;
    std::atomic<int> tail;
    std::atomic<int> ticket;
    volatile int served;
    
    public:
            // setup timer variables
        std::chrono::time_point<std::chrono::high_resolution_clock> start;
        std::chrono::time_point<std::chrono::high_resolution_clock> end;
        std::string file_name;
        double runtime = 0;
        double totruntime = 0;
        double HowFair = 0;
        double mean = 0;
        double var = 0;
        double tmp = 0;
        int tid;                                // thread ID
        long int counter = 0;                   // counter gets incremented in CS    

        std::atomic<bool> locked;
        LOCK* pred;

    LOCK()
    {
        state = false;
        ticket = 0;
        served = 0;
        locked = false;
        pred = nullptr;
    }
    
    void TAS_lock(){
        while (state.exchange(true))
        {}
    }

    void Ticket_lock() 
    {
        int next = ticket.fetch_add(1);
        while (served < next)
        {}
    }

    void TTAS_lock(){
        while (true) {
            while (state.load()) 
            {}
            if (!state.exchange(true))
                return;
        }
    }

    void unlock()
    {
        state.exchange(false);
    }

    void Ticket_unlock() 
    {
        served++;
    }

    

    void output(std::string Lockname, int iterations, int numthreads, int numofiter, double totruntime, double HowFair)
    {
        file_name = "data/Max"+std::to_string(iterations)+"_"+std::to_string(numthreads)+Lockname+".csv";
        std::ofstream myfile;
        myfile.open (file_name);
        myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair" << "; \n";
        myfile << iterations << ";" << numthreads << ";" << totruntime/numofiter << ";" << HowFair/numofiter << ";";	  
        myfile << std::endl;
        myfile.close();
    }

    double statistic(double tmp, int numthreads, double HowFair, double runtime)
    {   
        double temp;
        var = sqrt(tmp/numthreads);
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
        std::cout << "FunctionHowFair: " << HowFair << std::endl;
        temp = HowFair;
        return temp,totruntime;
    }

    void run_TAS_lock(int numthreads, int iterations, int numofiter) 
    {
        // setup timer variables
        std::vector<int> Verteilung (numthreads,0);
        double HowFair_funkt = 0;
        long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    
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
                    TAS_lock();
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
                    unlock();
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
            std::cout << "HowFair_bevor_funkt: " << HowFair << std::endl;
            HowFair,totruntime = statistic(tmp,numthreads,HowFair,runtime);
            std::cout << "HowFair_funkt: " << HowFair << std::endl;
            //std::cout << "Var: " << var << " HowFair: " << HowFair << " totruntime: " << totruntime << std::endl;
            //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
            //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
            
        }
        output("TAS_lock",iterations,numthreads,numofiter,totruntime,HowFair_funkt);
    }
    void run_TTAS_lock(int numthreads, int iterations, int numofiter) 
    {
        // setup timer variables                              // thread ID
        long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
        std::vector<int> Verteilung (numthreads,0);
        for (int i = 0; i < numthreads; ++i)
            turns[std::max(i*8-1,0)] = 0;
        

        omp_set_num_threads(numthreads);        // setting number of threads
        
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < numofiter; n++)
        {
            #pragma omp parallel private(tid) shared(counter)
            {		    
                tid = omp_get_thread_num();
                while(counter < iterations)
                {
                    TTAS_lock();
                    try 
                    {
                        if(counter < iterations)
                        {
                            counter++;
                            turns[std::max(tid*8-1,0)]++;
                            std::cout << "Thread " << tid << " is in critical section" << std::endl;
                        }
                    }
                    catch (int j)
                    {
                        std::cout << "Some error occured while in CS" << std::endl;
                    }
                    unlock();
                }       
            }
            end = std::chrono::high_resolution_clock::now();
            runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
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
            HowFair,totruntime = statistic(tmp,numthreads,HowFair,runtime);
                //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
        }       //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
        output("TTAS_lock",iterations,numthreads,numofiter,totruntime,HowFair);
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
    int numofiter = std::atoi(argv[4]);


		
    LOCK Locker;
    //for (int iterations = 0; iterations < maxiterations; iterations = iterations + 1)
    //{
        if (mode == 1) 
        {
            Locker.run_TAS_lock(numthreads, iterations, numofiter);
        }          
        if (mode == 2) 
        {
            Locker.run_TTAS_lock(numthreads, iterations, numofiter);
        }
}

