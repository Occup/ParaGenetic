#include <ga/GASStateGA.h>
#include <ga/std_stream.h>
#include <iostream>
#include <cmath>
#include <mpi.h> // Necessary for running the code in parallel
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <string>
#include "PolarGenome.hpp"
#include <filesystem>
#include <boost/program_options.hpp>

using namespace std;
namespace po=boost::program_options;
namespace fs=std::filesystem;
// const double INC = 0.01;
unsigned int NUM_POLs = 0;
// Dimensions of search space

int main(int argc, char **argv) {
    
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // for (int k = 0; k < argc; ++k) cout << argv[k] << endl;
    // Get the number of processes
    int mpi_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
    
    // Get the rank of the process
    int mpi_rank;
    string rank;
    string param;
    stringstream convert;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    convert << mpi_rank;
    rank = convert.str();
    convert.str("");
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    // Termination parameters:
    int nGen;   // Total number of generations for termination
    unsigned int popSize; // Generation size
    int nMigr; // Number of migrating individuals
    double pMu;
    double pCr;
    string outdir;
    po::options_description argu("帮助");
    argu.add_options()
        ("help,h","打印帮助信息")
        ("outputDirectory,o", po::value<string>(&outdir)->default_value("."),"运算结果输出目录")
        ("NumPolarDimension,D", po::value<unsigned int>(&NUM_POLs)->default_value(11),"设置偏振探头数量")
        ("populationSize,p", po::value<unsigned int>(&popSize)->default_value(3000), "设置每个线程的种群规模")
        ("NumGenerations,g", po::value<int>(&nGen)->default_value(30), "遗传迭代次数的上限")
        ("NumMigration,n", po::value<int>(&nMigr)->default_value(10), "种群之间交换的最优染色体数量")
        ("MutationPossibility,m", po::value<double>(&pMu)->default_value(0.8), "种群中染色体发生变异的概率")
        ("CrossPossibility,c", po::value<double>(&pCr)->default_value(0.7), "种群的染色体之间发生交叉的概率")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, argu), vm);
    po::notify(vm);
    if (vm.count("help")) {
        if (mpi_rank == 0) cout << argu << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        // Finalize the MPI environment.
        MPI_Finalize();
        return 0;
    }
    if (mpi_tasks < 3) {
        if (mpi_rank == 0) cout << "请至少分配3线程来执行此程序：" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        // Finalize the MPI environment.
        MPI_Finalize();
        return 0;
    }
    // if (vm.count("populationSize")) cout << vm["populationSize"].as<unsigned int>() << endl;
    // if (vm.count("NumGenerations")) cout << vm["NumGenerations"].as<int>() << endl;
    // if (vm.count("NumMigration")) cout << vm["NumMigration"].as<int>() << endl;
    // if (vm.count("MutationPossibility")) cout << vm["MutationPossibility"].as<double>() << endl;
    // if (vm.count("CrossPossibility")) cout << vm["CrossPossibility"].as<double>() << endl;
    convert << "nProbe=" << NUM_POLs << "_popSize=" << popSize << "_nGen=" << nGen << "_nMigr=" << nMigr << "_pMu=" << pMu << "_pCr=" << pCr;
    // cout << convert.str();
    // cout << "here !!!" << endl;
    fs::path storeDir(outdir + "/" + convert.str());
    
    if (mpi_rank == 0) cout << "使用以下配置参数运行:" << endl
                            << "偏振探头数量 = " << NUM_POLs << endl
                            << "每线程种群规模 = " << popSize << endl
                            << "线程数 = " << mpi_tasks << endl
                            << "遗传代数上限 = " << nGen << endl
                            << "种群每代交换染色体个数 = " << nMigr << endl
                            << "变异率 = " << pMu << endl
                            << "交叉率 = " << pCr << endl
                            << "运算结果将保存到 " << string(storeDir) << " 目录" << endl;
    
    if (mpi_rank == 0 && (!fs::exists(storeDir))) fs::create_directories(storeDir);
    MPI_Barrier(MPI_COMM_WORLD);
    ofstream outfile; // Initialising the output stream
    ofstream imm;     // Initialising the output stream for immigrants

    // At the moment, the master node is not used other than for communications
    PolarGenome genome;

    // Creating the genetic algorithm and setting its parameters:
    GASteadyStateGA ga(genome); // Steady state genetic algorithm

    // Setting the terminator conditions: either due to convergence or by no.
    // generations (default):
    // ga.terminator(GAGeneticAlgorithm::TerminateUponGeneration);    //default
    // ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);
    // //alternative
    // ga.terminator(GAGeneticAlgorithm::TerminateUponPopConvergence);

    // ga.minimize(); // Minimization problem
    ga.populationSize(popSize);
    ga.nGenerations(nGen); // When TerminateUponGeneration
    ga.minimize();
    // ga.nConvergence(ncon);    // When TerminateUponConvergence or
    // ga.pConvergence(pcon);    // TerminateUponPopConvergence

    // Setting the type of cross-over:
    // ga.crossover(GARealArithmeticCrossover);
    // ga.crossover(GARealBlendCrossover);
    // The following options are not advised:
    // ga.crossover(GARealUniformCrossover);
    // ga.crossover(GARealEvenOddCrossover);
    // ga.crossover(GARealOnePointCrossover);
    // ga.crossover(GARealTwoPointCrossover);
    // ga.crossover(GARealPartialMatchCrossover);
    // ga.crossover(GARealOrderCrossover);
    // ga.crossover(GARealCycleCrossover);

    // Setting the variables for the mutation and crossover probabilities:
    ga.pMutation(pMu); // likelihood of mutating new offspring
    ga.pCrossover(pCr); // likelihood of crossing over parents

    // Specifying settings for calculating the fitness function:
    string scorefile = string(storeDir) + "/score_of_generation" + rank + ".dat";
    ga.scoreFilename(scorefile.c_str()); // name of output file for scores
    ga.scoreFrequency(1);               // keep the scores of every generation
    ga.flushFrequency(1);               // specify how often to write the score
    // to disk
    ga.selectScores(GAStatistics::AllScores); /* writing out average, maximum,
                                                 minimum, standard deviation and
                                                 divergence values of the fitness
                                                 function for later plotting. */

    // Initializing the genetic algorithm:
    // cout << "Thread " + rank << " prepared for initialization..." << endl;
    ga.initialize();
    // cout << "Thread " + rank << " successfully initialized!" << endl;
    // Output the initial population to file:
    int counter = 0;
    string pop_evo = string(storeDir) +"/pop_evolution" + rank + ".dat";
    outfile.open(pop_evo.c_str());

    double emigrant[nMigr][2][NUM_POLs];
    double immigrant[nMigr][2][NUM_POLs];
    string immigrants = string(storeDir) + "/immigrants" + rank + ".dat";
    imm.open(immigrants.c_str());

    // Running the algorithm until the termination condition is met:
    while (!ga.done()) {
        counter++;

        outfile << endl << "Generation = " << counter << endl;
        for (int ii = 0; ii < popSize; ii++) {
            outfile << "No. " << ii+1 << " individual" << endl;
            genome = ga.population().individual(ii);
            outfile << genome;
        }
        // ga.step();
        // Printing the population to file for each generation for later plotting:
        GAPopulation swappop(ga.population());
            
        /* Migration - stepping stone approach, i.e. individuals are swapped
           between neighbouring populations. The best individuals of one
           population are copied across to the next one in stead of the worst
           performing ones:*/
        for (int i = 0; i < nMigr; i++) {
            genome = ga.population().individual(i);
            for (int j = 0; j < NUM_POLs; j++) {
                emigrant[i][0][j] = genome(0, j);
                emigrant[i][1][j] = genome(1, j);
            }
        }
        // The following if-loops are required to avoid dead-lock
        if (mpi_rank == (mpi_tasks - 1)) {
            cout << "Thread " + rank << " on Generation " << counter
                 << " emitted " <<(counter - 1) * (mpi_tasks - 1) << " to thread 0"
                 << " waiting for " << (counter - 1) * (mpi_tasks - 1) + mpi_rank << " from thread " << (mpi_rank - 1) << endl;
            MPI_Recv(&immigrant, (nMigr * 2 * NUM_POLs), MPI_DOUBLE, (mpi_rank - 1),
                     (counter - 1) * (mpi_tasks - 1) + mpi_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&emigrant, (nMigr * 2 * NUM_POLs), MPI_DOUBLE, 0, (counter - 1) * (mpi_tasks - 1), MPI_COMM_WORLD);
                
        } else if (mpi_rank == 0) {
            cout << "Thread " + rank << " on Generation " << counter
                 << " emitted " << (counter - 1) * (mpi_tasks - 1) + (mpi_rank + 1) << " to thread " << (mpi_rank + 1)
                 << " waiting for " << (counter - 1) * (mpi_tasks - 1) << " from thread " << (mpi_tasks - 1) << endl;
            MPI_Send(&emigrant, (nMigr * 2 * NUM_POLs), MPI_DOUBLE, (mpi_rank + 1),
                     (counter - 1) * (mpi_tasks - 1) + (mpi_rank + 1), MPI_COMM_WORLD);
            MPI_Recv(&immigrant, (nMigr * 2 * NUM_POLs), MPI_DOUBLE, (mpi_tasks - 1), (counter - 1) * (mpi_tasks - 1),
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            cout << "Thread " + rank << " on Generation " << counter
                 << " emitted " << (counter - 1) * (mpi_tasks - 1) + (mpi_rank + 1) << " to thread " << (mpi_rank + 1)
                 << " waiting for " << (counter - 1) * (mpi_tasks - 1) + mpi_rank << " from thread " << (mpi_rank - 1) << endl;
            MPI_Recv(&immigrant, (nMigr * 2 * NUM_POLs), MPI_DOUBLE, (mpi_rank - 1),
                     (counter - 1) * (mpi_tasks - 1) + mpi_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&emigrant, (nMigr * 2 * NUM_POLs), MPI_DOUBLE, (mpi_rank + 1),
                     (counter - 1) * (mpi_tasks - 1) + (mpi_rank + 1), MPI_COMM_WORLD);
        }
        // MPI_Barrier(MPI_COMM_WORLD);
        imm << endl <<"Counter = " << counter << endl;
        for (int i = 0; i < nMigr; i++) {
            imm << "No. " << i+1 << " immigrant" << endl;
            genome.read_params(immigrant[i][0], immigrant[i][1], NUM_POLs);
            swappop.remove();
            swappop.add(genome);
            // substituting the last (worst)
            // individuals with the best ones of
            // another population
                    
            // imm << "\t" << immigrant[i][j];
            imm << genome;
        }

        ga.population(swappop);
        ga.step();
    }
    imm.close();     // Closing the output file stream
    outfile.close(); // Closing the output file stream
    // Printing out the best genome that the GA found for the current population
    //  to file:
    string best_val = string(storeDir) + "/best_individual" + rank + ".dat";
    outfile.open(best_val.c_str());
    genome = ga.statistics().bestIndividual();
    outfile << genome;
    outfile.close(); // closing the output file stream
    // cout << "Thread " + rank
    //      << " emitted " << (mpi_rank) << " to thread 0" << endl;
    if (mpi_rank != 0) {
        double best_individual[2][NUM_POLs];
        for (int kk = 0; kk < NUM_POLs; ++kk) {
            best_individual[0][kk] = genome(0, kk);
            best_individual[1][kk] = genome(1, kk);
        }
        MPI_Send(&best_individual, 2*NUM_POLs, MPI_DOUBLE, 0, mpi_rank, MPI_COMM_WORLD);
    }
    else {
        double individual[mpi_tasks][2][NUM_POLs];
        int Map[mpi_tasks];
        vector<Polar> list;
        MPI_Status status;
        
        for (int i = 0; i < (mpi_tasks - 1); i++) {
            // cout << "Thread " + rank 
            //      << " waiting for No." << i+1 << " Msg from peer threads... " << endl;
            MPI_Recv(&individual[i], NUM_POLs*2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            Map[i] = status.MPI_SOURCE;
        }

        for (int kk = 0; kk < NUM_POLs; ++kk) {
            individual[mpi_tasks - 1][0][kk] = genome(0, kk);
            individual[mpi_tasks - 1][1][kk] = genome(1, kk);
        }
        Map[mpi_tasks - 1] = 0;

        double minimum = numeric_limits<double>::max();
        int pop = 0;
        for (int i = 0; i < mpi_tasks; i++) {
            Polar tmp(individual[i][0], individual[i][1], NUM_POLs);
            list.push_back(tmp);
            // Selecting the actual minimum:
            if (minimum > tmp.ErrorEstimation()) {
                minimum = tmp.ErrorEstimation();
                pop = i;
            }
        }
        /*Printing the data to a file as well so that it can be used as the
          starting point by the Nelder-Mead Simplex method in the hybrid
          optimization technique.*/
        Polar Best(list[pop]);
        Best.NeighborDownHill();
        outfile.open(string(storeDir) + "/best_individual.dat");
        if (Best.ErrorEstimation() > list[pop].ErrorEstimation()) outfile << list[pop];
        else outfile << Best;
        outfile << "Belonging to population " << Map[pop] << endl;
        outfile.close(); // Closing the output file stream*/
        cout << endl << "Best individual over all is:" << endl
             << list[pop] 
             << "Belonging to population " << Map[pop] << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
