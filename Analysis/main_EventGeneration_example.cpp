/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example EventGeneration.cpp
 * This is an example of how to perform a Bayesian Markov Chain Monte Carlo analysis with HEPfit and BAT.
 *
 */

#include <iostream>
#include <HEPfit.h>

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
#include <mpi.h>
#endif

int main(int argc, char** argv)
{

    /* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();
    MPI::Status status;
#else
    /* In our MPI implementation the process with rank 0 is the master. */
    int rank = 0;
#endif
    
    try {
        
        if(argc < 3){
            /* Print usage and exit. */
            if (rank == 0) std::cout << "\nusage: " << argv[0] << " ModelConf.conf nEvents [OutputFolder]" << std::endl;
            return EXIT_SUCCESS;
        }

        /* Define the model configuration file.                                */
        /* Here it is passed as the first argument to the executable. The      */
        /* model configuration file provides the values with errors for the    */
        /* mandatory model parameters, as well as the list of observables,     */
        /* observables2D, correlated Gaussian observables.                     */
        /* See documentation for details.                                      */
        std::string ModelConf = argv[1];
        
        /* Define the number of iterations.*/
        int nIterations = atoi(argv[2]);
        
        /* Define the output folder. Results will be put in ./GeneratedEvents/<outputFolder>*/
        std::string outputFolder = "";
        if(argc == 4) outputFolder = argv[3];
        
        /* Define the seed. 0 for using processor time.*/
        int seed = 1;
        
        /* Define whether weight will be calculated.*/
        bool weight = true;
        
        /* Define the optional job tag. */
        std::string JobTag = "";
        
        /* Create objects of the classes ModelFactory and ThObsFactory */
        ThObsFactory ThObsF;
        ModelFactory ModelF;

        /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
        /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */
        
        /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
        /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/
        
        /* Create an object of the class GenerateEvent. */
        GenerateEvent GE(ModelF, ThObsF, ModelConf, outputFolder, JobTag);
        
        /* Initiate the Mote Carlo run. */
        GE.generate(nIterations, seed, weight);

    /* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
        MPI::Finalize();
#endif
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
