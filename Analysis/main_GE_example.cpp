/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example GenerateEvents.cpp
 * This is an example of how to Generate Events using HEPfit.
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
        
        if(argc != 3){
            /* Print usage and exit. */
            if (rank == 0) std::cout << "\nusage: " << argv[0] << " ModelConf.conf no.of.iterations" << std::endl;
            return EXIT_SUCCESS;
        }

        /* Define the model configuration file.                                */
        /* Here it is passed as the first argument to the executable. The      */
        /* model configuration file provides the values with errors for the    */
        /* mandatory model parameters, as well as the list of observables,     */
        /* observables2D, correlated Gaussian observables.                     */
        /* See documentation for details.                                      */
        std::string ModelConf = argv[1];
        
        /* Define the number of iterations here.                               */
        int nIterations = atoi(argv[2]);
        
        /* Define the folder that will store the iterations                    */
        std::string FolderOut = "";
        
        /* Define the optional job tag. */
        std::string JobTag = "";
        
        /* Define the optional weight flag if you want weights. (Default: false)*/
        bool weight = false;
        
        /* Define the optional weight flag if you want weights. (Default 0)*/
        int seed = 1;
        
        /* Create objects of the classes ModelFactory and ThObsFactory */
        ThObsFactory ThObsF;
        ModelFactory ModelF;

        /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
        /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */
        
        /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
        /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/
        
        /* Create an object of the class GenerateEvents. */
        GenerateEvent GE(ModelF, ThObsF, ModelConf, FolderOut, JobTag);
        
        /* Initiate the Event Generation. */
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
