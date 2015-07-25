/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example myModel_MCMC.cpp
 * @example myModel.h
 * @example myModel.cpp
 * @example myObservables.h
 * @example myObservables.cpp
 * This is an example of how to add user-defined model and observables
 * and to perform a Bayesian analysis with the Markov Chain Monte Carlo.
 *
 */

#include <iostream>
#include <HEPfit.h>
#include <boost/bind.hpp>
#include "myModel.h"
#include "myObservables.h"

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
            if (rank == 0) std::cout << "\nusage: " << argv[0] << " ModelConf.conf MonteCarlo.conf\n" << std::endl;
            return EXIT_SUCCESS;
        }

        /* Define the model configuration file.                                */
        /* Here it is passed as the first argument to the executable. The      */
        /* model configuration file provides the values with errors for the    */
        /* mandatory model parameters, as well as the list of observables,     */
        /* observables2D, correlated Gaussian observables.                     */
        /* See documentation for details.                                      */
        std::string ModelConf = argv[1];
        
        /* Define the Monte Carlo configuration file.                         */
        /* Here it is passed as the second argument to the executable. The    */
        /* Monte Carlo configuration file provides the parameters used in the */
        /* Monte Carlo run. See documentation for details.                    */
        std::string MCMCConf = argv[2];
        
        /* Define the ROOT output file (w/o extension, empty string will set it to MCout) */
        std::string FileOut = "";        
        
        /* Define the optional job tag. */
        std::string JobTag = "";
        
        /* Create objects of the classes ModelFactory and ThObsFactory */
        ThObsFactory ThObsF;
        ModelFactory ModelF;
        myModel my_model;

        /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
        ModelF.addModelToFactory("myModel", boost::factory<myModel*>() );
        
        /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
        ThObsF.addObsToFactory("BIN1", boost::bind(boost::factory<yield*>(), _1, 1) );
        ThObsF.addObsToFactory("BIN2", boost::bind(boost::factory<yield*>(), _1, 2) );
        ThObsF.addObsToFactory("BIN3", boost::bind(boost::factory<yield*>(), _1, 3) );
        ThObsF.addObsToFactory("BIN4", boost::bind(boost::factory<yield*>(), _1, 4) );
        ThObsF.addObsToFactory("BIN5", boost::bind(boost::factory<yield*>(), _1, 5) );
        ThObsF.addObsToFactory("BIN6", boost::bind(boost::factory<yield*>(), _1, 6) );
        ThObsF.addObsToFactory("C_3", boost::factory<C_3*>() );
        ThObsF.addObsToFactory("C_4", boost::factory<C_4*>() );
        
        /* Create an object of the class MonteCarlo. */        
        MonteCarlo MC(ModelF, ThObsF, ModelConf, MCMCConf, FileOut, JobTag);
        
        /* Do a test run if you wish to see the values of the observables      */
        /* and the correlated Gaussian observables defined in the model        */
        /* configuration file computed with the central value of the mandatory */
        /* parameters defined in the same file.                                */
        // MC.TestRun(rank);
        
        /* Initiate the Mote Carlo run. */
        MC.Run(rank);

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
