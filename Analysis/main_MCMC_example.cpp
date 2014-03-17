/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <MonteCarlo.h>
#ifdef _MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) 
{

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
            //print usage and exit
            std::cout << "usage: " << argv[0] << " ModelConf MCMCConf" << std::endl;
            return EXIT_SUCCESS;
        }

        /* Define the model configuration file.                                */
        /* Here it is passed as the first argument to the executable. The      */
        /* model configuration file provides the values with errors for the    */
        /* mandatory model parameters, as well as the list of observables,     */
        /* observables2D, correlated Gaussian observables, model parameters vs.*/ 
        /* observables. See documentation for details.                         */
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
        
        /* Define the check theory range flag. */
        bool checkTheoryRange = false;

        /* Create an object of the class MonteCarlo. */        
        MonteCarlo MC(ModelConf, MCMCConf, FileOut, JobTag, checkTheoryRange);
        
        /* Do a test run if you wish to see the values of the observables      */
        /* and the correlated Gaussian observables defined in the model        */
        /* configuration file computed with the central value of the mandatory */
        /* parameters defined in the same file.                                */
        // MC.TestRun(rank);
        
        /* Initiate the Mote Carlo run. */
        MC.Run(rank);
        
#ifdef _MPI
        MPI::Finalize();
#endif
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
