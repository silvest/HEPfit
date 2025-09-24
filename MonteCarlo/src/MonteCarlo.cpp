/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MonteCarlo.h"
#include "ModelFactory.h"
#include "ThObsFactory.h"
#include <BAT/BCParameter.h>
#include <TSystem.h>
#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#ifdef _MPI
#include <mpi.h>
#endif
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <ctime>
#include <iostream>
#include <iomanip>


MonteCarlo::MonteCarlo(
        ModelFactory& ModelF, ThObsFactory& ThObsF,
        const std::string& ModelConf_i,
        const std::string& MonteCarloConf_i,
        const std::string& OutFile_i,
        const std::string& JobTag_i,
        const std::string& FileIn_i)
: myInputParser(ModelF, ThObsF), MCEngine(ModPars, Obs, Obs2D, CGO, CGP) {
    ModelConf = ModelConf_i;
    MCMCConf = MonteCarloConf_i;
    JobTag = JobTag_i;
    FileIn = FileIn_i;
    if (OutFile_i.compare("") == 0) OutFile = "MCout" + JobTag + ".root";
    else OutFile = OutFile_i + JobTag + ".root";
    ObsDirName = "Observables" + JobTag;
    FindModeWithMinuit = false;
    RunMinuitOnly = false;
    ComputeHessianOnly = false;
    CalculateNormalization = "false";
    NIterationNormalizationMC = 0;
    PrintAllMarginalized = false;
    PrintCorrelationMatrix = false;
    PrintKnowledgeUpdatePlots = false;
    PrintParameterPlot = false;
    PrintTrianglePlot = false;
    WritePreRunData = false;
    checkrun = false;
    writechains = false;
}

MonteCarlo::~MonteCarlo() 
{
    Obs.clear();
    boost::ptr_vector<Observable>().swap(Obs);
}

void MonteCarlo::TestRun(int rank) {
    if (checkrun == true) {
        if (rank == 0) throw std::runtime_error("ERROR: MonteCarlo::TestRun() cannot be called after calling MonteCarlo::Run().\n");
        else sleep (2);
    } else {
        checkrun = true;
    }

    if (rank == 0) {
        ModelName = myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO, CGP);
        std::map<std::string, double> DP;
        if (Obs.size() == 0 && CGO.size() == 0) throw std::runtime_error("\nMonteCarlo::TestRun(): No observables or correlated Gaussian observables defined in SomeModel.conf file\n");

        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            if (it->IsCorrelated()) {
                for (unsigned int i = 0; i < CGP.size(); i++) {
                    if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                        std::string index = it->getname().substr(CGP[i].getName().size());
                        long int lindex = strtol(index.c_str(), NULL, 10);
                        if (lindex > 0)
                            DP[CGP[i].getPar(lindex - 1).getname()] = CGP[i].getPar(lindex - 1).getave();
                        else {
                            std::stringstream out;
                            out << it->getname();
                            throw std::runtime_error("MonteCarlo::Run(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                        }
                    }
                }
            } else
                DP[it->getname()] = it->getave();
        }
        
        if (!myInputParser.getModel()->Init(DP)) {
            if (myInputParser.getModel()->getmissingModelParameters().size() > 0) {
                std::cout << "\nPlease set the following parameters in the model configuration files:\n" << std::endl;
                std::vector<std::string> missingParameters = myInputParser.getModel()->getmissingModelParameters(); // Store the result in a variable to avoid dangling pointer warning 
                for (std::vector<std::string>::iterator it = missingParameters.begin(); it != missingParameters.end(); it++) {
                    std::cout << "ModelParameter\t" << *it << std::endl;
                }
                std::cout << std::endl;
                myInputParser.getModel()->getmissingModelParameters().clear();
            }
            if (rank == 0) throw std::runtime_error("ERROR: Parameter(s) missing in model initialization. \n");
            else sleep (2);
        }
        std::vector<std::string> unknownParameters = myInputParser.getModel()->getUnknownParameters();
        if (unknownParameters.size() > 0 && rank == 0) {
            std::cout << "\n" << std::endl;
            for (std::vector<std::string>::iterator it = unknownParameters.begin(); it != unknownParameters.end(); it++)
                std::cout << "WARNING: unknown parameter " << *it << " not added to MCMC Test Run" << std::endl;
        }
        if (Obs.size() > 0) std::cout << "\nObservables: \n" << std::endl;
        for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
            double th = it->computeTheoryValue();
            std::cout << it->getName() << " = " << th << std::endl;
        }

        if (CGO.size() > 0) std::cout << "\nCorrelated Gaussian Observables: \n" << std::endl;
        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
            std::cout << it1->getName() << ":" << std::endl;
            std::vector<Observable> ObsInCGO = it1->getObs();
            for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                double th = it2->computeTheoryValue();
                std::cout << it2->getName() << " = " << th << std::endl;
            }
            std::cout << std::endl;
        }
    }
}

void MonteCarlo::Run(const int rank) {
    if (checkrun == true) {
        if (rank == 0) throw std::runtime_error("ERROR: MonteCarlo::Run() cannot be called after calling MonteCarlo::TestRun().\n");
        else sleep(2);
    } else {
        checkrun = true;
    }
    
    try {

        /* set model parameters */
        ModelName = myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO, CGP);
        int buffsize = 0;
        std::map<std::string, double> DP;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
                if (it->geterrg() > 0. || it->geterrf() > 0.)
                    buffsize++;
                if (it->IsCorrelated()) {
                    for (unsigned int i = 0; i < CGP.size(); i++) {
                        if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                            std::string index = it->getname().substr(CGP[i].getName().size());
                            long int lindex = strtol(index.c_str(), NULL, 10);
                            if (lindex > 0)
                                DP[CGP[i].getPar(lindex - 1).getname()] = CGP[i].getPar(lindex - 1).getave();
                            else {
                                std::stringstream out;
                                out << it->getname();
                                throw std::runtime_error("MonteCarlo::Run(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                            }
                        }
                    }
                } else
                    DP[it->getname()] = it->getave();
            }
        if (buffsize == 0){
            if (rank == 0) throw std::runtime_error("No parameters being varied. Aborting MCMC run.\n");
            else sleep(2);
        }
        buffsize++;
        if (!myInputParser.getModel()->Init(DP)) {
            if (myInputParser.getModel()->getmissingModelParameters().size() > 0) {
                if (rank == 0) std::cout << "\nPlease set the following parameters in the model configuration files:\n" << std::endl;
                std::vector<std::string> missingParameters = myInputParser.getModel()->getmissingModelParameters(); // Store the result in a variable to avoid dangling pointer warning
                for (std::vector<std::string>::iterator it = missingParameters.begin(); it != missingParameters.end(); it++) {
                    if (rank == 0) std::cout << "ModelParameter\t" << *it << std::endl;
                }
                std::cout << std::endl;
                myInputParser.getModel()->getmissingModelParameters().clear();
            }
        
            if (rank == 0) throw std::runtime_error("\nERROR: " + ModelName + " cannot be initialized.\n");
            else sleep(2);
        }
        int discardbuffer = 0;
        std::vector<std::string> unknownParameters = myInputParser.getModel()->getUnknownParameters();
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            if (std::find(unknownParameters.begin(), unknownParameters.end(), it->getname()) != unknownParameters.end()) {
                if (it->geterrg() > 0. || it->geterrf() > 0.) {
                    discardbuffer++;
                    if (rank == 0) std::cout << "\nWARNING: " << it->getname() << " has a flat or Gaussian error in the configuration file even though it is an unknown parameter!!!\n";
                }
            }
        }
        buffsize = buffsize - discardbuffer;
        if (buffsize == 1){
            if (rank == 0) throw std::runtime_error("No relevant parameters being varied. Aborting MCMC run.\n");
            else sleep(2);
        }

        /* create a directory for outputs */
        if (rank == 0) {
            FileStat_t info;
            if (gSystem->GetPathInfo(ObsDirName.c_str(), info) != 0) {
                if (gSystem->MakeDirectory(ObsDirName.c_str()) == 0) std::cout << ObsDirName << " directory has been created." << std::endl;
                else throw std::runtime_error("ERROR: " + ObsDirName + " directory cannot be created.\n");
            }
        }

        MCEngine.SetName(ModelName.c_str());
        MCEngine.Initialize(myInputParser.getModel());
        MCEngine.SetInitialPositionAttemptLimit(10000);   

#ifdef _MPI
        double *recvbuff = new double[buffsize];
#endif

        if (rank != 0) {
#ifdef _MPI
            double **sendbuff = new double*[1];
            sendbuff[0] = new double[1];

            std::vector<double> pars;
            int obsbuffsize = Obs.size() + 2 * Obs2D.size();
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++)
                obsbuffsize += it1->getObs().size();

            while (true) {
                MPI_Scatter(sendbuff[0], buffsize, MPI_DOUBLE, recvbuff, buffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                if (recvbuff[0] == 0.) { // do nothing and return ll
                    double ll = log(0.);
                    MPI_Gather(&ll, 1, MPI_DOUBLE, sendbuff[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                } else if (recvbuff[0] == 1.) { //compute log likelihood
                    pars.assign(recvbuff + 1, recvbuff + buffsize);
                    double ll = MCEngine.LogEval(pars);
                    MPI_Gather(&ll, 1, MPI_DOUBLE, sendbuff[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                } else if (recvbuff[0] == 2.) { // compute observables
                    double sbuff[obsbuffsize];
                    std::map<std::string, double> DPars;
                    pars.assign(recvbuff + 1, recvbuff + buffsize);
                    MCEngine.setDParsFromParameters(pars,DPars);
                    myInputParser.getModel()->Update(DPars);

                    int k = 0;
                    for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
                        sbuff[k++] = it->computeTheoryValue();
                    }
                    for (std::vector<Observable2D>::iterator it = Obs2D.begin(); it < Obs2D.end(); it++) {
                        sbuff[k++] = it->computeTheoryValue();
                        sbuff[k++] = it->computeTheoryValue2();
                    }

                    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                        std::vector<Observable> pino(it1->getObs());
                        for (std::vector<Observable>::iterator it = pino.begin(); it != pino.end(); ++it)
                            sbuff[k++] = it->computeTheoryValue();
                    }
                    MPI_Gather(sbuff, obsbuffsize, MPI_DOUBLE, sendbuff[0], obsbuffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                } else if (recvbuff[0] == 3.) { // do not compute observables, but gather the buffer
                    double sbuff[obsbuffsize];
                    MPI_Gather(sbuff, obsbuffsize, MPI_DOUBLE, sendbuff[0], obsbuffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                } else if (recvbuff[0] == -1.)
                    break;
                else {
                    std::cout << "recvbuff = " << recvbuff[0] << " rank " << rank << std::endl;
                    throw std::runtime_error("MonteCarlo::Run(): error in MPI message!\n");
                }
            }
            delete sendbuff[0];
            delete [] sendbuff;
#endif
        } else {
            ParseMCMCConfig(MCMCConf);
            if (!RunMinuitOnly && !ComputeHessianOnly) std::cout << std::endl << "\nRunning in MonteCarlo mode...\n" << std::endl;
            else if (RunMinuitOnly) std::cout << std::endl << "\nRunning Minuit Minimizer..." << std::endl;
            else std::cout << std::endl << "\nComputing the Hessian for central values of the input parameters..." << std::endl;
            std::cout << std::endl;
            if (ModPars.size() > 0) std::cout << ModPars.size() - unknownParameters.size() << " parameters defined." << std::endl;
            if (Obs.size() > 0) std::cout << Obs.size() << " observables defined." << std::endl;
            if (Obs2D.size() > 0) std::cout << Obs2D.size() << " 2D observables defined." << std::endl;
            if (CGO.size() > 0) std::cout << CGO.size() << " correlated gaussian observables defined:" << std::endl;
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                    it1 != CGO.end(); ++it1)
                std::cout << "  " << it1->getName() << " containing "
                << it1->getObs().size() << " observables." << std::endl;
            std::cout << std::endl;
                                  
            if (RunMinuitOnly) {
                BCLog::OpenLog(("MinuitMinimizationResults" + JobTag + ".txt").c_str(), BCLog::debug, BCLog::debug);
                RedirectHandle_t * rHandle = NULL;
                gSystem->RedirectOutput(("MinuitMinimizationResults" + JobTag + ".txt").c_str(), "a", rHandle);
                std::time_t ti = std::time(NULL);
#if __GNUC__ >= 5 || defined __clang__
                std::cout << "\nMinuit Minimization started at " << std::put_time(std::localtime(&ti), "%c %Z") << " (" << ti << "s since the Epoch)" << "\n" << std::endl;
#else
                char mitime[128];
                strftime(mitime,sizeof(mitime),"%c %Z", std::localtime(&ti));
                std::cout << "\nMinuit Minimization started at " << mitime << " (" << ti << "s since the Epoch)" << "\n" << std::endl;
#endif      
                std::cout << "          " <<std::endl;
                std::cout << "\t\t*** Minimizer Options ***" << std::endl;
                MCEngine.GetMinuit().Options().Print();
                std::cout << "          " <<std::endl;
                MCEngine.FindMode();
                std::time_t tf = std::time(NULL);
#if __GNUC__ >= 5 || defined __clang__
                std::cout << "\nMinuit Minimization ended at " << std::put_time(std::localtime(&tf), "%c %Z") << " (" << tf << "s since the Epoch)" << "\n" << std::endl;
#else
                char mftime[128];
                strftime(mitime,sizeof(mftime),"%c %Z", std::localtime(&tf));
                std::cout << "\nMinuit Minimization ended at " << mftime << " (" << tf << "s since the Epoch)" << "\n" << std::endl;
#endif                
                gSystem->RedirectOutput(0);
                BCLog::OutSummary(("Minuit results printed in MinuitMinimizationResults" + JobTag + ".txt").c_str());
                BCLog::SetPrefix(false);
                BCLog::OutSummary("");
                BCLog::SetPrefix(true);
                BCLog::CloseLog();               
                MCEngine.GetMinuit().PrintResults();
                std::cout << std::endl;
                return;
            }
            else if (ComputeHessianOnly) {
                BCLog::OpenLog(("HessianComputationResults" + JobTag + ".txt").c_str(), BCLog::debug, BCLog::debug);
                RedirectHandle_t * rHandle = NULL;
                gSystem->RedirectOutput(("HessianComputationResults" + JobTag + ".txt").c_str(), "a", rHandle);
                std::vector<double> point;
                unsigned int Npars = MCEngine.GetNParameters();
                for (unsigned int i = 0; i < Npars; i++) {
                    std::cout << "Parameter " << MCEngine.GetParameter(i).GetName() << " = " << MCEngine.GetParameter(i).GetPriorMean() << std::endl;
                    point.push_back(MCEngine.GetParameter(i).GetPriorMean());
                }
                gslpp::matrix<double> Hessian(Npars, Npars, 0.);
                std::vector<double> n;
                for (unsigned int i = 0; i < Npars; i++)
                    {
                        Hessian.assign(i, i, -MCEngine.SecondDerivative(MCEngine.GetParameter(i), MCEngine.GetParameter(i), point));
                        n.push_back(sqrt(fabs(Hessian(i, i))));
                        std::cout << "Sqrt(|Hessian(" << MCEngine.GetParameter(i).GetName() << "," << MCEngine.GetParameter(i).GetName() << ")|) = " << n[i] << std::endl;
                    }
                for (unsigned int i = 0; i < Npars; i++)
                        for (unsigned int j = i; j < Npars; j++) {
                        // calculate Hessian matrix element
                        Hessian.assign(i, j, -MCEngine.SecondDerivative(MCEngine.GetParameter(i), MCEngine.GetParameter(j), point));
                        if (fabs(Hessian(i, j))/n[i]/n[j] > .1)
                            std::cout << "Corr(" << MCEngine.GetParameter(i).GetName() << "," << MCEngine.GetParameter(j).GetName() << ") = " << Hessian(i, j)/n[i]/n[j] << std::endl;
                        }
                return;
            }

  
            MCEngine.CreateHistogramMaps();
            
            /* Open root file for storing data. */
            if (writechains) {
                MCEngine.WriteMarkovChain(OutFile, "RECREATE", true, false); /*Run: true, PreRun: false*/
                MCEngine.AddChains();
            } else {
                MCEngine.WriteMarkovChain(OutFile, "RECREATE", false, false);
                MCEngine.InitializeMarkovChainTree();
                MCEngine.WriteMarkovChainRun(false);
                MCEngine.WriteMarkovChainPreRun(false);
                if (MCEngine.getchainedObsSize() > 0 || MCEngine.getWriteLogLikelihoodChain() || MCEngine.getWriteParametersChain()) MCEngine.AddChains();
            }
            
            std::time_t ti = std::time(NULL);
#if __GNUC__ >= 5 || defined __clang__
            std::cout << "\nMCMC run started at " << std::put_time(std::localtime(&ti), "%c %Z") << " (" << ti << "s since the Epoch)" << "\n" << std::endl;
#else
            char mitime[128];
            strftime(mitime,sizeof(mitime),"%c %Z", std::localtime(&ti));
            std::cout << "\nMCMC run started at " << mitime << " (" << ti << "s since the Epoch)" << "\n" << std::endl;
#endif
            // open log file
            BCLog::OpenLog(("log" + JobTag + ".txt").c_str(), BCLog::debug, BCLog::debug);
            // run the MCMC and marginalize w.r.t. to all parameters
            MCEngine.BCIntegrate::SetNbins(NBINSMODELPARS);
            MCEngine.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
            if (FileIn != "")
                MCEngine.LoadMCMC(FileIn);
            MCEngine.MarginalizeAll();
            std::time_t tf = std::time(NULL);
#if __GNUC__ >= 5 || defined __clang__
            std::cout << "\nMCMC run ended at  " << std::put_time(std::localtime(&tf), "%c %Z") << " (" << tf << "s since the Epoch)" << std::endl;
#else
            char mftime[128];
            strftime(mftime,sizeof(mftime),"%c %Z", std::localtime(&tf));
            std::cout << "\nMCMC run ended at " << mftime << " (" << tf << "s since the Epoch)" << "\n" << std::endl;
#endif            
            int hour = (tf-ti)/3600;
            int min = ((tf-ti)%3600) / 60;
            int sec = ((tf-ti)%3600) % 60;
            std::cout << "Time taken by the MCMC run: ";
            std::cout << std::setfill('0') << std::setw(2) << hour << ":"; 
            std::cout << std::setfill('0') << std::setw(2) << min  << ":" ;
            std::cout << std::setfill('0') << std::setw(2) << sec   << "\n" << std::endl;

            // find mode using the best fit parameters as start values
            if (FindModeWithMinuit) MCEngine.FindMode(MCEngine.GetBestFitParameters());

            // draw all marginalized distributions into a pdf file
            if (PrintAllMarginalized) MCEngine.PrintAllMarginalized(("MonteCarlo_plots" + JobTag + ".pdf").c_str());
           
            // print histograms
            MCEngine.PrintHistogram(OutFile, ObsDirName);

            // draw the correlation matrix into a pdf file and make the tex files
            if (PrintCorrelationMatrix) {
                MCEngine.PrintCorrelationMatrix(("CorrelationMatrix" + JobTag + ".pdf").c_str());
                MCEngine.PrintCorrelationMatrix(("CorrelationMatrix" + JobTag + ".tex").c_str());
                MCEngine.PrintCorrelationMatrixToLaTeX(("ParameterCorrelations" + JobTag + ".tex").c_str());
            }
            
            // print comparisons of the prior knowledge to the posterior knowledge 
            // for all parameters into a pdf file
            if (PrintKnowledgeUpdatePlots) MCEngine.PrintKnowledgeUpdatePlots(("ParamUpdate" + JobTag + ".pdf").c_str());

            // draw an overview plot of the parameters into a pdf file
            if (PrintParameterPlot) MCEngine.PrintParameterPlot(("ParamSummary" + JobTag + ".pdf").c_str());

            MCEngine.WriteMarginalizedDistributions(OutFile, "UPDATE");

            // print logs for the histograms of the observables into a text file
            std::ofstream outHistoLog;
            outHistoLog.open((ObsDirName + "/HistoLog" + JobTag + ".txt").c_str(), std::ios::out);
            outHistoLog << MCEngine.getHistoLog();
            outHistoLog.close();

            // print statistics for the theory values of the observables into a text file
            std::ofstream outStatLog;
            outStatLog.open((ObsDirName + "/Statistics" + JobTag + ".txt").c_str(), std::ios::out);
            outStatLog << MCEngine.computeStatistics();
            outStatLog.close();
            
            // print global mode and scale factors for the 1st chain into a text file
            if (WritePreRunData) {
                std::ofstream outPreRun;
                outPreRun.open(("PreRun" + JobTag + ".txt").c_str(), std::ios::out);
                outPreRun << MCEngine.writePreRunData();
                outPreRun.close();
            }

            /* Number of events */
            std::stringstream ss;
            ss << ""; 
            BCLog::SetPrefix(false);
            BCLog::OutSummary(ss.str().c_str());
            BCLog::SetPrefix(true);
            ss << "Number of used events: " << MCEngine.getNumOfUsedEvents();
            BCLog::OutSummary(ss.str().c_str());
            ss.str("");
            ss << "Number of discarded events: " << MCEngine.getNumOfDiscardedEvents() << "\n";
            BCLog::OutSummary(ss.str().c_str());
            // close log file
            BCLog::CloseLog();
            
            // print results of the analysis into a text file
            BCLog::SetPrefix(false);
            BCLog::OpenLog(("MonteCarlo_results" + JobTag + ".txt").c_str(), BCLog::summary, BCLog::nothing);
            MCEngine.PrintSummary();
            BCLog::CloseLog();
            
            // Calculate and print normalization        
            if (CalculateNormalization.compare("false") != 0) {
                std::ofstream outStatLog;
                outStatLog.open((ObsDirName + "/Statistics" + JobTag + ".txt").c_str(), std::ios::app);
                
                if (CalculateNormalization.compare("LME") == 0) {
                    // Make sure mode is found by MINUIT for normalization computation.
                    if (!FindModeWithMinuit) MCEngine.FindMode(MCEngine.GetBestFitParameters());
                    outStatLog << "\nNormalization for " << ModelName.c_str() << ": " << MCEngine.computeNormalizationLME() << "\n" << std::endl;
                } 
                else if (CalculateNormalization.compare("MC") == 0) {
                    std::vector<double> normalizationMC = MCEngine.computeNormalizationMC(NIterationNormalizationMC);
                    outStatLog << "\nNormalization for " << ModelName.c_str() << ": " << normalizationMC[0] << " +- " << normalizationMC[1] << "\n" << std::endl;
                }
                else 
                    throw std::runtime_error(("\n ERROR: Normalization method" + CalculateNormalization + " not implemented.\n").c_str());
                
                outStatLog.close();
            }
            
            // Print the triangle plot (the size of the file is quite large for large number of parameters).
            // this is done last for now since BAT sometimes crashes while doing this. Did not look for the fix.
            if (PrintTrianglePlot) {
                BCLog::OpenLog(("log_CorrelationPlots" + JobTag + ".txt").c_str(), BCLog::debug, BCLog::debug);
                MCEngine.PrintCorrelationPlot("CorrelationPlots" + JobTag + ".pdf");
                BCLog::CloseLog();
            }
            
            
#ifdef _MPI
            double ** sendbuff = new double *[MCEngine.procnum];
            sendbuff[0] = new double[MCEngine.procnum * buffsize];
            for (int il = 1; il < MCEngine.procnum; il++) {
                sendbuff[il] = sendbuff[il - 1] + buffsize;
                sendbuff[il][0] = -1.; //Exit command
            }
            MPI_Scatter(sendbuff[0], buffsize, MPI_DOUBLE, recvbuff, buffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            delete sendbuff[0];
            delete [] sendbuff;
#endif
        }
#ifdef _MPI
        delete [] recvbuff;
#endif
    } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
    }
}

void MonteCarlo::ParseMCMCConfig(std::string file)
{
    std::ifstream ifile(file.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nERROR: " + MCMCConf + " does not exist. Make sure to specify a valid Monte Carlo configuration file.\n");
    std::string line;
    bool IsEOF = false;
    do {
        IsEOF = getline(ifile, line).eof();
        if(IsEOF) continue;
        if (!line.empty() && *line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
        if (line.empty() || line.at(0) == '#')
            continue;
        boost::char_separator<char> sep(" \t");
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
        if (beg->compare("NChains") == 0) {
            ++beg;
            if (isdigit(beg->at(0)) && atoi((*beg).c_str()) > 0) MCEngine.setNChains(atoi((*beg).c_str()));
            else if (isdigit(beg->at(0)) && atoi((*beg).c_str()) == 0) 
                if (MCEngine.getMPIWorldSize() > 1) MCEngine.setNChains(MCEngine.getMPIWorldSize());
                else throw std::runtime_error("\nERROR: NChains in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer > 0.\n");
            else
                throw std::runtime_error("\nERROR: NChains in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer >= 0.\n");
        } else if (beg->compare("PrerunMaxIter") == 0) {
            ++beg;
            if (isdigit(beg->at(0)) && atoi((*beg).c_str()) > 0) MCEngine.SetNIterationsPreRunMax(atoi((*beg).c_str()));
            else
                throw std::runtime_error("\nERROR: PrerunMaxIter in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer > 0.\n");
        } else if (beg->compare("NIterationsUpdateMax") == 0) {
            ++beg;
            if (isdigit(beg->at(0)) && atoi((*beg).c_str()) > 0) MCEngine.SetNIterationsPreRunCheck(atoi((*beg).c_str()));
            else
                throw std::runtime_error("\nERROR: NIterationsUpdateMax in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer > 0.\n");
        } else if (beg->compare("NIterationsPreRunFactorized") == 0) {
            ++beg;
            if (isdigit(beg->at(0)) && atoi((*beg).c_str()) >= 0) MCEngine.SetNIterationsPreRunFactorized(atoi((*beg).c_str()));
            else
                throw std::runtime_error("\nERROR: NIterationsPrerunFactorized in the MonteCarlo configuration file: " + MCMCConf + " can only be a non-negative integer.\n");
        } else if (beg->compare("Seed") == 0) {
            ++beg;
            if (!isdigit(beg->at(0))) throw std::runtime_error("\nERROR: Seed in the MonteCarlo configuration file: " + MCMCConf + " can only be a number.\n");
            int seed = atoi((*beg).c_str());
            if (seed != 0)
                MCEngine.SetRandomSeed(seed);
        } else if (beg->compare("Iterations") == 0) {
            ++beg;
            if (isdigit(beg->at(0))) MCEngine.SetNIterationsRun(atoi((*beg).c_str()));
            else
                throw std::runtime_error("\nERROR: Iterations in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer.\n");
        } else if (beg->compare("MinimumEfficiency") == 0) {
            ++beg;
            double efficiency = atof((*beg).c_str());
            if (efficiency > 0. && efficiency <= 1.) MCEngine.SetMinimumEfficiency(efficiency);
            else
                throw std::runtime_error("\nERROR: MinimumEfficiency in the MonteCarlo configuration file: " + MCMCConf + " can only be a real number greater than 0.0 and less than or equal to 1.0.\n");
        } else if (beg->compare("MaximumEfficiency") == 0) {
            ++beg;
            double efficiency = atof((*beg).c_str());
            if (efficiency > 0. && efficiency <= 1.) MCEngine.SetMaximumEfficiency(efficiency);
            else
                throw std::runtime_error("\nERROR: MaximumEfficiency in the MonteCarlo configuration file: " + MCMCConf + " can only be a real number greater than 0.0 and less than or equal to 1.0.\n");
        } else if (beg->compare("RValueForConvergence") == 0) {
            ++beg;
            double rvalue = atof((*beg).c_str());
            if (rvalue > 1. ) MCEngine.SetRValueParametersCriterion(rvalue);
            else
                throw std::runtime_error("\nERROR: RValueForConvergence in the MonteCarlo configuration file: " + MCMCConf + " can only be a real number greater than 1.0.\n");
        } else if (beg->compare("WriteChain") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) writechains = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: WriteChain in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("FindModeWithMinuit") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) FindModeWithMinuit = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: FindModeWithMinuit in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("RunMinuitOnly") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) RunMinuitOnly = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: RunMinuitOnly in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("ComputeHessianOnly") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) ComputeHessianOnly = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: ComputeHessianOnly in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("CalculateNormalization") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) CalculateNormalization = *beg;
            else
                throw std::runtime_error("\nERROR: CalculateNormalization in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("NIterationNormalizationMC") == 0) {
            ++beg;
            if (isdigit(beg->at(0))) NIterationNormalizationMC = atoi((*beg).c_str());
            else
                throw std::runtime_error("\nERROR: NIterationNormalizationMC in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer.\n");
        } else if (beg->compare("PrintAllMarginalized") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) PrintAllMarginalized = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: PrintAllMarginalized in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("PrintCorrelationMatrix") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) PrintCorrelationMatrix = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: PrintCorrelationMatrix in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("PrintKnowledgeUpdatePlots") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) PrintKnowledgeUpdatePlots = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: PrintKnowledgeUpdatePlots in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("PrintParameterPlot") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) PrintParameterPlot = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: PrintParameterPlot in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("PrintTrianglePlot") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) PrintTrianglePlot = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: PrintTrianglePlot in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("WritePreRunData") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) WritePreRunData = (beg->compare("true") == 0);
            else
                throw std::runtime_error("\nERROR: WritePreRunData in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("ReadPreRunData") == 0) {
            ++beg;
            ReadPreRunData(*beg);
        } else if (beg->compare("MultivariateProposal") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.SetProposeMultivariate((beg->compare("true") == 0));
            else
                throw std::runtime_error("\nERROR: MultivariateProposal in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("Histogram1DSmooth") == 0) {
            ++beg;
            if (beg->compare("true") == 0) {
                MCEngine.setSmooth(1);
            } else if (beg->compare("false") == 0) {
                MCEngine.setSmooth(0); /* Default */
            } else if (isdigit(beg->at(0))) {
                if (atoi((*beg).c_str()) >= 0 && atoi((*beg).c_str()) <= 5) MCEngine.setSmooth(atoi((*beg).c_str()));
            } else
                throw std::runtime_error("\nERROR: Histogram1DSmooth in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true', 'false' or an integer from 0 to 5.\n");
        } else if (beg->compare("Histogram2DType") == 0) {
            ++beg;
            if (!isdigit(beg->at(0)))
                throw std::runtime_error("\nERROR: Histogram2DType in the MonteCarlo configuration file : " + MCMCConf + "can only be an integer amongst 1001 -> Lego, 101 -> Filled, 1 -> Contour.\n");
            int type = atoi((*beg).c_str());
            if (type == 1 || type == 101 || type == 1001) {
                MCEngine.setHistogram2DType(type);
            } else
                throw std::runtime_error("\nERROR: Histogram2DType in the MonteCarlo configuration file : " + MCMCConf + "can only be an integer amongst 1001 -> Lego, 101 -> Filled, 1 -> Contour.\n");
        } else if (beg->compare("MCMCInitialPosition") == 0) {
            ++beg;
            if (beg->compare("Center") == 0) {
                MCEngine.SetInitialPositionScheme(BCEngineMCMC::kInitCenter); /* Default */
            } else if (beg->compare("RandomUniform") == 0) {
                MCEngine.SetInitialPositionScheme(BCEngineMCMC::kInitRandomUniform);
            } else if (beg->compare("RandomPrior") == 0) {
                MCEngine.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
            } else
                throw std::runtime_error("\nERROR: MCMCInitialPosition in MonteCarlo config file: " + MCMCConf + " can only be 'Center', 'RandomUniform' or 'RandomPrior'.\n");
        } else if (beg->compare("PrintLogo") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.setPrintLogo((beg->compare("true") == 0));
            else
                throw std::runtime_error("\nERROR: PrintLogo in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("NoHistogramLegend") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.setNoLegend((beg->compare("true") == 0));
            else
                throw std::runtime_error("\nERROR: PrintLogo in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("PrintLoglikelihoodPlots") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.setPrintLoglikelihoodPlots((beg->compare("true") == 0));
            else
                throw std::runtime_error("\nERROR: PrintLoglikelihoodPlots in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
            if (beg->compare("true") == 0) gSystem->MakeDirectory((ObsDirName + "/LogLikelihoodPlots").c_str());
        } else if (beg->compare("WriteLogLikelihoodChain") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.setWriteLogLikelihoodChain((beg->compare("true") == 0));
            else
                throw std::runtime_error("\nERROR: WriteLogLikelihoodChain in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("WriteParametersChains") == 0) {
            ++beg;
            if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.setWriteParametersChain((beg->compare("true") == 0));
            else
                throw std::runtime_error("\nERROR: WriteParametersChain in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
        } else if (beg->compare("Histogram2DAlpha") == 0) {
            ++beg;
            double alpha = atof((*beg).c_str());
            if (alpha > 0. && alpha <= 1.) MCEngine.setAlpha2D(alpha);
            else
                throw std::runtime_error("\nERROR: Histogram2DAlpha in the MonteCarlo configuration file: " + MCMCConf + " can only be a real number greater than 0.0 and less than or equal to 1.0.\n");
        } else if (beg->compare("NBinsHistogram1D") == 0) {
            ++beg;
            unsigned int nBins1D = atoi((*beg).c_str());
            if (nBins1D >= 0) MCEngine.setNBins1D(nBins1D);
            else throw std::runtime_error("\nERROR: NBinsHistogram1D in the MonteCarlo configuration file: " + MCMCConf + " can only be a integer greater than 0 or 0 to set to default value (100).\n");
        } else if (beg->compare("NBinsHistogram2D") == 0) {
            ++beg;
            unsigned int nBins2D = atoi((*beg).c_str());
            if (nBins2D >= 0) MCEngine.setNBins2D(nBins2D);
            else throw std::runtime_error("\nERROR: NBinsHistogram2D in the MonteCarlo configuration file: " + MCMCConf + " can only be a integer greater than 0 or 0 to set to default value (100).\n");
        } else if (beg->compare("InitialPositionAttemptLimit") == 0) {
            ++beg;
            unsigned int max_tries = atoi((*beg).c_str());
            if (max_tries > 0) MCEngine.SetInitialPositionAttemptLimit(max_tries);
            else if (max_tries == 0) MCEngine.SetInitialPositionAttemptLimit(10000);
            else throw std::runtime_error("\nERROR: InitialPositionAttemptLimit in the MonteCarlo configuration file: " + MCMCConf + " can only be a integer greater than 0 or 0 to set to default value (100).\n");
        } else if (beg->compare("SignificantDigits") == 0) {
            ++beg;
            unsigned int significants = atoi((*beg).c_str());
            if (significants >= 0) MCEngine.setSignificants(significants);
            else throw std::runtime_error("\nERROR: SignificantDigits in the MonteCarlo configuration file: " + MCMCConf + " can only be a integer greater than 0 or 0 to set to default values.\n");
        } else if (beg->compare("HistogramBufferSize") == 0) {
            ++beg;
            unsigned int histogramBufferSize = atoi((*beg).c_str());
            if (histogramBufferSize >= 0) MCEngine.setHistogramBufferSize(histogramBufferSize);
            else throw std::runtime_error("\nERROR: HistogramBufferSize in the MonteCarlo configuration file: " + MCMCConf + " can only be a integer greater than 0 or 0 to set to default values.\n");
        } else if (beg->compare("ScaleFactorLowerLimit") == 0) {
            ++beg;
            double ScaleFactorLowerLimit = atof((*beg).c_str());
            if (ScaleFactorLowerLimit >= 0) MCEngine.SetScaleFactorLowerLimit(ScaleFactorLowerLimit);
            else throw std::runtime_error("\nERROR: ScaleFactorLowerLimit in the MonteCarlo configuration file: " + MCMCConf + " can only be a double greater or equal to 0.\n");
        } else
            throw std::runtime_error("\nERROR: Wrong keyword in MonteCarlo config file: " + MCMCConf + "\n Make sure to specify a valid Monte Carlo configuration file.\n");
    } while (!IsEOF);
    if (FindModeWithMinuit || RunMinuitOnly) 
        MCEngine.SetOptimizationMethod(BCIntegrate::kOptMinuit);
    if (MCEngine.GetMaximumEfficiency() <= MCEngine.GetMinimumEfficiency()) 
        throw std::runtime_error("\nERROR: MaximumEfficiency (default 0.5) must be greater than MaximumEfficiency (default 0.15) in the MonteCarlo configuration file: " + MCMCConf + ".\n");
     if (CalculateNormalization.compare("MC") == 0 && NIterationNormalizationMC <= 0) 
        throw std::runtime_error(("\nMonteCarlo ERROR: CalculateNormalization cannot be set to MC without setting NIterationNormalizationMC > 0 in " + MCMCConf + " .\n").c_str());
}

void MonteCarlo::ReadPreRunData(std::string file)
{
    std::ifstream ifile(file.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nMonteCarlo::ReadPreRunData ERROR: " + file + " does not exist.\n");
    std::string line;
    bool IsEOF = false;
    std::vector<double> mode;
    std::vector<double> scale;
    do {
        IsEOF = getline(ifile, line).eof();
        if (line.empty())
            continue;   
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
        ++beg;
        mode.push_back(atof((*beg).c_str()));
        ++beg;
        scale.push_back(atof((*beg).c_str()));
    } while (!IsEOF);
    if (mode.size() != MCEngine.GetNParameters())
        throw std::runtime_error("\nMonteCarlo::ReadPreRunData ERROR: wrong data size.\n");
    std::vector<double> mode_all;
    std::vector<double> scale_all;
    for (unsigned int i = 0; i < MCEngine.GetNChains(); i++){
        mode_all.insert(mode_all.end(), mode.begin(), mode.end());
        scale_all.insert(scale_all.end(), scale.begin(), scale.end());
    }
    MCEngine.SetInitialPositions(mode_all);
    MCEngine.SetInitialScaleFactors(scale_all);
}

void MonteCarlo::addCustomObservableType(const std::string name, boost::function<Observable*() > funct) {
    myInputParser.addCustomObservableType(name, funct);
}
