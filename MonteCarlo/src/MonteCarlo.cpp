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


MonteCarlo::MonteCarlo(
        ModelFactory& ModelF, ThObsFactory& ThObsF,
        const std::string& ModelConf_i,
        const std::string& MonteCarloConf_i,
        const std::string& OutFile_i,
        const std::string& JobTag_i)
: myInputParser(ModelF, ThObsF), MCEngine(ModPars, Obs, Obs2D, CGO, CGP) {
    ModelConf = ModelConf_i;
    MCMCConf = MonteCarloConf_i;
    JobTag = JobTag_i;
    if (OutFile_i.compare("") == 0) OutFile = "MCout" + JobTag + ".root";
    else OutFile = OutFile_i + JobTag + ".root";
    ObsDirName = "Observables" + JobTag;
    FindModeWithMinuit = false;
    CalculateNormalization = "false";
    NIterationNormalizationMC = 0;
    PrintAllMarginalized = false;
    PrintCorrelationMatrix = false;
    PrintKnowledgeUpdatePlots = false;
    PrintParameterPlot = false;
    WritePreRunData = false;
    checkrun = false;
}

MonteCarlo::~MonteCarlo() {}

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
            if (rank == 0) throw std::runtime_error("ERROR: Parameter(s) missing in model initialization. \n");
            else sleep (2);
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
        MCEngine.SetName(ModelName);
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
        if (!myInputParser.getModel()->Init(DP)){
            if (rank == 0) throw std::runtime_error("\nERROR: Model cannot be initialization.\n");
            else sleep(2);
        }

        if (rank == 0) std::cout << std::endl << "Running in MonteCarlo mode...\n" << std::endl;

        /* create a directory for outputs */
        if (rank == 0) {
            FileStat_t info;
            if (gSystem->GetPathInfo(ObsDirName.c_str(), info) != 0) {
                if (gSystem->MakeDirectory(ObsDirName.c_str()) == 0)
                    std::cout << ObsDirName << " directory has been created." << std::endl;
                else
                    throw std::runtime_error("ERROR: " + ObsDirName + " directory cannot be created.\n");
            }
        }

        MCEngine.SetName(ModelName.c_str());
        MCEngine.Initialize(myInputParser.getModel());

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
                MPI_Scatter(sendbuff[0], buffsize, MPI_DOUBLE,
                        recvbuff, buffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
                    throw "MonteCarlo::Run(): error in MPI message!\n";
                }
            }
            delete sendbuff[0];
            delete [] sendbuff;
#endif
        } else {

            bool writechains = false;
            std::cout << std::endl;
            std::cout << std::endl;
            if (ModPars.size() > 0) std::cout << ModPars.size() << " parameters defined." << std::endl;
            if (Obs.size() > 0) std::cout << Obs.size() << " observables defined." << std::endl;
            if (Obs2D.size() > 0) std::cout << Obs2D.size() << " 2D observables defined." << std::endl;
            if (CGO.size() > 0) std::cout << CGO.size() << " correlated gaussian observables defined:" << std::endl;
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                    it1 != CGO.end(); ++it1)
                std::cout << "  " << it1->getName() << " containing "
                << it1->getObs().size() << " observables." << std::endl;
            //MonteCarlo configuration parser
            std::ifstream ifile(MCMCConf.c_str());
            if (!ifile.is_open())
                throw std::runtime_error("\nERROR: " + MCMCConf + " does not exist. Make sure to specify a valid Monte Carlo configuration file.\n");
            std::string line;
            bool IsEOF = false;
            do {
                IsEOF = getline(ifile, line).eof();
                if (*line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
                if (line.empty() || line.at(0) == '#')
                    continue;
                boost::char_separator<char> sep(" \t");
                boost::tokenizer<boost::char_separator<char> > tok(line, sep);
                boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
                if (beg->compare("NChains") == 0) {
                    ++beg;
                    if (isdigit(beg->at(0)) && atoi((*beg).c_str()) > 0) MCEngine.setNChains(atoi((*beg).c_str()));
                    else
                        throw std::runtime_error("\nERROR: NChains in the MonteCarlo configuration file: " + MCMCConf + " can only be an integer > 0.\n");
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
                    if(efficiency > 0. && efficiency <= 1.) MCEngine.SetMinimumEfficiency(efficiency);
                    else
                        throw std::runtime_error("\nERROR: MinimumEfficiency in the MonteCarlo configuration file: " + MCMCConf + " can only be an real number greater than 0.0 and less than or equal to 1.0.\n");
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
                    if (beg->compare("true") == 0 || beg->compare("false") == 0) PrintParameterPlot = beg->compare("true") == 0;
                    else
                        throw std::runtime_error("\nERROR: PrintParameterPlot in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
                } else if (beg->compare("WritePreRunData") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0 || beg->compare("false") == 0) WritePreRunData = beg->compare("true") == 0;
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
                    } else if (beg->compare("RandomUniform") == 0){
                        MCEngine.SetInitialPositionScheme(BCEngineMCMC::kInitRandomUniform);
                    } else if (beg->compare("RandomPrior") == 0){
                        MCEngine.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
                    } else 
                        throw std::runtime_error("\nERROR: MCMCInitialPosition in MonteCarlo config file: " + MCMCConf + " can only be 'Center', 'RandomUniform' or 'RandomPrior'.\n");
                } else if (beg->compare("PrintLogo") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0 || beg->compare("false") == 0) MCEngine.setPrintLogo(true);
                    else 
                        throw std::runtime_error("\nERROR: PrintLogo in the MonteCarlo configuration file: " + MCMCConf + " can only be 'true' or 'false'.\n");
                } else
                    throw std::runtime_error("\nERROR: Wrong keyword in MonteCarlo config file: " + MCMCConf + "\n Make sure to specify a valid Monte Carlo configuration file.\n");
            } while (!IsEOF);

            if (CalculateNormalization.compare("MC") == 0 && NIterationNormalizationMC <= 0) 
                throw std::runtime_error(("\nMonteCarlo ERROR: CalculateNormalization cannot be set to MC without setting NIterationNormalizationMC > 0 in " + MCMCConf + " .\n").c_str());
            
            /* Open root file for storing data. */
            if (writechains) {
                MCEngine.WriteMarkovChain(OutFile, "RECREATE", true, false); /*Run: true, PreRun: false*/
                MCEngine.AddChains();
            } else {
                MCEngine.WriteMarkovChain(OutFile, "RECREATE", false, false);
                MCEngine.InitializeMarkovChainTree();
                MCEngine.WriteMarkovChainRun(false);
                MCEngine.WriteMarkovChainPreRun(false);
            }
            // set nicer style for drawing than the ROOT default
            BCAux::SetStyle();

            // open log file
            BCLog::OpenLog(("log" + JobTag + ".txt").c_str(), BCLog::debug, BCLog::debug);
//            BCLog::SetLogLevel(BCLog::debug);

            // run the MCMC and marginalize w.r.t. to all parameters
            MCEngine.BCIntegrate::SetNbins(NBINSMODELPARS);
            MCEngine.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
            std::time_t ti = std::time(NULL);
            char mbstr[100];
            if (std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S %d/%m/%y ", std::localtime(&ti))) 
                std::cout << "MCMC Run started at " <<  mbstr << std::endl; 
            MCEngine.MarginalizeAll();
            std::time_t tf = std::time(NULL);
            if (std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S %d/%m/%y ", std::localtime(&tf))) 
                std::cout << "MCMC Run ended at " <<  mbstr << std::endl; 

            // find mode using the best fit parameters as start values
            if (FindModeWithMinuit)
                MCEngine.FindMode(MCEngine.GetBestFitParameters());

            // draw all marginalized distributions into a pdf file
            if (PrintAllMarginalized)
                MCEngine.PrintAllMarginalized(("MonteCarlo_plots" + JobTag + ".pdf").c_str());
           
            // print histograms
            MCEngine.PrintHistogram(OutFile, ObsDirName);

            // draw the correlation matrix into a pdf file and make the tex files
            if (PrintCorrelationMatrix) {
                MCEngine.PrintCorrelationMatrix(("CorrelationMatrix" + JobTag + ".pdf").c_str());
                MCEngine.PrintCorrelationMatrix(("CorrelationMatrix" + JobTag + ".tex").c_str());
                MCEngine.PrintCorrelationPlot("CorrelationPlots" + JobTag + ".pdf");
                MCEngine.PrintCorrelationMatrixToLaTeX(("ParameterCorrelations" + JobTag + ".tex").c_str());
            }
               
            // print comparisons of the prior knowledge to the posterior knowledge 
            // for all parameters into a pdf file
            if (PrintKnowledgeUpdatePlots)
                MCEngine.PrintKnowledgeUpdatePlots(("ParamUpdate" + JobTag + ".pdf").c_str());

            // draw an overview plot of the parameters into a pdf file
            if (PrintParameterPlot)
                MCEngine.PrintParameterPlot(("ParamSummary" + JobTag + ".pdf").c_str());
            
            // print a LaTeX table of the parameters into a tex file
            //MCEngine.PrintParameterLatex(("ParamSummary" + JobTag + ".tex").c_str());

            MCEngine.WriteMarginalizedDistributions(OutFile, "UPDATE");
            //out.FillAnalysisTree();
            //out.Close();

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
            ss << "Number of used events: " << MCEngine.getNumOfUsedEvents();
            BCLog::OutSummary(ss.str().c_str());
            ss.str("");
            ss << "Number of discarded events: " << MCEngine.getNumOfDiscardedEvents();
            BCLog::OutSummary(ss.str().c_str());
            // close log file
            BCLog::CloseLog();
            
            // print results of the analysis into a text file
            BCLog::OpenLog(("MonteCarlo_results" + JobTag + ".txt").c_str(), BCLog::results, BCLog::nothing);
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
            
            
#ifdef _MPI
            double ** sendbuff = new double *[MCEngine.procnum];
            sendbuff[0] = new double[MCEngine.procnum * buffsize];
            for (int il = 1; il < MCEngine.procnum; il++) {
                sendbuff[il] = sendbuff[il - 1] + buffsize;
                sendbuff[il][0] = -1.; //Exit command
            }
            MPI_Scatter(sendbuff[0], buffsize, MPI_DOUBLE,
                    recvbuff, buffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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