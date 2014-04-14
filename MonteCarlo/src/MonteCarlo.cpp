/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MonteCarlo.h"
#include <BAT/BCParameter.h>
#include <TSystem.h>
#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCSummaryTool.h>
#ifdef _MPI
#include <mpi.h>
#endif
#include <fstream>
#include <sstream>
#include <stdexcept>

MonteCarlo::MonteCarlo(const std::string& ModelConf_i,
                       const std::string& MonteCarloConf_i,
                       const std::string& OutFile_i,
                       const std::string& JobTag_i,
                       const bool checkTheoryRange_i)
: myInputParser(), MCEngine(ModPars, Obs, Obs2D, CGO, ParaObs, checkTheoryRange_i)
{
    ModelConf = ModelConf_i;
    MCMCConf = MonteCarloConf_i;
    JobTag = JobTag_i;
    if (OutFile_i.compare("") == 0) OutFile = "MCout" + JobTag + ".root";
    else OutFile = OutFile_i + JobTag + ".root";
    ObsDirName = "Observables" + JobTag;
    FindModeWithMinuit = false;
    PrintAllMarginalized = false;
    PrintCorrelationMatrix = false;
    PrintKnowledgeUpdatePlots = false;
    PrintParameterPlot = false;
    checkrun = false;
}

//MonteCarlo::~MonteCarlo() {}

void MonteCarlo::TestRun(int rank) {
    if (checkrun == true){
        if (rank == 0) throw std::runtime_error("ERROR: MonteCarlo::TestRun() cannot be called after calling MonteCarlo::Run().\n");
    } else {
        checkrun = true;
    }
    
    if (rank == 0){
        std::string ModelName = myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO, ParaObs);
        std::map<std::string, double> DP;
        if (Obs.size() == 0 && CGO.size() == 0) throw std::runtime_error("\nMonteCarlo::TestRun(): No observables or correlated Gaussian observables defined in SomeModel.conf file\n");

        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            DP[it->name] = it->ave;
        }

        if (!myInputParser.getMyModel()->Init(DP)) {
            throw std::runtime_error("ERROR: Parameter(s) missing in model initialization. \n");
        }

        if (Obs.size() > 0) std::cout << "\nOservables: \n" << std::endl;
        for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
            double th = it->computeTheoryValue();
            std::cout << it->getName() << " = " << th << std::endl;
        }

        if (CGO.size() > 0) std::cout << "\nCorrelated Gaussian Oservables: \n" << std::endl;
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

void MonteCarlo::Run(const int rank)
{
    if (checkrun == true){
        if (rank == 0) throw std::runtime_error("ERROR: MonteCarlo::Run() cannot be called after calling MonteCarlo::TestRun().\n");
    } else {
        checkrun = true;
    }
    
    try {

        /* set model parameters */
        std::string ModelName = myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO, ParaObs);
        int buffsize = 0;
        std::map<std::string, double> DP;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            DP[it->name] = it->ave;
            if (it->errg > 0. || it->errf > 0.)
                buffsize++;
        }
        buffsize++;
        if (!myInputParser.getMyModel()->Init(DP))
            throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");

        if (rank == 0) std::cout << std::endl << "Running in MonteCarlo mode...\n" << std::endl;

        /* create a directory for outputs */
        if (rank == 0){
        FileStat_t info;
            if (gSystem->GetPathInfo(ObsDirName.c_str(), info) != 0) {
                if (gSystem->MakeDirectory(ObsDirName.c_str()) == 0)
                    std::cout << ObsDirName << " directory has been created." << std::endl;
                else
                    throw std::runtime_error("ERROR: " + ObsDirName + " director cannot be created.\n");
            }
        }

        MCEngine.SetName(ModelName.c_str());
        MCEngine.Initialize(myInputParser.getMyModel());

#ifdef _MPI
        double *recvbuff = new double[buffsize];
#endif

        if (rank != 0) {
#ifdef _MPI
            double **sendbuff = new double*[1];
            sendbuff[0] = new double[1];

            std::vector<double> pars;
            int obsbuffsize = Obs.size() + 2 * Obs2D.size() + 2 * ParaObs.size();
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++)
                obsbuffsize += it1->getObs().size();
            
            while (true) {
                MPI::COMM_WORLD.Scatter(sendbuff[0], buffsize, MPI::DOUBLE,
                                        recvbuff, buffsize, MPI::DOUBLE, 0);
                
                if (recvbuff[0] == 0.){ // do nothing and return ll
                    double ll = log(0.);
                    MPI::COMM_WORLD.Gather(&ll, 1, MPI::DOUBLE, sendbuff[0], 1, MPI::DOUBLE, 0);
                }
                else if (recvbuff[0] == 1.) { //compute log likelihood
                    pars.assign(recvbuff + 1, recvbuff + buffsize);
                    double ll = MCEngine.LogEval(pars);
                    MPI::COMM_WORLD.Gather(&ll, 1, MPI::DOUBLE, sendbuff[0], 1, MPI::DOUBLE, 0);
                } 
                else if (recvbuff[0] == 2.) { // compute observables
                    double sbuff[obsbuffsize];
                    std::map<std::string, double> DPars;
                    pars.assign(recvbuff + 1, recvbuff + buffsize);
                    for (unsigned int k = 0; k < pars.size(); k++) {
                        DPars[MCEngine.GetParameter(k)->GetName()] = pars[k];
                    }
                    myInputParser.getMyModel()->Update(DPars);

                    int k = 0;
                    for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++){
                        sbuff[k++] = it->computeTheoryValue();
                    }
                    for (std::vector<Observable2D>::iterator it = Obs2D.begin(); it < Obs2D.end(); it++) {
                        sbuff[k++] = it->computeTheoryValue();
                        sbuff[k++] = it->computeTheoryValue2();
                    }

                    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++){
                        std::vector<Observable> pino(it1->getObs());
                        for (std::vector<Observable>::iterator it = pino.begin(); it != pino.end(); ++it)
                            sbuff[k++] = it->computeTheoryValue();
                    }
                    for (std::vector<ModelParaVsObs>::iterator it = ParaObs.begin(); it < ParaObs.end(); it++) {
                        sbuff[k++] = DPars[it->getParaName()];
                        sbuff[k++] = it->computeTheoryValue();
                    }
                    MPI::COMM_WORLD.Gather(sbuff, obsbuffsize, MPI::DOUBLE, sendbuff[0], obsbuffsize, MPI::DOUBLE, 0);
                }
                else if (recvbuff[0] == 3.) { // do not compute observables, but gather the buffer
                    double sbuff[obsbuffsize];
                    MPI::COMM_WORLD.Gather(sbuff, obsbuffsize, MPI::DOUBLE, sendbuff[0], obsbuffsize, MPI::DOUBLE, 0);
                }
                else if (recvbuff[0] == -1.)
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
            std::cout << ModPars.size() << " parameters defined." << std::endl;
            std::cout << Obs.size() << " observables defined." << std::endl;
            std::cout << CGO.size() << " correlated gaussian observables defined";
            if (CGO.size() == 0)
                std::cout << "." << std::endl;
            else
                std::cout << ":" << std::endl;
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                    it1 != CGO.end(); ++it1)
                std::cout << "  " << it1->getName() << " containing "
                          << it1->getObs().size() << " observables." << std::endl;
            std::cout << ParaObs.size() << " ModelParaVsObs defined." << std::endl;
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
                boost::char_separator<char> sep(" ");
                boost::tokenizer<boost::char_separator<char> > tok(line, sep);
                boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
                if (beg->compare("NChains") == 0) {
                    ++beg;
                    MCEngine.setNChains(atoi((*beg).c_str()));
                } else if (beg->compare("PrerunMaxIter") == 0) {
                    ++beg;
                    MCEngine.MCMCSetNIterationsMax(atoi((*beg).c_str()));
                } else if (beg->compare("Seed") == 0) {
                    ++beg;
                    int seed = atoi((*beg).c_str());
                    if (seed != 0)
                        MCEngine.MCMCSetRandomSeed(seed);
                } else if (beg->compare("Iterations") == 0) {
                    ++beg;
                    MCEngine.MCMCSetNIterationsRun(atoi((*beg).c_str()));
                } else if (beg->compare("WriteChain") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0)
                        writechains = true;
                } else if (beg->compare("FindModeWithMinuit") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        FindModeWithMinuit = true;
                    }
                } else if (beg->compare("PrintAllMarginalized") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        PrintAllMarginalized = true;
                    }
                } else if (beg->compare("PrintCorrelationMatrix") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        PrintCorrelationMatrix = true;
                    }
                } else if (beg->compare("PrintKnowledgeUpdatePlots") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        PrintKnowledgeUpdatePlots = true;
                    }
                } else if (beg->compare("PrintParameterPlot") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        PrintParameterPlot = true;
                    }

                } else
                    throw std::runtime_error("\nERROR: Wrong keyword in MonteCarlo config file: " + *beg + "\n Make sure to specify a valid Monte Carlo configuration file.\n");
            } while (!IsEOF);

            BCModelOutput out(&MCEngine, OutFile.c_str());
            if (writechains) {
                out.WriteMarkovChain(true);
                MCEngine.AddChains();
            }

            // set nicer style for drawing than the ROOT default
            BCAux::SetStyle();

            // open log file
            BCLog::OpenLog(("log" + JobTag + ".txt").c_str());
            BCLog::SetLogLevel(BCLog::debug);

            // run the MCMC and marginalize w.r.t. to all parameters
            MCEngine.BCIntegrate::SetNbins(NBINSMODELPARS);
            MCEngine.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
            MCEngine.MarginalizeAll();

            // find mode using the best fit parameters as start values
            if (FindModeWithMinuit)
                MCEngine.FindMode(MCEngine.GetBestFitParameters());

            // draw all marginalized distributions into a PostScript file
            if (PrintAllMarginalized)
                MCEngine.PrintAllMarginalized(("MonteCarlo_plots" + JobTag + ".ps").c_str());

            // print results of the analysis into a text file
            MCEngine.PrintResults(("MonteCarlo_results" + JobTag + ".txt").c_str());

            // print histograms
            MCEngine.PrintHistogram(out, ObsDirName);
            
            BCSummaryTool myBCSummaryTool(&MCEngine);

            // draw the correlation matrix into an eps file
            if (PrintCorrelationMatrix)
                myBCSummaryTool.PrintCorrelationMatrix(("ParamCorrelations" + JobTag + ".eps").c_str());

            // print the correlation matrix into a tex file
            if (PrintCorrelationMatrix)
                MCEngine.PrintCorrelationMatrix(("ParamCorrelations" + JobTag + ".tex").c_str());

            // print comparisons of the prior knowledge to the posterior knowledge 
            // for all parameters into a PostScript file
            if (PrintKnowledgeUpdatePlots)
                myBCSummaryTool.PrintKnowledgeUpdatePlots(("ParamUpdate" + JobTag + ".ps").c_str());

            // draw an overview plot of the parameters into an eps file
            if (PrintParameterPlot)
                myBCSummaryTool.PrintParameterPlot(("ParamSummary" + JobTag + ".eps").c_str());
            
            // print a LaTeX table of the parameters into a tex file
            //myBCSummaryTool.PrintParameterLatex(("ParamSummary" + JobTag + ".tex").c_str());
        
            out.WriteMarginalizedDistributions();
            out.FillAnalysisTree();
            out.Close();

            // print logs for the histograms of the observables into a text file
            std::ofstream outHistoLog;
            outHistoLog.open((ObsDirName + "/HistoLog.txt").c_str(), std::ios::out);
            outHistoLog << MCEngine.getHistoLog();
            outHistoLog.close();

            /* Number of events */
            std::stringstream ss;
            ss << "Number of used events: " << MCEngine.getNumOfUsedEvents();
            BCLog::OutSummary(ss.str().c_str());
            ss.str("");
            ss << "Number of discarded events: " << MCEngine.getNumOfDiscardedEvents();
            BCLog::OutSummary(ss.str().c_str());

            // close log file
            BCLog::CloseLog();

#ifdef _MPI
            double ** sendbuff = new double *[MCEngine.procnum];
            sendbuff[0] = new double[MCEngine.procnum * buffsize];
            for (int il = 1; il < MCEngine.procnum; il++) {
                sendbuff[il] = sendbuff[il - 1] + buffsize;
                sendbuff[il][0] = -1.; //Exit command
            }
            MPI::COMM_WORLD.Scatter(sendbuff[0], buffsize, MPI::DOUBLE,
                                    recvbuff, buffsize, MPI::DOUBLE, 0);
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
