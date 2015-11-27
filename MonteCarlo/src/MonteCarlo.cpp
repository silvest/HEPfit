/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
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
#include <BAT/BCSummaryTool.h>
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
    CalculateNormalization = false;
    PrintAllMarginalized = false;
    PrintCorrelationMatrix = false;
    PrintKnowledgeUpdatePlots = false;
    PrintParameterPlot = false;
    WritePreRunData = false;
    checkrun = false;
    normalization = 0.;
}

//MonteCarlo::~MonteCarlo() {}

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
                for (int i = 0; i < CGP.size(); i++) {
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
                for (int i = 0; i < CGP.size(); i++) {
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
                MPI::COMM_WORLD.Scatter(sendbuff[0], buffsize, MPI::DOUBLE,
                        recvbuff, buffsize, MPI::DOUBLE, 0);

                if (recvbuff[0] == 0.) { // do nothing and return ll
                    double ll = log(0.);
                    MPI::COMM_WORLD.Gather(&ll, 1, MPI::DOUBLE, sendbuff[0], 1, MPI::DOUBLE, 0);
                } else if (recvbuff[0] == 1.) { //compute log likelihood
                    pars.assign(recvbuff + 1, recvbuff + buffsize);
                    double ll = MCEngine.LogEval(pars);
                    MPI::COMM_WORLD.Gather(&ll, 1, MPI::DOUBLE, sendbuff[0], 1, MPI::DOUBLE, 0);
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
                    MPI::COMM_WORLD.Gather(sbuff, obsbuffsize, MPI::DOUBLE, sendbuff[0], obsbuffsize, MPI::DOUBLE, 0);
                } else if (recvbuff[0] == 3.) { // do not compute observables, but gather the buffer
                    double sbuff[obsbuffsize];
                    MPI::COMM_WORLD.Gather(sbuff, obsbuffsize, MPI::DOUBLE, sendbuff[0], obsbuffsize, MPI::DOUBLE, 0);
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
                } else if (beg->compare("NIterationsUpdateMax") == 0) {
                    ++beg;
                    MCEngine.MCMCSetNIterationsUpdateMax(atoi((*beg).c_str()));
                } else if (beg->compare("Seed") == 0) {
                    ++beg;
                    int seed = atoi((*beg).c_str());
                    if (seed != 0)
                        MCEngine.MCMCSetRandomSeed(seed);
                } else if (beg->compare("Iterations") == 0) {
                    ++beg;
                    MCEngine.MCMCSetNIterationsRun(atoi((*beg).c_str()));
                } else if (beg->compare("MinimumEfficiency") == 0) {
                    ++beg;
                    MCEngine.MCMCSetMinimumEfficiency(atof((*beg).c_str()));
                } else if (beg->compare("WriteChain") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0)
                        writechains = true;
                } else if (beg->compare("FindModeWithMinuit") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        FindModeWithMinuit = true;
                    }
                } else if (beg->compare("CalculateNormalization") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        CalculateNormalization = true;
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
                } else if (beg->compare("WritePreRunData") == 0) {
                    ++beg;
                    if (beg->compare("true") == 0) {
                        WritePreRunData = true;
                    }
                } else if (beg->compare("ReadPreRunData") == 0) {
                    ++beg;
                    ReadPreRunData(*beg);
                } else if (beg->compare("OrderParameters") == 0) {
                    ++beg;
                    if (beg->compare("false") == 0) {
                        MCEngine.MCMCSetFlagOrderParameters(false);
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

            if (CalculateNormalization) normalization = MCEngine.computeNormalization();

            // draw all marginalized distributions into a pdf file
            if (PrintAllMarginalized)
                MCEngine.PrintAllMarginalized(("MonteCarlo_plots" + JobTag + ".pdf").c_str());

            // print results of the analysis into a text file
            MCEngine.PrintResults(("MonteCarlo_results" + JobTag + ".txt").c_str());

            // print histograms
            MCEngine.PrintHistogram(out, ObsDirName);

            BCSummaryTool myBCSummaryTool(&MCEngine);

            // draw the correlation matrix into a pdf file
            if (PrintCorrelationMatrix)
                myBCSummaryTool.PrintCorrelationMatrix(("ParamCorrelations" + JobTag + ".pdf").c_str());

            // print the correlation matrix into a tex file
            if (PrintCorrelationMatrix)
                MCEngine.PrintCorrelationMatrix(("ParamCorrelations" + JobTag + ".tex").c_str());

            // print comparisons of the prior knowledge to the posterior knowledge 
            // for all parameters into a pdf file
            if (PrintKnowledgeUpdatePlots)
                myBCSummaryTool.PrintKnowledgeUpdatePlots(("ParamUpdate" + JobTag + ".pdf").c_str());

            // draw an overview plot of the parameters into a pdf file
            if (PrintParameterPlot)
                myBCSummaryTool.PrintParameterPlot(("ParamSummary" + JobTag + ".pdf").c_str());

            // print a LaTeX table of the parameters into a tex file
            //myBCSummaryTool.PrintParameterLatex(("ParamSummary" + JobTag + ".tex").c_str());

            out.WriteMarginalizedDistributions();
            out.FillAnalysisTree();
            out.Close();

            // print logs for the histograms of the observables into a text file
            std::ofstream outHistoLog;
            outHistoLog.open((ObsDirName + "/HistoLog" + JobTag + ".txt").c_str(), std::ios::out);
            outHistoLog << MCEngine.getHistoLog();
            outHistoLog.close();

            // print statistics for the theory values of the observables into a text file
            std::ofstream outStatLog;
            outStatLog.open((ObsDirName + "/Statistics" + JobTag + ".txt").c_str(), std::ios::out);
            if (CalculateNormalization) outStatLog << "Normalization for " << ModelName.c_str() << ": " << normalization << "\n" << std::endl;
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
    for (int i = 0; i < MCEngine.MCMCGetNChains(); i++){
        mode_all.insert(mode_all.end(), mode.begin(), mode.end());
        scale_all.insert(scale_all.end(), scale.begin(), scale.end());
    }
    MCEngine.MCMCSetInitialPositions(mode_all);
    MCEngine.MCMCSetTrialFunctionScaleFactor(scale_all);
}

void MonteCarlo::addCustomObservableType(const std::string name, boost::function<Observable*() > funct) {
    myInputParser.addCustomObservableType(name, funct);
}