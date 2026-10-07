/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MonteCarloEngine.h"
#include "StandardModel.h"
#include <BAT/BCParameter.h>
#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCTF1Prior.h>
#include <BAT/BCCombinedPrior.h>
#include <BAT/BCLog.h>
#ifdef _MPI
#include <mpi.h>
#endif
#include <TF1.h>
#include <TTree.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <limits>
#include <gsl/gsl_eigen.h>
#include <memory>
#include <sstream>

MonteCarloEngine::MonteCarloEngine(
        const std::vector<ModelParameter>& ModPars_i,
        boost::ptr_vector<Observable>& Obs_i,
        std::vector<Observable2D>& Obs2D_i,
        std::vector<CorrelatedGaussianObservables>& CGO_i,
        std::vector<CorrelatedGaussianParameters>& CGP_i)
: BCModel(""), ModPars(ModPars_i), CGP(CGP_i), Obs_ALL(Obs_i), Obs2D_ALL(Obs2D_i),
  CGO(CGO_i), NumOfUsedEvents(0), NumOfDiscardedEvents(0) {
    SetMultivariateCovarianceUpdateLambda(0.5);
    obval = NULL;
    obweight = NULL;
    Mod = NULL;
    hessianNoise = 0.;
    cindex = 0;
    printLogo = false;
    nSmooth = 0;
    histogram2Dtype = 1001;
    noLegend = true;
    PrintLoglikelihoodPlots = false;
    WriteLogLikelihoodChain = false;
    WriteParametersChain = false;
    WriteMCMCweights = false;
    alpha2D = 1.;
    kchainedObs = 0;
    kwmcmc = 0;
    nBins1D = NBINS1D;
    nBins2D = NBINS2D;
    significants = 0;
    histogramBufferSize = 0;
    LogLikelihood_max = std::numeric_limits<double>::lowest();

#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif
    if (rank == 0) {
        TH1::StatOverflows(kTRUE);
        TH1::SetDefaultBufferSize(100000);

#if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
        gIdx = TColor::GetFreeColorIndex();
        rIdx = TColor::GetFreeColorIndex() + 1;
#else
        gIdx = 1000;
        rIdx = 1001;
#endif

        HEPfit_green = new TColor(gIdx, 0.0, 0.56, 0.57, "HEPfit_green");
        HEPfit_red = new TColor(rIdx, 0.57, 0.01, 0.00, "HEPfit_red");
    }
};


void MonteCarloEngine::Initialize(StandardModel* Mod_i)
{
    Mod = Mod_i;
    int k = 0, kweight = 0, kmcmcweight = 0;

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            k++;
            if (it->getDistr().compare("noweight") != 0) kweight++;
            if (it->isWriteChain()) kchainedObs++;
        } else {
            kmcmcweight++;
        }
        thMin[it->getName()] = std::numeric_limits<double>::max();
        thMax[it->getName()] = -std::numeric_limits<double>::max();
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            if (!it->isTMCMC())
                throw std::runtime_error("ERROR: cannot handle noMCMC for Observable2D file yet!");
        } else if (it->getDistr().compare("weight") == 0)
            throw std::runtime_error("ERROR: do not use Observable2D for analytic 2D weights!");
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 != CGO.end(); ++it1) {
        std::vector<Observable> ObsV(it1->getObs());
        for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it) {
            if ((it->getDistr()).compare("file") == 0)
                throw std::runtime_error("Cannot use file in CorrelatedGaussianObservables!");
            if (!(it->isTMCMC())) {
                k++;
                if (it->getDistr().compare("noweight") != 0)
                    throw std::runtime_error("Cannot use weight in CorrelatedGaussianObservables!");
            }
            thMin[it->getName()] = std::numeric_limits<double>::max();
            thMax[it->getName()] = -std::numeric_limits<double>::max();
        }
        if (it1->isPrediction()) {
            CorrelationMap[it1->getName()] = new TPrincipal(it1->getObs().size(), "N");
        } else {
            kmcmcweight += 1 + ObsV.size(); // full weight + leave-one-out weights
        }
    }
    kmax = k;
    kwmax = kweight;
    kwmcmc = kmcmcweight;

    unknownParameters = Mod->getUnknownParameters();
    DefineParameters();
};

void MonteCarloEngine::CreateHistogramMaps() 
{
    if (histogramBufferSize != 0) TH1::SetDefaultBufferSize(histogramBufferSize);
    
    TH1D * lhisto = new TH1D("LogLikelihood", "LogLikelihood", nBins1D, 1., -1.);
    lhisto->GetXaxis()->SetTitle("LogLikelihood");
    BCH1D bclhisto = BCH1D(lhisto);
    Histo1D["LogLikelihood"] = bclhisto;

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it != Obs_ALL.end(); it++) {
        std::string HistName = it->getName();
        if (Histo1D.find(HistName) == Histo1D.end()) {
            TH1D * histo = new TH1D(HistName.c_str(), it->getLabel().c_str(), nBins1D, it->getMin(), it->getMax());
            histo->GetXaxis()->SetTitle(it->getLabel().c_str());
            BCH1D bchisto = BCH1D(histo);
            Histo1D[HistName] = bchisto;
        }
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it != Obs2D_ALL.end(); it++) {
        std::string HistName = it->getName();
        if (Histo2D.find(HistName) == Histo2D.end()) {
            TH2D * histo2 = new TH2D(HistName.c_str(), (it->getLabel() + " vs. " + it->getLabel2()).c_str(), nBins2D, it->getMin(), it->getMax(), nBins2D, it->getMin2(), it->getMax2());
            histo2->GetXaxis()->SetTitle(it->getLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel2().c_str());
            BCH2D bchisto2 = BCH2D(histo2);
            Histo2D[HistName] = bchisto2;
        }
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 != CGO.end(); ++it1) {
        std::vector<Observable> ObsV(it1->getObs());
        for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it) {
            std::string HistName = it->getName();
            if (Histo1D.find(HistName) == Histo1D.end()) {
                TH1D * histo = new TH1D(HistName.c_str(), it->getLabel().c_str(), nBins1D, it->getMin(), it->getMax());
                histo->GetXaxis()->SetTitle(it->getLabel().c_str());
                BCH1D bchisto = BCH1D(histo);
                Histo1D[HistName] = bchisto;
            }
        }
    }

    if (PrintLoglikelihoodPlots) {
        for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it != ModPars.end(); it++) {
            if (it->IsFixed()) continue;
            if (std::find(unknownParameters.begin(), unknownParameters.end(), it->getname()) != unknownParameters.end()) continue;
            std::string HistName = it->getname() + "_vs_LogLikelihood";
            if (Histo2D.find(HistName) == Histo2D.end()) {
                TH2D * histo2 = new TH2D(HistName.c_str(), (it->getname() + " vs. LogLikelihood").c_str(), nBins2D, 1., -1., nBins2D, 1., -1.);
                histo2->GetXaxis()->SetTitle(it->getname().c_str());
                histo2->GetYaxis()->SetTitle("LogLikelihood");
                BCH2D bchisto2 = BCH2D(histo2);
                Histo2D[HistName] = bchisto2;
            }
        }
    }
};

void MonteCarloEngine::setNChains(unsigned int i) {
    SetNChains(i);
    obval = new double[fMCMCNChains * kmax];
    obweight = new double[fMCMCNChains * kwmax];
}

// ---------------------------------------------------------

MonteCarloEngine::~MonteCarloEngine()
{
    if (rank == 0) { // These are created only by the master
        delete [] obval;
        delete [] obweight;
        delete HEPfit_red;
        delete HEPfit_green;
        HEPfit_red = NULL;
        HEPfit_green = NULL;
        if (CorrelationMap.size() > 0) {
            for (std::map<std::string, TPrincipal *>::iterator it = CorrelationMap.begin(); it != CorrelationMap.end(); it++) {
                delete it->second;
                it->second = NULL;
            }
        }
    }
};

// ---------------------------------------------------------

void MonteCarloEngine::DefineParameters() {
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    if (rank == 0) std::cout << "\nParameters varied in this run:" << std::endl;
    unsigned int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin();
            it < ModPars.end(); it++) {
        if (std::find(unknownParameters.begin(), unknownParameters.end(), it->getname()) == unknownParameters.end()) {
            if (it->geterrf() == 0. && it->geterrg() == 0.)
                continue;

            AddParameter(it->getname().c_str(), it->getmin(), it->getmax());
            if (rank == 0) std::cout << k << ": " << it->getname() << ", ";

            if (it->IsCorrelated()) {
                for (unsigned int i = 0; i < CGP.size(); i++) {
                    if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                        std::string index = it->getname().substr(CGP[i].getName().size());
                        long int lindex = strtol(index.c_str(), NULL, 10);
                        if (lindex > 0)
                            DPars[CGP[i].getPar(lindex - 1).getname()] = 0.;
                        else {
                            std::stringstream out;
                            out << it->getname();
                            throw std::runtime_error("MonteCarloEngine::DefineParameters(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                        }
                    }
                }
            } else
                DPars[it->getname()] = 0.;
            if (it->geterrf() == 0.) GetParameter(k).SetPrior(std::make_shared<BCGaussianPrior>(it->getave(), it->geterrg())); //SetPriorGauss(k, it->getave(), it->geterrg());
            else if (it->geterrg() == 0.) GetParameter(k).SetPriorConstant(); //SetPriorConstant(k);
            else {
	      GetParameter(k).SetPrior(std::make_shared<BCCombinedPrior>(it->getave(), it->geterrg(), it->geterrf())); //SetPrior(k, combined);
            }
            k++;
        }

    }
    if (unknownParameters.size() > 0 && rank == 0) {
        std::cout << "\n" << std::endl;
        for (std::vector<std::string>::iterator it = unknownParameters.begin(); it != unknownParameters.end(); it++)
            std::cout << "WARNING: unknown parameter " << *it << " not added to MCMC" << std::endl;
    }
}

void MonteCarloEngine::setDParsFromParameters(const std::vector<double>& parameters, 
        std::map<std::string,double>& DPars_i) 
{
    std::map<std::string, std::vector<double> > cgpmap;

    unsigned int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it != ModPars.end(); it++){
        if(it->IsFixed())
            continue;
        if (std::find(unknownParameters.begin(), unknownParameters.end(), it->getname()) != unknownParameters.end())
            continue;
        if(it->getname().compare(GetParameter(k).GetName()) != 0)
            {
                        std::stringstream out;
                        out << it->getname();
                        throw std::runtime_error("MonteCarloEngine::setDParsFromParameters(): " + out.str() + "is sitting at the wrong position in the BAT parameters vector");
                    }
        if (it->IsCorrelated()) {
            std::string index = it->getname().substr(it->getCgp_name().size());
            unsigned long int lindex = strtol(index.c_str(),NULL,10);
            if (lindex - 1 == cgpmap[it->getCgp_name()].size())
                cgpmap[it->getCgp_name()].push_back(parameters[k]);
            else {
                std::stringstream out;
                out << it->getname() << " " << lindex;
                throw std::runtime_error("MonteCarloEngine::setDParsFromParameters(): " + out.str() + "seems to be a CorrelatedGaussianParameters object but the corresponding parameters are missing or not in the right order");
            }

        } else
            DPars_i[it->getname()] = parameters[k];
        k++;
    }

    for (unsigned int j = 0; j < CGP.size(); j++) {
        std::vector<double> current = cgpmap.at(CGP[j].getName());
        if (current.size() != CGP[j].getPars().size()) {
            std::stringstream out;
            out << CGP[j].getName();
            throw std::runtime_error("MonteCarloEngine::setDParsFromParameters(): " + out.str() + " appears to be represented in cgpmap with a wrong size");
        }
        
        std::vector<double> porig = CGP[j].getOrigParsValue(current);

        for(unsigned int l = 0; l < porig.size(); l++) {
            DPars_i[CGP[j].getPar(l).getname()] = porig[l];
        }
    }
}

// ---------------------------------------------------------

double MonteCarloEngine::LogLikelihood(const std::vector<double>& parameters) {
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.

    double logprob = 0.;
    
    setDParsFromParameters(parameters, DPars);

    // if update false set probability equal zero
    if (!Mod->Update(DPars)) {
#ifdef _MCDEBUG
        std::cout << "event discarded" << std::endl;

        /* Debug */
        //for (int k = 0; k < parameters.size(); k++)
        //    std::cout << "  " << GetParameter(k)->GetName() << " = "
        //              << DPars[GetParameter(k)->GetName()] << std::endl;
#endif
        NumOfDiscardedEvents++;
        return (log(0.));
    }
#ifdef _MCDEBUG
    //std::cout << "event used in MC" << std::endl;
#endif

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it != Obs_ALL.end(); it++) {
        if (it->isTMCMC()) logprob += it->computeWeight();
    }

    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it != Obs2D_ALL.end(); it++) {
        if (it->isTMCMC()) logprob += it->computeWeight();
    }

    for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++) {
        if(!(it->isPrediction())) logprob += it->computeWeight();
    }
    if (!std::isfinite(logprob) || !Mod->isQCDsuccess() || !Mod->isSMSuccess()) {
        NumOfDiscardedEvents++;
#ifdef _MCDEBUG
//        std::cout << "Event discarded since logprob evaluated to: " << logprob << std::endl ;
#endif
        return (log(0.));
    }
    NumOfUsedEvents++;
    return logprob;
}

void MonteCarloEngine::MCMCUserIterationInterface() {  
#ifdef _MPI
    unsigned mychain = 0;
    int iproc = 0;
    unsigned npars = GetNParameters();
    int buffsize = npars + 1;
    int index_chain[procnum];
    double *recvbuff = new double[buffsize];
    std::vector<double> pars;
    double **buff;

    buff = new double*[procnum];
    int obsbuffsize = Obs_ALL.size() + 2 * Obs2D_ALL.size();
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++)
        obsbuffsize += it1->getObs().size();
    buff[0] = new double[procnum * obsbuffsize];
    for (int i = 1; i < procnum; i++) {
        buff[i] = buff[i - 1] + obsbuffsize;
        index_chain[i] = -1;
    }

    double ** sendbuff = new double *[procnum];
    sendbuff[0] = new double[procnum * buffsize];
    for (int il = 1; il < procnum; il++)
        sendbuff[il] = sendbuff[il - 1] + buffsize;

    while (mychain < fMCMCNChains) {
        pars.clear();
        pars = fMCMCStates.at(mychain).parameters;

        if (PrintLoglikelihoodPlots || WriteParametersChain) {
            std::map<std::string, double> tmpDPars;
            setDParsFromParameters(pars, tmpDPars);
            if (PrintLoglikelihoodPlots) DPars_allChains.push_back(tmpDPars);
            if (WriteParametersChain) {
                int k = 0;
                for (std::map<std::string, double>::iterator it = tmpDPars.begin(); it != tmpDPars.end(); it++) hMCMCParameters[mychain][k++] = it->second;
            }
        }

        index_chain[iproc] = mychain;
        iproc++;
        mychain++;
        if (iproc < procnum && mychain < fMCMCNChains)
            continue;

        for (int il = 0; il < iproc; il++) {
            //The first entry of the array specifies the task to be executed.

            sendbuff[il][0] = 2.; // 2 = observables calculation
            for (int im = 1; im < buffsize; im++) sendbuff[il][im] = fMCMCStates.at(index_chain[il]).parameters.at(im - 1);
        }
        for (int il = iproc; il < procnum; il++) {
            sendbuff[il][0] = 3.; // 3 = nothing to execute, but return a buffer of observables
            index_chain[il] = -1;
        }
        //       double inittime = MPI::Wtime();
        MPI_Scatter(sendbuff[0], buffsize, MPI_DOUBLE,
                recvbuff, buffsize, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

        if (recvbuff[0] == 2.) { // compute observables
            double sbuff[obsbuffsize];
            pars.assign(recvbuff + 1, recvbuff + buffsize);
            setDParsFromParameters(pars,DPars);
            Mod->Update(DPars);
                
            int k = 0;
            for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
                sbuff[k++] = it->computeTheoryValue();
            }
            for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end(); it++) {
                sbuff[k++] = it->computeTheoryValue();
                sbuff[k++] = it->computeTheoryValue2();
            }

            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                std::vector<Observable> ObsV(it1->getObs());
                for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it)
                    sbuff[k++] = it->computeTheoryValue();
            }
            MPI_Gather(sbuff, obsbuffsize, MPI_DOUBLE, buff[0], obsbuffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else if (recvbuff[0] == 3.) { // do not compute observables, but gather the buffer
            double sbuff[obsbuffsize];
            MPI_Gather(sbuff, obsbuffsize, MPI_DOUBLE, buff[0], obsbuffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        for (int il = 0; il < procnum; il++) {
            if (index_chain[il] >= 0) {
                int k = 0;
                // fill the histograms for observables
                int ko = 0, kweight = 0, kweight_ob = 0, k_all = 0, k_cObs = 0;
                for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                        it < Obs_ALL.end(); it++) {
                    double th = buff[il][k++];
                    /* set the min and max of theory values */
                    if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                    if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                    Histo1D[it->getName()].GetHistogram()->Fill(th);
                    if (fMCMCFlagWriteChainToFile) hMCMCObservables[index_chain[il]][k_all++] = th;
                    else if (getchainedObsSize() > 0 && it->isWriteChain()) hMCMCObservables[index_chain[il]][k_cObs++] = th;
                    if (!it->isTMCMC()) {
                        obval[index_chain[il] * kmax + ko] = th;
                        ko++;
                        if (it->getDistr().compare("noweight") != 0 && it->getDistr().compare("writeChain") != 0) {
                            double weight = it->computeWeight(th);
                            obweight[index_chain[il] * kwmax + kweight_ob] = weight;
                            kweight_ob++;
                            if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[index_chain[il]][kweight] = weight;
                            kweight++;
                        }
                    } else if (WriteMCMCweights) {
                        double weight = it->computeWeight(th);
                        if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[index_chain[il]][kweight] = weight;
                        kweight++;
                    }
                }

                // fill the 2D histograms for observables
                for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
                        it < Obs2D_ALL.end(); it++) {
                    double th1 = buff[il][k++];
                    double th2 = buff[il][k++];
                    Histo2D[it->getName()].GetHistogram()->Fill(th1, th2);
                    if (fMCMCFlagWriteChainToFile) {
                        hMCMCObservables[index_chain[il]][k_all++] = th1;
                        hMCMCObservables[index_chain[il]][k_all++] = th2;
                    }
                }

                // fill the histograms for correlated observables
                for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                    std::vector<Observable> ObsV(it1->getObs());
                    Double_t * COdata = new Double_t[ObsV.size()];
                    int nObs = 0;
                    for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it) {
                        double th = buff[il][k++];
                        /* set the min and max of theory values */
                        if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                        if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                        Histo1D[it->getName()].GetHistogram()->Fill(th);
                        if (fMCMCFlagWriteChainToFile) hMCMCObservables[index_chain[il]][k_all++] = th;                        
                        if (it1->isPrediction()) COdata[nObs++] = th;
                    }
                    if (it1->isPrediction()) CorrelationMap[it1->getName()]->AddRow(COdata);
                    delete [] COdata;
                    if (WriteMCMCweights && !it1->isPrediction()) {
                        std::vector<double> weights = it1->computeLeaveOneOutWeights();
                        for (unsigned int iw = 0; iw < weights.size(); iw++) {
                            if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[index_chain[il]][kweight] = weights[iw];
                            kweight++;
                        }
                    }
                }
            }
        }
        iproc = 0;
    }
    if (fMCMCFlagWriteChainToFile || getchainedObsSize() > 0 || WriteLogLikelihoodChain) InChainFillObservablesTree();
    if (WriteParametersChain) InChainFillParametersTree();
    delete sendbuff[0];
    delete [] sendbuff;
    delete [] recvbuff;
    delete buff[0];
    delete [] buff;
#else
    for (unsigned int i = 0; i < fMCMCNChains; ++i) {
        // NOTE: BAT syncs fMCMCThreadLocalStorage with fMCMCStates before calling MCMCUserIterationInterface.
        std::vector<double>::const_iterator first = fMCMCStates.at(i).parameters.begin(); 
        std::vector<double>::const_iterator last = first + GetNParameters();
        std::vector<double> currvec(first, last);
        setDParsFromParameters(currvec,DPars);
        if (PrintLoglikelihoodPlots) DPars_allChains.push_back(DPars);
        if (WriteParametersChain) {
            int k = 0;
            for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++) hMCMCParameters[i][k++] = it->second;
        }

        Mod->Update(DPars);
        // fill the histograms for observables
        int k = 0, kweight = 0, kweight_ob = 0, k_all = 0, k_cObs = 0;
        for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) {
            double th = it->computeTheoryValue();
            /* set the min and max of theory values */
            if (th < thMin[it->getName()]) thMin[it->getName()] = th;
            if (th > thMax[it->getName()]) thMax[it->getName()] = th;
            Histo1D[it->getName()].GetHistogram()->Fill(th);
            if (fMCMCFlagWriteChainToFile) hMCMCObservables[i][k_all++] = th;
            else if (getchainedObsSize() > 0 && it->isWriteChain()) hMCMCObservables[i][k_cObs++] = th;
            if (!it->isTMCMC()) {
                obval[i * kmax + k] = th;
                k++;
                if (it->getDistr().compare("noweight") != 0 && it->getDistr().compare("writeChain") != 0) {
                    double weight = it->computeWeight(th);
                    obweight[i * kwmax + kweight_ob] = weight;
                    kweight_ob++;
                    if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[i][kweight] = weight;
                    kweight++;
                }
            } else if (WriteMCMCweights) {
                double weight = it->computeWeight(th);
                if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[i][kweight] = weight;
                kweight++;
            }
        }
        
        // fill the 2D histograms for observables
        for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
                it < Obs2D_ALL.end(); it++) {
            double th1 = it->computeTheoryValue();
            double th2 = it->computeTheoryValue2();
            Histo2D[it->getName()].GetHistogram()->Fill(th1, th2);
            if (fMCMCFlagWriteChainToFile) {
                hMCMCObservables[i][k_all++] = th1;
                hMCMCObservables[i][k_all++] = th2;
            }
        }

        // fill the histograms for correlated observables
        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                it1 < CGO.end(); it1++) {
            std::vector<Observable> ObsV(it1->getObs());
            Double_t * COdata = new Double_t[ObsV.size()];
            int nObs = 0;
            for (std::vector<Observable>::iterator it = ObsV.begin();
                    it != ObsV.end(); ++it) {
                double th = it->computeTheoryValue();
                /* set the min and max of theory values */
                if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                Histo1D[it->getName()].GetHistogram()->Fill(th);
                if (fMCMCFlagWriteChainToFile) hMCMCObservables[i][k_all++] = th;                        
                if (it1->isPrediction()) COdata[nObs++] = th;
            }
            if (it1->isPrediction()) CorrelationMap[it1->getName()]->AddRow(COdata);
            delete [] COdata;
            if (WriteMCMCweights && !it1->isPrediction()) {
                std::vector<double> weights = it1->computeLeaveOneOutWeights();
                for (unsigned int iw = 0; iw < weights.size(); iw++) {
                    if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[i][kweight] = weights[iw];
                    kweight++;
                }
            }
        }
    }
    
    if (fMCMCFlagWriteChainToFile || getchainedObsSize() > 0  || WriteLogLikelihoodChain) InChainFillObservablesTree();
    if (WriteParametersChain) InChainFillParametersTree();
#endif
    for (unsigned int i = 0; i < fMCMCNChains; i++) {
        double LogLikelihood = fMCMCStates.at(i).log_likelihood;
        if (LogLikelihood > LogLikelihood_max) {
            LogLikelihood_max = std::max(LogLikelihood, LogLikelihood_max);
            par_at_LL_max = Getx(i); // NOTE: Unrotated, to be rotated later.
        }
        Histo1D["LogLikelihood"].GetHistogram()->Fill(LogLikelihood);
        if (PrintLoglikelihoodPlots) {
            for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it != ModPars.end(); it++) {
                if (it->IsFixed()) continue;
                if (std::find(unknownParameters.begin(), unknownParameters.end(), it->getname()) != unknownParameters.end()) continue;
                std::string HistName = it->getname() + "_vs_LogLikelihood";
                Histo2D[HistName].GetHistogram()->Fill(DPars_allChains.at(i).at(it->getname()), LogLikelihood);
            }
        }
    }
    if (PrintLoglikelihoodPlots) DPars_allChains.clear();
}

void MonteCarloEngine::CheckHistogram(TH1& hist, const std::string name) {
    double UnderFlowContent = hist.GetBinContent(0);
    double OverFlowContent = hist.GetBinContent(nBins1D + 1);
    double Integral = hist.Integral();
    double TotalContent = 0.0;
    for (unsigned int n = 0; n <= nBins1D + 1; n++)
        TotalContent += hist.GetBinContent(n);
    HistoLog << name << ": "
            << Integral / TotalContent * 100. << "% within the range, "
            << UnderFlowContent / TotalContent * 100. << "% underflow, "
            << OverFlowContent / TotalContent * 100. << "% overflow"
            << std::endl;
}

void MonteCarloEngine::CheckHistogram(TH2& hist, const std::string name) {
    double Integral = hist.Integral();
    double TotalContent = 0.0;
    for (unsigned int m = 0; m <= nBins2D + 1; m++)
        for (unsigned int n = 0; n <= nBins2D + 1; n++)
            TotalContent += hist.GetBinContent(m, n);
    HistoLog << name << ": "
            << Integral / TotalContent * 100. << "% within the ranges"
            << std::endl;
}

void MonteCarloEngine::Print1D(BCH1D bch1d, const char* filename, int ww, int wh) {
    TCanvas * c;
    cindex++;
    if(ww > 0 && wh > 0)
        c = new TCanvas(TString::Format("c_bch1d_%d",cindex), TString::Format("c_bch1d_%d",cindex), ww, wh);
    else
        c = new TCanvas(TString::Format("c_bch1d_%d",cindex));

    bch1d.GetHistogram()->Scale(1./bch1d.GetHistogram()->Integral("width"));
    
    bch1d.SetBandType(BCH1D::kSmallestInterval);
    bch1d.SetBandColor(0, gIdx);
    bch1d.SetBandColor(1, rIdx);
    bch1d.SetBandColor(2, kOrange - 3);
    bch1d.SetNBands(3);
    bch1d.SetNSmooth(nSmooth);
    bch1d.SetDrawGlobalMode(true);
    bch1d.SetDrawMean(true, true);
    bch1d.SetDrawLegend(!noLegend);
    if (noLegend) gStyle->SetOptStat("emr");
    bch1d.SetNLegendColumns(1);
    bch1d.SetStats(true);
    
    bch1d.Draw();
    
    if (printLogo) {
        double xRange = (bch1d.GetHistogram()->GetXaxis()->GetXmax() - bch1d.GetHistogram()->GetXaxis()->GetXmin())*3./4.;
        double yRange = (bch1d.GetHistogram()->GetMaximum() - bch1d.GetHistogram()->GetMinimum());
        
        double xL;
        if (noLegend) xL = bch1d.GetHistogram()->GetXaxis()->GetXmin()+0.0475*xRange;
        else xL = bch1d.GetHistogram()->GetXaxis()->GetXmin() + 0.0375 * xRange;
        double yL = bch1d.GetHistogram()->GetYaxis()->GetXmin() + 0.89 * yRange;

        double xR;
        if (noLegend) xR = xL + 0.21*xRange;
        else xR = xL + 0.18 * xRange;
        double yR = yL + 0.09 * yRange;

        TBox b1 = TBox(xL, yL, xR, yR);
        b1.SetFillColor(gIdx);
        
        TBox b2; 
        b2 = TBox(xL+0.008*xRange, yL+0.008*yRange, xR-0.008*xRange, yR-0.008*yRange);
        b2.SetFillColor(kWhite);
        
        TPaveText b3 = TPaveText(xL+0.014*xRange, yL+0.013*yRange, xL+0.70*(xR-xL), yR-0.013*yRange);
        if (noLegend) b3.SetTextSize(0.056);
        else b3.SetTextSize(0.051);
        b3.SetTextAlign(22);
        b3.SetTextColor(kWhite);
        b3.AddText("HEP");
        b3.SetFillColor(rIdx);
        
        TPaveText * b4; 
        if (noLegend) {
            b4 = new TPaveText(xL+0.72*(xR-xL), yL+0.030*yRange, xR-0.008*xRange, yR-0.013*yRange);
            b4->SetTextSize(0.048);
        } else {
            b4 = new TPaveText(xL + 0.75 * (xR - xL), yL + 0.024 * yRange, xR - 0.008 * xRange, yR - 0.013 * yRange);
            b4->SetTextSize(0.039);
        }
        b4->SetTextAlign(33);
        b4->SetTextColor(rIdx);
        b4->AddText("fit");
        b4->SetFillColor(kWhite);

        b1.Draw("SAME");
        b2.Draw("SAME");
        b3.Draw("SAME");
        b4->Draw("SAME");
        
        c->Print(filename);
        delete b4;
        b4 = NULL;
    } else c->Print(filename);
    
    delete c;
    c = NULL;
}

void MonteCarloEngine::Print2D(BCH2D bch2d, const char * filename, int ww, int wh)
{
    TCanvas * c;
    cindex++;
    if(ww > 0 && wh > 0)
        c = new TCanvas(TString::Format("c_bch2d_%d",cindex), TString::Format("c_bch2d_%d",cindex), ww, wh);
    else
        c = new TCanvas(TString::Format("c_bch2d_%d",cindex));
    
    bch2d.GetHistogram()->Scale(1./bch2d.GetHistogram()->Integral("width"));
    bch2d.GetHistogram()->GetYaxis()->SetTitleOffset(1.45);

    bch2d.SetBandType(BCH2D::kSmallestInterval);
    bch2d.SetBandColor(0, TColor::GetColorTransparent(kOrange - 3, alpha2D)); 
    bch2d.SetBandColor(1, TColor::GetColorTransparent(rIdx, alpha2D));
    bch2d.SetBandColor(2, TColor::GetColorTransparent(gIdx, alpha2D));
    bch2d.SetNBands(3);
    bch2d.SetBandFillStyle(histogram2Dtype);// Type of 2D Histogram 1001 -> box pixel, 101 -> filled, 1 -> contour.
    if (histogram2Dtype == 1 || histogram2Dtype == 101) bch2d.SetNSmooth(1);
    else bch2d.SetNSmooth(0);
    if (histogram2Dtype == 1) bch2d.GetHistogram()->SetLineWidth(3);
    bch2d.SetDrawLocalMode(false);
    bch2d.SetDrawGlobalMode(true);
    bch2d.SetDrawMean(true, true);
    bch2d.SetDrawLegend(!noLegend);
    if (noLegend) gStyle->SetOptStat("emr");
    bch2d.SetNLegendColumns(1);
    bch2d.SetStats(true);
    
    bch2d.Draw();

    if (printLogo) {
        double xRange = (bch2d.GetHistogram()->GetXaxis()->GetXmax() - bch2d.GetHistogram()->GetXaxis()->GetXmin())*3./4.;
        double yRange = bch2d.GetHistogram()->GetYaxis()->GetXmax() - bch2d.GetHistogram()->GetYaxis()->GetXmin();

        double xL;
        if (noLegend) xL = bch2d.GetHistogram()->GetXaxis()->GetXmin()+0.0475*xRange;
        else xL = bch2d.GetHistogram()->GetXaxis()->GetXmin() + 0.0375 * xRange;
        double yL = bch2d.GetHistogram()->GetYaxis()->GetXmin() + 0.89 * yRange;

        double xR;
        if (noLegend) xR = xL + 0.21*xRange;
        else xR = xL + 0.18 * xRange;
        double yR = yL + 0.09 * yRange;

        TBox b1 = TBox(xL, yL, xR, yR);
        b1.SetFillColor(gIdx);

        TBox b2 = TBox(xL + 0.008 * xRange, yL + 0.008 * yRange, xR - 0.008 * xRange, yR - 0.008 * yRange);
        b2.SetFillColor(kWhite);

        TPaveText b3 = TPaveText(xL + 0.014 * xRange, yL + 0.013 * yRange, xL + 0.70 * (xR - xL), yR - 0.013 * yRange);
        if (noLegend) b3.SetTextSize(0.056);
        else b3.SetTextSize(0.051);
        b3.SetTextAlign(22);
        b3.SetTextColor(kWhite);
        b3.AddText("HEP");
        b3.SetFillColor(rIdx);

        TPaveText * b4;
        if (noLegend) {
            b4 = new TPaveText(xL+0.72*(xR-xL), yL+0.030*yRange, xR-0.008*xRange, yR-0.013*yRange);
            b4->SetTextSize(0.048);
        } else {
            b4 = new TPaveText(xL + 0.75 * (xR - xL), yL + 0.024 * yRange, xR - 0.008 * xRange, yR - 0.013 * yRange);
            b4->SetTextSize(0.039);
        }
        b4->SetTextAlign(33);
        b4->SetTextColor(rIdx);
        b4->AddText("fit");
        b4->SetFillColor(kWhite);
        
        b1.Draw("SAME");
        b2.Draw("SAME");
        b3.Draw("SAME");
        b4->Draw("SAME");
        
        c->Print(filename);
        delete b4;
        b4 = NULL;
    } else c->Print(filename);
    
    delete c;
    c = NULL;
}

void MonteCarloEngine::PrintHistogram(std::string& OutFile, Observable& it, const std::string OutputDir) 
{
    
    std::string HistName = it.getName();
    double min = thMin[it.getName()];
    double max = thMax[it.getName()];
    if (Histo1D[HistName].GetHistogram()->Integral() > 0.0) {
        std::string fname = OutputDir + "/" + HistName + ".pdf";
        Histo1D[HistName].SetGlobalMode(it.computeTheoryValue());
        Print1D(Histo1D[HistName], fname.c_str());
        std::cout << " " + HistName + ".pdf" << std::endl;
        if(OutFile.compare("") == 0) {
            throw std::runtime_error("\nMonteCarloEngine::PrintHistogram ERROR: No root file specified for writing histograms.");
        }
        TDirectory * dir = gDirectory;
        GetOutputFile()->cd();
        Histo1D[HistName].GetHistogram()->Write();
        gDirectory = dir;
        CheckHistogram(*Histo1D[HistName].GetHistogram(), it.getName());
    } else
        HistoLog << "WARNING: The histogram of "
            << it.getName() << " is empty!" << std::endl;

    HistoLog.precision(10);
    HistoLog << "  [min, max]=[" << min << ", " << max << "]" << std::endl;
    HistoLog.precision(6);
}

void MonteCarloEngine::PrintHistogram(std::string& OutFile, const std::string OutputDir) 
{
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables histogram cannot be generated.\n");
    setDParsFromParameters(mode,DPars);

    Mod->Update(DPars);

    if (Obs_ALL.size() != 0 || CGO.size() != 0) std::cout << "\nPrinting 1D histograms in the Observables directory: " << std::endl;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) PrintHistogram(OutFile, *it, OutputDir);
    
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
        std::vector<Observable> ObsV(it1->getObs());
        for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it) PrintHistogram(OutFile, *it, OutputDir);
    }
    
    if (Obs2D_ALL.size() != 0) std::cout << "\nPrinting 2D histograms in the Observables directory: " << std::endl;
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end(); it++) {
        std::string HistName = it->getName();
        if (Histo2D[HistName].GetHistogram()->Integral() > 0.0) {
            std::string fname = OutputDir + "/" + HistName + ".pdf";
            std::vector<double> th;
            th.push_back(it->computeTheoryValue());
            th.push_back(it->computeTheoryValue2());
            Histo2D[HistName].SetGlobalMode(th);
            Print2D(Histo2D[HistName], fname.c_str());
            std::cout << " " + HistName + ".pdf" << std::endl;
            if(OutFile.compare("") == 0) throw std::runtime_error("\nMonteCarloEngine::PrintHistogram ERROR: No root file specified for writing histograms.");
            TDirectory * dir = gDirectory;
            GetOutputFile()->cd();
            Histo2D[HistName].GetHistogram()->Write();
            gDirectory = dir;
            CheckHistogram(*Histo2D[HistName].GetHistogram(), HistName);
        } else HistoLog << "WARNING: The histogram of " << HistName << " is empty!" << std::endl;
    }
    
    std::cout << "\nPrinting LogLikelihood histogram in the Observables directory: " << std::endl;
    if (Histo1D["LogLikelihood"].GetHistogram()->Integral() > 0.0) {
        std::string fname = OutputDir + "/LogLikelihood.pdf";
        Print1D(Histo1D["LogLikelihood"], fname.c_str());
        std::cout << " LogLikelihood.pdf" << std::endl;
        TDirectory * dir = gDirectory;
        GetOutputFile()->cd();
        Histo1D["LogLikelihood"].GetHistogram()->Write();
        gDirectory = dir;
        CheckHistogram(*Histo1D["LogLikelihood"].GetHistogram(), "LogLikelihood");
    }
    
    if (PrintLoglikelihoodPlots) {
        std::cout << "\nPrinting LogLikelihood vs. parameter 2D histograms in the Observables directory: " << std::endl;
        for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it != ModPars.end(); it++) {
            if (it->IsFixed()) continue;
            if (std::find(unknownParameters.begin(), unknownParameters.end(), it->getname()) != unknownParameters.end()) continue;
            std::string HistName = it->getname() + "_vs_LogLikelihood";
            if (Histo2D[HistName].GetHistogram()->Integral() > 0.0) {
                std::string fname = OutputDir + "/LogLikelihoodPlots/" + HistName + ".pdf";
                Print2D(Histo2D[HistName], fname.c_str());
                std::cout << " " + HistName + ".pdf" << std::endl;
                if (OutFile.compare("") == 0) throw std::runtime_error("\nMonteCarloEngine::PrintHistogram ERROR: No root file specified for writing histograms.");
                TDirectory * dir = gDirectory;
                GetOutputFile()->cd();
                Histo2D[HistName].GetHistogram()->Write();
                gDirectory = dir;
                CheckHistogram(*Histo2D[HistName].GetHistogram(), HistName);
            } else HistoLog << "WARNING: The histogram of " << HistName << " is empty!" << std::endl;
        }
    }
    if (noLegend) gStyle->SetOptStat("emrn");
}

void MonteCarloEngine::AddChains() {
    if (fMCMCFlagWriteChainToFile) InitializeMarkovChainTree();
    TDirectory* dir = gDirectory;
    GetOutputFile()->cd();
    
    hMCMCObservableTree = new TTree(TString::Format("%s_Observables", GetSafeName().data()), TString::Format("%s_Observables", GetSafeName().data()));
    hMCMCObservableTree->Branch("Chain", &fMCMCTree_Chain, "chain/i");
    hMCMCObservableTree->Branch("Iteration", &fMCMCCurrentIteration, "iteration/i");
    if (WriteLogLikelihoodChain) {
        hMCMCObservableTree->Branch("LogLikelihood", &hMCMCLogLikelihood, "loglikelihood/D");
        hMCMCObservableTree->Branch("LogProbability", &hMCMCLogProbability, "logprobability/D");
        hMCMCObservableTree->Branch("LogPriorProbability", &hMCMCLogPriorProbability, "logpriorprobability/D");
    }
    if (fMCMCFlagWriteChainToFile) {
        //compute size of observables to be written, including Observables in Obs_ALL, Obs2D_ALL and CGO
        int ObsNum = Obs_ALL.size() + Obs2D_ALL.size() * 2;
        for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++) {
            ObsNum += it->getObs().size();
        }
        hMCMCObservables.assign(fMCMCNChains, std::vector<double>(ObsNum, 0.));
        hMCMCTree_Observables.assign(ObsNum, 0.);
        unsigned int kwtotal = kwmax + (WriteMCMCweights ? kwmcmc : 0);
        if (kwtotal > 0) {
            hMCMCObservables_weight.assign(fMCMCNChains, std::vector<double>(kwtotal, 0.));
            hMCMCTree_Observables_weight.assign(kwtotal, 0.);
        }
        int k = 0, kweight = 0;
        for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
            hMCMCObservableTree->Branch(it->getName().data(), &hMCMCTree_Observables[k], (it->getName() + "/D").data());
            hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables%i", k), it->getName().data());
            k++;
            if (!it->isTMCMC() && it->getDistr().compare("weight") == 0) {
                hMCMCObservableTree->Branch((it->getName() + "_weight").data(), &hMCMCTree_Observables_weight[kweight], (it->getName() + "_weight/D").data());
                hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables_weight%i", kweight), (it->getName() + "_weight").data());
                kweight++;
            } else if (WriteMCMCweights && it->isTMCMC()) {
                hMCMCObservableTree->Branch((it->getName() + "_weight").data(), &hMCMCTree_Observables_weight[kweight], (it->getName() + "_weight/D").data());
                hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables_weight%i", kweight), (it->getName() + "_weight").data());
                kweight++;
            }
        }
        for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end(); it++) {
            hMCMCObservableTree->Branch(it->getThname().data(), &hMCMCTree_Observables[k], (it->getThname() + "/D").data());
            hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables%i", k), it->getThname().data());
            k++;           
            hMCMCObservableTree->Branch(it->getThname2().data(), &hMCMCTree_Observables[k], (it->getThname2() + "/D").data());
            hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables%i", k), it->getThname2().data());
            k++;           
        }
        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
            std::vector<Observable> ObsV(it1->getObs());
            for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it) {
                hMCMCObservableTree->Branch(it->getName().data(), &hMCMCTree_Observables[k], (it->getName() + "/D").data());
                hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables%i", k), it->getName().data());
                k++;
            }
            if (WriteMCMCweights && !it1->isPrediction()) {
                hMCMCObservableTree->Branch((it1->getName() + "_weight").data(), &hMCMCTree_Observables_weight[kweight], (it1->getName() + "_weight/D").data());
                hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables_weight%i", kweight), (it1->getName() + "_weight").data());
                kweight++;
                for (std::vector<Observable>::iterator it = ObsV.begin(); it != ObsV.end(); ++it) {
                    hMCMCObservableTree->Branch((it->getName() + "_weight").data(), &hMCMCTree_Observables_weight[kweight], (it->getName() + "_weight/D").data());
                    hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables_weight%i", kweight), (it->getName() + "_weight").data());
                    kweight++;
                }
            }
        }
    } else if (getchainedObsSize() > 0) {
        hMCMCObservables.assign(fMCMCNChains, std::vector<double>(getchainedObsSize(), 0.));
        hMCMCTree_Observables.assign(getchainedObsSize(), 0.);
        int k = 0;
        for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
            if (it->isWriteChain()) {
                hMCMCObservableTree->Branch(it->getName().data(), &hMCMCTree_Observables[k], (it->getName() + "/D").data());
                hMCMCObservableTree->SetAlias(TString::Format("HEPfit_Observables%i", k), it->getName().data());
                k++;
            }
        }
    }
    hMCMCObservableTree->SetAutoSave(10 * fMCMCNIterationsPreRunCheck);
    hMCMCObservableTree->AutoSave("SelfSave");
    
    hMCMCParameterTree = new TTree(TString::Format("%s_Parameters", GetSafeName().data()), TString::Format("%s_Parameters", GetSafeName().data()));
    hMCMCParameterTree->Branch("Chain", &fMCMCTree_Chain, "chain/i");
    hMCMCParameterTree->Branch("Iteration", &fMCMCCurrentIteration, "iteration/i");
    if (WriteParametersChain) {
        hMCMCParameters.assign(fMCMCNChains, std::vector<double>(DPars.size(), 0.));
        hMCMCTree_Parameters.assign(DPars.size(), 0.);
        int k = 0;
        for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++) {
            hMCMCParameterTree->Branch(it->first.data(), &hMCMCTree_Parameters[k], (it->first + "/D").data());
            hMCMCParameterTree->SetAlias(TString::Format("HEPfit_Parameters%i", k), it->first.data());
            k++;
        }
    }
    hMCMCParameterTree->SetAutoSave(10 * fMCMCNIterationsPreRunCheck);
    hMCMCParameterTree->AutoSave("SelfSave");
    
    gDirectory = dir;
}

void MonteCarloEngine::InChainFillObservablesTree()
{
    if (!hMCMCObservableTree) return;
    for (fMCMCTree_Chain = 0; fMCMCTree_Chain < fMCMCNChains; ++fMCMCTree_Chain) {
        if (getchainedObsSize() > 0 || fMCMCFlagWriteChainToFile) hMCMCTree_Observables = hMCMCObservables[fMCMCTree_Chain];
        if (WriteLogLikelihoodChain) {
            hMCMCLogLikelihood = fMCMCStates.at(fMCMCTree_Chain).log_likelihood;
            hMCMCLogProbability = fMCMCStates.at(fMCMCTree_Chain).log_probability;
            hMCMCLogPriorProbability = fMCMCStates.at(fMCMCTree_Chain).log_prior;
        }
        if ((kwmax > 0 || (WriteMCMCweights && kwmcmc > 0)) && fMCMCFlagWriteChainToFile) hMCMCTree_Observables_weight = hMCMCObservables_weight[fMCMCTree_Chain];
        hMCMCObservableTree->Fill();
    }
}

void MonteCarloEngine::InChainFillParametersTree()
{
    if (!hMCMCParameterTree) return;
    for (fMCMCTree_Chain = 0; fMCMCTree_Chain < fMCMCNChains; ++fMCMCTree_Chain) {
        hMCMCTree_Parameters = hMCMCParameters[fMCMCTree_Chain];
        hMCMCParameterTree->Fill();
    }
}

void MonteCarloEngine::PrintCorrelationMatrixToLaTeX(const std::string filename) {
    std::ofstream out;
    out.open(filename.c_str(), std::ios::out);

    int npar = GetNParameters();

    for (int i = 0; i < npar; ++i)
        out << " & " << GetParameter(i).GetName();
    out << " \\\\" << std::endl;

    for (int i = 0; i < npar; ++i) {
        out << GetParameter(i).GetName() << " & $";
        for (int j = 0; j < npar; ++j) {
            if (i != j) {
                BCH2D* bch2d_temp = new BCH2D(GetMarginalized(GetParameter(i).GetName(), GetParameter(j).GetName()));
                if (bch2d_temp != NULL)
                    out << bch2d_temp->GetHistogram()->GetCorrelationFactor();
                else
                    out << 0.;
                delete bch2d_temp;
                bch2d_temp = NULL;
            } else
                out << 1.;
            if (j == npar - 1) out << "$ \\\\" << std::endl;
            else out << "$ & $";
        }
    }

    out.close();
}

int MonteCarloEngine::getPrecision(double value, double rms, int rmsPrecision) {
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    if (significants == 0) return rmsPrecision + ceil(log10(fabs(value)))-ceil(log10(rms));   
    else return significants;
}

std::string MonteCarloEngine::goodFormat(double value, double rms, int rmsPrecision) {
    
    std::ostringstream out;
    // We want to give value with significant digits up to the precision of the rms. This requires 
    // a bit of manipulation if value is smaller than rms.
    int precision = getPrecision(value, rms, rmsPrecision);
    if (precision > 0)  
        out.precision(precision);
    else if (precision == 0)
    {
        double factor=pow(10,-ceil(log10(fabs(value))));
        value=round(value*factor)/factor;
        out.precision(1);
    }
    else
    {
        value = 0.;
        if(ceil(log10(rms))<=0)
            out << std::fixed << std::setprecision(ceil(log10(rms))+rmsPrecision);
        else
            out.precision(1);
    }
    out << value; 
    return out.str();
}

std::string MonteCarloEngine::computeStatistics() {
    
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables statistics cannot be generated.\n");
    std::streamsize ss_prec = std::cout.precision();
    int rmsPrecision = 2;
    if (significants > 0) rmsPrecision = significants;
    std::ostringstream StatsLog;
    int i = 0;
    StatsLog << "Statistics file for Observables, Binned Observables and Correlated Gaussian Observables.\n" << std::endl;
    if (Obs_ALL.size() > 0) StatsLog << "Observables:\n" << std::endl;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
        StatsLog.precision(ss_prec); /* resets precision*/
        if (it->getObsType().compare("BinnedObservable") == 0) {
            StatsLog << "  (" << ++i << ") Binned Observable \"";
            StatsLog << it->getName() << "[" << it->getTho()->getBinMin() << ", " << it->getTho()->getBinMax() << "]" << "\":";
        } else if (it->getObsType().compare("FunctionObservable") == 0) {
            StatsLog << "  (" << ++i << ") Function Observable \"";
            StatsLog << it->getName() << "[" << it->getTho()->getBinMin() << "]" << "\":";
        } else if (it->getObsType().compare("HiggsObservable") == 0) {
            StatsLog << "  (" << ++i << ") Higgs Observable \"";
            StatsLog << it->getName() << "\":";
        } else {
            StatsLog << "  (" << ++i << ") " << it->getObsType() << " \"";
            StatsLog << it->getName() << "\":";
        }
        StatsLog << std::endl;

        BCH1D bch1d = Histo1D[it->getName()];
        
        if (bch1d.GetHistogram()->Integral() > 0.0) {
            double rms = bch1d.GetHistogram()->GetRMS();
            StatsLog << "      Mean +- sqrt(V):                " << goodFormat(bch1d.GetHistogram()->GetMean(),rms,rmsPrecision)
                    << " +- " << std::setprecision(rmsPrecision)
                    << rms << std::endl
                    << "      (Marginalized) mode:            " << goodFormat(bch1d.GetLocalMode(0),rms) << std::endl;
            std::vector<double> intervals;
            intervals.push_back(0.682689492137);
            intervals.push_back(0.954499736104);
            intervals.push_back(0.997300203937);
            std::vector<BCH1D::BCH1DSmallestInterval> v = bch1d.GetSmallestIntervals(intervals);
            for (unsigned int i = 0; i < v.size(); i++) {
                StatsLog << "      Smallest interval(s) containing at least " << std::setprecision(ss_prec) << v[i].total_mass * 100 << "% and local mode(s):"
                        << std::endl;
                for (unsigned j = 0; j < v[i].intervals.size(); j++) {
                    double interval_xmin = v[i].intervals[j].xmin;
                    double interval_xmax = v[i].intervals[j].xmax;
                    double interval_mode = v[i].intervals[j].mode;
                    double interval_heignt = v[i].intervals[j].relative_height;
                    double interval_relative_mass = v[i].intervals[j].relative_mass;
                    StatsLog << "       (" << goodFormat(interval_xmin, rms) << ", " << goodFormat(interval_xmax, rms) 
                            << ") corresponding to " <<  goodFormat((interval_xmin + interval_xmax)/2.,rms) << " +- " << goodFormat((-interval_xmin + interval_xmax)/2./(i+1),rms) << " (local mode at " << goodFormat(interval_mode, rms) << " with rel. height "
                            << std::setprecision(getPrecision(interval_heignt, ss_prec)) << interval_heignt << "; rel. area " << std::setprecision(getPrecision(interval_relative_mass, ss_prec)) << interval_relative_mass << ")"
                            << std::endl;
                    StatsLog << std::endl;
                }
            }
        } else {
            StatsLog << "\nWARNING: The histogram of " << it->getName() << " is empty! Statistics cannot be generated\n" << std::endl;
        }
    }
    
    if (CGO.size() > 0) StatsLog << "\nCorrelated (Gaussian) Observables:\n" << std::endl;
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
        StatsLog << "\n" << it1->getName() << ":\n" << std::endl;
        i = 0;
        std::vector<Observable> CGObs(it1->getObs());
        for (std::vector<Observable>::iterator it2 = CGObs.begin(); it2 < CGObs.end(); it2++) {
            StatsLog.precision(ss_prec); /* resets precision*/
            if (it2->getObsType().compare("BinnedObservable") == 0) {
                StatsLog << "  (" << ++i << ") Binned Observable \"";
                StatsLog << it2->getName() << "[" << it2->getTho()->getBinMin() << ", " << it2->getTho()->getBinMax() << "]" << "\":";
            } else if (it2->getObsType().compare("FunctionObservable") == 0) {
                StatsLog << "  (" << ++i << ") Function Observable \"";
                StatsLog << it2->getName() << "[" << it2->getTho()->getBinMin() << "]" << "\":";
            } else if (it2->getObsType().compare("HiggsObservable") == 0) {
                StatsLog << "  (" << ++i << ") Higgs Observable \"";
                StatsLog << it2->getName() << "\":";
            } else {
                StatsLog << "  (" << ++i << ") " << it2->getObsType() << " \"";
                StatsLog << it2->getName() << "\":";
            }

            StatsLog << std::endl;
            BCH1D bch1d = Histo1D[it2->getName()];
            if (bch1d.GetHistogram()->Integral() > 0.0) {
                double rms = bch1d.GetHistogram()->GetRMS();
                StatsLog << "      Mean +- sqrt(V):                " << goodFormat(bch1d.GetHistogram()->GetMean(), rms)
                        << " +- " << std::setprecision(rmsPrecision)
                        << rms << std::endl
                        << "      (Marginalized) mode:            " << goodFormat(bch1d.GetLocalMode(0), rms)  << std::endl;

                std::vector<double> intervals;
                intervals.push_back(0.682689492137);
                intervals.push_back(0.954499736104);
                intervals.push_back(0.997300203937);

                std::vector<BCH1D::BCH1DSmallestInterval> v = bch1d.GetSmallestIntervals(intervals);
                for (unsigned int i = 0; i < v.size(); i++) {
                    StatsLog << "      Smallest interval(s) containing at least " << std::setprecision(ss_prec) << v[i].total_mass * 100 << "% and local mode(s):" << std::endl;
                    for (unsigned j = 0; j < v[i].intervals.size(); j++) {
                        double interval_xmin = v[i].intervals[j].xmin;
                        double interval_xmax = v[i].intervals[j].xmax;
                        double interval_mode = v[i].intervals[j].mode;
                        double interval_heignt = v[i].intervals[j].relative_height;
                        double interval_relative_mass = v[i].intervals[j].relative_mass;
                        StatsLog << "       (" <<goodFormat(interval_xmin, rms) << ", " << goodFormat(interval_xmax, rms)
                                << ") (local mode at " << goodFormat(interval_mode, rms) << " with rel. height "
                                << std::setprecision(getPrecision(interval_heignt, ss_prec)) << interval_heignt << "; rel. area " << std::setprecision(getPrecision(interval_relative_mass, ss_prec)) << interval_relative_mass << ")"
                                << std::endl;
                        StatsLog << std::endl;
                    }
                }
            } else {
                StatsLog << "\nWARNING: The histogram of " << it2->getName() << " is empty! Statistics cannot be generated\n" << std::endl;
            }
        }
        if (it1->isPrediction()) {
            int size = it1->getObs().size();
            CorrelationMap[it1->getName()]->MakePrincipals();
            //CorrelationMap[it1->getName()]->Print();
            TMatrixD * corr = const_cast<TMatrixD*>(CorrelationMap[it1->getName()]->GetCovarianceMatrix()); // This returns the normalized correlation matrix, i.e. the correlation matrix
            TVectorD * mean = const_cast<TVectorD*>(CorrelationMap[it1->getName()]->GetMeanValues()); // This returns the vector of mean values.
            TVectorD * sigma = const_cast<TVectorD*>(CorrelationMap[it1->getName()]->GetSigmas()); // This returns the vector of standard deviations.
            *corr *= (double)size; // Get rid of the normalization which is just the size of the matrix.
            gslpp::matrix<double> inverseCovariance(size, size);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    inverseCovariance(i, j) = (*corr)(i, j) * (*sigma)(i) * (*sigma)(j);
                    inverseCovariance(j, i) = inverseCovariance(i, j);
                }
            }
            bool SingularCovariance = inverseCovariance.isSingular();
            if (!SingularCovariance) inverseCovariance = inverseCovariance.inverse(); // Invert to finally produce the inverse covariance (the name is misleading).
            StatsLog << "\nThe correlation matrix for " << it1->getName() << " is given by the " << size << "x"<< size << " matrix:\n" << std::endl;

            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(4) << "" << " | ";
                else StatsLog << std::setw(6) << i << std::setw(6) << "     |";
            }
            StatsLog << std::endl;
            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(8) << "--------";
                else StatsLog << std::setw(12) << "------------";
            }
            StatsLog << std::endl;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size + 1; j++) {
                    if (j == 0) StatsLog << std::setw(4) << i+1 << " |";
                    else StatsLog << std::setprecision(5) << std::setw(12) << (*corr)(i, j - 1);
                }
            StatsLog << std::endl;
            }            
            StatsLog << std::endl;
            
            StatsLog << " The corresponding means and sqrt(V) in a list:\n" << std::endl;
            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(4) << "Mean" << "|";
                else StatsLog << std::setprecision(5) << std::setw(12) << (*mean)(i - 1);
            }
            StatsLog << std::endl;
            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(4) << "sqrt(V)" << "|";
                else StatsLog << std::setprecision(5) << std::setw(12) << (*sigma)(i - 1);
            }
            StatsLog << std::endl;
            
            StatsLog << std::endl;
            if (!SingularCovariance) {
                StatsLog << " The inverse of the square root of the diagonal elements of the inverse covariance matrix are:\n" << std::endl;
                for (int i = 0; i < size + 1; i++) {
                    if (i == 0) StatsLog << std::setw(4) << "sigma" << "|";
                    else StatsLog << std::setprecision(5) << std::setw(12) << 1. / sqrt(inverseCovariance(i - 1, i - 1));
                }
            } else StatsLog << " The covariance matrix cannot be inverted.\n" << std::endl;
            StatsLog << std::endl;
            
            StatsLog << "\nThe correlation matrix for " << it1->getName() << " in Latex form:\n" << std::endl;

            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << " " << " & ";
                else if (i < size) StatsLog << "$" << it1->getObs(i-1).getLabel() << "$" << " & ";
                else StatsLog << it1->getObs(i-1).getLabel() << " \\\\ \\hline";
            }
            StatsLog << std::endl;

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size + 1; j++) {
                    if (j == 0) StatsLog << "$" << it1->getObs(i).getLabel() << "$ & ";
                    else if (j < size) StatsLog << std::setprecision(5) << (*corr)(i, j - 1) << " & ";
                    else StatsLog << std::setprecision(5) << (*corr)(i, j - 1) << " \\\\ ";
                }
            StatsLog << std::endl;
            }
            StatsLog << "\\hline" << std::endl;

            
            TMatrixD * EigVec = const_cast<TMatrixD*>(CorrelationMap[it1->getName()]->GetEigenVectors()); // This returns a matrix with the eigenvectors of the PCA.
            TVectorD * EigVal = const_cast<TVectorD*>(CorrelationMap[it1->getName()]->GetEigenValues()); // This returns the vector of eigenvalues.

            StatsLog << "\nThe matrix of the PCA eigenvectors (columns) for " << it1->getName() << " is given by the " << size << "x"<< size << " matrix:\n" << std::endl;

            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(4) << "" << " | ";
                else StatsLog << std::setw(6) << i << std::setw(6) << "     |";
            }
            StatsLog << std::endl;
            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(8) << "--------";
                else StatsLog << std::setw(12) << "------------";
            }
            StatsLog << std::endl;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size + 1; j++) {
                    if (j == 0) StatsLog << std::setw(4) << i+1 << " |";
                    else StatsLog << std::setprecision(5) << std::setw(12) << (*EigVec)(i, j - 1);
                }
            StatsLog << std::endl;
            }
            StatsLog << std::endl;
            
            StatsLog << " The corresponding PCA eigenvalues are:\n" << std::endl;
            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(4) << "Eigenvalues" << "|";
                else StatsLog << std::setprecision(5) << std::setw(12) << (*EigVal)(i - 1);
            }
            StatsLog << std::endl;
        }
    }

    setDParsFromParameters(mode,DPars);
    Mod->Update(DPars);
    
    // values for controlling format
    int par_width = 0;
    int obs_width = 0;
    int value_width = 15;
    for (std::map<std::string,double>::iterator it = DPars.begin(); it != DPars.end(); it++) par_width = std::max(par_width, (int)(it->first).size());
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) obs_width = std::max(obs_width, (int)(it->getName()).size());
    par_width = par_width + 1;
    obs_width = obs_width + 1;
    const std::string sep = " |";
    const std::string par_line = sep + std::string(par_width + value_width + sep.size() * 2 - 1, '-') + '|';
    const std::string obs_line = sep + std::string(obs_width + value_width + sep.size() * 2 - 1, '-') + '|';
    StatsLog << std::setprecision(5);
    
    StatsLog << "\n*** Statistical details using global mode ***\n" << std::endl;

    StatsLog << "\nValue of the parameters at the global mode:" << std::endl;
    StatsLog << std::endl;
    StatsLog << par_line << '\n' << sep
                 << std::left << std::setw(par_width) << "parameter" << sep << std::right << std::setw(value_width) << "value at mode" << sep << '\n' << par_line << '\n';

    for (std::map<std::string,double>::iterator it = DPars.begin(); it != DPars.end(); it++)
        StatsLog << sep << std::left << std::setw(par_width) << it->first << sep << std::right << std::setw(value_width) << it->second << sep << '\n';
    
    StatsLog << par_line << '\n';
    StatsLog << std::endl;
    
    StatsLog << "\nValue of the observables at the global mode:" << std::endl;
    StatsLog << std::endl;
    StatsLog << obs_line << '\n' << sep
                 << std::left << std::setw(obs_width) << "observable" << sep << std::right << std::setw(value_width) << "value at mode" << sep << '\n' << obs_line << '\n';

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++)
        StatsLog << sep << std::left << std::setw(obs_width) << it->getName() << sep << std::right << std::setw(value_width) << it->computeTheoryValue() << sep << '\n';
    
    StatsLog << obs_line << '\n';
    StatsLog << std::endl;
    
    StatsLog << "LogProbability at mode: " << LogLikelihood(mode) + LogAPrioriProbability(mode) << std::endl;
    StatsLog << "LogLikelihood at mode: " << LogLikelihood(mode) << std::endl;
    StatsLog << "LogAPrioriProbability at mode: " << LogAPrioriProbability(mode) << "\n\n" << std::endl;
    
    double llika = Histo1D["LogLikelihood"].GetHistogram()->GetMean();
    StatsLog << "LogLikelihood mean value: " << llika << std::endl;
    double llikv = Histo1D["LogLikelihood"].GetHistogram()->GetRMS();
    llikv *= llikv;
    StatsLog << "LogLikelihood variance: " << llikv << std::endl;
    double dbar = -2.*llika; //Wikipedia notation... 
    double pd = 2.*llikv; //Wikipedia notation...
    StatsLog << "IC value: " << dbar + 2.*pd << std::endl; 
    StatsLog << "DIC value: " << dbar + pd << std::endl; 
    StatsLog << std::endl;
    StatsLog << std::endl;
    
    
    //For testing purposes:
    const BCEngineMCMC::Statistics& st = GetStatistics();
    //get mean value of parameters from BAT
    std::vector<double> parmeans = st.mean;
    
    setDParsFromParameters(parmeans,DPars);
    Mod->Update(DPars);
    
    StatsLog << "*** Statistical details using mean values of parameters ***\n" << std::endl;
    
    StatsLog << "\nMean value of the parameters:" << std::endl;
    StatsLog << std::endl;
    StatsLog << par_line << '\n' << sep
                 << std::left << std::setw(par_width) << "parameter" << sep << std::right << std::setw(value_width) << "mean value" << sep << '\n' << par_line << '\n';

    for (std::map<std::string,double>::iterator it = DPars.begin(); it != DPars.end(); it++)
        StatsLog << sep << std::left << std::setw(par_width) << it->first << sep << std::right << std::setw(value_width) << it->second << sep << '\n';
    
    StatsLog << par_line << '\n';
    StatsLog << std::endl;
    
    StatsLog << "Mean of LogProbability: " << st.probability_mean << std::endl; 
    StatsLog << "Variance of LogProbability: " << st.probability_variance << std::endl; 
    StatsLog << "LogProbability at mode: " << st.probability_at_mode << std::endl; 
    StatsLog << std::endl;
    
    double llonmean = LogLikelihood(parmeans);
    StatsLog << "LogLikelihood on mean value of parameters: " << llonmean << std::endl;
    StatsLog << "pD computed using variance: " << pd << std::endl; 
    pd = 2.*llonmean-2.*llika;
    StatsLog << "pD computed using 2LL(thetabar) - 2LLbar: " << pd << std::endl; 
    StatsLog << "IC value computed from BAT with alternate pD definition: " << dbar + 2.*pd << std::endl; 
    StatsLog << "DIC value computed from BAT with alternate pD definition: " << dbar + pd << std::endl;
    StatsLog << std::endl;
    StatsLog << std::endl;
   
    
    setDParsFromParameters(par_at_LL_max,DPars);
    Mod->Update(DPars);
    
    StatsLog << "*** Statistical details using parameter values at maximum LogLikelihood ***\n" << std::endl;
    
    StatsLog << "\nValue of the parameters at maximum LogLikelihood:" << std::endl;
    StatsLog << std::endl;
    StatsLog << par_line << '\n' << sep
                 << std::left << std::setw(par_width) << "parameter" << sep << std::right << std::setw(value_width) << "value at max." << sep << '\n' << par_line << '\n';

    for (std::map<std::string,double>::iterator it = DPars.begin(); it != DPars.end(); it++)
        StatsLog << sep << std::left << std::setw(par_width) << it->first << sep << std::right << std::setw(value_width) << it->second << sep << '\n';
    
    StatsLog << par_line << '\n';
    StatsLog << std::endl;
    
    StatsLog << "\nValue of the observables at the maximum LogLikelihood:" << std::endl;
    StatsLog << std::endl;
    StatsLog << obs_line << '\n' << sep
                 << std::left << std::setw(obs_width) << "observable" << sep << std::right << std::setw(value_width) << "value at max." << sep << '\n' << obs_line << '\n';

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++)
        StatsLog << sep << std::left << std::setw(obs_width) << it->getName() << sep << std::right << std::setw(value_width) << it->computeTheoryValue() << sep << '\n';
    
    StatsLog << obs_line << '\n';
    StatsLog << std::endl;
    
    StatsLog << "Maximum LogLikelihood: " << LogLikelihood_max << std::endl;
    
    StatsLog << std::endl;

    return StatsLog.str().c_str();
}

std::string MonteCarloEngine::writePreRunData() 
{
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) {
        throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. PreRun Data cannot be stored.\n");
    }
    std::vector<double> scales(GetScaleFactors().at(0));
    std::ostringstream StatsLog;
    for (unsigned int i = 0; i < mode.size(); i++)
        StatsLog << GetParameter(i).GetName() << " " << mode.at(i) << " " << scales.at(i) << std::endl;
    return StatsLog.str().c_str();
}

std::vector<double> MonteCarloEngine::computeNormalizationMC(int NIterationNormalizationMC) {
    // Number of MC iterations
    SetNIterationsMin(NIterationNormalizationMC);
    SetIntegrationMethod(BCIntegrate::kIntMonteCarlo);
    Integrate();
    std::vector<double> norm;
    norm.clear();
    norm.push_back(GetIntegral());
    norm.push_back(GetError());
    if (norm[0] < 0.) {
        throw std::runtime_error("\n ERROR: Normalization computation cannot be completed since integral is negative.\n");
    }
    
    return norm;
}

/**
 * @brief A duration, in whichever unit reads best.
 * @param[in] seconds the duration
 * @return the duration as a string
 */
static std::string duration(double seconds) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(1);
    if (seconds < 90.) out << seconds << " s";
    else if (seconds < 5400.) out << seconds / 60. << " min";
    else out << seconds / 3600. << " h";
    return out.str();
}

/**
 * @brief Report progress on a computation too long to sit through in silence.
 * @details Written to stdout, which MonteCarlo::Run() redirects into the results
 * file, and flushed so that it can be followed while the job runs.
 * @param[in] line the message
 */
static void progress(const std::string& line) {
    std::cout << line << std::endl;
    fflush(stdout);
}

double MonteCarloEngine::computeNormalizationLME() {
/* PENDING REVIEW FOR USE WITH BAT v1.0. */
    unsigned int Npars = GetNParameters();
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) {
        throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Normalization computation cannot be completed.\n");
    }
    gslpp::matrix<double> Hessian(Npars, Npars, 0.);
    
    for (unsigned int i = 0; i < Npars; i++)
        for (unsigned int j = 0; j < Npars; j++) {
            // calculate Hessian matrix element
            Hessian.assign(i, j, -SecondDerivative(GetParameter(i), GetParameter(j), mode));
        }
    double det_Hessian = Hessian.determinant();

    return exp(Npars / 2. * log(2. * M_PI) + 0.5 * log(1. / det_Hessian) + LogLikelihood(mode) + LogAPrioriProbability(mode));
}

std::vector<double> MonteCarloEngine::computeFunction_h(const std::vector<std::vector<double> >& points) {

    unsigned int npoints = points.size();
    std::vector<double> values(npoints, 0.);

#ifdef _MPI
    if (procnum > 1) {
        unsigned int npars = GetNParameters();
        int buffsize = npars + 1;
        std::vector<double> sendbuff(procnum * buffsize, 0.);
        std::vector<double> recvbuff(buffsize, 0.);
        std::vector<double> gathbuff(procnum, 0.);
        std::vector<double> pars(npars, 0.);

        for (unsigned int first = 0; first < npoints; first += procnum) {
            unsigned int nchunk = std::min((unsigned int) procnum, npoints - first);

            for (unsigned int il = 0; il < nchunk; il++) {
                // The first entry of the array specifies the task to be executed.
                sendbuff[il * buffsize] = 4.; // 4 = evaluate Function_h
                for (unsigned int im = 0; im < npars; im++)
                    sendbuff[il * buffsize + im + 1] = points.at(first + il).at(im);
            }
            for (unsigned int il = nchunk; il < (unsigned int) procnum; il++)
                sendbuff[il * buffsize] = 0.; // 0 = nothing to execute

            MPI_Scatter(&sendbuff[0], buffsize, MPI_DOUBLE, &recvbuff[0], buffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            double fh;
            if (recvbuff[0] == 4.) {
                pars.assign(recvbuff.begin() + 1, recvbuff.end());
                fh = Function_h(pars);
            } else
                fh = log(0.);

            MPI_Gather(&fh, 1, MPI_DOUBLE, &gathbuff[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            for (unsigned int il = 0; il < nchunk; il++)
                values[first + il] = gathbuff[il];
        }
        return values;
    }
#endif

    for (unsigned int i = 0; i < npoints; i++)
        values[i] = Function_h(points[i]);

    return values;
}

double MonteCarloEngine::getDerivativeStep(unsigned int i, const std::vector<double>& point, double relStep) const {

    const BCParameter& par = GetParameter(i);

    /* The width of the prior of the parameter, which is the scale the curvature of
       the posterior is measured against. HEPfit builds the prior out of a Gaussian
       of width errg and a flat part of half width errf, so its variance is known
       analytically and needs no integration. */
    double scale = 0.;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it < ModPars.end(); it++)
        if (it->getname().compare(par.GetName()) == 0) {
            scale = sqrt(it->geterrg() * it->geterrg() + it->geterrf() * it->geterrf() / 3.);
            break;
        }
    if (!std::isfinite(scale) || scale <= 0.) scale = par.GetRangeWidth() / sqrt(12.);

    double step = relStep * scale;

    /* Keep point +- 2 steps inside the range of the parameter: outside it the priors
       of BAT are still finite, so leaving the range is silently allowed and the
       model is evaluated at unphysical values. */
    double room = std::min(par.GetUpperLimit() - point.at(i), point.at(i) - par.GetLowerLimit());
    if (room <= 0.)
        BCLog::OutWarning(("MonteCarloEngine::getDerivativeStep(): " + par.GetName() + " is on or outside the boundary of its range; its derivatives are unreliable.").c_str());
    else if (step > room / 2.)
        step = room / 2.;

    return step;
}

std::vector<double> MonteCarloEngine::calibrateDerivativeSteps(const std::vector<double>& point,
        double relStep, double target) {

    unsigned int Npars = GetNParameters();
    if (point.size() != Npars) {
        throw std::runtime_error("MonteCarloEngine::calibrateDerivativeSteps(): Invalid number of entries in the vector.");
    }

    std::vector<double> step(Npars, 0.);
    std::vector<double> maxstep(Npars, 0.);
    for (unsigned int i = 0; i < Npars; i++) {
        step[i] = getDerivativeStep(i, point, relStep);
        /* the largest step that keeps point +- 2 steps inside the range */
        const BCParameter& par = GetParameter(i);
        double room = std::min(par.GetUpperLimit() - point.at(i), point.at(i) - par.GetLowerLimit());
        maxstep[i] = (room > 0.) ? room / 2. : step[i];
    }

    derivativeStepStatus.assign(Npars, 0); // 0 = still searching
    std::vector<bool> failing(Npars, false);

    std::vector<std::vector<double> > centre(1, point);
    double f0 = computeFunction_h(centre).at(0);
    if (!std::isfinite(f0))
        throw std::runtime_error("MonteCarloEngine::calibrateDerivativeSteps(): the model cannot be evaluated at the point itself.");

    std::time_t tstart = std::time(NULL);
    unsigned long nevaluations = 1;
    progress("Calibrating the finite-difference steps...");

    for (unsigned int iter = 0; iter < HESSIAN_MAXITER; iter++) {

        std::vector<unsigned int> active;
        for (unsigned int i = 0; i < Npars; i++)
            if (derivativeStepStatus[i] == 0) active.push_back(i);
        if (active.empty()) break;

        /* every parameter still being calibrated goes into one batch, so a round
           costs 2 x active points spread over the ranks */
        std::vector<std::vector<double> > points;
        for (unsigned int k = 0; k < active.size(); k++) {
            std::vector<double> p(point);
            p[active[k]] += step[active[k]];
            points.push_back(p);
            p = point;
            p[active[k]] -= step[active[k]];
            points.push_back(p);
        }
        std::vector<double> f(computeFunction_h(points));
        nevaluations += points.size();
        {
            std::ostringstream line;
            line << "  round " << iter + 1 << ": " << active.size()
                 << " parameters still being calibrated, " << points.size() << " evaluations";
            progress(line.str());
        }

        for (unsigned int k = 0; k < active.size(); k++) {
            unsigned int i = active[k];

            if (!std::isfinite(f.at(2 * k)) || !std::isfinite(f.at(2 * k + 1))) {
                /* the model does not survive this step: back off */
                step[i] *= 0.25;
                failing[i] = true;
                continue;
            }
            failing[i] = false;

            double change = fabs(f.at(2 * k) + f.at(2 * k + 1) - 2. * f0);

            if (change == 0.) {
                /* nothing resolved yet: grow, unless the whole range has been used */
                if (step[i] >= maxstep[i]) {
                    step[i] = maxstep[i];
                    derivativeStepStatus[i] = 2; // flat over the range of the parameter
                } else
                    step[i] = std::min(step[i] * 4., maxstep[i]);
                continue;
            }

            if (change >= target / HESSIAN_WINDOW && change <= target * HESSIAN_WINDOW) {
                derivativeStepStatus[i] = 1; // accepted
                continue;
            }

            /* the change grows as the square of the step, so this lands on the
               target in one round whenever the posterior is quadratic */
            double factor = sqrt(target / change);
            if (factor > 10.) factor = 10.;
            if (factor < 0.1) factor = 0.1;
            step[i] *= factor;
            if (step[i] > maxstep[i]) step[i] = maxstep[i];
        }
    }

    unsigned int naccept = 0, nflat = 0, nfail = 0;
    std::ostringstream flatnames, failnames;
    for (unsigned int i = 0; i < Npars; i++) {
        if (derivativeStepStatus[i] == 0)
            derivativeStepStatus[i] = failing[i] ? 3 : 1;
        if (derivativeStepStatus[i] == 1) naccept++;
        else if (derivativeStepStatus[i] == 2) {
            nflat++;
            flatnames << " " << GetParameter(i).GetName();
        } else {
            nfail++;
            failnames << " " << GetParameter(i).GetName();
        }
    }

    {
        std::ostringstream line;
        line << "  calibration complete: " << naccept << " accepted, " << nflat << " flat, "
             << nfail << " unusable (" << nevaluations << " evaluations, "
             << duration(std::difftime(std::time(NULL), tstart)) << ")";
        progress(line.str());
    }

    if (nflat > 0) {
        std::ostringstream message;
        message << "MonteCarloEngine::calibrateDerivativeSteps(): the log posterior does not change "
                << "anywhere inside the range of " << nflat << " parameters, whose curvature is therefore zero:"
                << flatnames.str();
        BCLog::OutWarning(message.str().c_str());
    }
    if (nfail > 0) {
        std::ostringstream message;
        message << "MonteCarloEngine::calibrateDerivativeSteps(): the model could not be evaluated at any step "
                << "along " << nfail << " parameters, whose derivatives are unusable:" << failnames.str();
        BCLog::OutWarning(message.str().c_str());
    }

    return step;
}

gslpp::matrix<double> MonteCarloEngine::computeHessian(const std::vector<double>& point, double relStep, bool adaptive,
        double target) {

    unsigned int Npars = GetNParameters();
    if (point.size() != Npars) {
        throw std::runtime_error("MonteCarloEngine::computeHessian(): Invalid number of entries in the vector.");
    }

    gslpp::matrix<double> Hessian(Npars, Npars, 0.);

    std::time_t tstart = std::time(NULL);
    unsigned long long ntotal = 2ULL * Npars * Npars + 2ULL * Npars + 1ULL;
    unsigned long long ndone = 0;
    {
        std::ostringstream line;
        line << "Computing the Hessian: " << Npars << " parameters, " << ntotal
             << " evaluations of the model";
        progress(line.str());
    }

    std::vector<double> step(Npars, 0.);
    if (adaptive)
        step = calibrateDerivativeSteps(point, relStep, target);
    else {
        for (unsigned int i = 0; i < Npars; i++)
            step[i] = getDerivativeStep(i, point, relStep);
        derivativeStepStatus.assign(Npars, 1);
    }
    derivativeSteps = step;

    /* the stencil begins here: the calibration is a one-off, and folding it into
       the rate would make the first estimates of the time left far too pessimistic */
    std::time_t tstencil = std::time(NULL);

    unsigned int nonfinite = 0;
    std::vector<bool> badpar(Npars, false);

    /* The centre and the points at +- one and +- two steps along each parameter.
       The doubled step gives a second estimate of each diagonal element, whose
       disagreement with the first measures how far the step has sunk into the
       numerical noise of the model. */
    static const double axis[4] = {1., -1., 2., -2.};
    std::vector<std::vector<double> > points;
    points.push_back(point);
    for (unsigned int i = 0; i < Npars; i++)
        for (unsigned int k = 0; k < 4; k++) {
            std::vector<double> p(point);
            p[i] += axis[k] * step[i];
            points.push_back(p);
        }
    std::vector<double> f(computeFunction_h(points));
    ndone += points.size();
    {
        std::ostringstream line;
        line << "  diagonal complete: " << ndone << " evaluations, "
             << duration(std::difftime(std::time(NULL), tstart));
        progress(line.str());
    }

    double worstDiscrepancy = 0.;
    std::string worstParameter;
    for (unsigned int i = 0; i < Npars; i++) {
        double d2 = (f.at(4 * i + 1) - 2. * f.at(0) + f.at(4 * i + 2)) / step[i] / step[i];
        double d2double = (f.at(4 * i + 3) - 2. * f.at(0) + f.at(4 * i + 4)) / 4. / step[i] / step[i];
        if (!std::isfinite(d2)) {
            nonfinite++;
            badpar[i] = true;
            d2 = 0.;
        } else if (std::isfinite(d2double) && (fabs(d2) > 0. || fabs(d2double) > 0.)) {
            double discrepancy = fabs(d2 - d2double) / std::max(fabs(d2), fabs(d2double));
            if (discrepancy > worstDiscrepancy) {
                worstDiscrepancy = discrepancy;
                worstParameter = GetParameter(i).GetName();
            }
        }
        Hessian.assign(i, i, d2);
    }

    /* The mixed derivatives, four evaluations per pair of parameters, sent out one
       row at a time so that the buffer stays small and every rank stays busy. */
    static const double si[4] = {1., 1., -1., -1.};
    static const double sj[4] = {1., -1., 1., -1.};
    for (unsigned int i = 0; i < Npars; i++) {
        std::vector<std::vector<double> > row;
        for (unsigned int j = i + 1; j < Npars; j++)
            for (unsigned int k = 0; k < 4; k++) {
                std::vector<double> p(point);
                p[i] += si[k] * step[i];
                p[j] += sj[k] * step[j];
                row.push_back(p);
            }
        if (row.empty()) continue;

        std::vector<double> fr(computeFunction_h(row));
        ndone += row.size();
        {
            std::time_t now = std::time(NULL);
            std::ostringstream line;
            line << "  row " << i + 1 << "/" << Npars << ": " << ndone << "/" << ntotal
                 << " evaluations (" << (100 * ndone) / ntotal << "%), "
                 << duration(std::difftime(now, tstart));
            /* the remaining rows shrink, so scale what is left by the evaluations
               still to do rather than by the rows still to do */
            double stencil = std::difftime(now, tstencil);
            if (ndone > 0 && ndone < ntotal)
                line << " elapsed, about " << duration(stencil * (ntotal - ndone) / ndone) << " left";
            progress(line.str());
        }
        for (unsigned int j = i + 1; j < Npars; j++) {
            unsigned int b = 4 * (j - i - 1);
            double d2 = (fr.at(b) - fr.at(b + 1) - fr.at(b + 2) + fr.at(b + 3)) / 4. / step[i] / step[j];
            if (!std::isfinite(d2)) {
                nonfinite++;
                badpar[i] = true;
                badpar[j] = true;
                d2 = 0.;
            }
            Hessian.assign(i, j, d2);
            Hessian.assign(j, i, d2);
        }
    }

    if (nonfinite > 0) {
        std::ostringstream message;
        message << "MonteCarloEngine::computeHessian(): " << nonfinite
                << " elements could not be evaluated and were set to zero. The model fails along:";
        for (unsigned int i = 0; i < Npars; i++)
            if (badpar[i]) message << " " << GetParameter(i).GetName();
        BCLog::OutWarning(message.str().c_str());
    }

    {
        std::ostringstream line;
        line << "Hessian complete: " << ndone << " evaluations in "
             << duration(std::difftime(std::time(NULL), tstart));
        progress(line.str());
    }

    if (worstDiscrepancy > 1.e-2) {
        std::ostringstream message;
        message << "MonteCarloEngine::computeHessian(): the curvature along " << worstParameter
                << " changes by " << worstDiscrepancy * 100. << "% when the step is doubled."
                << " Either the posterior is not quadratic over this step, or the step is small enough"
                << " for the numerical noise of the model to show; raise HessianTarget to widen it.";
        BCLog::OutWarning(message.str().c_str());
    }

    return Hessian;
}

/**
 * @brief The eigenvalues, ascending, and eigenvectors of a symmetric matrix.
 * @param[in] a the matrix, row major, N x N
 * @param[in] N the dimension
 * @param[out] eval the eigenvalues
 * @param[out] evec the eigenvectors, evec[k] being the one of eval[k]
 */
static void symmetricEigensystem(const std::vector<double>& a, unsigned int N,
        std::vector<double>& eval, std::vector<std::vector<double> >& evec) {
    gsl_matrix * m = gsl_matrix_alloc(N, N);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++)
            gsl_matrix_set(m, i, j, a[i * N + j]);
    gsl_vector * val = gsl_vector_alloc(N);
    gsl_matrix * vec = gsl_matrix_alloc(N, N);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(N);
    gsl_eigen_symmv(m, val, vec, w);
    gsl_eigen_symmv_sort(val, vec, GSL_EIGEN_SORT_VAL_ASC);
    eval.assign(N, 0.);
    evec.assign(N, std::vector<double>(N, 0.));
    for (unsigned int k = 0; k < N; k++) {
        eval[k] = gsl_vector_get(val, k);
        for (unsigned int i = 0; i < N; i++)
            evec[k][i] = gsl_matrix_get(vec, i, k);
    }
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(vec);
    gsl_vector_free(val);
    gsl_matrix_free(m);
}

gslpp::matrix<double> MonteCarloEngine::refineHessianSoftDirections(const std::vector<double>& point,
        const gslpp::matrix<double>& curvature, double threshold, double target, double noise) {

    unsigned int N = GetNParameters();
    if (point.size() != N || curvature.size_i() != N || curvature.size_j() != N) {
        throw std::runtime_error("MonteCarloEngine::refineHessianSoftDirections(): Invalid dimensions.");
    }

    std::time_t tstart = std::time(NULL);

    /* Unit-diagonal coordinates y, with x = point + d y. A unit step in y changes the
       log posterior by about one half along every parameter, so the eigenvalues of
       the normalised matrix S = d K d compare directions on an equal footing. */
    std::vector<double> d(N, 0.);
    for (unsigned int i = 0; i < N; i++) {
        if (curvature(i, i) > 0.)
            d[i] = 1. / sqrt(curvature(i, i));
        else if (derivativeSteps.size() == N && derivativeSteps[i] > 0.)
            d[i] = derivativeSteps[i];
        else
            d[i] = getDerivativeStep(i, point);
    }
    std::vector<double> S(N * N, 0.);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++)
            S[i * N + j] = 0.5 * (curvature(i, j) + curvature(j, i)) * d[i] * d[j];

    std::vector<double> lambda;
    std::vector<std::vector<double> > v;
    symmetricEigensystem(S, N, lambda, v);

    unsigned int k = 0;
    while (k < N && lambda[k] < threshold) k++;

    {
        std::ostringstream line;
        line << "Refining the Hessian along its soft directions: " << k
             << " eigenvalues of the normalised matrix below " << threshold
             << " (lowest " << lambda[0] << ", highest " << lambda[N - 1] << ")";
        progress(line.str());
    }
    if (k == 0) return curvature;

    /* The change in the log posterior the steps aim for. It sets the balance between
       the two errors of a second difference: the higher derivatives of a posterior that
       is not quadratic, which grow with the step (as the target), and the numerical
       noise deltaf of Function_h(), which enters the normalised curvature as
       deltaf / target and is amplified many times over in the covariance, since the
       soft eigenvalues span orders of magnitude. So the target stays large for all
       directions (HESSIAN_TARGET unless set), and only the directions shown to be
       non-quadratic get a smaller step, below.
       The noise deltaf is measured on the first pass, which the steps of this one
       cannot contaminate: each of its normalised elements carries a noise of about
       deltaf / (2 c), c being the change its steps produced (the median of
       K_ii h_i^2), and its n eigenvalues within the band that noise spreads,
       |lambda| <= |lambda_min|, form the edge of a random matrix, so that
       lambda_min = -2 sqrt(n) deltaf / (2 c) (an upper bound on deltaf if lambda_min
       is positive). It limits how far a step may be shrunk, and decides which of two
       estimates of an element to trust. A later pass, whose input has no such band,
       is given the estimate of the first. */
    bool given = noise > 0.;
    if (!given) {
        std::vector<double> c;
        if (derivativeSteps.size() == N)
            for (unsigned int i = 0; i < N; i++)
                if (curvature(i, i) > 0. && derivativeSteps[i] > 0.)
                    c.push_back(curvature(i, i) * derivativeSteps[i] * derivativeSteps[i]);
        double cfirst = HESSIAN_TARGET;
        if (!c.empty()) {
            std::nth_element(c.begin(), c.begin() + c.size() / 2, c.end());
            cfirst = c[c.size() / 2];
        }
        unsigned int nband = 0;
        for (unsigned int a = 0; a < N; a++)
            if (fabs(lambda[a]) <= fabs(lambda[0])) nband++;
        noise = fabs(lambda[0]) * cfirst / sqrt((double) nband);
    }
    hessianNoise = noise;
    bool automatic = !(target > 0.);
    if (automatic) target = HESSIAN_TARGET;
    {
        std::ostringstream line;
        line << "  noise of the log posterior " << (given || lambda[0] < 0. ? "about " : "at most ") << noise
             << (given ? " (from the first pass)" : "") << ", target " << target
             << (automatic ? " (default)" : " (set in the configuration)");
        progress(line.str());
    }

    /* The direction in parameter space of a unit step along each eigenvector, and
       the largest step that keeps point +- 2 steps inside every range. All of them
       are needed: the soft x stiff elements must be recomputed too, since in the
       eigenbasis of the first pass they are of the order of its noise, and that,
       squared over a stiff eigenvalue, can exceed the soft curvatures. */
    std::vector<double> room(N, 0.);
    for (unsigned int i = 0; i < N; i++) {
        const BCParameter& par = GetParameter(i);
        room[i] = std::max(0., std::min(par.GetUpperLimit() - point.at(i), point.at(i) - par.GetLowerLimit()));
    }
    std::vector<std::vector<double> > u(N, std::vector<double>(N, 0.));
    std::vector<double> maxstep(N, std::numeric_limits<double>::max());
    for (unsigned int a = 0; a < N; a++)
        for (unsigned int i = 0; i < N; i++) {
            u[a][i] = d[i] * v[a][i];
            if (u[a][i] != 0.) maxstep[a] = std::min(maxstep[a], room[i] / 2. / fabs(u[a][i]));
        }

    /* point + sa ta ua + sb tb ub */
    std::vector<double> t(N, 0.);
    struct Shift {
        static std::vector<double> at(const std::vector<double>& x0, const std::vector<double>& ua, double ha,
                const std::vector<double>* ub = NULL, double hb = 0.) {
            std::vector<double> p(x0);
            for (unsigned int i = 0; i < p.size(); i++) {
                p[i] += ha * ua[i];
                if (ub) p[i] += hb * (*ub)[i];
            }
            return p;
        }
    };

    /* Calibrate the step along each direction, as calibrateDerivativeSteps() does
       along each parameter. The first guess comes from the first-pass eigenvalue,
       so the stiff directions are usually accepted at once; a direction that stays
       flat ends at the largest step allowed. An accepted step is then set on the
       target itself, since the change grows as its square: nothing is evaluated at
       it before the stencil, and a change left anywhere in the window would carry up
       to HESSIAN_WINDOW times the noise, or the truncation error, that the target is
       meant to allow. The soft directions need it most: their first guesses come from
       eigenvalues inflated by the noise, so they tend to be accepted low in the window. */
    std::vector<int> status(N, 0); // 0 searching, 1 accepted, 2 flat, 3 model fails
    std::vector<bool> failing(N, false);
    std::vector<std::vector<double> > centre(1, point);
    double f0 = computeFunction_h(centre).at(0);
    if (!std::isfinite(f0))
        throw std::runtime_error("MonteCarloEngine::refineHessianSoftDirections(): the model cannot be evaluated at the point itself.");
    unsigned long nevaluations = 1;
    for (unsigned int a = 0; a < N; a++)
        t[a] = std::min(maxstep[a], sqrt(target / std::max(fabs(lambda[a]), 1.e-12)));

    for (unsigned int iter = 0; iter < HESSIAN_MAXITER; iter++) {
        std::vector<unsigned int> active;
        for (unsigned int a = 0; a < N; a++)
            if (status[a] == 0) active.push_back(a);
        if (active.empty()) break;

        std::vector<std::vector<double> > points;
        for (unsigned int m = 0; m < active.size(); m++) {
            unsigned int a = active[m];
            points.push_back(Shift::at(point, u[a], t[a]));
            points.push_back(Shift::at(point, u[a], -t[a]));
        }
        std::vector<double> f(computeFunction_h(points));
        nevaluations += points.size();
        {
            std::ostringstream line;
            line << "  round " << iter + 1 << ": " << active.size()
                 << " directions still being calibrated, " << points.size() << " evaluations";
            progress(line.str());
        }

        for (unsigned int m = 0; m < active.size(); m++) {
            unsigned int a = active[m];
            if (!std::isfinite(f.at(2 * m)) || !std::isfinite(f.at(2 * m + 1))) {
                t[a] *= 0.25;
                failing[a] = true;
                continue;
            }
            failing[a] = false;
            double change = fabs(f.at(2 * m) + f.at(2 * m + 1) - 2. * f0);
            if (change == 0.) {
                if (t[a] >= maxstep[a]) {
                    t[a] = maxstep[a];
                    status[a] = 2;
                } else
                    t[a] = std::min(t[a] * 4., maxstep[a]);
                continue;
            }
            if (change >= target / HESSIAN_WINDOW && change <= target * HESSIAN_WINDOW) {
                status[a] = 1;
                t[a] = std::min(t[a] * sqrt(target / change), maxstep[a]);
                continue;
            }
            if (change < target / HESSIAN_WINDOW && t[a] >= maxstep[a]) {
                status[a] = 2; // as far as the ranges allow: take what it gives
                continue;
            }
            double factor = sqrt(target / change);
            if (factor > 10.) factor = 10.;
            if (factor < 0.1) factor = 0.1;
            t[a] = std::min(t[a] * factor, maxstep[a]);
        }
    }
    for (unsigned int a = 0; a < N; a++)
        if (status[a] == 0) status[a] = failing[a] ? 3 : 1;

    /* The signs of the four points of a mixed difference, and the points of a
       diagonal one at one and two steps. */
    static const double sa[4] = {1., 1., -1., -1.};
    static const double sb[4] = {1., -1., 1., -1.};
    static const double axis[4] = {1., -1., 2., -2.};

    /* Every element involving a soft direction is computed with the step h and with
       2h, and extrapolated as (4 M(h) - M(2h)) / 3, which cancels the h^2 error of a
       posterior that is not quadratic over the step (Richardson) and changes the
       noise very little. It matters: the soft eigenvalues span orders of magnitude, so
       the covariance amplifies errors of these elements many times over, and even a
       per-cent truncation error of an element that passes any sensible check is too
       much. The relative difference between M(h) and M(2h) of each soft diagonal is
       kept as the measure of how far the direction is from quadratic.
       A mixed element at the steps of its two directions reuses their points on the
       axes, -(f(+a+b) + f(-a-b) - f(+a) - f(-a) - f(+b) - f(-b) + 2 f0) / (2 ha hb),
       which is exact for a quadratic posterior and needs two new points per step
       instead of four; a pair at any other steps (see below) takes its own four,
       -(f(+a+b) - f(+a-b) - f(-a+b) + f(-a-b)) / (4 ha hb).
       Points of one stencil, in order: for each direction flagged in diag, +- h and
       +- 2h; for each pair with a soft index flagged in pair, its four points at h
       and its four at 2h. */
    struct Stencil {
        static void points(const std::vector<double>& x0, const std::vector<std::vector<double> >& u,
                unsigned int k, unsigned int N, const std::vector<double>& hd,
                const std::vector<bool>& diag, const std::vector<std::vector<double> >& hp,
                const std::vector<bool>& pair, std::vector<std::vector<double> >& out) {
            for (unsigned int a = 0; a < diag.size(); a++)
                if (diag[a])
                    for (unsigned int s = 0; s < 4; s++)
                        out.push_back(Shift::at(x0, u[a], axis[s] * hd[a]));
            unsigned int p = 0;
            for (unsigned int a = 0; a < k; a++)
                for (unsigned int b = a + 1; b < N; b++, p++)
                    if (pair[p])
                        for (unsigned int m = 1; m <= 2; m++)
                            for (unsigned int s = 0; s < 4; s++)
                                out.push_back(Shift::at(x0, u[a], m * sa[s] * hp[p][0], &u[b], m * sb[s] * hp[p][1]));
        }
    };
    unsigned int npairs = 0;
    for (unsigned int a = 0; a < k; a++) npairs += N - a - 1;
    const std::vector<bool> nopair(npairs, false);

    /* Second differences: the diagonal at h and 2h from the values of one stencil,
       the pairs from their own four points or from two and those on the axes. */
    struct Differences {
        static double diag(const std::vector<double>& f, unsigned int i, double f0, double h, double m) {
            return -(f.at(i + (m == 1. ? 0 : 2)) + f.at(i + (m == 1. ? 1 : 3)) - 2. * f0) / (m * h) / (m * h);
        }
        static double mixed(const std::vector<double>& f, unsigned int i, double ha, double hb) {
            return -(f.at(i) - f.at(i + 1) - f.at(i + 2) + f.at(i + 3)) / 4. / ha / hb;
        }
        static double mixedaxes(double fpp, double fmm, double fpa, double fma, double fpb, double fmb,
                double f0, double ha, double hb) {
            return -(fpp + fmm - fpa - fma - fpb - fmb + 2. * f0) / 2. / ha / hb;
        }
    };

    /* M = minus the second derivatives of Function_h along the eigenvectors, in y
       units: the first-pass eigenvalues, with every element involving a soft
       direction recomputed; a failed element keeps its first-pass value (zero off the
       diagonal). */
    std::vector<double> M(N * N, 0.);
    for (unsigned int a = 0; a < N; a++) M[a * N + a] = lambda[a];
    std::vector<double> d2single(k, 0.), d2double(k, 0.); // the soft diagonal at h and 2h
    std::vector<double> discrepancy(k, 0.);                // their relative difference
    unsigned int nonfinite = 0;

    /* The points on the axes of every direction come first. Those of the soft ones
       give the diagonal, which alone decides which directions are retried with a
       smaller step (below), so that every pair is computed only once, with its final
       steps; all of them serve the mixed elements. */
    std::vector<std::vector<double> > hp(npairs, std::vector<double>(2, 0.));
    std::vector<std::vector<double> > points;
    Stencil::points(point, u, k, N, t, std::vector<bool>(N, true), hp, nopair, points);
    {
        std::ostringstream line;
        line << "  axes (steps h and 2h): " << points.size() << " evaluations";
        progress(line.str());
    }
    const std::vector<double> faxis(computeFunction_h(points)); // at 4 a: +h, -h, +2h, -2h
    nevaluations += points.size();
    for (unsigned int a = 0; a < k; a++) {
        double m1 = Differences::diag(faxis, 4 * a, f0, t[a], 1.);
        double m2 = Differences::diag(faxis, 4 * a, f0, t[a], 2.);
        d2single[a] = m1;
        d2double[a] = m2;
        if (!std::isfinite(m1) || !std::isfinite(m2)) {
            nonfinite++;
            continue;
        }
        if (fabs(m1) > 0. || fabs(m2) > 0.)
            discrepancy[a] = fabs(m1 - m2) / std::max(fabs(m1), fabs(m2));
        M[a * N + a] = (4. * m1 - m2) / 3.;
    }

    /* One retry, never more, for the soft directions whose curvature still changes by
       more than HESSIAN_SOFTTOL between h and 2h: there the terms beyond h^2, which the
       extrapolation leaves, need not be small. For a posterior that is not quadratic
       that change grows as the square of the step, so the step is shrunk by
       s = sqrt(HESSIAN_SOFTTOL / (2 x change)), within [1/16, 1/2], but not below the
       step at which the noise would reach 1 / HESSIAN_SOFTNOISEFACTOR of the change in
       the log posterior. The extrapolation from s h and 2 s h replaces the first only
       if the two differ by more than three times its noise (3.13 deltaf / (s h)^2):
       otherwise the first was right, as it is whenever the terms beyond h^2 vanish (a
       log posterior quartic in the parameters, such as a Gaussian likelihood of
       predictions quadratic in them), and it has less noise. The direction is then
       confirmed rather than retried. */
    std::vector<bool> retried(k, false), kept(k, false), confirmed(k, false);
    std::vector<double> th(t), shrinkf(k, 1.);
    unsigned int nretry = 0;
    double smin = std::min(0.5, std::max(1. / 16., sqrt(HESSIAN_SOFTNOISEFACTOR * noise / target)));
    for (unsigned int a = 0; a < k; a++)
        if (discrepancy[a] > HESSIAN_SOFTTOL && smin < 0.5) {
            double shrink = sqrt(HESSIAN_SOFTTOL / (2. * discrepancy[a]));
            shrink = std::min(0.5, std::max(smin, shrink));
            retried[a] = true;
            shrinkf[a] = shrink;
            th[a] = shrink * t[a];
            nretry++;
        }
    if (nretry > 0) {
        std::vector<std::vector<double> > rpoints;
        Stencil::points(point, u, k, N, th, retried, hp, nopair, rpoints);
        {
            std::ostringstream line;
            line << "  retry of " << nretry << " soft directions with a smaller step: "
                 << rpoints.size() << " evaluations";
            progress(line.str());
        }
        std::vector<double> fr(computeFunction_h(rpoints));
        nevaluations += rpoints.size();

        unsigned int i = 0;
        for (unsigned int a = 0; a < k; a++)
            if (retried[a]) {
                double m1 = Differences::diag(fr, i, f0, th[a], 1.);
                double m2 = Differences::diag(fr, i, f0, th[a], 2.);
                i += 4;
                if (!std::isfinite(m1) || !std::isfinite(m2) || !(fabs(m1) > 0. || fabs(m2) > 0.)) continue;
                double newdiscrepancy = fabs(m1 - m2) / std::max(fabs(m1), fabs(m2));
                double r = (4. * m1 - m2) / 3.;
                if (fabs(r - M[a * N + a]) > 3. * 3.13 * noise / (th[a] * th[a])) {
                    kept[a] = true;
                    M[a * N + a] = r;
                    d2single[a] = m1;
                    d2double[a] = m2;
                    discrepancy[a] = newdiscrepancy;
                } else
                    confirmed[a] = true;
            }
    }

    /* Every pair with a soft index, once. A direction that is not quadratic generally
       leans on the others too (the eigenvectors of the first pass carry its noise), so
       a pair involving a direction whose smaller step was kept has both of its steps
       shrunk, by the smaller factor of its two directions if both were retried. Each
       direction alone stays inside the ranges of the parameters, but two together need
       not: a pair that would leave them has both steps shrunk until it does not. Both
       kinds of pair take their own four points. */
    std::vector<bool> ownpoints(npairs, false);
    unsigned int nlimited = 0;
    {
        unsigned int p = 0;
        for (unsigned int a = 0; a < k; a++)
            for (unsigned int b = a + 1; b < N; b++, p++) {
                bool shrunk = kept[a] || (b < k && kept[b]);
                double s = shrunk ? std::min(shrinkf[a], b < k ? shrinkf[b] : 1.) : 1.;
                bool limited = false;
                for (unsigned int i = 0; i < N; i++) {
                    double reach = 2. * s * (t[a] * fabs(u[a][i]) + t[b] * fabs(u[b][i]));
                    if (reach > room[i]) {
                        s *= room[i] / reach;
                        limited = true;
                    }
                }
                if (limited) nlimited++;
                ownpoints[p] = shrunk || limited;
                hp[p][0] = s * t[a];
                hp[p][1] = s * t[b];
            }
    }
    points.clear();
    Stencil::points(point, u, k, N, t, std::vector<bool>(), hp, ownpoints, points);
    unsigned int nown = points.size();
    {
        unsigned int p = 0;
        for (unsigned int a = 0; a < k; a++)
            for (unsigned int b = a + 1; b < N; b++, p++)
                if (!ownpoints[p])
                    for (unsigned int m = 1; m <= 2; m++) {
                        points.push_back(Shift::at(point, u[a], m * hp[p][0], &u[b], m * hp[p][1]));
                        points.push_back(Shift::at(point, u[a], -(m * hp[p][0]), &u[b], -(m * hp[p][1])));
                    }
    }
    {
        std::ostringstream line;
        line << "  pairs with a soft direction (steps h and 2h): " << points.size() << " evaluations";
        if (nlimited > 0) line << " (" << nlimited << " pairs with smaller steps, to stay inside the ranges)";
        progress(line.str());
    }
    std::vector<double> f(computeFunction_h(points));
    nevaluations += points.size();
    std::vector<double> pdiscrepancy(npairs, 0.);
    {
        unsigned int p = 0, i = 0, j = nown;
        for (unsigned int a = 0; a < k; a++)
            for (unsigned int b = a + 1; b < N; b++, p++) {
                double m1, m2;
                if (ownpoints[p]) {
                    m1 = Differences::mixed(f, i, hp[p][0], hp[p][1]);
                    m2 = Differences::mixed(f, i + 4, 2. * hp[p][0], 2. * hp[p][1]);
                    i += 8;
                } else {
                    m1 = Differences::mixedaxes(f.at(j), f.at(j + 1), faxis.at(4 * a), faxis.at(4 * a + 1),
                            faxis.at(4 * b), faxis.at(4 * b + 1), f0, hp[p][0], hp[p][1]);
                    m2 = Differences::mixedaxes(f.at(j + 2), f.at(j + 3), faxis.at(4 * a + 2), faxis.at(4 * a + 3),
                            faxis.at(4 * b + 2), faxis.at(4 * b + 3), f0, 2. * hp[p][0], 2. * hp[p][1]);
                    j += 4;
                }
                double m = (4. * m1 - m2) / 3.;
                if (!std::isfinite(m)) {
                    nonfinite++;
                    m = 0.;
                } else
                    pdiscrepancy[p] = fabs(m1 - m2) / sqrt(fabs(M[a * N + a] * M[b * N + b]));
                M[a * N + b] = m;
                M[b * N + a] = m;
            }
    }

    /* A pair can be far from quadratic even where neither of its directions is: in the
       eigenbasis the elements off the diagonal are small, so a term shared by the two
       directions that is not quadratic can dominate their element while it is lost in
       each diagonal. So the change between h and 2h of each pair is measured too, on
       the scale of its element, sqrt(M_aa M_bb), which the noise barely reaches; a pair
       changing by more than HESSIAN_SOFTTOL is recomputed once with its own four
       points and both steps shrunk as that change requires, within the same limits as
       a direction, and the new extrapolation replaces the first only if the two differ
       by more than three times its noise (0.668 deltaf / (s^2 ha hb)). */
    std::vector<unsigned int> rpairs;
    std::vector<double> pshrink;
    unsigned int npairnoisy = 0; // changing too much, but too noisy to be retried
    for (unsigned int p = 0; p < npairs; p++)
        if (pdiscrepancy[p] > HESSIAN_SOFTTOL) {
            if (smin < 0.5) {
                rpairs.push_back(p);
                pshrink.push_back(std::min(0.5, std::max(smin, sqrt(HESSIAN_SOFTTOL / (2. * pdiscrepancy[p])))));
            } else
                npairnoisy++;
        }
    unsigned int nrpairkept = 0, nrpairconfirmed = 0, nrpairopen = 0;
    if (!rpairs.empty()) {
        std::vector<bool> flag(npairs, false);
        std::vector<std::vector<double> > rhp(hp);
        for (unsigned int m = 0; m < rpairs.size(); m++) {
            flag[rpairs[m]] = true;
            rhp[rpairs[m]][0] *= pshrink[m];
            rhp[rpairs[m]][1] *= pshrink[m];
        }
        points.clear();
        Stencil::points(point, u, k, N, t, std::vector<bool>(), rhp, flag, points);
        {
            std::ostringstream line;
            line << "  retry of " << rpairs.size() << " pairs with smaller steps: " << points.size() << " evaluations";
            progress(line.str());
        }
        std::vector<double> fr(computeFunction_h(points));
        nevaluations += points.size();
        unsigned int p = 0, i = 0;
        for (unsigned int a = 0; a < k; a++)
            for (unsigned int b = a + 1; b < N; b++, p++)
                if (flag[p]) {
                    double m1 = Differences::mixed(fr, i, rhp[p][0], rhp[p][1]);
                    double m2 = Differences::mixed(fr, i + 4, 2. * rhp[p][0], 2. * rhp[p][1]);
                    i += 8;
                    double m = (4. * m1 - m2) / 3.;
                    if (!std::isfinite(m)) continue;
                    double newdiscrepancy = fabs(m1 - m2) / sqrt(fabs(M[a * N + a] * M[b * N + b]));
                    if (fabs(m - M[a * N + b]) > 3. * 0.668 * noise / (rhp[p][0] * rhp[p][1])) {
                        M[a * N + b] = m;
                        M[b * N + a] = m;
                        pdiscrepancy[p] = newdiscrepancy;
                        nrpairkept++;
                        if (newdiscrepancy > HESSIAN_SOFTTOL) nrpairopen++;
                    } else
                        nrpairconfirmed++;
                }
    }
    for (unsigned int a = 0; a < k; a++)
        if (kept[a]) t[a] = th[a];

    /* S' = V M V^T, then back to parameter units */
    std::vector<double> W(N * N, 0.); // V M
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int b = 0; b < N; b++) {
            double w = 0.;
            for (unsigned int a = 0; a < N; a++) w += v[a][i] * M[a * N + b];
            W[i * N + b] = w;
        }
    std::vector<double> Sp(N * N, 0.);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++) {
            double w = 0.;
            for (unsigned int b = 0; b < N; b++) w += W[i * N + b] * v[b][j];
            Sp[i * N + j] = w;
        }
    gslpp::matrix<double> refined(N, N, 0.);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++)
            refined.assign(i, j, 0.5 * (Sp[i * N + j] + Sp[j * N + i]) / d[i] / d[j]);

    std::vector<double> lambdap;
    std::vector<std::vector<double> > vp;
    symmetricEigensystem(Sp, N, lambdap, vp);
    unsigned int nneg = 0, nnegp = 0;
    for (unsigned int i = 0; i < N; i++) {
        if (lambda[i] <= 0.) nneg++;
        if (lambdap[i] <= 0.) nnegp++;
    }

    /* Report: every soft direction, its step and its curvature before and after. */
    std::streamsize prec = std::cout.precision();
    std::cout << std::setprecision(6);
    std::cout << std::endl << "Soft directions of the Hessian (normalised to unit diagonal; target "
              << target << "):" << std::endl;
    std::vector<unsigned int> nonquadratic, nonquadraticconfirmed;
    for (unsigned int a = 0; a < k; a++) {
        std::vector<std::pair<double, unsigned int> > comp;
        for (unsigned int i = 0; i < N; i++) comp.push_back(std::make_pair(-fabs(v[a][i]), i));
        std::sort(comp.begin(), comp.end());
        std::cout << "  " << a + 1 << ": first pass " << lambda[a] << ", recomputed " << M[a * N + a]
                  << " (step h " << d2single[a] << ", 2h " << d2double[a] << ", " << 100. * discrepancy[a]
                  << "%), h = " << t[a];
        if (status[a] == 2) std::cout << " (largest step allowed by the ranges)";
        else if (status[a] == 3) std::cout << " (the model could not be evaluated)";
        if (retried[a]) std::cout << (kept[a] ? " [smaller step kept]" : confirmed[a] ? " [confirmed with a smaller step]" : " [smaller step unusable]");
        std::cout << ";";
        for (unsigned int m = 0; m < std::min(4u, N); m++)
            std::cout << " " << GetParameter(comp[m].second).GetName() << ":" << v[a][comp[m].second];
        std::cout << std::endl;
        if (discrepancy[a] > HESSIAN_SOFTTOL) (confirmed[a] ? nonquadraticconfirmed : nonquadratic).push_back(a);
    }
    std::cout << "Normalised eigenvalues <= 0: " << nneg << " before, " << nnegp << " after refinement;"
              << " lowest " << lambda[0] << " before, " << lambdap[0] << " after." << std::endl;
    std::cout << "Soft directions retried with a smaller step: " << nretry << " (kept "
              << std::count(kept.begin(), kept.end(), true) << ", confirmed "
              << std::count(confirmed.begin(), confirmed.end(), true) << "); pairs retried: " << rpairs.size()
              << " (kept " << nrpairkept << ", confirmed " << nrpairconfirmed << ")." << std::endl;
    if (nonquadratic.empty() && nonquadraticconfirmed.empty())
        std::cout << "Every soft curvature changes by less than " << 100. * HESSIAN_SOFTTOL
                  << "% when the step is doubled." << std::endl;
    if (!nonquadraticconfirmed.empty()) {
        std::cout << "Soft directions not quadratic over their step, their curvature at the point confirmed"
                  << " with a smaller step (the posterior is not Gaussian along them over its width):";
        for (unsigned int m = 0; m < nonquadraticconfirmed.size(); m++)
            std::cout << " " << nonquadraticconfirmed[m] + 1 << " (" << 100. * discrepancy[nonquadraticconfirmed[m]] << "%)";
        std::cout << std::endl;
    }
    if (!nonquadratic.empty()) {
        std::cout << "Soft directions still not quadratic over their step, their curvature not confirmed"
                  << " (the covariance along them is unreliable):";
        for (unsigned int m = 0; m < nonquadratic.size(); m++)
            std::cout << " " << nonquadratic[m] + 1 << " (" << 100. * discrepancy[nonquadratic[m]] << "%)";
        std::cout << std::endl;
    }
    std::cout << std::setprecision(prec) << std::endl;

    {
        std::ostringstream line;
        line << "Refinement complete: " << nevaluations << " evaluations in "
             << duration(std::difftime(std::time(NULL), tstart));
        progress(line.str());
    }

    if (nonfinite > 0) {
        std::ostringstream message;
        message << "MonteCarloEngine::refineHessianSoftDirections(): " << nonfinite
                << " elements could not be evaluated and kept their first-pass values (off-diagonal: zero).";
        BCLog::OutWarning(message.str().c_str());
    }
    if (nnegp > 0) {
        std::ostringstream message;
        message << "MonteCarloEngine::refineHessianSoftDirections(): " << nnegp
                << " eigenvalues of the refined matrix are still not positive (lowest " << lambdap[0]
                << " in normalised units). Either the noise of the log posterior (about " << noise
                << ") is too large for curvatures this small, or the posterior is not quadratic"
                << " along them; see the list of soft directions.";
        BCLog::OutWarning(message.str().c_str());
    }
    if (!nonquadratic.empty() || nrpairopen + npairnoisy > 0) {
        std::ostringstream message;
        message << "MonteCarloEngine::refineHessianSoftDirections(): " << nonquadratic.size()
                << " soft directions and " << nrpairopen + npairnoisy << " pairs are not quadratic over their"
                << " steps (curvature changes by more than " << 100. * HESSIAN_SOFTTOL
                << "% when the step is doubled) and could not be confirmed with smaller steps; see the list"
                << " of soft directions.";
        BCLog::OutWarning(message.str().c_str());
    }

    return refined;
}

gslpp::matrix<double> MonteCarloEngine::computeHessianLegacy(const std::vector<double>& point) {

    unsigned int Npars = GetNParameters();
    if (point.size() != Npars) {
        throw std::runtime_error("MonteCarloEngine::computeHessianLegacy(): Invalid number of entries in the vector.");
    }

    gslpp::matrix<double> Hessian(Npars, Npars, 0.);
    for (unsigned int i = 0; i < Npars; i++)
        for (unsigned int j = i; j < Npars; j++) {
            double d2 = SecondDerivative(GetParameter(i), GetParameter(j), point);
            Hessian.assign(i, j, d2);
            Hessian.assign(j, i, d2);
        }

    return Hessian;
}

double MonteCarloEngine::SecondDerivative(const BCParameter& par1, const BCParameter& par2, const std::vector<double>& point) {

    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MonteCarloEngine::SecondDerivative : Invalid number of entries in the vector.");
    }

    // define steps
    const double dy1 = par2.GetRangeWidth() / NSTEPS;
    const double dy2 = dy1 * 2.;
    const double dy3 = dy1 * 3.;

    // define points at which to evaluate
    std::vector<double> y1p = point;
    std::vector<double> y1m = point;
    std::vector<double> y2p = point;
    std::vector<double> y2m = point;
    std::vector<double> y3p = point;
    std::vector<double> y3m = point;

    unsigned idy = GetParameters().Index(par2.GetName());

    y1p[idy] += dy1;
    y1m[idy] -= dy1;
    y2p[idy] += dy2;
    y2m[idy] -= dy2;
    y3p[idy] += dy3;
    y3m[idy] -= dy3;

    const double m1 = (FirstDerivative(par1, y1p) - FirstDerivative(par1, y1m)) / 2. / dy1;
    const double m2 = (FirstDerivative(par1, y2p) - FirstDerivative(par1, y2m)) / 4. / dy1;
    const double m3 = (FirstDerivative(par1, y3p) - FirstDerivative(par1, y3m)) / 6. / dy1;

    return 3. / 2. * m1 - 3. / 5. * m2 + 1. / 10. * m3;
}

double MonteCarloEngine::FirstDerivative(const BCParameter& par, const std::vector<double>& point) {

    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MonteCarloEngine::FirstDerivative : Invalid number of entries in the vector.");
    }

    // define steps
    const double dx1 = par.GetRangeWidth() / NSTEPS;
    const double dx2 = dx1 * 2.;
    const double dx3 = dx1 * 3.;

    // define points at which to evaluate
    std::vector<double> x1p = point;
    std::vector<double> x1m = point;
    std::vector<double> x2p = point;
    std::vector<double> x2m = point;
    std::vector<double> x3p = point;
    std::vector<double> x3m = point;

    unsigned idx = GetParameters().Index(par.GetName());

    x1p[idx] += dx1;
    x1m[idx] -= dx1;
    x2p[idx] += dx2;
    x2m[idx] -= dx2;
    x3p[idx] += dx3;
    x3m[idx] -= dx3;

    const double m1 = (Function_h(x1p) - Function_h(x1m)) / 2. / dx1;
    const double m2 = (Function_h(x2p) - Function_h(x2m)) / 4. / dx1;
    const double m3 = (Function_h(x3p) - Function_h(x3m)) / 6. / dx1;

    return 3. / 2. * m1 - 3. / 5. * m2 + 1. / 10. * m3;
}

double MonteCarloEngine::Function_h(const std::vector<double>& point) {
    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MonteCarloEngine::Function_h : Invalid number of entries in the vector.");
    }
    return LogLikelihood(point) + LogAPrioriProbability(point);
}
