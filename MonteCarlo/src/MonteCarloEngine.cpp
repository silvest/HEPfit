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
#ifdef _MPI
#include <mpi.h>
#endif
#include <TF1.h>
#include <TTree.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <limits>
#include <memory>

MonteCarloEngine::MonteCarloEngine(
        const std::vector<ModelParameter>& ModPars_i,
        boost::ptr_vector<Observable>& Obs_i,
        std::vector<Observable2D>& Obs2D_i,
        std::vector<CorrelatedGaussianObservables>& CGO_i,
        std::vector<CorrelatedGaussianParameters>& CGP_i)
: BCModel(""), ModPars(ModPars_i), CGP(CGP_i), Obs_ALL(Obs_i), Obs2D_ALL(Obs2D_i),
  CGO(CGO_i), NumOfUsedEvents(0), NumOfDiscardedEvents(0) {
    obval = NULL;
    obweight = NULL;
    Mod = NULL;
    cindex = 0;
    printLogo = false;
    nSmooth = 0;
    histogram2Dtype = 1001;
    noLegend = true;
    PrintLoglikelihoodPlots = false;
    WriteLogLikelihoodChain = false;
    WriteParametersChain = false;
    alpha2D = 1.;
    kchainedObs = 0;
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
    int k = 0, kweight = 0;

    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            k++;
            if (it->getDistr().compare("noweight") != 0) kweight++;
            if (it->isWriteChain()) kchainedObs++;
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
        }
    }
    kmax = k;
    kwmax = kweight;

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
    if (!std::isfinite(logprob)) {
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
                int ko = 0, kweight = 0, k_all = 0, k_cObs = 0;
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
                            obweight[index_chain[il] * kwmax + kweight] = weight;
                            if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[index_chain[il]][kweight] = weight;
                            kweight++;
                        }
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
        int k = 0, kweight = 0, k_all = 0, k_cObs = 0;
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
                    obweight[i * kwmax + kweight] = weight;
                    if (fMCMCFlagWriteChainToFile) hMCMCObservables_weight[i][kweight] = weight;
                    kweight++;
                }
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
        if (kwmax > 0) {
            hMCMCObservables_weight.assign(fMCMCNChains, std::vector<double>(kwmax, 0.));
            hMCMCTree_Observables_weight.assign(kwmax, 0.);
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
        if (kwmax > 0 && fMCMCFlagWriteChainToFile) hMCMCTree_Observables_weight = hMCMCObservables_weight[fMCMCTree_Chain];
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

int MonteCarloEngine::getPrecision(double value, double rms) {
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    if (significants == 0) return 2 + ceil(log10(fabs(value)))-ceil(log10(rms));   
    else return significants;
}

std::string MonteCarloEngine::computeStatistics() {
    
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables statistics cannot be generated.\n");
    std::streamsize ss_prec = std::cout.precision();
    unsigned int rmsPrecision = 2;
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
            StatsLog << "      Mean +- sqrt(V):                " << std::setprecision(getPrecision(bch1d.GetHistogram()->GetMean(),rms))
                    << bch1d.GetHistogram()->GetMean() << " +- " << std::setprecision(rmsPrecision)
                    << rms << std::endl

                    << "      (Marginalized) mode:            " << std::setprecision(getPrecision(bch1d.GetLocalMode(0),rms)) << bch1d.GetLocalMode(0) << std::endl;
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
                    StatsLog << "       (" << std::setprecision(getPrecision(interval_xmin, rms)) << interval_xmin << ", " << std::setprecision(getPrecision(interval_xmax, rms)) << interval_xmax
                            << ") corresponding to " << (interval_xmin + interval_xmax)/2. << " +- " << (-interval_xmin + interval_xmax)/2./(i+1) << " (local mode at " << std::setprecision(getPrecision(interval_mode, rms)) << interval_mode << " with rel. height "
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
                StatsLog << "      Mean +- sqrt(V):                " << std::setprecision(getPrecision(bch1d.GetHistogram()->GetMean(), rms))
                        << bch1d.GetHistogram()->GetMean() << " +- " << std::setprecision(rmsPrecision)
                        << rms << std::endl

                        << "      (Marginalized) mode:            " << std::setprecision(getPrecision(bch1d.GetLocalMode(0), rms)) << bch1d.GetLocalMode(0) << std::endl;

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
                        StatsLog << "       (" << std::setprecision(getPrecision(interval_xmin, rms)) << interval_xmin << ", " << std::setprecision(getPrecision(interval_xmax, rms)) << interval_xmax
                                << ") (local mode at " << std::setprecision(getPrecision(interval_mode, rms)) << interval_mode << " with rel. height "
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

double MonteCarloEngine::SecondDerivative(BCParameter par1, BCParameter par2, std::vector<double> point) {

    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MCMCENgine::SecondDerivative : Invalid number of entries in the vector.");
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

double MonteCarloEngine::FirstDerivative(BCParameter par, std::vector<double> point) {

    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MCMCENgine::FirstDerivative : Invalid number of entries in the vector.");
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

double MonteCarloEngine::Function_h(std::vector<double> point) {
    if (point.size() != GetNParameters()) {
        throw std::runtime_error("MCMCENgine::Function_h : Invalid number of entries in the vector.");
    }
    return LogLikelihood(point) + LogAPrioriProbability(point);
}
