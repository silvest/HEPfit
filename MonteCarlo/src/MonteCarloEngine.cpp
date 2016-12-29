/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MonteCarloEngine.h"
#include <BAT/BCParameter.h>
#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCTF1Prior.h>
#ifdef _MPI
#include <mpi.h>
#endif
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TROOT.h>
#include <TColor.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <fstream>
#include <stdexcept>
#include <iomanip>

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
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif
    if (rank == 0) {
        TH1::StatOverflows(kTRUE);
        TH1::SetDefaultBufferSize(100000);
    }
};

void MonteCarloEngine::Initialize(StandardModel* Mod_i) 
{   
    Mod = Mod_i;
    int k = 0, kweight = 0;
    
    if (rank == 0) {
        TH1D * lhisto = new TH1D("LogLikelihood", "LogLikelihood", NBINS1D, 1., -1.);
        lhisto->GetXaxis()->SetTitle("LogLikelihood");
        BCH1D bclhisto = BCH1D(lhisto);
        Histo1D["LogLikelihood"] = bclhisto;
    }
    
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            k++;
            if (it->getDistr().compare("noweight") != 0)
                kweight++;
        }
        std::string HistName = it->getName();
        thMin[it->getName()] = std::numeric_limits<double>::max();
        thMax[it->getName()] = -std::numeric_limits<double>::max();
        if (rank == 0 && Histo1D.find(HistName) == Histo1D.end()) {
            TH1D * histo = new TH1D(HistName.c_str(), it->getLabel().c_str(), NBINS1D, it->getMin(), it->getMax());
            histo->GetXaxis()->SetTitle(it->getLabel().c_str());
            BCH1D bchisto = BCH1D(histo);
            Histo1D[HistName] = bchisto;
        }
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            if (!it->isTMCMC())
                throw std::runtime_error("ERROR: cannot handle noMCMC for Observable2D file yet!");
        } else if (it->getDistr().compare("weight") == 0)
            throw std::runtime_error("ERROR: do not use Observable2D for analytic 2D weights!");
        std::string HistName = it->getName();
        if (rank == 0 && Histo2D.find(HistName) == Histo2D.end()) {
            TH2D * histo2 = new TH2D(HistName.c_str(),
                    (it->getLabel() + " vs " + it->getLabel2()).c_str(),
                    NBINS2D, it->getMin(), it->getMax(),
                    NBINS2D, it->getMin2(), it->getMax2());
            histo2->GetXaxis()->SetTitle(it->getLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel2().c_str());
            BCH2D bchisto2 = BCH2D(histo2);
            Histo2D[HistName] = bchisto2;
        }
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 != CGO.end(); ++it1) {
        std::vector<Observable> ObsV(it1->getObs());
        for (std::vector<Observable>::iterator it = ObsV.begin();
                it != ObsV.end(); ++it) {
            if ((it->getDistr()).compare("file") == 0)
                throw std::runtime_error("Cannot use file in CorrelatedGaussianObservables!");
            if (!(it->isTMCMC())) {
                k++;
                if (it->getDistr().compare("noweight") != 0)
                    throw std::runtime_error("Cannot use weight in CorrelatedGaussianObservables!");
            }
            std::string HistName = it->getName();
            thMin[HistName] = std::numeric_limits<double>::max();
            thMax[HistName] = -std::numeric_limits<double>::max();
            if (rank == 0 && Histo1D.find(HistName) == Histo1D.end()) {
                TH1D * histo = new TH1D(HistName.c_str(), it->getLabel().c_str(),
                        NBINS1D, it->getMin(), it->getMax());
                histo->GetXaxis()->SetTitle(it->getLabel().c_str());
                BCH1D bchisto = BCH1D(histo);
                Histo1D[HistName] = bchisto;
            }
        }
        if (it1->isPrediction()) {
            CorrelationMap[it1->getName()] = new TPrincipal(it1->getObs().size());
        }
    }
    kmax = k;
    kwmax = kweight;

    unknownParameters = Mod->getUnknownParameters();
    DefineParameters();
    SetMaximumEfficiency(0.5);

};

void MonteCarloEngine::setNChains(unsigned int i) {
    SetNChains(i);
    obval = new double[fMCMCNChains * kmax];
    obweight = new double[fMCMCNChains * kwmax];
}

// ---------------------------------------------------------

MonteCarloEngine::~MonteCarloEngine()
// default destructor
{
    delete [] obval;
    delete [] obweight;
    /* The following code has been commented out pending further review.
       It is causing crashes at the termination of the code if the histograms
       are accessed from the main program.*/
    //    for (std::map<std::string, BCH1D *>::iterator it = Histo1D.begin();
    //            it != Histo1D.end(); it++)
    //        delete it->second;
    //    for (std::map<std::string, BCH2D *>::iterator it = Histo2D.begin();
    //            it != Histo2D.end(); it++)
    //        delete it->second;
    for (std::map<std::string, TPrincipal *>::iterator it = CorrelationMap.begin();
            it != CorrelationMap.end(); it++)
        delete it->second;
//    for (unsigned int i = 0; i < fMCMCNChains; ++i) 
//        if (fMCMCObservablesTrees[i]) delete fMCMCObservablesTrees[i];
    
};

// ---------------------------------------------------------

void MonteCarloEngine::DefineParameters() {
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    if (rank == 0) std::cout << "\nParameters varied in Monte Carlo:" << std::endl;
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
            if (it->geterrf() == 0.) GetParameter(k).SetPrior(new BCGaussianPrior(it->getave(), it->geterrg())); //SetPriorGauss(k, it->getave(), it->geterrg());
            else if (it->geterrg() == 0.) GetParameter(k).SetPriorConstant(); //SetPriorConstant(k);
            else {
                TF1 combined = TF1(it->getname().c_str(),
                        "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                        it->getmin(), it->getmax());
                combined.SetParameter(0, it->getave());
                combined.SetParameter(1, it->geterrg());
                combined.SetParameter(2, it->geterrf());
                GetParameter(k).SetPrior(new BCTF1Prior(combined)); //SetPrior(k, combined);
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
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it < ModPars.end(); it++){
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
    if (std::isnan(logprob)) {
        NumOfDiscardedEvents++;
#ifdef _MCDEBUG
        std::cout << "Event discarded since logprob evaluated to NAN.\n" << logprob << std::endl ;
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
    std::vector<std::vector<double> > fMCMCxvect;
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
        pars = fMCMCx.at(mychain);

        fMCMCxvect.push_back(pars);
        index_chain[iproc] = mychain;
        iproc++;
        mychain++;
        if (iproc < procnum && mychain < fMCMCNChains)
            continue;

        for (int il = 0; il < iproc; il++) {
            //The first entry of the array specifies the task to be executed.

            sendbuff[il][0] = 2.; // 2 = observables calculation
            for (int im = 1; im < buffsize; im++)
                sendbuff[il][im] = fMCMCxvect[il][im - 1];
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
//            std::map<std::string, double> DPars;
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
                int ko = 0, kweight = 0, k_all = 0;
                for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                        it < Obs_ALL.end(); it++) {
                    double th = buff[il][k++];
                    /* set the min and max of theory values */
                    if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                    if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                    Histo1D[it->getName()].GetHistogram()->Fill(th);
                    if (fMCMCFlagWriteChainToFile) hMCMCObservables[index_chain[il]][k_all++] = th;
                    if (!it->isTMCMC()) {
                        obval[index_chain[il] * kmax + ko] = th;
                        
                        ko++;
                        if (it->getDistr().compare("noweight") != 0) {
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
                }

                // fill the histograms for correlated observables
                for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                        it1 < CGO.end(); it1++) {
                    std::vector<Observable> ObsV(it1->getObs());
                    Double_t * COdata = new Double_t[ObsV.size()];
                    int nObs = 0;
                    for (std::vector<Observable>::iterator it = ObsV.begin();
                            it != ObsV.end(); ++it) {
                        double th = buff[il][k++];
                        /* set the min and max of theory values */
                        if (th < thMin[it->getName()]) thMin[it->getName()] = th;
                        if (th > thMax[it->getName()]) thMax[it->getName()] = th;
                        Histo1D[it->getName()].GetHistogram()->Fill(th);
                        if (it1->isPrediction()) COdata[nObs++] = th;
                    }
                    if (it1->isPrediction()) CorrelationMap[it1->getName()]->AddRow(COdata);
                    delete [] COdata;
                }
            }
        }
        iproc = 0;
        fMCMCxvect.clear();
    }
    if (fMCMCFlagWriteChainToFile) InChainFillObservablesTree();
    delete sendbuff[0];
    delete [] sendbuff;
    delete [] recvbuff;
    delete buff[0];
    delete [] buff;
#else
    for (unsigned int i = 0; i < fMCMCNChains; ++i) {
        std::vector<double>::const_iterator first = fMCMCx.at(i).begin();
        std::vector<double>::const_iterator last = first + GetNParameters();
        std::vector<double> currvec(first, last);
        setDParsFromParameters(currvec,DPars);

        Mod->Update(DPars);
        // fill the histograms for observables
        int k = 0, kweight = 0, k_all = 0;
        for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) {
            double th = it->computeTheoryValue();
            /* set the min and max of theory values */
            if (th < thMin[it->getName()]) thMin[it->getName()] = th;
            if (th > thMax[it->getName()]) thMax[it->getName()] = th;
            Histo1D[it->getName()].GetHistogram()->Fill(th);
            if (fMCMCFlagWriteChainToFile) hMCMCObservables[i][k_all++] = th;
            if (!it->isTMCMC()) {
                obval[i * kmax + k] = th;
                k++;
                if (it->getDistr().compare("noweight") != 0) {
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
                if (it1->isPrediction()) COdata[nObs++] = th;
            }
            if (it1->isPrediction()) CorrelationMap[it1->getName()]->AddRow(COdata);
            delete [] COdata;
        }
    }
    
    if (fMCMCFlagWriteChainToFile) InChainFillObservablesTree();
#endif
    for (unsigned int i = 0; i < fMCMCNChains; i++)
    {
        Histo1D["LogLikelihood"].GetHistogram()->Fill(GetLogProbx(i)-LogAPrioriProbability(Getx(i)));
    }
}

void MonteCarloEngine::CheckHistogram(TH1& hist, const std::string name) {
    // output the portions of underflow and overflow bins
    double UnderFlowContent = hist.GetBinContent(0);
    double OverFlowContent = hist.GetBinContent(NBINS1D + 1);
    double Integral = hist.Integral();
    double TotalContent = 0.0;
    for (int n = 0; n <= NBINS1D + 1; n++)
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
    for (int m = 0; m <= NBINS2D + 1; m++)
        for (int n = 0; n <= NBINS2D + 1; n++)
            TotalContent += hist.GetBinContent(m, n);
    HistoLog << name << ": "
            << Integral / TotalContent * 100. << "% within the ranges"
            << std::endl;
}

void MonteCarloEngine::Print1D(BCH1D bch1d, const char* filename, int ww, int wh) {
    // create temporary canvas
    TCanvas * c;
    cindex++;
    if(ww > 0 && wh > 0)
        c = new TCanvas(TString::Format("c_bch1d_%d",cindex), TString::Format("c_bch1d_%d",cindex), ww, wh);
    else
        c = new TCanvas(TString::Format("c_bch1d_%d",cindex));

    bch1d.GetHistogram()->Scale(1./bch1d.GetHistogram()->Integral("width"));
    
    double xRange = bch1d.GetHistogram()->GetXaxis()->GetXmax() - bch1d.GetHistogram()->GetXaxis()->GetXmin();
    double yRange = bch1d.GetHistogram()->GetMaximum() - bch1d.GetHistogram()->GetMinimum();

    double xL = bch1d.GetHistogram()->GetXaxis()->GetXmin()+0.035*xRange;
    double yL = bch1d.GetHistogram()->GetMinimum()+1.25*yRange;
    
    double xR = xL+0.18*xRange;
    double yR = yL+0.13*yRange;
    
    int gIdx = 1000;
    int rIdx = 1001;

    TColor green = TColor(gIdx, 0.0, 0.56, 0.57);
    TColor red = TColor(rIdx, 0.57, 0.01, 0.00);

    TBox b1 = TBox(xL, yL, xR, yR);
    b1.SetFillColor(gIdx);
    
    TBox b2 = TBox(xL+0.008*xRange, yL+0.013*yRange, xR-0.008*xRange, yR-0.013*yRange);
    b2.SetFillColor(kWhite);
    
    TPaveText b3 = TPaveText(xL+0.014*xRange, yL+0.02*yRange, xL+0.70*(xR-xL), yR-0.02*yRange);
    b3.SetTextAlign(22);
    b3.SetTextSize(0.044);
    b3.SetTextColor(kWhite);
    b3.AddText("HEP");
    b3.SetFillColor(rIdx);
    
    TPaveText b4 = TPaveText(xL+0.7*(xR-xL), yL+0.027*yRange, xR-0.008*xRange, yR-0.015*yRange);
    b4.SetTextAlign(33);
    b4.SetTextSize(0.037);
    b4.SetTextColor(rIdx);
    b4.AddText("fit");
    b4.SetFillColor(kWhite);
    
    bch1d.SetBandType(BCH1D::kSmallestInterval);
    bch1d.SetBandColor(0, gIdx);
    bch1d.SetBandColor(1, rIdx);
    bch1d.SetBandColor(2, kOrange - 3);
    bch1d.SetNBands(3);
    bch1d.SetNSmooth(0);
    bch1d.SetDrawGlobalMode(true);
    bch1d.SetDrawMean(true, true);
    bch1d.SetDrawLegend(true);
    bch1d.SetNLegendColumns(1);
    bch1d.SetStats(true);
    
    bch1d.Draw();
    
    // draw logo
    if (printLogo) {
        b1.Draw("SAME");
        b2.Draw("SAME");
        b3.Draw("SAME");
        b4.Draw("SAME");
    }
    
    // print to file.
    c->Print(filename);
}

void MonteCarloEngine::Print2D(BCH2D bch2d, const char * filename, int ww, int wh)
{
    // create temporary canvas
    TCanvas * c;
    cindex++;
    if(ww > 0 && wh > 0)
        c = new TCanvas(TString::Format("c_bch2d_%d",cindex), TString::Format("c_bch2d_%d",cindex), ww, wh);
    else
        c = new TCanvas(TString::Format("c_bch2d_%d",cindex));
    
    bch2d.GetHistogram()->Scale(1./bch2d.GetHistogram()->Integral("width"));
    
    double xRange = bch2d.GetHistogram()->GetXaxis()->GetXmax() - bch2d.GetHistogram()->GetXaxis()->GetXmin();
    double yRange = bch2d.GetHistogram()->GetYaxis()->GetXmax() - bch2d.GetHistogram()->GetYaxis()->GetXmin();

    double xL = bch2d.GetHistogram()->GetXaxis()->GetXmin()+0.035*xRange;
    double yL = bch2d.GetHistogram()->GetYaxis()->GetXmin()+0.89*yRange;
    
    double xR = xL+0.18*xRange;
    double yR = yL+0.09*yRange;
    
    int gIdx = 1000;
    int rIdx = 1001;

    TColor green = TColor(gIdx, 0.0, 0.56, 0.57);
    TColor red = TColor(rIdx, 0.57, 0.01, 0.00);

    TBox b1 = TBox(xL, yL, xR, yR);
    b1.SetFillColor(gIdx);
    
    TBox b2 = TBox(xL+0.008*xRange, yL+0.008*yRange, xR-0.008*xRange, yR-0.008*yRange);
    b2.SetFillColor(kWhite);
    
    TPaveText b3 = TPaveText(xL+0.014*xRange, yL+0.013*yRange, xL+0.70*(xR-xL), yR-0.013*yRange);
    b3.SetTextAlign(22);
    b3.SetTextSize(0.044);
    b3.SetTextColor(kWhite);
    b3.AddText("HEP");
    b3.SetFillColor(rIdx);
    
    TPaveText b4 = TPaveText(xL+0.75*(xR-xL), yL+0.024*yRange, xR-0.008*xRange, yR-0.013*yRange);
    b4.SetTextAlign(33);
    b4.SetTextSize(0.038);
    b4.SetTextColor(rIdx);
    b4.AddText("fit");
    b4.SetFillColor(kWhite);
    
    bch2d.SetBandType(BCH2D::kSmallestInterval);
    bch2d.SetBandColor(0, kOrange - 3); 
    bch2d.SetBandColor(1, rIdx);
    bch2d.SetBandColor(2, gIdx);
    bch2d.SetNBands(3);
    bch2d.SetBandFillStyle(1001);
    bch2d.SetNSmooth(0);
    bch2d.SetDrawLocalMode(false);
    bch2d.SetDrawGlobalMode(true);
    bch2d.SetDrawMean(true, true);
    bch2d.SetDrawLegend(true);
    bch2d.SetNLegendColumns(1);
    bch2d.SetStats(true);
    
    bch2d.Draw();

    //draw logo
    if (printLogo) {
        b1.Draw("SAME");
        b2.Draw("SAME");
        b3.Draw("SAME");
        b4.Draw("SAME");
    }
    
    // print to file
    c->Print(filename);
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
        std::cout << fname << " has been created." << std::endl;
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
    if (mode.size() == 0) {
        if (rank == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables histogram cannot be generated.\n");
        else sleep(2);
    }
    setDParsFromParameters(mode,DPars);

    Mod->Update(DPars);

    // print the histograms to pdf files
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end();
            it++) {
        PrintHistogram(OutFile, *it, OutputDir);
    }
    
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 < CGO.end(); it1++) {
        std::vector<Observable> ObsV(it1->getObs());
        for (std::vector<Observable>::iterator it = ObsV.begin();
                it != ObsV.end(); ++it)
            PrintHistogram(OutFile, *it, OutputDir);
    }
    
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        std::string HistName = it->getName();
        if (Histo2D[HistName].GetHistogram()->Integral() > 0.0) {
            std::string fname = OutputDir + "/" + HistName + ".pdf";
            std::vector<double> th;
            th.push_back(it->computeTheoryValue());
            th.push_back(it->computeTheoryValue2());
            Histo2D[HistName].SetGlobalMode(th);
            Print2D(Histo2D[HistName], fname.c_str());
            std::cout << fname << " has been created." << std::endl;
            if(OutFile.compare("") == 0) {
                throw std::runtime_error("\nMonteCarloEngine::PrintHistogram ERROR: No root file specified for writing histograms.");
            }
            TDirectory * dir = gDirectory;
            GetOutputFile()->cd();
            Histo2D[HistName].GetHistogram()->Write();
            gDirectory = dir;
            CheckHistogram(*Histo2D[HistName].GetHistogram(), HistName);
        } else
            HistoLog << "WARNING: The histogram of "
                << HistName << " is empty!" << std::endl;
    }
        
    if (Histo1D["LogLikelihood"].GetHistogram()->Integral() > 0.0) {
        std::string fname = OutputDir + "/LogLikelihood.pdf";
        Print1D(Histo1D["LogLikelihood"], fname.c_str());
        std::cout << fname << " has been created." << std::endl;
        TDirectory * dir = gDirectory;
        GetOutputFile()->cd();
        Histo1D["LogLikelihood"].GetHistogram()->Write();
        gDirectory = dir;
        CheckHistogram(*Histo1D["LogLikelihood"].GetHistogram(), "LogLikelihood");
    }
}

void MonteCarloEngine::AddChains() {
    InitializeMarkovChainTree();
    TDirectory* dir = gDirectory;
    GetOutputFile()->cd();
    
    hMCMCObservableTree = new TTree(TString::Format("%s_Observables", GetSafeName().data()), TString::Format("%s_Observables", GetSafeName().data()));
    hMCMCObservableTree->Branch("Chain", &fMCMCTree_Chain, "chain/i");
    hMCMCObservableTree->Branch("Iteration", &fMCMCCurrentIteration, "iteration/i");
    // hMCMCObservableTree->Branch("Phase", &fMCMCPhase, "phase/I");
    // hMCMCObservableTree->Branch("LogProbability", &fMCMCTree_Prob, "log(probability)/D");
    hMCMCObservables.assign(fMCMCNChains, std::vector<double>(Obs_ALL.size(), 0.));
    hMCMCTree_Observables.assign(Obs_ALL.size(), 0.);  
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
    hMCMCObservableTree->SetAutoSave(10 * fMCMCNIterationsPreRunCheck);
    hMCMCObservableTree->AutoSave("SelfSave");
    gDirectory = dir;
}

void MonteCarloEngine::InChainFillObservablesTree()
{
    if (!hMCMCObservableTree)
        return;
    // loop over all chains
    for (fMCMCTree_Chain = 0; fMCMCTree_Chain < fMCMCNChains; ++fMCMCTree_Chain) {
        hMCMCTree_Observables = hMCMCObservables[fMCMCTree_Chain];
        if (kwmax > 0) hMCMCTree_Observables_weight = hMCMCObservables_weight[fMCMCTree_Chain];
        hMCMCObservableTree->Fill();
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

    return 2 + ceil(log10(fabs(value)))-ceil(log10(rms));   
}

std::string MonteCarloEngine::computeStatistics() {
    
    std::vector<double> mode(GetBestFitParameters());
    if (mode.size() == 0) {
        if(rank == 0) throw std::runtime_error("\n ERROR: Global Mode could not be determined possibly because of infinite loglikelihood. Observables statistics cannot be generated.\n");
        else sleep (2);
    }
    std::ostringstream StatsLog;
    int i = 0;
    StatsLog << "Statistics file for Observables, Binned Observables and Correlated Gaussian Observables.\n" << std::endl;
    if (Obs_ALL.size() > 0) StatsLog << "Observables:\n" << std::endl;
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end(); it++) {

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
                    << bch1d.GetHistogram()->GetMean() << " +- " << std::setprecision(2)
                    << rms << std::endl

                    << "      (Marginalized) mode:            " << std::setprecision(getPrecision(bch1d.GetLocalMode(0),rms)) << bch1d.GetLocalMode(0) << std::endl;
            std::vector<double> intervals;
            intervals.push_back(0.682689492137);
            intervals.push_back(0.954499736104);
            intervals.push_back(0.997300203937);
            std::vector<BCH1D::BCH1DSmallestInterval> v = bch1d.GetSmallestIntervals(intervals);
            for (unsigned int i = 0; i < v.size(); i++) {
                StatsLog << "      Smallest interval(s) containing at least " << v[i].total_mass * 100 << "% and local mode(s):"
                        << std::endl;
                for (unsigned j = 0; j < v[i].intervals.size(); j++) {
                    double interval_xmin = v[i].intervals[j].xmin;
                    double interval_xmax = v[i].intervals[j].xmax;
                    double interval_mode = v[i].intervals[j].mode;
                    double interval_heignt = v[i].intervals[j].relative_height;
                    double interval_relative_mass = v[i].intervals[j].relative_mass;
                    StatsLog << "       (" << std::setprecision(getPrecision(interval_xmin, rms)) << interval_xmin << ", " << std::setprecision(getPrecision(interval_xmax, rms)) << interval_xmax
                            << ") (local mode at " << std::setprecision(getPrecision(interval_mode, rms)) << interval_mode << " with rel. height "
                            << std::setprecision(getPrecision(interval_heignt, rms)) << interval_heignt << "; rel. area " << std::setprecision(getPrecision(interval_relative_mass, rms)) << interval_relative_mass << ")"
                            << std::endl;
                    StatsLog << std::endl;
                }
            }
        } else {
            StatsLog << "\nWARNING: The histogram of " << it->getName() << " is empty! Statistics cannot be generated\n" << std::endl;
        }
    }
    
    if (CGO.size() > 0) StatsLog << "\nCorrelated (Gaussian) Observables:\n" << std::endl;
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 < CGO.end(); it1++) {
        StatsLog << "\n" << it1->getName() << ":\n" << std::endl;
        i = 0;
        std::vector<Observable> CGObs(it1->getObs());
        for (std::vector<Observable>::iterator it2 = CGObs.begin(); it2 < CGObs.end(); it2++) {

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
                        << bch1d.GetHistogram()->GetMean() << " +- " << std::setprecision(2)
                        << rms << std::endl

                        << "      (Marginalized) mode:            " << std::setprecision(getPrecision(bch1d.GetLocalMode(0), rms)) << bch1d.GetLocalMode(0) << std::endl;

                std::vector<double> intervals;
                intervals.push_back(0.682689492137);
                intervals.push_back(0.954499736104);
                intervals.push_back(0.997300203937);

                std::vector<BCH1D::BCH1DSmallestInterval> v = bch1d.GetSmallestIntervals(intervals);
                for (unsigned int i = 0; i < v.size(); i++) {
                    StatsLog << "      Smallest interval(s) containing at least " << v[0].total_mass * 100 << "% and local mode(s):"
                            << std::endl;
                    for (unsigned j = 0; j < v[0].intervals.size(); j++) {
                        double interval_xmin = v[0].intervals[j].xmin;
                        double interval_xmax = v[0].intervals[j].xmax;
                        double interval_mode = v[0].intervals[j].mode;
                        double interval_heignt = v[0].intervals[j].relative_height;
                        double interval_relative_mass = v[0].intervals[j].relative_mass;
                        StatsLog << "       (" << std::setprecision(getPrecision(interval_xmin, rms)) << interval_xmin << ", " << std::setprecision(getPrecision(interval_xmax, rms)) << interval_xmax
                                << ") (local mode at " << std::setprecision(getPrecision(interval_mode, rms)) << interval_mode << " with rel. height "
                                << std::setprecision(getPrecision(interval_heignt, rms)) << interval_heignt << "; rel. area " << std::setprecision(getPrecision(interval_relative_mass, rms)) << interval_relative_mass << ")"
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
            TMatrixD * corr = const_cast<TMatrixD*>(CorrelationMap[it1->getName()]->GetCovarianceMatrix());
            *corr *= (double)size;
            StatsLog << "\nThe correlation matrix for " << it1->getName() << " is given by the " << size << "x"<< size << " matrix:\n" << std::endl;

            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(4) << "" << " | ";
                else StatsLog << std::setw(5) << i << std::setw(5) << "    |";
            }
            StatsLog << std::endl;
            for (int i = 0; i < size + 1; i++) {
                if (i == 0) StatsLog << std::setw(8) << "--------";
                else StatsLog << std::setw(10) << "----------";
            }
            StatsLog << std::endl;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size + 1; j++) {
                    if (j == 0) StatsLog << std::setw(4) << i+1 << " |";
                    else StatsLog << std::setprecision(5) << std::setw(10) << (*corr)(i, j - 1);
                }
            StatsLog << std::endl;
            }
        }
    }

    double llika = Histo1D["LogLikelihood"].GetHistogram()->GetMean();
    StatsLog << "LogLikelihood mean value: " << llika << std::endl;
    double llikv = Histo1D["LogLikelihood"].GetHistogram()->GetRMS();
    llikv *= llikv;
    StatsLog << "LogLikelihood variance: " << llikv << std::endl;
    double dbar = -2.*llika; //Wikipedia notation... 
    double pd = 2.*llikv; //Wikipedia notation...
    StatsLog << "IC value (don't ask me what it means...): " << dbar + 2.*pd << std::endl; 
    StatsLog << "DIC value (same as above...): " << dbar + pd << std::endl; 
    StatsLog << std::endl;
    StatsLog << "Value of the Parameters and Observables at the global mode:" << std::endl;
    StatsLog << std::endl;
    
    setDParsFromParameters(mode,DPars);

    Mod->Update(DPars);

    StatsLog << std::setprecision(5);
    
    for (std::map<std::string,double>::iterator it = DPars.begin(); it != DPars.end(); it++)
        StatsLog << it->first << ": " << it->second << std::endl;
       
    StatsLog << std::endl;
    
    for (boost::ptr_vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) 
        StatsLog << it->getName() << ": " << it->computeTheoryValue() << std::endl;

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

double MonteCarloEngine::computeNormalizationMC(int NIterationNormalizationMC) {
    // Number of MC iterations
    SetNIterationsMin(NIterationNormalizationMC);
    SetIntegrationMethod(BCIntegrate::kIntMonteCarlo);
    Integrate();
    double norm = GetIntegral();
    if (norm < 0.) {
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
