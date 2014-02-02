/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MonteCarloEngine.h"
#include <BAT/BCParameter.h>
#include <BAT/BCMath.h>
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <fstream>
#include <stdexcept>

MonteCarloEngine::MonteCarloEngine(
        const std::vector<ModelParameter>& ModPars_i,
        std::vector<Observable>& Obs_i,
        std::vector<Observable2D>& Obs2D_i,
        std::vector<CorrelatedGaussianObservables>& CGO_i,
        std::vector<ModelParaVsObs>& ParaObs_i, const bool checkHistRange_i)
: BCModel(""), ModPars(ModPars_i), Obs_ALL(Obs_i), Obs2D_ALL(Obs2D_i),
        CGO(CGO_i), ParaObs(ParaObs_i), NumOfUsedEvents(0), NumOfDiscardedEvents(0),
        checkTheoryRange(checkHistRange_i)
{
    obval = NULL;
    obweight = NULL;
    Mod = NULL;
};

void MonteCarloEngine::Initialize(Model* Mod_i) 
{
    Mod = Mod_i;
    int k = 0, kweight = 0;
    for (std::vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            TFile *lik = new TFile((it->getFilename() + ".root").c_str(), "read");
            TH1D *htmp = (TH1D*) (lik->Get(it->getHistoname().c_str()));
            if (htmp == NULL)
                throw std::runtime_error("ERROR: nonexistent histogram called "
                                         + it->getHistoname() + " in "
                                         + it->getFilename() + ".root");
            TH1D *inhisto = (TH1D *) htmp->Clone((it->getFilename() + "_" + it->getHistoname()).c_str());
            inhisto->SetDirectory(gROOT);
            InHisto1D[it->getFilename() + it->getHistoname()] = inhisto;
            std::cout << "added input histogram "
                      << InHisto1D[it->getFilename() + it->getHistoname()]->GetName()
                      << std::endl;
            it->setMin(inhisto->GetXaxis()->GetXmin());
            it->setMax(inhisto->GetXaxis()->GetXmax());
            lik->Close();
            delete lik;
        }
        if (it->isTMCMC()) {
            Obs_MCMC.push_back(*it);
        } else {
            k++;
            if (it->getDistr().compare("noweight") != 0)
                kweight++;
        }
        if (Histo1D.find(it->getThname()) == Histo1D.end()) {
            TH1D * histo = new TH1D(it->getThname().c_str(), it->getLabel().c_str(), NBINS1D,
                    it->getMin(), it->getMax());
            histo->GetXaxis()->SetTitle(it->getLabel().c_str());
            BCH1D * bchisto = new BCH1D(histo);
            Histo1D[it->getThname()] = bchisto;
            thMin[it->getThname()] = std::numeric_limits<double>::max();
            thMax[it->getThname()] = - std::numeric_limits<double>::max();
        }
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            TFile *lik2 = new TFile((it->getFilename() + ".root").c_str(), "read");
            TH2D *htmp2 = (TH2D*) (lik2->Get(it->getHistoname().c_str()));
            if (htmp2 == NULL)
                throw std::runtime_error("ERROR: nonexistent histogram called "
                                         + it->getHistoname()
                                         + " in " + it->getFilename() + ".root");
            TH2D *inhisto2 = (TH2D *) htmp2->Clone((it->getFilename() + "_" + it->getHistoname()).c_str());
            inhisto2->SetDirectory(gROOT);
            InHisto2D[it->getFilename() + it->getHistoname()] = inhisto2;
            std::cout << "added 2D input histogram "
                      << InHisto2D[it->getFilename() + it->getHistoname()]->GetName()
                      << std::endl;
            it->setMin(inhisto2->GetXaxis()->GetXmin());
            it->setMax(inhisto2->GetXaxis()->GetXmax());
            it->setMin2(inhisto2->GetYaxis()->GetXmin());
            it->setMax2(inhisto2->GetYaxis()->GetXmax());
            lik2->Close();
            delete lik2;
            if (it->isTMCMC()) {
                Obs2D_MCMC.push_back(*it);
            } else
                throw std::runtime_error("ERROR: cannot handle noMCMC for Observable2D file yet!");
        } else if (it->getDistr().compare("weight") == 0)
            throw std::runtime_error("ERROR: do not use Observable2D for analytic 2D weights!");
        if (Histo2D.find(it->getThname() + "_vs_" + it->getThname2()) == Histo2D.end()) {
            TH2D * histo2 = new TH2D((it->getThname() + "_vs_" + it->getThname2()).c_str(),
                                     (it->getLabel() + " vs " + it->getLabel2()).c_str(), 
                                     NBINS2D, it->getMin(), it->getMax(),
                                     NBINS2D, it->getMin2(), it->getMax2());
            histo2->GetXaxis()->SetTitle(it->getLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel2().c_str());
            BCH2D * bchisto2 = new BCH2D(histo2);
            Histo2D[it->getThname() + "_vs_" + it->getThname2()] = bchisto2;
        }
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); 
            it1 != CGO.end(); ++it1) {
        std::vector<Observable> pino(it1->getObs());
        for (std::vector<Observable>::iterator it = pino.begin();
                it != pino.end(); ++it) {
            //for (int i = 0; i < it1->getObs().size(); i++) {
            //Observable * it = &(it1->getObs().at(i));
            if ((it->getDistr()).compare("file") == 0)
                throw std::runtime_error("Cannot use file in CorrelatedGaussianObservables!");
            if (!(it->isTMCMC())) {
                k++;
                if (it->getDistr().compare("noweight") != 0)
                    throw std::runtime_error("Cannot use weight in CorrelatedGaussianObservables!");
            }
            if (Histo1D.find(it->getThname()) == Histo1D.end()) {
                TH1D * histo = new TH1D(it->getThname().c_str(), it->getLabel().c_str(),
                                        NBINS1D, it->getMin(), it->getMax());
                histo->GetXaxis()->SetTitle(it->getLabel().c_str());
                BCH1D * bchisto = new BCH1D(histo);	
                Histo1D[it->getThname()] = bchisto;
                thMin[it->getThname()] = std::numeric_limits<double>::max();
                thMax[it->getThname()] = - std::numeric_limits<double>::max();
            }
        }
    }
    kmax = k;
    kwmax = kweight;

    DefineParameters();

    for (std::vector<ModelParaVsObs>::iterator it = ParaObs.begin();
            it < ParaObs.end(); it++) {

        /* check if the parameter in ModelParaVsObs exists in MCMCparameters */
        bool checkParam = false;
        for (int k = 0; k < GetNParameters(); k++)
            if (it->getParaName().compare(GetParameter(k)->GetName())==0)
                checkParam = true;
        if(!checkParam)
            throw std::runtime_error(it->getParaName() + " cannot be used in ModelParaVsObs!");

        if (Histo2D.find(it->getParaName() + "_vs_" + it->getThname()) == Histo2D.end()) {
            TH2D * histo2 = new TH2D((it->getParaName() + "_vs_" + it->getThname()).c_str(),
                                     (it->getParaLabel() + " vs " + it->getLabel()).c_str(),
                                     NBINS2D, it->getParaMin(), it->getParaMax(),
                                     NBINS2D, it->getMin(), it->getMax());
            histo2->GetXaxis()->SetTitle(it->getParaLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel().c_str());
            BCH2D * bchisto2 = new BCH2D(histo2);
            Histo2D[it->getParaName() + "_vs_" + it->getThname()] = bchisto2;
        }
    }
};

void MonteCarloEngine::setNChains(unsigned int i) 
{
    MCMCSetNChains(i);
    obval = new double[fMCMCNChains * kmax];
    obweight = new double[fMCMCNChains * kwmax];
}

// ---------------------------------------------------------

MonteCarloEngine::~MonteCarloEngine()
// default destructor
{
    delete [] obval;
    delete [] obweight;
    for (std::map<std::string, BCH1D *>::iterator it = Histo1D.begin();
            it != Histo1D.end(); it++)
        delete it->second;
    for (std::map<std::string, BCH2D *>::iterator it = Histo2D.begin();
            it != Histo2D.end(); it++)
        delete it->second;
    for (std::map<std::string, TH1D *>::iterator it = InHisto1D.begin();
            it != InHisto1D.end(); it++)
        delete it->second;
    for (std::map<std::string, TH2D *>::iterator it = InHisto2D.begin();
            it != InHisto2D.end(); it++)
        delete it->second;
};

// ---------------------------------------------------------

void MonteCarloEngine::DefineParameters() 
{
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    std::cout << "Parameters varied in Monte Carlo:" << std::endl;
    int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin();
            it < ModPars.end(); it++) {
        if (it->errf == 0. && it->errg == 0.)
            continue;
        AddParameter(it->name.c_str(), it->min, it->max);
        std::cout << "  " << it->name << " " << k << std::endl;
        DPars[it->name] = 0.;
        if (it->errf == 0.) SetPriorGauss(k, it->ave, it->errg);
        else if (it->errg == 0.) SetPriorConstant(k);
        else {
            TF1 * combined = new TF1(it->name.c_str(),
                    "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                    it->min, it->max);
            combined->SetParameter(0, it->ave);
            combined->SetParameter(1, it->errg);
            combined->SetParameter(2, it->errf);
            SetPrior(k, combined);
            delete combined;
        }
        k++;
    }
}

// ---------------------------------------------------------

double MonteCarloEngine::Weight(const Observable& obs, const double& th) 
{
    double logprob;
    if (obs.getDistr().compare("weight") == 0) {
        if (obs.getErrf() == 0.)
            logprob = BCMath::LogGaus(th, obs.getAve(), obs.getErrg());
        else if (obs.getErrg() == 0.) {
            if (th < obs.getAve() + obs.getErrf() && th > obs.getAve() - obs.getErrf())
                logprob = 1.;
            else
                logprob = log(0.);
        } else
            logprob = log(TMath::Erf((th - obs.getAve() + obs.getErrf()) / sqrt(2.) / obs.getErrg())
                      - TMath::Erf((th - obs.getAve() - obs.getErrf()) / sqrt(2.) / obs.getErrg()));
    } else if (obs.getDistr().compare("file") == 0) {
        TH1D * h = InHisto1D[obs.getFilename() + obs.getHistoname()];
        int i = h->FindBin(th);
        if (h->IsBinOverflow(i) || h->IsBinUnderflow(i))
            logprob = log(0.);
        else
            logprob = log(h->GetBinContent(i));
        //logprob = log(h->GetBinContent(h->FindBin(th)));
    } else
        throw std::runtime_error("ERROR: MonteCarloEngine::Weight() called without weight for "
                                 + obs.getName());
    return (logprob);
}

double MonteCarloEngine::Weight(const Observable2D& obs, const double& th1, const double& th2) 
{
    double logprob;
    if (obs.getDistr().compare("file") == 0) {
        TH2D * h = InHisto2D[obs.getFilename() + obs.getHistoname()];
        int i = h->FindBin(th1, th2);
        if (h->IsBinOverflow(i) || h->IsBinUnderflow(i))
            logprob = log(0.);
        else
            logprob = log(h->GetBinContent(i));
        //logprob = log(h->GetBinContent(h->GetXaxis()->FindBin(th1),h->GetYaxis()->FindBin(th2)));
    } else
        throw std::runtime_error("ERROR: 2D MonteCarloEngine::Weight() called without file for "
                                 + obs.getName());
    return (logprob);
}

double MonteCarloEngine::Weight(const CorrelatedGaussianObservables& obs) 
{

    int size = obs.getObs().size();
    gslpp::vector<double> x(size);

    for (int i = 0; i < size; i++) {
        x(i) = obs.getObs().at(i).computeTheoryValue() - obs.getObs().at(i).getAve();
    }

    return (-0.5 * x * (obs.getCov() * x));
}

double MonteCarloEngine::LogLikelihood(const std::vector<double>& parameters) 
{
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.

    double logprob = 0.;

    for (int k = 0; k < parameters.size(); k++) {
        //        std::string pippo = GetParameter(k)->GetName();
        //        double pluto = parameters[k];
        //        DPars[pippo]=pluto;
        DPars[GetParameter(k)->GetName()] = parameters[k];
    }

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
    NumOfUsedEvents++;
#ifdef _MCDEBUG
    //std::cout << "event used in MC" << std::endl;
#endif

    for (std::vector<Observable>::iterator it = Obs_MCMC.begin(); it < Obs_MCMC.end(); it++) {
        double th = it->computeTheoryValue();
        logprob += Weight(*it, th);
    }

    for (std::vector<Observable2D>::iterator it = Obs2D_MCMC.begin(); it < Obs2D_MCMC.end(); it++) {
        double th1 = it->computeTheoryValue();
        double th2 = it->computeTheoryValue2();
        logprob += Weight(*it, th1, th2);
    }

    for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++) {
        logprob += Weight(*it);
    }
    //std::cout << "logprob " << logprob <<std::endl;    
    return logprob;
}

void MonteCarloEngine::MCMCIterationInterface() 
{

    for (int i = 0; i < fMCMCNChains; ++i) {
        for (int k = 0; k < GetNParameters(); k++) {
            //        std::string pippo = GetParameter(k)->GetName();
            //        double pluto = parameters[k];
            //        DPars[pippo]=pluto;
            DPars[GetParameter(k)->GetName()] = fMCMCx.at(i * GetNParameters() + k);
        }

        Mod->Update(DPars);

        // fill the histograms for observables
        int k = 0, kweight = 0;
        for (std::vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) {
            double th = it->computeTheoryValue();

            /* set the min and max of theory values */
            if (checkTheoryRange) {
                if (th < thMin[it->getThname()]) thMin[it->getThname()] = th;
                if (th > thMax[it->getThname()]) thMax[it->getThname()] = th;
            }
            
            Histo1D[it->getThname()]->GetHistogram()->Fill(th);
            if (!it->isTMCMC()) {
                obval[i * kmax + k] = th;
                k++;
                if (it->getDistr().compare("noweight") != 0) {
                    obweight[i * kwmax + kweight] = Weight(*it, th);
                    kweight++;
                }
            }
        }

        // fill the 2D histograms for observables
        for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
                it < Obs2D_ALL.end(); it++) {
            double th1 = it->computeTheoryValue();
            double th2 = it->computeTheoryValue2();
            Histo2D[it->getThname() + "_vs_" + it->getThname2()]->GetHistogram()->Fill(th1, th2);
        }

        // fill the histograms for correlated observables
        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
                it1 < CGO.end(); it1++){
            std::vector<Observable> pino(it1->getObs());
            for (std::vector<Observable>::iterator it = pino.begin();
                    it != pino.end(); ++it) {
                double th = it->computeTheoryValue();

                /* set the min and max of theory values */
                if (checkTheoryRange) {
                    if (th < thMin[it->getThname()]) thMin[it->getThname()] = th;
                    if (th > thMax[it->getThname()]) thMax[it->getThname()] = th;
                }

                Histo1D[it->getThname()]->GetHistogram()->Fill(th);
            }
        }

        // fill the 2D histograms for ModelParaVsObs
        for (std::vector<ModelParaVsObs>::iterator it = ParaObs.begin();
                it < ParaObs.end(); it++) {
            double par = DPars[it->getParaName()];
            double th = it->computeTheoryValue();
            Histo2D[it->getParaName() + "_vs_" + it->getThname()]->GetHistogram()->Fill(par, th);
        }
    }
}

void MonteCarloEngine::CheckHistogram(const TH1D& hist, const std::string name)
{
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

void MonteCarloEngine::CheckHistogram(const TH2D& hist, const std::string name)
{
    double Integral = hist.Integral();
    double TotalContent = 0.0;
    for (int m = 0; m <= NBINS2D + 1; m++)
        for (int n = 0; n <= NBINS2D + 1; n++)
            TotalContent += hist.GetBinContent(m,n);
    HistoLog << name << ": "
             << Integral / TotalContent * 100. << "% within the ranges"
             << std::endl;
}

void MonteCarloEngine::PrintHistogram(BCModelOutput & out, 
                                      std::vector<Observable>::iterator it,
                                      const std::string OutputDir)
{
    if (Histo1D[it->getThname()]->GetHistogram()->Integral() > 0.0) {
        std::string fname = OutputDir + "/" + it->getThname() + ".pdf";
        //        BCH1D* pippo =  Histo1D[it->getThname()];
        //        double x = pippo->GetMean();
        //        pippo->Print("Dmd1.pdf");
        Histo1D[it->getThname()]->SetGlobalMode(it->computeTheoryValue());
        Histo1D[it->getThname()]->Print(fname.c_str());
        std::cout << fname << " has been created." << std::endl;
        out.Write(Histo1D[it->getThname()]->GetHistogram());
        CheckHistogram(*Histo1D[it->getThname()]->GetHistogram(), it->getThname());
    } else
        HistoLog << "WARNING: The histogram of "
                 << it->getThname() << " is empty!" << std::endl;

    if (checkTheoryRange) {
        double min = thMin[it->getThname()];
        double max = thMax[it->getThname()];
        double range = max - min;
        HistoLog << "  [" << min << ", " << max << "] --> suggested range: "
                 << min - range/7.0 << " " << max + range/7.0 << std::endl;
    }
}

void MonteCarloEngine::PrintHistogram(BCModelOutput & out, const std::string OutputDir)
{
    std::vector<double> mode(GetBestFitParameters());
    for (int k = 0; k < GetNParameters(); k++)
        DPars[GetParameter(k)->GetName()] = mode[k];
    Mod->Update(DPars);

    // print the histograms to pdf files
    for (std::vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end();
            it++) {
        PrintHistogram(out, it, OutputDir);
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin();
            it1 < CGO.end(); it1++){
        std::vector<Observable> pino(it1->getObs());
        for (std::vector<Observable>::iterator it = pino.begin();
                it != pino.end(); ++it)
            PrintHistogram(out, it, OutputDir);
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        std::string HistName = it->getThname() + "_vs_" + it->getThname2();
        if (Histo2D[HistName]->GetHistogram()->Integral() > 0.0) {
            std::string fname = OutputDir + "/" + HistName + ".pdf";
            double th[2];
            th[0] = it->computeTheoryValue();
            th[1] = it->computeTheoryValue2();
            Histo2D[HistName]->SetGlobalMode(th);
            Histo2D[HistName]->Print(fname.c_str());
            std::cout << fname << " has been created." << std::endl;
            out.Write(Histo2D[HistName]->GetHistogram());
            CheckHistogram(*Histo2D[HistName]->GetHistogram(), HistName);
        } else
            HistoLog << "WARNING: The histogram of "
                     << HistName << " is empty!" << std::endl;
    }
    for (std::vector<ModelParaVsObs>::iterator it = ParaObs.begin();
            it < ParaObs.end(); it++) {
        std::string HistName = it->getParaName() + "_vs_" + it->getThname();
        if (Histo2D[HistName]->GetHistogram()->Integral() > 0.0) {
            std::string fname = OutputDir + "/" + HistName + ".pdf";
            double th[2];
            th[0] = DPars[it->getParaName()];
            th[1] = it->computeTheoryValue();
            Histo2D[HistName]->SetGlobalMode(th);
            Histo2D[HistName]->Print(fname.c_str());
            std::cout << fname << " has been created." << std::endl;
            out.Write(Histo2D[HistName]->GetHistogram());
            CheckHistogram(*Histo2D[HistName]->GetHistogram(), HistName);
        } else
            HistoLog << "WARNING: The histogram of "
                     << HistName << " is empty!" << std::endl;
    }
}

void MonteCarloEngine::AddChains() 
{
    fMCMCFlagWritePreRunToFile = false;
    int k = 0, kweight = 0;
    for (std::vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            for (int i = 0; i < fMCMCNChains; ++i)
                fMCMCTrees[i]->Branch(it->getName().c_str(), &obval[i * kmax + k],
                    (it->getName() + "/D").c_str());
            k++;
            if (it->getDistr().compare("noweight") != 0) {
                for (int i = 0; i < fMCMCNChains; ++i)
                    fMCMCTrees[i]->Branch((it->getName() + "_weight").c_str(),
                                          &obweight[i * kwmax + kweight],
                                          (it->getName() + "_weight/D").c_str());
                kweight++;
            }
        }
    }
}

void MonteCarloEngine::PrintCorrelationMatrix(const std::string filename) 
{
    std::ofstream out;
    out.open(filename.c_str(), std::ios::out);

    int npar = GetNParameters();

    for (int i = 0; i < npar; ++i) 
        out << " & " << GetParameter(i)->GetName();
    out << " \\\\"  << std::endl;

    for (int i = 0; i < npar; ++i) {
        out << GetParameter(i)->GetName() << " & $";
        for (int j = 0; j < npar; ++j) {
            if (i!=j) {
                BCH2D* bch2d_temp = GetMarginalized(GetParameter(i), GetParameter(j));
                if (bch2d_temp != NULL)
                    out << bch2d_temp->GetHistogram()->GetCorrelationFactor();
                else
                    out << 0.;
            }
            else
                out << 1.;            
            if (j==npar-1) out<< "$ \\\\" << std::endl;
            else out << "$ & $";
        }
    }

    out.close();
}
 
