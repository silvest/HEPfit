/* 
 * File:   MonteCarloEngine.cpp
 * Author: silvest
 * 
 * Created on March 8, 2011, 3:18 PM
 */

#include "MonteCarloEngine.h"
#include <BAT/BCMath.h>
#include <root/TF1.h>
#include <root/TMath.h>
#include <root/TTree.h>
#include <root/TROOT.h>

MonteCarloEngine::MonteCarloEngine(
        const std::vector<ModelParameter>& ModPars_i,
        std::vector<Observable>& Obs_i,
        std::vector<Observable2D>& Obs2D_i) : BCModel(""),
ModPars(ModPars_i), Obs_ALL(Obs_i), Obs2D_ALL(Obs2D_i) {
    Mod = NULL;
};

void MonteCarloEngine::Initialize(Model* Mod_i) {
    Mod = Mod_i;
    int k = 0, kweight = 0;
    for (std::vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            TFile *lik = new TFile((it->getFilename() + ".root").c_str(), "read");
            TH1D *htmp = (TH1D*) (lik->Get(it->getHistoname().c_str()));
            if (htmp == NULL) {
                std::cout << "nonexistent histogram called " + it->getHistoname()
                        + " in " + it->getFilename() + ".root\n";
                exit(EXIT_FAILURE);
            }
            TH1D *inhisto = (TH1D *) htmp->Clone((it->getFilename() + "_" + it->getHistoname()).c_str());
            inhisto->SetDirectory(gROOT);
            InHisto1D[it->getFilename() + it->getHistoname()] = inhisto;
            std::cout << "added input histogram " << InHisto1D[it->getFilename() + it->getHistoname()]->GetName() << std::endl;
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
        }
    }
    obval = new double[fMCMCNChains * k];
    obweight = new double[fMCMCNChains * kweight];
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
            it < Obs2D_ALL.end(); it++) {
        if ((it->getDistr()).compare("file") == 0) {
            TFile *lik2 = new TFile((it->getFilename() + ".root").c_str(), "read");
            TH2D *htmp2 = (TH2D*) (lik2->Get(it->getHistoname().c_str()));
            if (htmp2 == NULL) {
                std::cout << "nonexistent histogram called " + it->getHistoname()
                        + " in " + it->getFilename() + ".root\n";
                exit(EXIT_FAILURE);
            }
            TH2D *inhisto2 = (TH2D *) htmp2->Clone((it->getFilename() + "_" + it->getHistoname()).c_str());
            inhisto2->SetDirectory(gROOT);
            InHisto2D[it->getFilename() + it->getHistoname()] = inhisto2;
            std::cout << "added 2D input histogram " << InHisto2D[it->getFilename() + it->getHistoname()]->GetName() << std::endl;
            it->setMin(inhisto2->GetXaxis()->GetXmin());
            it->setMax(inhisto2->GetXaxis()->GetXmax());
            it->setMin2(inhisto2->GetYaxis()->GetXmin());
            it->setMax2(inhisto2->GetYaxis()->GetXmax());
            lik2->Close();
            delete lik2;
            if (it->isTMCMC()) {
                Obs2D_MCMC.push_back(*it);
            } else {
                std::cout << "cannot handle noMCMC for Observable2D file yet!\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (it->getDistr().compare("weight") == 0) {
            std::cout << "do not use Observable2D for analytic 2D weights!\n";
            exit(EXIT_FAILURE);
        }
        if (Histo2D.find(it->getThname() + "_vs_" + it->getThname2()) == Histo2D.end()) {
            TH2D * histo2 = new TH2D((it->getThname() + "_vs_" + it->getThname2()).c_str(), (it->getLabel() + " vs " + it->getLabel2()).c_str(), NBINS2D,
                    it->getMin(), it->getMax(), NBINS2D, it->getMin2(), it->getMax2());
            histo2->GetXaxis()->SetTitle(it->getLabel().c_str());
            histo2->GetYaxis()->SetTitle(it->getLabel2().c_str());
            BCH2D * bchisto2 = new BCH2D(histo2);
            Histo2D[it->getThname() + "_vs_" + it->getThname2()] = bchisto2;
        }
    }
    DefineParameters();
};

// ---------------------------------------------------------

MonteCarloEngine::~MonteCarloEngine()
// default destructor
{
    delete obval, obweight;
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

void MonteCarloEngine::DefineParameters() {
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        if(it->errf == 0. && it->errg == 0.)
            continue;
        AddParameter(it->name.c_str(), it->min, it->max);
        std::cout << it->name << " " << k << std::endl;
        DPars[it->name]=0.; 
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

double MonteCarloEngine::Weight(const Observable& obs, const double& th) {
    double logprob;
    if (obs.getDistr().compare("weight") == 0) {
        if (obs.getErrf() == 0.)
            logprob = BCMath::LogGaus(th, obs.getAve(), obs.getErrg());
        else if (obs.getErrg() == 0.) {
            std::cout << "cannot use purely flat error in MCMC " << obs.getName() << std::endl;
            exit(EXIT_FAILURE);
        } else
            logprob = log(TMath::Erf((th - obs.getAve() + obs.getErrf()) / sqrt(2.) / obs.getErrg()) -
                TMath::Erf((th - obs.getAve() - obs.getErrf()) / sqrt(2.) / obs.getErrg()));
    } 
    else if (obs.getDistr().compare("file") == 0) {
        TH1D * h = InHisto1D[obs.getFilename()+obs.getHistoname()];
       logprob = log(h->GetBinContent(h->FindBin(th)));
    }
    else {
        std::cout << "Weight called without weight! " << obs.getName() << std::endl;
        exit(EXIT_FAILURE);
    }
    return (logprob);
}

double MonteCarloEngine::Weight(const Observable2D& obs, const double& th1, const double& th2) {
    double logprob;
    if (obs.getDistr().compare("file") == 0) {
        TH2D * h = InHisto2D[obs.getFilename()+obs.getHistoname()];
        logprob = log(h->GetBinContent(h->GetXaxis()->FindBin(th1),h->GetYaxis()->FindBin(th2)));
    }
    else {
        std::cout << "2D Weight called without file! " << obs.getName() << std::endl;
        exit(EXIT_FAILURE);
    }
    return (logprob);
}

double MonteCarloEngine::LogLikelihood(std::vector <double> parameters) {
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.

    double logprob = 0.;

    for(int k=0; k<parameters.size(); k++) {
//        std::string pippo = GetParameter(k)->GetName();
//        double pluto = parameters[k];
//        DPars[pippo]=pluto;
        DPars[GetParameter(k)->GetName()]=parameters[k];
    }

    Mod->update(DPars);
    // std::cout << "loglike" << parameters[0];
    for (std::vector<Observable>::iterator it = Obs_MCMC.begin(); it < Obs_MCMC.end(); it++) {
        double th = it->getTheoryValue();
        logprob += Weight(*it, th);
    }

    for (std::vector<Observable2D>::iterator it = Obs2D_MCMC.begin(); it < Obs2D_MCMC.end(); it++) {
        double th1 = it->getTheoryValue();
        double th2 = it->getTheoryValue2();
        logprob += Weight(*it, th1, th2);
    }

    return logprob;
}

void MonteCarloEngine::MCMCIterationInterface() {
    // get number of chains
    int nchains = MCMCGetNChains();

    // get number of parameters
    int npar = GetNParameters();

    for (int i = 0; i < nchains; ++i) {
        for (int k = 0; k < npar; k++) {
            //        std::string pippo = GetParameter(k)->GetName();
            //        double pluto = parameters[k];
            //        DPars[pippo]=pluto;
            DPars[GetParameter(k)->GetName()] = fMCMCx.at(i * npar + k);
        }

        Mod->update(DPars);

        // fill the histograms for all observables
        int k = 0, kweight = 0;
        for (std::vector<Observable>::iterator it = Obs_ALL.begin();
                it < Obs_ALL.end(); it++) {
            double th = it->getTheoryValue();
            Histo1D[it->getThname()]->GetHistogram()->Fill(th);
            if (!it->isTMCMC()) {
                obval[(i + 1)*(k + 1) - 1] = th;
                k++;
                if (it->getDistr().compare("noweight") != 0) {
                    obweight[(i + 1)*(kweight + 1) - 1] = Weight(*it, th);
                    kweight++;
                }
            }
        }
        for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin();
                it < Obs2D_ALL.end(); it++) {
            double th1 = it->getTheoryValue();
            double th2 = it->getTheoryValue2();
            Histo2D[it->getThname() + "_vs_" + it->getThname2()]->GetHistogram()->Fill(th1,th2);
        }
    }
}


void MonteCarloEngine::PrintHistogram(BCModelOutput& out) {
    //print the BAT histograms to an eps file
    for (std::vector<Observable>::iterator it = Obs_ALL.begin(); it < Obs_ALL.end();
            it++) {
//        std::string fname = "Observables/" + it->name + ".pdf";
//        Histo1D[it->name]->Print(fname.c_str());
        out.Write(Histo1D[it->getThname()]->GetHistogram());
    }
    for (std::vector<Observable2D>::iterator it = Obs2D_ALL.begin(); it < Obs2D_ALL.end();
            it++) {
//        std::string fname = "Observables/" + it->name + ".pdf";
//        Histo1D[it->name]->Print(fname.c_str());
        out.Write(Histo2D[it->getThname() + "_vs_" + it->getThname2()]->GetHistogram());
    }
}

void MonteCarloEngine::AddChains() {
    fMCMCFlagWritePreRunToFile = false;
    int k = 0, kweight = 0;
    for (std::vector<Observable>::iterator it = Obs_ALL.begin();
            it < Obs_ALL.end(); it++) {
        if (!it->isTMCMC()) {
            for (int i = 0; i < fMCMCNChains; ++i)
                fMCMCTrees[i]->Branch(it->getName().c_str(), &obval[(i + 1)*(k + 1) - 1],
                    (it->getName() + "/D").c_str());
            k++;
            if (it->getDistr().compare("noweight") == 0) {
                for (int i = 0; i < fMCMCNChains; ++i)
                    fMCMCTrees[i]->Branch((it->getName() + "_weight").c_str(),
                        &obweight[(i + 1)*(kweight + 1) - 1], (it->getName() + "_weight/D").c_str());
                kweight++;
            }
        }
    }
}
