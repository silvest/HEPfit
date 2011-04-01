/* 
 * File:   MonteCarlo.cpp
 * Author: silvest
 * 
 * Created on March 8, 2011, 3:18 PM
 */

#include "MonteCarlo.h"

#include <BCMath.h>
#include <root/TF1.h>
#include <root/TMath.h>

// ---------------------------------------------------------

MonteCarlo::MonteCarlo(const char * name, Model * Mod_i, 
        std::vector<ModelParameter> ModPars_i, std::vector<Observable> Obs_i) : BCModel(name) {
    // default constructor
    Mod = Mod_i;
    ModPars = ModPars_i;
    for (std::vector<Observable>::iterator it = Obs_i.begin(); it < Obs_i.end(); it++) {
        if (it->tMCMC) {
            if((it->distr).compare("file")==0) {
                TFile *lik = new TFile((it->filename+".root").c_str(),"read");
                if(lik->Get(it->histoname.c_str())==NULL){
                    std::cout << "nonexistent histogram called " + it->histoname + " in " + it->filename+".root\n";
                    exit(EXIT_FAILURE);
                }
                TH1D * inhisto = (TH1D*) lik->Get(it->histoname.c_str())->Clone();
                InHisto1D[it->histoname]=inhisto;
                it->min = inhisto->GetXaxis()->GetXmin();
                it->max = inhisto->GetXaxis()->GetXmax();
            }
            Obs_MCMC.push_back(*it);
        }
        else {
            Obs_noMCMC.push_back(*it);
        }
        TH1D * histo = new TH1D(it->name.c_str(), it->name.c_str(), NBINS1D, it->min, it->max);
        BCH1D * bchisto = new BCH1D();
        bchisto->SetHistogram(histo);
        Histo1D[it->name] = bchisto;
    }
    DefineParameters();
};

// ---------------------------------------------------------

MonteCarlo::~MonteCarlo()
// default destructor
{
};

// ---------------------------------------------------------

void MonteCarlo::DefineParameters() {
    // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    int k = 0;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        if(it->errf == 0. && it->errg == 0.)
            continue;
        AddParameter(it->name.c_str(), it->min, it->max);
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
        }
        k++;
    }
}

// ---------------------------------------------------------

double MonteCarlo::LogLikelihood(std::vector <double> parameters) {
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
        if (it->distr.compare("weight")==0)
            if (it->errf == 0.)
                logprob += BCMath::LogGaus(th, it->ave, it->errg);
            else if (it->errg == 0.) {
                std::cout << "cannot use purely flat error in MCMC" << it->name << std::endl;
                exit(EXIT_FAILURE);
            } else
                logprob += log(TMath::Erf((th - it->ave + it->errf) / sqrt(2.) / it->errg) -
                    TMath::Erf((th - it->ave - it->errf) / sqrt(2.) / it->errg));
        else if(it->distr.compare("file")==0){
            TH1D * h = InHisto1D[it->histoname];
            logprob += log(h->GetBinContent(h->FindBin(th)));
        }
    }

    return logprob;
}

//void MonteCarlo::MCMCIterationInterface()
//{
//  // get number of chains
//  int nchains = MCMCGetNChains();
//
//  // get number of parameters
//  int npar = GetNParameters();
//
//  // loop over all chains and fill histogram
//  for (int i = 0; i < nchains; ++i) {
//    // get the current values of the parameters x and y. These are
//    // stored in fMCMCx. 
//    double x = fMCMCx.at(i * npar + 0); 
//    double y = fMCMCx.at(i * npar + 1); 
//
//    // fill the ratio histogram
//    fHistRatio->GetHistogram()->Fill(x/y);
//  }
//}
void MonteCarlo::MCMCCurrentPointInterface(std::vector <double> & point, int ichain, bool accepted)
{
    //std::cout << " currentpoint " << point[0] << std::endl;
}

void MonteCarlo::PrintHistogram()
{
  // print the BAT histogram to an eps file
  //fHistRatio->Print("ratio.eps");
}
