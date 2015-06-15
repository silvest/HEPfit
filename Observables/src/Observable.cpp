/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Observable.h"
#include <TNamed.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>


Observable::Observable (const std::string name_i,
                        const std::string thname_i,
                        const std::string label_i,
                        const bool tMCMC_i,
                        const double min_i,
                        const double max_i,
                        ThObservable * tho_i)
{
    name = name_i;
    thname = thname_i;
    label = label_i;
    min = min_i;
    max = max_i;
    tMCMC = tMCMC_i;
    tho = tho_i;
    distr = "";
    filename = "";
    histoname = "";
    ave = 0.;
    errg = 0.;
    errf = 0.;
    obsType = 0;
}

Observable::Observable(const Observable& orig) 
{
    name = orig.name;
    thname = orig.thname;
    label = orig.label;
    min = orig.min;
    max = orig.max;
    tMCMC = orig.tMCMC;
    tho = orig.tho;
    distr = orig.distr; 
    filename = orig.filename; 
    histoname = orig.histoname;
    ave = orig.ave; 
    errg = orig.errg; 
    errf = orig.errf;
}

Observable::~Observable() {}

std::ostream& operator<<(std::ostream& output, const Observable& o)
{
    output << "Observable name, tMCMC, min, max, distribution, distribution parameters" << std::endl;
    output << o.name << " " << o.tMCMC << " " << o.min << " " << o.max << " " 
           << o.distr << " " << o.filename << " " << o.histoname << " " << o.ave 
           << " " << o.errg << " " << o.errf << std::endl;
    return output;
}

void Observable::setLikelihoodFromHisto(std::string filename, std::string histoname)
    {
        this->filename = filename;
        this->histoname = histoname;
        TFile *lik = new TFile((filename + ".root").c_str(), "read");
        TH1D *htmp = (TH1D*) (lik->Get(histoname.c_str()));
        if (htmp == NULL)
            throw std::runtime_error("ERROR: nonexistent histogram called "
                    + histoname + " in "
                    + filename + ".root");
        inhisto = (TH1D *) htmp->Clone((filename + "/" + histoname).c_str());
        inhisto->SetDirectory(gROOT);
        std::cout << "added input histogram " << inhisto->GetName() << std::endl;
        setMin(inhisto->GetXaxis()->GetXmin());
        setMax(inhisto->GetXaxis()->GetXmax());
        lik->Close();
        delete lik;
    }


double Observable::computeTheoryValue()
{
    return tho->computeThValue();
}

double Observable::LogSplitGaussian(double x, double ave, double errl, double errr)
{
    double sigma = (x > ave ? errr : errl);
    return -0.5 * (x-ave) * (x-ave) / sigma / sigma;
}

double Observable::LogGaussian(double x, double ave, double sigma)
{
    return -0.5 * (x-ave) * (x-ave) / sigma / sigma;
}

double Observable::computeWeight(double th)
{
    double logprob;
    if (distr.compare("weight") == 0) {
        if (errf == 0.)
            logprob = LogGaussian(th, ave, errg);
        else if (errg == 0.) {
            if (th < ave + errf && th > ave - errf)
                logprob = 1.;
            else
                logprob = log(0.);
        } else
            logprob = log(TMath::Erf((th - ave + errf) / sqrt(2.) / errg)
                - TMath::Erf((th - ave - errf) / sqrt(2.) / errg));
    } else if (distr.compare("file") == 0) {
        int i = inhisto->FindBin(th);
        if (inhisto->IsBinOverflow(i) || inhisto->IsBinUnderflow(i))
            logprob = log(0.);
        else
            logprob = log(inhisto->GetBinContent(i));
        //logprob = log(h->GetBinContent(h->FindBin(th)));
    } else
        throw std::runtime_error("ERROR: MonteCarloEngine::Weight() called without weight for "
            + name);
    return (logprob);
}

double Observable::computeWeight(double th, double ave_i, double errg_i, double errf_i)
{
    double logprob;
    if (distr.compare("weight") == 0) {
        if (errf_i == 0.)
            logprob = LogGaussian(th, ave_i, errg_i);
        else if (errg_i == 0.) {
            if (th < ave_i + errf_i && th > ave_i - errf_i)
                logprob = 1.;
            else
                logprob = log(0.);
        } else
            logprob = log(TMath::Erf((th - ave_i + errf_i) / sqrt(2.) / errg_i)
                - TMath::Erf((th - ave_i - errf_i) / sqrt(2.) / errg_i));
    } else
        throw std::runtime_error("ERROR: MonteCarloEngine::Weight() called without weight for "
            + name);
    return (logprob);
}
