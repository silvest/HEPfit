/* 
 * Copyright (C) 2012 HEPfit Collaboration
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
    obsType = "";
    bin_min = 0.;
    bin_max = 0.;
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
    obsType = orig.obsType;
    bin_min = orig.bin_min;
    bin_min = orig.bin_max;
}

Observable::Observable() 
{
    name = "";
    thname = "";
    label = "";
    min = 0.;
    max = 0.;
    tMCMC = false;
    tho = NULL;
    distr = ""; 
    filename = ""; 
    histoname = "";
    ave = 0.; 
    errg = 0.; 
    errf = 0.;
    obsType = "";
    bin_min = 0.;
    bin_max = 0.;
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

boost::tokenizer<boost::char_separator<char> >::iterator & Observable::ParseObservable(std::string& type, 
                                                                                       boost::tokenizer<boost::char_separator<char> >* tok, 
                                                                                       boost::tokenizer<boost::char_separator<char> >::iterator & beg, 
                                                                                       std::string& filepath, 
                                                                                       int rank) 
{
    obsType = type;
    name = *beg;
    ++beg;
    thname = *beg;
    ++beg;
    label = *beg;
    size_t pos = 0;
    while ((pos = label.find("~", pos)) != std::string::npos)
        label.replace(pos++, 1, " ");
    ++beg;
    min = atof((*beg).c_str());
    ++beg;
    max = atof((*beg).c_str());
    ++beg;
    std::string toMCMC = *beg;
    if (toMCMC.compare("MCMC") == 0)
        tMCMC = true;
    else if (toMCMC.compare("noMCMC") == 0)
        tMCMC = false;
    else
        throw std::runtime_error("ERROR: wrong MCMC flag in " + name);

    if (obsType.compare("Observable") == 0 || obsType.compare("BinnedObservable") == 0) {
        ++beg;
        distr = *beg;
        if (distr.compare("file") == 0) {
            if (std::distance(tok->begin(), tok->end()) < 10)
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                        + *beg + " in " + filename);
            filename = filepath + *(++beg);
            histoname = *(++beg);
            setLikelihoodFromHisto(filename, histoname);
            if (rank == 0) std::cout << "added input histogram " << filename << "/" << histoname << std::endl;
        } else if (distr.compare("weight") == 0) {
            if (std::distance(tok->begin(), tok->end()) < 11 && rank == 0) throw std::runtime_error("ERROR: lack of information on " + *beg + " in " + filename);
            ++beg;
            ave = atof((*beg).c_str());
            ++beg;
            errg = atof((*beg).c_str());
            ++beg;
            errf = atof((*beg).c_str());
            if (errf == 0. && errg == 0.) {
                if (rank == 0) throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + name + " cannot both be 0. in the " + filename + " file.\n");
            }
        } else if (distr.compare("noweight") == 0) {
            if (obsType.compare("BinnedObservable") == 0) {
                ++beg;
                ++beg;
                ++beg;
            }
        } else if (rank == 0)
            throw std::runtime_error("ERROR: wrong distribution flag in " + name);
        ++beg;
        if (obsType.compare("BinnedObservable") == 0) {
            bin_min = atof((*beg).c_str());
            ++beg;
            bin_max = atof((*beg).c_str());
            ++beg;
        }
        if (beg != tok->end() && rank == 0) std::cout << "WARNING: unread information in observable " << name << std::endl;
    }
    return beg;
}
