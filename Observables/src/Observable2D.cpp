/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Observable2D.h"
#include <TFile.h>
#include <TROOT.h>

Observable2D::Observable2D(const std::string name_i,
        const std::string thname_i,
        const std::string thname2_i,
        const std::string label_i,
        const std::string label2_i,
        const bool tMCMC_i,
        const double min_i,
        const double max_i,
        const double min2_i,
        const double max2_i,
        ThObservable * tho_i,
        ThObservable * tho2_i)
: Observable(name_i, thname_i, label_i, tMCMC_i, min_i, max_i, tho_i)
{
    thname2 = thname2_i;
    label2 = label2_i;
    min2 = min2_i;
    max2 = max2_i;
    tho2 = tho2_i;
}

Observable2D::Observable2D(const Observable& o1d)
: Observable(o1d)
{
    thname2 = "";
    label2 = "";
    min2 = 0.;
    max2 = 0.;
    tho2 = NULL;
}

Observable2D::Observable2D(const Observable2D& orig)
: Observable(orig.name, orig.thname, orig.label, orig.tMCMC, orig.min, orig.max, orig.tho)
{

    distr = orig.distr;
    filename = orig.filename;
    histoname = orig.histoname;
    ave = orig.ave;
    errg = orig.errg;
    errf = orig.errf;

    thname2 = orig.thname2;
    label2 = orig.label2;
    min2 = orig.min2;
    max2 = orig.max2;
    tho2 = orig.tho2;
}

Observable2D::~Observable2D()
{
}

double Observable2D::computeTheoryValue2()
{
    return tho2->computeThValue();
}

void Observable2D::setLikelihoodFromHisto(std::string filename, std::string histoname)
{
    this->filename = filename;
    this->histoname = histoname;
    TFile *lik2 = new TFile((filename + ".root").c_str(), "read");
    TH2D *htmp2 = (TH2D*) (lik2->Get(histoname.c_str()));
    if (htmp2 == NULL)
        throw std::runtime_error("ERROR: nonexistent histogram called "
            + histoname
            + " in " + filename + ".root");
    inhisto2d = (TH2D *) htmp2->Clone((filename + "/" + histoname).c_str());
    inhisto2d->SetDirectory(gROOT);
    std::cout << "added 2D input histogram " << name << std::endl;
    setMin(inhisto2d->GetXaxis()->GetXmin());
    setMax(inhisto2d->GetXaxis()->GetXmax());
    setMin2(inhisto2d->GetYaxis()->GetXmin());
    setMax2(inhisto2d->GetYaxis()->GetXmax());
    lik2->Close();
    delete lik2;
}

double Observable2D::computeWeight(double th1, double th2)
{
    double logprob;
    if (distr.compare("file") == 0) {
        int i = inhisto2d->FindBin(th1, th2);
        if (inhisto2d->IsBinOverflow(i) || inhisto2d->IsBinUnderflow(i))
            logprob = log(0.);
        else
            logprob = log(inhisto2d->GetBinContent(i));
        //logprob = log(h->GetBinContent(h->GetXaxis()->FindBin(th1),h->GetYaxis()->FindBin(th2)));
    } else if (distr.compare("weight") == 0) {
        logprob = Observable::computeWeight(th1) + Observable::computeWeight(th2, ave2, errg2, errf2);
    } else
        throw std::runtime_error("ERROR: 2D MonteCarloEngine::Weight() called without file for "
            + name);
    return (logprob);
}
