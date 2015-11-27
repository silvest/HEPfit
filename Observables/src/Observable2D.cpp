/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Observable2D.h"
#include <TFile.h>
#include <fstream>
#include <TROOT.h>
#include <limits>

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
: Observable(name_i, thname_i, label_i, tMCMC_i, min_i, max_i, tho_i), 
  bin_min(2, 0.),
  bin_max(2, 0.)
{
    thname2 = thname2_i;
    label2 = label2_i;
    min2 = min2_i;
    max2 = max2_i;
    tho2 = tho2_i;
    obsType2 = "";
    filepath = "";
    iterationNo2 = std::numeric_limits<int>::max();
}

Observable2D::Observable2D(const Observable& o1d)
: Observable(o1d),
  bin_min(2, 0.),
  bin_max(2, 0.)
{
    thname2 = "";
    label2 = "";
    min2 = 0.;
    max2 = 0.;
    tho2 = NULL;
    obsType2 = "";
    filepath = "";
    iterationNo2 = std::numeric_limits<int>::max();
}

Observable2D::Observable2D()
: Observable(),
  bin_min(2, 0.),
  bin_max(2, 0.)
{
    thname2 = "";
    label2 = "";
    min2 = 0.;
    max2 = 0.;
    tho2 = NULL;
    obsType2 = "";
    filepath = "";
    iterationNo2 = std::numeric_limits<int>::max();
}

Observable2D::Observable2D(const Observable2D& orig)
: Observable(orig.name, orig.thname, orig.label, orig.tMCMC, orig.min, orig.max, orig.tho),
  bin_min(orig.bin_min),
  bin_max(orig.bin_max)
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
    filepath = orig.filepath;
    iterationNo2 = orig.iterationNo2;
}

Observable2D::~Observable2D()
{}

double Observable2D::computeTheoryValue2()
{
   if (tho2->getModel().getIterationNo() == iterationNo2) {
        return thValue2;
    } else {
       iterationNo2 = tho2->getModel().getIterationNo();
        thValue2 = tho2->computeThValue();
        return thValue2;
    }
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

int Observable2D::ParseObservable2D(std::string& type, 
                                     boost::tokenizer<boost::char_separator<char> >* tok, 
                                     boost::tokenizer<boost::char_separator<char> >::iterator& beg, 
                                     std::string& infilename, 
                                     std::ifstream& ifile,
                                     int lineNo,
                                     int rank)
{
    if (infilename.find("\\/") == std::string::npos) filepath = infilename.substr(0, infilename.find_last_of("\\/") + 1);
    if (std::distance(tok->begin(), tok->end()) < 12) {
        setName(*beg);
        ++beg;
        if (std::distance(tok->begin(), tok->end()) < 4) {
            if(rank == 0) throw std::runtime_error("ERROR: lack of information on " + name + " in " + infilename + " at line number" + boost::lexical_cast<std::string>(lineNo));
            else sleep (2);
        }
        std::string toMCMC = *beg;
        if (toMCMC.compare("MCMC") == 0)
            setTMCMC(true);
        else if (toMCMC.compare("noMCMC") == 0)
            setTMCMC(false);
        else {
            if (rank == 0) throw std::runtime_error("ERROR: wrong MCMC flag in Observable2D" + name + " at line number:" + boost::lexical_cast<std::string>(lineNo) + " of file " + infilename + ".\n");
            else sleep(2);
        }
        
        ++beg;
        setDistr(*beg);
        if (distr.compare("file") == 0) {
            if (std::distance(tok->begin(), tok->end()) < 6) {
                if(rank == 0) throw std::runtime_error("ERROR: lack of information on "+ *beg + " in " + infilename);
                else sleep (2);
            }
            setFilename(filepath + *(++beg));
            setHistoname(*(++beg));
        }

        std::vector<double> min(2, 0.);
        std::vector<double> max(2, 0.);
        std::vector<double> ave(2, 0.);
        std::vector<double> errg(2, 0.);
        std::vector<double> errf(2, 0.);
        std::vector<std::string> thname(2, "");
        std::vector<std::string> label(2, "");
        std::vector<std::string> type2D(2, "");
        std::string line;
        size_t pos = 0;
        boost::char_separator<char> sep(" \t");
        for (int i = 0; i < 2; i++) {
            IsEOF = getline(ifile, line).eof();
            if (line.empty() || line.at(0) == '#') {
                if (rank == 0) throw std::runtime_error("ERROR: no comments or empty lines in Observable2D please! In file " + infilename + " at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n");
                else sleep(2);
            }
            lineNo++;
            boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
            beg = mytok.begin();
            type2D[i] = *beg;
            if (type2D[i].compare("Observable") != 0 && type2D[i].compare("BinnedObservable") != 0 && type2D[i].compare("FunctionObservable") != 0) {
                if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + infilename + ", expecting an Observable or BinnedObservable or FunctionObservable type here...\n");
                else sleep(2);
            }
            ++beg;
            thname[i] = *beg;
            ++beg;
            label[i] = *beg;
            while ((pos = label[i].find("~", pos)) != std::string::npos)
                label[i].replace(pos++, 1, " ");
            ++beg;
            min[i] = atof((*beg).c_str());
            ++beg;
            max[i] = atof((*beg).c_str());
            if (distr.compare("weight") == 0) {
                ++beg;
                ave[i] = atof((*beg).c_str());
                ++beg;
                errg[i] = atof((*beg).c_str());
                ++beg;
                errf[i] = atof((*beg).c_str());
                if (errg[i] == 0. && errg[i] == 0.) {
                    if (rank == 0) throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + name + " cannot both be 0. in the " + infilename + " file, line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n");
                    else sleep(2);
                }
            } else if (distr.compare("noweight") == 0 || distr.compare("file") == 0) {
                if (type2D[i].compare("BinnedObservable") == 0 || type2D[i].compare("FunctionObservable") == 0) {
                    ++beg;
                    ++beg;
                    ++beg;
                }
            } else {
                if (rank == 0) throw std::runtime_error("ERROR: wrong distribution flag in " + name + " in file " + infilename + ".\n");
                else sleep(2);
            }
            if (type2D[i].compare("BinnedObservable") == 0) {
                ++beg;
                bin_min[i] = atof((*beg).c_str());
                ++beg;
                bin_max[i] = atof((*beg).c_str());
            } else if (type2D[i].compare("FunctionObservable") == 0) {
                ++beg;
                bin_min[i] = atof((*beg).c_str());
                ++beg;
            }
        }
        setObsType(type2D[0]);
        obsType2 = type2D[1];
        setThname(thname[0]);
        thname2 = thname[1];
        setLabel(label[0]);
        label2 = label[1];
        setMin(min[0]);
        min2 = min[1];
        setMax(max[0]);
        max2= max[1];
        setAve(ave[0]);
        ave2 = ave[1];
        setErrg(errg[0]);
        errg2 = errg[1];
        setErrf(errf[0]);
        errf2 = errf[1];
        if (distr.compare("file") == 0) {
            setLikelihoodFromHisto(filename, histoname);
            if (rank == 0) std::cout << "added input histogram " << filename << "/" << histoname << std::endl;
        }
        return lineNo;
    } else {
        beg = ParseObservable(type, tok, beg, filepath, filename, rank);
        ++beg;
        std::string distr = *beg;
        if (distr.compare("file") == 0) {
            if (std::distance(tok->begin(), tok->end()) < 14) {
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on " + *beg + " in " + infilename + ".\n");
                else sleep(2);
            }
            setFilename(filepath + *(++beg));
            setHistoname(*(++beg));
            setLikelihoodFromHisto(filename, histoname);
            if (rank == 0) std::cout << "added input histogram " << filename << "/" << histoname << std::endl;
        } else if (distr.compare("noweight") == 0) {
        } else { 
            if (rank == 0) throw std::runtime_error("ERROR: wrong distribution flag in " + name);
            else sleep(2);
        }
        setDistr(distr);
        ++beg;
        thname2 = *beg;
        ++beg;
        std::string label = *beg;
        size_t pos = 0;
        while ((pos = label.find("~", pos)) != std::string::npos)
            label.replace(pos, 1, " ");
        label2 = label;
        ++beg;
        min2 = atof((*beg).c_str());
        ++beg;
        max2 = atof((*beg).c_str());
        ++beg;
        if (beg != tok->end())
            if (rank == 0) std::cout << "WARNING: unread information in observable2D " << name << std::endl;
        return lineNo;
    }
}
