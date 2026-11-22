/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsObservable.h"
#include "ThObsFactory.h"
#include <TNamed.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <string>
#include <sstream>

namespace {
// Relative file names are relative to the directory of the config file (filepath), as for histogram files.
std::string configRelativePath(const std::string& filepath, const std::string& name)
{
    return (!name.empty() && name[0] == '/') ? name : filepath + name;
}
}

HiggsObservable::HiggsObservable(const Observable& Obs) :
Observable(Obs)
{
    isnew = true;
    thobsvsize = 0.;
    channelsize = 0.;
    isCorrelated = false;
    covarianceFromConfig = false;
    InvCov = NULL;
    rank = 0;
};

HiggsObservable::HiggsObservable(const HiggsObservable& orig) :
Observable(orig), channels(orig.channels), isnew(orig.isnew), thObsV(orig.thObsV), theoryValues(orig.theoryValues),
thobsvsize(orig.thobsvsize), channelsize(orig.channelsize), isCorrelated(orig.isCorrelated),
covarianceFromConfig(orig.covarianceFromConfig),
InvCov(orig.InvCov != NULL ? new gslpp::matrix<double>(*orig.InvCov) : NULL), // each copy owns its matrix
rank(orig.rank)
{
}

HiggsObservable::~HiggsObservable()
{
    delete InvCov;
}

void HiggsObservable::setParametricLikelihood(std::string filename, std::vector<ThObservable*> thObsV)
{
    this->filename = filename;
    this->thObsV = thObsV;
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid Higgs parameters configuration file.\n");
    std::string line;
    bool IsEOF;
    bool readCorrelations = false;
    unsigned int i = 0; 
    unsigned int nrows = 0;
    do {
        IsEOF = getline(ifile, line).eof();
        if (line.compare(0, 11, "MEASUREMENT") == 0) nrows++;
    } while (!IsEOF);
    if (nrows == 0) {
        std::stringstream ss;
        ss << "\nERROR: " << filename << " should contain at least 1 measurement.\n";
        throw std::runtime_error(ss.str());
    }
    const unsigned int meanColumn = thObsV.size() - (isnew ? 0 : 1);
    const unsigned int errorColumn = meanColumn + 2;
    channels.ResizeTo(nrows, meanColumn + 3);
    ifile.clear();
    ifile.seekg(0, ifile.beg);
    int lineNo = 0;
    do {
        IsEOF = getline(ifile, line).eof();
        lineNo++;
        if (!line.empty() && line.back() == '\r') line.pop_back(); // for CR+LF
        // Only MEASUREMENT lines are measurements: comments and blank lines between them are skipped also
        // for a correlated observable (whose correlation rows are read right after the last measurement).
        if (line.compare(0, 11, "MEASUREMENT") != 0) {
            continue;
        } else {
            boost::char_separator<char> sep(" |\t");
            boost::tokenizer<boost::char_separator<char> > tok(line, sep);
            boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
            if (isCorrelated) readCorrelations = true;
            // Read the necessary information from the config file. Each row contains:
            // ggH fraction
            // VBF fraction
            // VH fraction (ttH is computed as 1-ggH-VBF-VH)
            // average value of mu
            // left-side error
            // right-side error
            beg++; // skip label
            for (int j = 0; j < channels.GetNcols(); j++)
                channels(i, j) = atof((*(++beg)).c_str());
            if (isCorrelated && (std::fabs(channels(i, meanColumn + 1)) != std::fabs(channels(i, errorColumn)))) {
                if (rank == 0) throw std::runtime_error("ERROR: Gaussian errors must be symmetric for CorrelatedHiggsObservables in " + filename + " at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n");
                else sleep(2);
            }
            i++;
            if (readCorrelations && (i == nrows)) {
                gslpp::matrix<double> myCorr(gslpp::matrix<double>::Id(nrows));
                int ni = 0;
                for (unsigned int irow = 0; irow < nrows; irow++) {
                    IsEOF = getline(ifile, line).eof();
                    if (line.empty() || line.at(0) == '#') {
                        if (rank == 0) throw std::runtime_error("ERROR: no comments or empty lines in CorrelatedHiggsObservables please! In file " + filename + " at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n");
                        else sleep(2);
                    }
                    lineNo++;
                    boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                    beg = mytok.begin();
                    int nj = 0;
                    for (unsigned int jcol = 0; jcol < nrows; jcol++) {
                        if ((*beg).compare(0, 1, "0") == 0
                                || (*beg).compare(0, 1, "1") == 0
                                || (*beg).compare(0, 1, "-") == 0 || (covarianceFromConfig == true)) {
                            if (std::distance(mytok.begin(), mytok.end()) < nrows) {
                                if (rank == 0) throw std::runtime_error(("ERROR: Correlation matrix is of wrong size in CorrelatedHiggsObservables: " + name + +" at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n").c_str());
                                else sleep(2);
                            }
                            myCorr(ni, nj) = atof((*beg).c_str());
                            nj++;
                            beg++;
                        } else {
                            if (rank == 0) throw std::runtime_error("ERROR: invalid correlation matrix for " + name + ". Check element (" + boost::lexical_cast<std::string>(ni + 1) + "," + boost::lexical_cast<std::string>(nj + 1) + ") in line number " + boost::lexical_cast<std::string>(lineNo) + " in file " + filename + ".\n");
                            else sleep(2);
                        }
                    }
                    ni++;
                }
                if (myCorr.size_i() != nrows || myCorr.size_j() != nrows)
                    throw std::runtime_error("The size of the correlated observables in " + name + " does not match the size of the correlation matrix!");
                delete InvCov;
                InvCov = new gslpp::matrix<double>(nrows, nrows, 0.);
                if (covarianceFromConfig) {
                    for (unsigned int i = 0; i < nrows; i++)
                        for (unsigned int j = 0; j < nrows; j++)
                            (*InvCov)(i, j) = myCorr(i, j);
                } else {
                    for (unsigned int i = 0; i < nrows; i++)
                        for (unsigned int j = 0; j < nrows; j++)
                            (*InvCov)(i, j) = std::fabs(channels(i, errorColumn)) * myCorr(i, j) * std::fabs(channels(j, errorColumn));
                    *InvCov = InvCov->inverse();
                }
                readCorrelations = false;
            }
        }
    } while (!IsEOF);
    if (i != nrows) {
        std::stringstream ss;
        ss << "\nERROR: " << filename << " should contain " << nrows << " measurements, but I have read " << i << " ones instead.\n";
        throw std::runtime_error(ss.str());
    }

    thobsvsize = thObsV.size();
    channelsize = channels.GetNrows();
    theoryValues.resize(thobsvsize + 1 + channelsize);
}

double HiggsObservable::computeWeight()
{
    double logprob = 0;
    const unsigned int meanColumn = thobsvsize - (isnew ? 0 : 1);

    if (isnew) {
        for (int i = 0; i < channels.GetNrows(); i++) {
            double mu = 0;
            for (unsigned int j = 0; j < thobsvsize; j++) {
                theoryValues.at(j) = thObsV.at(j)->computeThValue();
                mu += channels(i, j) * theoryValues.at(j);
            }
            theoryValues.at(thobsvsize) = tho->computeThValue();
            if (thObsV.at(0)->getModel().isModelLinearized())
                theoryValues.at(thobsvsize + i + 1) = -1. + mu + theoryValues.at(thobsvsize);
            else
                theoryValues.at(thobsvsize + i + 1) = mu * theoryValues.at(thobsvsize);
        }
    } else {
        for (int i = 0; i < channels.GetNrows(); i++) {
            double mu = 0, sum = 0.;
            for (unsigned int j = 0; j < thobsvsize - 1; j++) {
                theoryValues.at(j) = thObsV.at(j)->computeThValue();
                mu += channels(i, j) * theoryValues.at(j);
                sum += channels(i, j);
            }
            mu += (1. - sum) * thObsV.at(thobsvsize - 1)->computeThValue();
            theoryValues.at(thobsvsize) = tho->computeThValue();
            theoryValues.at(thobsvsize + i + 1) = mu * theoryValues.at(thobsvsize);
        }
    }

    if (isCorrelated) {
        gslpp::vector<double> residuals(channels.GetNrows(), 0.);
        for (int i = 0; i < channels.GetNrows(); i++)
            residuals(i) = theoryValues.at(thobsvsize + i + 1) - channels(i, meanColumn);
        return -0.5 * residuals * ((*InvCov) * residuals);
    }
    for (int i = 0; i < channels.GetNrows(); i++)
        logprob += LogSplitGaussian(theoryValues.at(thobsvsize + i + 1), channels(i, meanColumn),
                                    channels(i, meanColumn + 1), channels(i, meanColumn + 2));
    return (logprob);
}

boost::tokenizer<boost::char_separator<char> >::iterator & HiggsObservable::ParseHiggsObservable(
        boost::tokenizer<boost::char_separator<char> >::iterator & beg,
        ThObsFactory& myObsFactory,
        StandardModel * myModel,
        int rank,
        const std::string& filepath)
{
    this->rank = rank;
    std::string type = "HiggsObservable";
    setObsType(type);
    std::vector<ThObservable*> hthobs;
    ++beg;
    distr = *beg;
    if (distr.compare("parametric") == 0) {
        setIsnew(false);
        ++beg;
        distr = *beg;
        if (distr.compare("LHC7") == 0) {
            hthobs.push_back(myObsFactory.CreateThMethod("ggH7", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("VBF7", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("VH7", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("ttH7", *myModel));
        } else if (distr.compare("LHC8") == 0) {
            hthobs.push_back(myObsFactory.CreateThMethod("ggH8", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("VBF8", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("VH8", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("ttH8", *myModel));
        } else if (distr.compare("TeV196") == 0) {
            hthobs.push_back(myObsFactory.CreateThMethod("ggH196", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("VBF196", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("VH196", *myModel));
            hthobs.push_back(myObsFactory.CreateThMethod("ttH196", *myModel));
        } else if (rank == 0)
            throw std::runtime_error("ERROR: wrong keyword " + distr + " in " + getName());
        ++beg;
        setParametricLikelihood(configRelativePath(filepath, *beg), hthobs);
    } else if (distr.compare("new_parametric") == 0) {
        ++beg;
        std::string filename_h = configRelativePath(filepath, *beg);
        std::ifstream ifile(filename_h.c_str());
        if (!ifile.is_open()) {
            if (rank == 0) throw std::runtime_error("\nERROR: " + filename_h + " does not exist. Make sure to specify a valid Higgs parameters configuration file.\n");
            else sleep(2);
        }
        std::string line_h;
        bool IsEOF_h;
        do {
            IsEOF_h = getline(ifile, line_h).eof();
            if (line_h.compare(0, 10, "CATEGORIES") == 0) {
                boost::char_separator<char> sep(" \t");
                boost::tokenizer<boost::char_separator<char> > mytok(line_h, sep);
                boost::tokenizer<boost::char_separator<char> >::iterator beg2 = mytok.begin();
                beg2++;
                while (beg2 != mytok.end()) {
                    std::string cat = *beg2;
                    hthobs.push_back(myObsFactory.CreateThMethod(cat, *myModel));
                    beg2++;
                }
            }
        } while (!IsEOF_h);
        if (hthobs.size() > 0)
            setParametricLikelihood(filename_h, hthobs);
        else {
            if (rank == 0) throw std::runtime_error("\nERROR: " + getName() + " does not provide at least one category\n");
            else sleep(2);
        }
    } else {
        if (rank == 0) throw std::runtime_error("ERROR: wrong distribution flag " + distr + " in " + getName());
        else sleep(2);
    }

    return beg;
}

