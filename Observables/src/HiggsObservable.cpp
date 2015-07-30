/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsObservable.h"
#include <TNamed.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <string>
#include <sstream>

//HiggsObservable::~HiggsObservable() {}

void HiggsObservable::setParametricLikelihood(std::string filename, std::vector<ThObservable*> thObsV)
{
    this->filename = filename;
    this->thObsV = thObsV;
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid Higgs parameters configuration file.\n");
    std::string line;
    bool IsEOF;
    int i = 0, nrows = 0;
    do {
        IsEOF = getline(ifile, line).eof();
        if (line.compare(0, 11, "MEASUREMENT") == 0) nrows++;
    } while (!IsEOF);
    if (nrows == 0) {
        std::stringstream ss;
        ss << "\nERROR: " << filename << " should contain at least 1 measurement.\n";
        throw std::runtime_error(ss.str());
    }
    int implicit_tth = (isnew ? 0 : -1);
    channels.ResizeTo(nrows, thObsV.size() + 3 + implicit_tth);
    ifile.clear();
    ifile.seekg(0, ifile.beg);
    do {
        IsEOF = getline(ifile, line).eof();
        if (*line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
        if (line.compare(0, 11, "MEASUREMENT") != 0)
            continue;
        boost::char_separator<char> sep(" |");
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
        // Read the necessary information from the config file. Each row contains:
        // ggH fraction
        // VBF fraction
        // VH fraction (ttH is computed as 1-ggH-VBF-VH)
        // average value of mu
        // left-side error
        // right-side error
        beg++; // skip label
        for (unsigned int j = 0; j < channels.GetNcols(); j++) 
            channels(i, j) = atof((*(++beg)).c_str());
        i++;

    } while (!IsEOF);
    if (i != nrows) {
        std::stringstream ss;
        ss << "\nERROR: " << filename << " should contain " << nrows << " measurements, but I have read " << i << " ones instead.\n";
        throw std::runtime_error(ss.str());
    }

}

double HiggsObservable::computeWeight()
{
    double logprob = 0;

    double thobsvsize = thObsV.size();
    
    if (isnew) {
        for (int i = 0; i < channels.GetNrows(); i++) {
            double mu = 0;
            for (unsigned int j = 0; j < thobsvsize; j++)
                mu += channels(i, j) * thObsV.at(j)->computeThValue();
            logprob += LogSplitGaussian(mu * tho->computeThValue(), channels(i, thobsvsize), channels(i, thobsvsize + 1), channels(i, thobsvsize + 2));
        }
      

    } else {
        for (int i = 0; i < channels.GetNrows(); i++) {
            double mu = 0, sum = 0.;
            for (unsigned int j = 0; j < thobsvsize - 1; j++) {
                mu += channels(i, j) * thObsV.at(j)->computeThValue();
                sum += channels(i, j);
            }
            mu += (1. - sum) * thObsV.at(thobsvsize - 1)->computeThValue();
            logprob += LogSplitGaussian(mu * tho->computeThValue(), channels(i, 3), channels(i, 4), channels(i, 5));
        }
    }

    return (logprob);
}


