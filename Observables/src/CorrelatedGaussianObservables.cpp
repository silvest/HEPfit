/* 
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "CorrelatedGaussianObservables.h"
#include "ThObsFactory.h"
#include <fstream>

CorrelatedGaussianObservables::CorrelatedGaussianObservables(std::string name_i)
{
    name = name_i;
    IsPrediction = false;
    IsEOF = false;
    covarianceFromConfig = false;
    InvCov = NULL;
}

CorrelatedGaussianObservables::CorrelatedGaussianObservables()
{
    name = "";
    IsPrediction = false;
    IsEOF = false;
    covarianceFromConfig = false;
    InvCov = NULL;
}

CorrelatedGaussianObservables::CorrelatedGaussianObservables(const CorrelatedGaussianObservables& orig)
{
    Obs = orig.Obs;
    name = orig.name;
    if (orig.InvCov) InvCov = new gslpp::matrix<double>(*orig.InvCov);
    else InvCov = NULL;
    IsPrediction = orig.IsPrediction;
    IsEOF = orig.IsEOF;
    covarianceFromConfig = false;
}

CorrelatedGaussianObservables::~CorrelatedGaussianObservables()
{   
//    Cov = new gslpp::matrix<double>(100, 100, 0.); /** Put in to prevent seg fault during error handling in InputParser **/
    if (InvCov) delete(InvCov);
}

void CorrelatedGaussianObservables::AddObs(Observable& Obs_i)
{
    Obs.push_back(Obs_i);
}

void CorrelatedGaussianObservables::ComputeCov(gslpp::matrix<double> Corr)
{
    unsigned int size = Obs.size();
    if (Corr.size_i() != size || Corr.size_j() != size)
        throw std::runtime_error("The size of the correlated observables in " + name + " does not match the size of the correlation matrix!");
    InvCov = new gslpp::matrix<double>(size, size, 0.);
    if (covarianceFromConfig) {
        for (unsigned int i = 0; i < size; i++)
            for (unsigned int j = 0; j < size; j++)
                (*InvCov)(i, j) = Corr(i, j);
    } else {
        for (unsigned int i = 0; i < size; i++) {
            for (unsigned int j = 0; j < size; j++) {
                if (Corr(i, j) != Corr(j, i)) throw std::runtime_error("Invalid correlation matrix for " + name + ". The matrix is not symmetric.");
                (*InvCov)(i, j) = Obs.at(i).getErrg() * Corr(i, j) * Obs.at(j).getErrg();
            }
        }
        *InvCov = InvCov->inverse();
    }
}

double CorrelatedGaussianObservables::computeWeight()
{
    gslpp::vector<double> x(Obs.size());

    for (unsigned int i = 0; i < Obs.size(); i++)
        x(i) = Obs.at(i).computeTheoryValue() - Obs.at(i).getAve();

    return (-0.5 * x * ((*InvCov) * x));
}

int CorrelatedGaussianObservables::ParseCGO(boost::ptr_vector<Observable>& Observables, 
                                            std::ifstream& ifile, 
                                            boost::tokenizer<boost::char_separator<char> >::iterator& beg, 
                                            std::string& infilename, 
                                            ThObsFactory& myObsFactory,
                                            StandardModel * myModel,
                                            int lineNo,
                                            int rank) {
    if (infilename.find("\\/") == std::string::npos) filepath = infilename.substr(0, infilename.find_last_of("\\/") + 1);
    boost::char_separator<char> sep(" \t");
    name = *beg;
    ++beg;
    int size = atoi((*beg).c_str());
    
    int nlines = 0;
    std::vector<bool> lines;
    std::string line;
    for (int i = 0; i < size; i++) {
        IsEOF = getline(ifile, line).eof();
        if (line.empty() || line.at(0) == '#') {
            if (rank == 0) throw std::runtime_error("ERROR: no comments or empty lines in CorrelatedGaussianObservables please! In file " + infilename + " at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n");
            else sleep(2);
        }
        lineNo++;
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        beg = tok.begin();
        std::string type = *beg;
        ++beg;
        if (type.compare("Observable") != 0 && type.compare("BinnedObservable") != 0 && type.compare("FunctionObservable") != 0) {
            if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + infilename + ", expecting an Observable or BinnedObservable or FunctionObservable type here...\n");
            else sleep(2);
        }
        Observable * tmpObs = new Observable();
        tmpObs->setHasInverseCovariance(covarianceFromConfig);
        beg = tmpObs->ParseObservable(type, &tok, beg, filepath, infilename, rank);
        tmpObs->setTho(myObsFactory.CreateThMethod(tmpObs->getThname(), *myModel));
        if (!IsPrediction) {
            if (tmpObs->isTMCMC()) {
                AddObs(*tmpObs);
                lines.push_back(true);
                nlines++;
            } else {
                Observables.push_back(tmpObs);
                lines.push_back(false);
            }
        } else {
            if (tmpObs->isTMCMC()) {
                if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + infilename + ", Observable/BinnedObservable cannot be set to MCMC in CorrelatedObservable. Use CorrelatedGaussianObservable instead.\n");
                else sleep(2);
            } else {
                AddObs(*tmpObs);
                nlines++;
            }
            if (tmpObs->getDistr().compare("weight") == 0 || tmpObs->getDistr().compare("file") == 0) {
                if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + infilename + ", Observable/BinnedObservable cannot be set to weight or file in CorrelatedObservable. Use CorrelatedGaussianObservable instead.\n");
                else sleep(2);
            }
        }
    }

    if (nlines > 1) {
        if (!IsPrediction) {
            gslpp::matrix<double> myCorr(gslpp::matrix<double>::Id(nlines));
            int ni = 0;
            for (int i = 0; i < size; i++) {
                IsEOF = getline(ifile, line).eof();
                if (line.empty() || line.at(0) == '#') {
                    if (rank == 0) throw std::runtime_error("ERROR: no comments or empty lines in CorrelatedGaussianObservables please! In file " + infilename + " at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n");
                    else sleep(2);
                }
                lineNo++;
                if (lines.at(i)) {
                    boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                    beg = mytok.begin();
                    int nj = 0;
                    for (int j = 0; j < size; j++) {
                        if ((*beg).compare(0, 1, "0") == 0
                                || (*beg).compare(0, 1, "1") == 0
                                || (*beg).compare(0, 1, "-") == 0 || (covarianceFromConfig == true)) {
                            if (std::distance(mytok.begin(), mytok.end()) < size) {
                                if (rank == 0) throw std::runtime_error(("ERROR: Correlation matrix is of wrong size in Correlated Gaussian Observables: " + name + +" at line number:" + boost::lexical_cast<std::string>(lineNo) + ".\n").c_str());
                                else sleep(2);
                            }
                            if (lines.at(j)) {
                                myCorr(ni, nj) = atof((*beg).c_str());
                                nj++;
                            }
                            beg++;
                        } else {
                            if (rank == 0) throw std::runtime_error("ERROR: invalid correlation matrix for " + name + ". Check element (" + boost::lexical_cast<std::string>(ni + 1) + "," + boost::lexical_cast<std::string>(nj + 1) + ") in line number " + boost::lexical_cast<std::string>(lineNo) + " in file " + infilename + ".\n");
                            else sleep(2);
                        }
                    }
                    ni++;
                }
            }
            ComputeCov(myCorr);
        } else {
            InvCov = new gslpp::matrix<double>(size, size, 0.);
        }
    } else {
        if (rank == 0) std::cout << "\nWARNING: Correlated (Gaussian) Observable " << name.c_str() << " defined with less than two observables" << " in file " << infilename << ". The set is being marked as normal Observables." << std::endl;
        if (getObs().size() == 1) Observables.push_back(new Observable(getObs(0)));
        for (int i = 0; i < size; i++) {
            IsEOF = getline(ifile, line).eof();
            lineNo++;
        }
    }
    return lineNo;
}