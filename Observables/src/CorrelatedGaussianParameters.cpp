/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include "CorrelatedGaussianParameters.h"

CorrelatedGaussianParameters::CorrelatedGaussianParameters(std::string name_i)
{
    name = name_i;
    Cov = NULL; 
    v = NULL;
    e = NULL;
}

CorrelatedGaussianParameters::CorrelatedGaussianParameters()
{
    Cov = NULL; 
    v = NULL;
    e = NULL;
}

CorrelatedGaussianParameters::CorrelatedGaussianParameters(const CorrelatedGaussianParameters& orig)
{
    Pars = orig.Pars;
    name = orig.name;
    Cov = new gslpp::matrix<double>(*orig.Cov);
    v = new gslpp::matrix<double>(*orig.v);
    e = new gslpp::vector<double>(*orig.e);
    DiagPars = orig.DiagPars;
}

CorrelatedGaussianParameters::~CorrelatedGaussianParameters()
{   
    //Cov = new gslpp::matrix<double>(10, 10, 0.); /** Put in to prevent seg fault during error handling in InputParser **/
    if(Cov != NULL)
        delete(Cov);
    if(v != NULL)
        delete(v);
    if(e != NULL)
        delete(e);
}

void CorrelatedGaussianParameters::AddPar(ModelParameter& Par_i)
{
    Pars.push_back(Par_i);
}

void CorrelatedGaussianParameters::DiagonalizePars(gslpp::matrix<double> Corr)
{
    unsigned int size = Pars.size();
    if (Corr.size_i() != size || Corr.size_j() != size)
        throw std::runtime_error("The size of the correlated parameters in " + name + " does not match the size of the correlation matrix!");
    Cov = new gslpp::matrix<double>(size, size, 0.);
    for (unsigned int i = 0; i < size; i++)
        for (unsigned int j = 0; j < size; j++)
            (*Cov)(i, j) = Pars.at(i).geterrg() * Corr(i, j) * Pars.at(j).geterrg();

    *Cov = Cov->inverse();
    
    gslpp::matrix<gslpp::complex>vv(size, size, 0.);
    gslpp::vector<gslpp::complex>ee(size, 0.);

    Cov->eigensystem(vv, ee);

    v = new gslpp::matrix<double>(size, size, 0.);
    e = new gslpp::vector<double>(size, 0.);
    *v = vv.real();
    *e = ee.real();

    gslpp::vector<double> ave_in(size, 0.);
    
    int ind = 0;
    for(std::vector<ModelParameter>::iterator it = Pars.begin(); it != Pars.end(); it++)
    {
        ave_in(ind) = it->getave();
        ind++;
    }
    
    gslpp::vector<double> ave = v->transpose() * ave_in;
    
    for (int i = 0; i < size; i++)
    {
        std::stringstream ss;
        ss << (i+1);
        std::string namei = name + ss.str();
        ModelParameter current(namei,ave(i),1./sqrt((*e)(i)),0.);
        current.setCgp_name(name);
        DiagPars.push_back(current);
    }
}

std::vector<double> CorrelatedGaussianParameters::getOrigParsValue(const std::vector<double>& DiagPars_i) const 
{
    if (DiagPars_i.size() != DiagPars.size()) {
        std::stringstream out;
        out << DiagPars_i.size();
        throw std::runtime_error("CorrelatedGaussianParameters::getOrigParsValue(DiagPars_i): DiagPars_i.size() = " + out.str() + " does not match the size of DiagPars");
    }
    gslpp::vector<double> pars_in(DiagPars_i.size(), 0.);

    int ind = 0;
    for (std::vector<double>::const_iterator it = DiagPars_i.begin(); it != DiagPars_i.end(); it++) {
        pars_in(ind) = *it;
        ind++;
    }

    gslpp::vector<double> val = (*v) * pars_in;

    std::vector<double> res;

    for (unsigned int i = 0; i < DiagPars_i.size(); i++) {
        res.push_back(val(i));
    }
    return (res);
}

int CorrelatedGaussianParameters::ParseCGP(std::vector<ModelParameter>& ModPars, 
                                            std::string& filename,
                                            std::ifstream& ifile, 
                                            boost::tokenizer<boost::char_separator<char> >::iterator & beg,
                                            int lineNo,
                                            int rank)
{
    name = *beg;
    ++beg;
    int size = atoi((*beg).c_str());
    int nlines = 0;
    std::string line;
    boost::char_separator<char>sep(" \t");
    for (int i = 0; i < size; i++) {
        IsEOF = getline(ifile, line).eof();
        if (line.empty() || line.at(0) == '#') {
            if (rank == 0) throw std::runtime_error("ERROR: no comments or empty lines in CorrelatedGaussianParameters please! In line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + filename);
            else sleep(2);
        }
        lineNo++;
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        beg = tok.begin();
        std::string type = *beg;
        ++beg;
        if (type.compare("ModelParameter") != 0){
            if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + filename + ", expecting a ModelParameter type here...\n");
            else sleep(2);
        }
        ModelParameter tmpMP;
        beg = tmpMP.ParseModelParameter(beg);
        if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in parameter " << tmpMP.getname() << std::endl;
        tmpMP.setCgp_name(name);
        AddPar(tmpMP);
        nlines++;
    }
    if (nlines > 1) {
        gslpp::matrix<double> myCorr(gslpp::matrix<double>::Id(nlines));
        int ni = 0;
        for (int i = 0; i < size; i++) {
            IsEOF = getline(ifile, line).eof();
            if (line.empty() || line.at(0) == '#') {
                if (rank == 0) throw std::runtime_error("ERROR: no comments or empty lines in CorrelatedGaussianParameters please! In line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + filename);
                else sleep(2);
            }
            lineNo++;
            boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
            beg = mytok.begin();
            int nj = 0;
            for (int j = 0; j < size; j++) {
                if ((*beg).compare(0, 1, "0") == 0
                        || (*beg).compare(0, 1, "1") == 0
                        || (*beg).compare(0, 1, "-") == 0) {
                    if (std::distance(mytok.begin(), mytok.end()) < size) { 
                        if (rank == 0) throw std::runtime_error("ERROR: invalid correlation matrix for " + name + ". Check element (" + boost::lexical_cast<std::string>(ni + 1) + "," + boost::lexical_cast<std::string>(nj + 1) + ") in line number " + boost::lexical_cast<std::string>(lineNo) + " in file " + filename + ".\n");
                        else sleep(2);
                    }
                    myCorr(ni, nj) = atof((*beg).c_str());
                    nj++;
                    beg++;
                } else {
                    if (rank == 0) std::cout << "ERROR: invalid correlation matrix for "
                            << name << ". Check element (" << ni + 1 << "," << nj + 1 << ") in line number " + boost::lexical_cast<std::string>(lineNo) << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            ni++;
        }
        DiagonalizePars(myCorr);
        ModPars.insert(ModPars.end(), getDiagPars().begin(), getDiagPars().end());

    } else {
        if (rank == 0) std::cout << "\nWARNING: Correlated Gaussian Parameters " << name.c_str() << " defined with less than two correlated parameters. The set is being marked as normal Parameters." << std::endl;
        if (getPars().size() == 1) ModPars.push_back(ModelParameter(getPar(0)));
        for (int i = 0; i < size; i++) {
            IsEOF = getline(ifile, line).eof();
            lineNo++;
        }
    }
    return lineNo;
}