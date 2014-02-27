 /*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "InputParser.h"
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>

InputParser::InputParser() 
{
    myModel = NULL;
    thf = NULL;
}

InputParser::InputParser(const InputParser& orig) 
{
    myModel = new StandardModel(*orig.myModel);
    thf = new ThFactory(*orig.thf);
}

InputParser::~InputParser() 
{
    if (myModel != NULL)
        delete myModel;
    if (thf != NULL)
        delete thf;
}

Observable InputParser::ParseObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg) 
{
    std::string name = *beg;
    ++beg;
    std::string thname = *beg;
    ++beg;
    std::string label = *beg;
    size_t pos = -1;
    while ((pos = label.find("~", pos + 1)) != std::string::npos)
        label.replace(pos, 1, " ");
    ++beg;
    double min = atof((*beg).c_str());
    ++beg;
    double max = atof((*beg).c_str());
    ++beg;
    std::string toMCMC = *beg;
    bool tMCMC;
    if (toMCMC.compare("MCMC") == 0)
        tMCMC = true;
    else if (toMCMC.compare("noMCMC") == 0)
        tMCMC = false;
    else
        throw std::runtime_error("ERROR: wrong MCMC flag in " + name);
    Observable o(name, thname, label, tMCMC, min, max, thf->getThMethod(thname));
    ++beg;
    std::string distr = *beg;
    if (distr.compare("file") == 0) {
        ++beg;
        o.setFilename(*beg);
        ++beg;
        o.setHistoname(*beg);
    } else if (distr.compare("weight") == 0) {
        ++beg;
        o.setAve(atof((*beg).c_str()));
        ++beg;
        o.setErrg(atof((*beg).c_str()));
        ++beg;
        o.setErrf(atof((*beg).c_str()));
    } else if (distr.compare("noweight") == 0) {
    } else
        throw std::runtime_error("ERROR: wrong distribution flag in " + name);
    o.setDistr(distr);
    return (o);
}

std::string InputParser::ReadParameters(const std::string filename,
                                        std::vector<ModelParameter>& ModelPars,
                                        std::vector<Observable>& Observables,
                                        std::vector<Observable2D>& Observables2D,
                                        std::vector<CorrelatedGaussianObservables>& CGO,
                                        std::vector<ModelParaVsObs>& ParaObs)
{
    std::string modname = "";
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid model configuration file.\n");
    std::string line;
    bool IsEOF = false;
    do {
        IsEOF = getline(ifile, line).eof();
        if (*line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
        if (line.empty() || line.at(0) == '#')
            continue;
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
        if (beg->compare("StandardModel") == 0) {
            modname = *beg;
            myModel = new StandardModel();
            continue;
        } else if (beg->compare("NPSTU") == 0) {
            modname = *beg;
            myModel = new NPSTU();
            continue;
        } else if (beg->compare("NPSTUVWXY") == 0) {
            modname = *beg;
            myModel = new NPSTUVWXY();
            continue;
        } else if (beg->compare("NPEpsilons") == 0) {
            modname = *beg;
            myModel = new NPEpsilons();
            continue;
        } else if (beg->compare("NPEpsilons_pureNP") == 0) {
            modname = *beg;
            myModel = new NPEpsilons_pureNP();
            continue;
        } else if (beg->compare("NPHiggsST") == 0) {
            modname = *beg;
            myModel = new NPHiggsST();
            continue;
        } else if (beg->compare("NPZbbbar") == 0) {
            modname = *beg;
            myModel = new NPZbbbar();
            continue;
        } else if (beg->compare("NPEffective1") == 0) {
            modname = *beg;
            myModel = new NPEffective1();
            continue;
        } else if (beg->compare("NPEffective2") == 0) {
            modname = *beg;
            myModel = new NPEffective2();
            continue;
        } else if (beg->compare("MFV") == 0) {
            modname = *beg;
            myModel = new MFV();
            continue;
        } else if (beg->compare("GeneralSUSY") == 0) {
            modname = *beg;
            myModel = new GeneralSUSY();
            continue;
        } else if (beg->compare("pMSSM") == 0) {
            modname = *beg;
            myModel = new pMSSM();
            continue;
        } else if (beg->compare("SusyMI") == 0) {
            modname = *beg;
            myModel = new SUSYMassInsertion();
            continue;
        } else if (beg->compare("THDM") == 0) {
            modname = *beg;
            myModel = new THDM();
            continue;
        }
        if (!myModel->IsModelInitialized()) {
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
        }
        
        std::string type = *beg;
        ++beg;
        if (type.compare("ModelParameter") == 0) {
            if (std::distance(tok.begin(),tok.end()) < 5)
                throw std::runtime_error("ERROR: lack of information on "
                                         + *beg + " in " + filename);
            std::string name = *beg;
            ++beg;
            double mean = atof((*beg).c_str());
            ++beg;
            double errg = atof((*beg).c_str());
            ++beg;
            double errf = atof((*beg).c_str());
            ++beg;
            ModelParameter m(name, mean, errg, errf);
            ModelPars.push_back(m);
            if (beg != tok.end())
                std::cout << "WARNING: unread information in parameter " << name << std::endl;
        } else if (type.compare("Observable") == 0) {
            if (std::distance(tok.begin(),tok.end()) < 8)
                throw std::runtime_error("ERROR: lack of information on "
                                         + *beg + " in " + filename);
            Observables.push_back(ParseObservable(beg));
            ++beg;
            if (beg != tok.end())
                std::cout << "WARNING: unread information in observable "
                          << Observables.back().getName() << std::endl;
        } else if (type.compare("Observable2D") == 0) {
            if (std::distance(tok.begin(),tok.end()) < 12)
                throw std::runtime_error("ERROR: lack of information on "
                                         + *beg + " in " + filename);
            Observable2D o2(ParseObservable(beg));
            ++beg;
            o2.setThname2(*beg);
            o2.setTho2(thf->getThMethod(*beg));
            ++beg;
            std::string label = *beg;
            size_t pos = 0;
            while ((pos = label.find("~", pos)) != std::string::npos)
                label.replace(pos, 1, " ");
            o2.setLabel2(label);
            ++beg;
            o2.setMin2(atof((*beg).c_str()));
            ++beg;
            o2.setMax2(atof((*beg).c_str()));
            Observables2D.push_back(o2);
            ++beg;
            if (beg != tok.end())
                std::cout << "WARNING: unread information in observable2D "
                          << Observables2D.back().getName() << std::endl;
        } else if (type.compare("CorrelatedGaussianObservables") == 0) {
            std::string name = *beg;
            ++beg;
            int size = atoi((*beg).c_str());
            CorrelatedGaussianObservables o3(name);
            int nlines = 0;
            std::vector<bool> lines;
            for (int i = 0; i < size; i++) {
                IsEOF = getline(ifile, line).eof();
                if (line.empty() || line.at(0) == '#') {
                    std::cout << "ERROR: no comments or empty lines in CorrelatedGaussianObservables please!"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                beg = mytok.begin();
                std::string type = *beg;
                ++beg;
                if (type.compare("Observable") != 0)
                    throw std::runtime_error("ERROR: expecting an Observable type here...");
                Observable tmp = ParseObservable(beg);
                if (tmp.isTMCMC()) {
                    o3.AddObs(tmp);
                    lines.push_back(true);
                    nlines++;
                } else {
                    Observables.push_back(tmp);
                    lines.push_back(false);
                }
            }
            gslpp::matrix<double> myCorr(gslpp::matrix<double>::Id(nlines));
            int ni = 0;
            for (int i = 0; i < size; i++) {
                IsEOF = getline(ifile, line).eof();
                if (line.empty() || line.at(0) == '#') {
                    std::cout << "ERROR: no comments or empty lines in CorrelatedGaussianObservables please!"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (lines.at(i)) {
                    boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                    beg = mytok.begin();
                    int nj = 0;
                    for (int j = 0; j < size; j++) {
                        if ((*beg).compare(0,1,"0") == 0
                                || (*beg).compare(0,1,"1") == 0
                                || (*beg).compare(0,1,"-") == 0 ) {
                            if (lines.at(j)) {
                                myCorr(ni, nj) = atof((*beg).c_str());
                                nj++;
                            }
                            beg++;
                        } else {
                            std::cout << "ERROR: invalid correlation matrix for "
                                      << name << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                    ni++;
                }
            }
            o3.ComputeCov(myCorr);
            CGO.push_back(o3);
        } else if (type.compare("ModelParaVsObs") == 0) {
            if (std::distance(tok.begin(),tok.end()) < 10)
                throw std::runtime_error("ERROR: lack of information on "
                                         + *beg + " in " + filename);
            std::string name = *beg;
            ++beg;
            std::string ParaName = *beg;
            ++beg;
            std::string ParaLabel = *beg;
            size_t pos = -1;
            while ((pos = ParaLabel.find("~", pos + 1)) != std::string::npos)
                ParaLabel.replace(pos, 1, " ");
            ++beg;
            double ParaMin = atof((*beg).c_str());
            ++beg;
            double ParaMax = atof((*beg).c_str());
            ++beg;
            std::string ObsName = *beg;
            ++beg;
            std::string ObsLabel = *beg;
            pos = -1;
            while ((pos = ObsLabel.find("~", pos + 1)) != std::string::npos)
                ObsLabel.replace(pos, 1, " ");
            ++beg;
            double ObsMin = atof((*beg).c_str());
            ++beg;
            double ObsMax = atof((*beg).c_str());

            ModelParaVsObs pm(name, ParaName, ParaLabel, ParaMin, ParaMax,
                              ObsName, ObsLabel, ObsMin, ObsMax,
                              thf->getThMethod(ObsName));
            ParaObs.push_back(pm);
            ++beg;
            if (beg != tok.end()) std::cout << "WARNING: unread information in ModelParaVsObs "
                    << ParaObs.back().getName() << std::endl;
        } else if (type.compare("ModelFlag") == 0) {
            if (std::distance(tok.begin(),tok.end()) < 3)
                throw std::runtime_error("ERROR: lack of information on "
                                         + *beg + " in " + filename);
            std::string flagname = *beg;
            ++beg;
            if ( boost::iequals(*beg, "true") || boost::iequals(*beg, "false") ) {
                /* Boolean flags */
                bool value_bool;
                if ( boost::iequals(*beg, "true") )
                    value_bool = 1;
                else
                    value_bool = 0;
                if (!myModel->setFlag(flagname, value_bool))
                    throw std::runtime_error("ERROR: setFlag error for " + flagname);
                else
                    std::cout << "set flag " << flagname << "=" << *beg << std::endl;
            } else {
                /* String flags */
                std::string value_str = *beg;
                if (!myModel->setFlagStr(flagname, value_str))
                    throw std::runtime_error("ERROR: setFlag error for " + flagname);
                else
                    std::cout << "set flag " << flagname << "=" << value_str << std::endl;
            }
            ++beg;
            if (beg != tok.end())
                std::cout << "WARNING: unread information in Flag " << flagname << std::endl;
        } else
            throw std::runtime_error("\nERROR: wrong keyword " + type + " in config file. Make sure to specify a valid model configuration file.\n" );
    } while (!IsEOF);

    if (!myModel->CheckFlags())
        throw std::runtime_error("ERROR: incompatible flag(s)");

    return (modname);
}


