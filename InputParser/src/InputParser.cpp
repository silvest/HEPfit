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
    modelnotset = 0;
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

StandardModel* InputParser::ModelFactory(std::string& ModelName){
    
    std::map<std::string, StandardModel* > DModel;
    
    if (ModelName.compare("StandardModel") == 0) DModel["StandardModel"] = new StandardModel();
    else if (ModelName.compare("NPSTU") == 0) DModel["NPSTU"] = new NPSTU();
    else if (ModelName.compare("NPSTUVWXY") == 0) DModel["NPSTUVWXY"] = new NPSTUVWXY();
    else if (ModelName.compare("NPEpsilons") == 0) DModel["NPEpsilons"] =  new NPEpsilons();
    else if (ModelName.compare("NPEpsilons_pureNP") == 0) DModel["NPEpsilons_pureNP"] = new NPEpsilons_pureNP();
    else if (ModelName.compare("NPHiggs") == 0) DModel["NPHiggs"] = new NPHiggs();
    else if (ModelName.compare("NPZbbbar") == 0) DModel["NPZbbbar"] = new NPZbbbar();
    else if (ModelName.compare("NPEffective1") == 0) DModel["NPEffective1"] = new NPEffective1();
    else if (ModelName.compare("NPEffective2") == 0) DModel["NPEffective2"] = new NPEffective2();
    else if (ModelName.compare("MFV") == 0) DModel["MFV"] = new MFV();
    else if (ModelName.compare("GeneralSUSY") == 0) DModel["GeneralSUSY"] =  new GeneralSUSY();
    else if (ModelName.compare("pMSSM") == 0) DModel["pMSSM"] = new pMSSM();
    else if (ModelName.compare("SusyMassInsertion") == 0) DModel["SusyMassInsertion"] = new SUSYMassInsertion();
    else if (ModelName.compare("THDM") == 0) DModel["THDM"] = new THDM();
    else return NULL;
    
    return DModel[ModelName];
}

std::string InputParser::ReadParameters(const std::string filename,
                                        std::vector<ModelParameter>& ModelPars,
                                        std::vector<Observable>& Observables,
                                        std::vector<Observable2D>& Observables2D,
                                        std::vector<CorrelatedGaussianObservables>& CGO,
                                        std::vector<ModelParaVsObs>& ParaObs)
{
    modname = "";
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
        
        if (modelnotset == 0) {
            modname = *beg;
            myModel = ModelFactory(modname);
            if (myModel == NULL) continue;
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            modelnotset = 1;
            continue;
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

    if (modelnotset == 0)
        throw std::runtime_error("ERROR: Incorrect or missing model name in model configuration file.\n");
    if (!myModel->CheckFlags())
        throw std::runtime_error("ERROR: incompatible flag(s)\n");

    return (modname);
}


