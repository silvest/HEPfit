/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "InputParser.h"
#include "ModelFactory.h"
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>

InputParser::InputParser(ModelFactory& ModF, ThObsFactory& ObsF) : myModelFactory(ModF), myObsFactory(ObsF)
{
    modelset = 0;
}

InputParser::InputParser(const InputParser& orig) : myModelFactory(orig.myModelFactory), myObsFactory(orig.myObsFactory)
{
    myModel = new StandardModel(*orig.myModel);
}

InputParser::~InputParser()
{
}

Observable InputParser::ParseObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg)
{
    double min;
    double max;
    std::string name = *beg;
    ++beg;
    std::string thname = *beg;
    ++beg;
    std::string label = *beg;
    size_t pos = -1;
    while ((pos = label.find("~", pos + 1)) != std::string::npos)
        label.replace(pos, 1, " ");
    ++beg;
    min = atof((*beg).c_str());
    ++beg;
    max = atof((*beg).c_str());
    ++beg;
    std::string toMCMC = *beg;
    bool tMCMC;
    if (toMCMC.compare("MCMC") == 0)
        tMCMC = true;
    else if (toMCMC.compare("noMCMC") == 0)
        tMCMC = false;
    else
        throw std::runtime_error("ERROR: wrong MCMC flag in " + name);
    Observable o(name, thname, label, tMCMC, min, max, myObsFactory.CreateThMethod(thname, *myModel));
    return (o);
}

std::string InputParser::ReadParameters(const std::string filename,
        const int rank,
        std::vector<ModelParameter>& ModelPars,
        boost::ptr_vector<Observable>& Observables,
        std::vector<Observable2D>& Observables2D,
        std::vector<CorrelatedGaussianObservables>& CGO)
{
    modname = "";
    int lineNo = 0;
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid model configuration file.\n");
    std::string filepath; 
    if (filename.find("\\/") == std::string::npos) filepath = filename.substr(0, filename.find_last_of("\\/")+1);
    std::cout << filepath << std::endl;
    std::string line;
    bool IsEOF = false;
    do {
        IsEOF = getline(ifile, line).eof();
        lineNo += 1;
        if (*line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
        if (line.empty() || line.at(0) == '#')
            continue;
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();

        if (modelset == 0) {
            modname = *beg;
            myModel = myModelFactory.CreateModel(modname);
            myModel->setModelName(modname);
            myModel->InitializeModel();
            if (myModel->IsModelInitialized()) {
                if (rank == 0) std::cout << "\nModel Initialized: " << modname << std::endl;
                modeldefinedinfile = filename;
            } else {
                throw std::runtime_error("\nERROR: " + modname + " not initialized successfully.\n");
            }
            modelset = 1;
            continue;
        } else if (modelset == 1 && beg->compare(myModel->ModelName()) == 0) {
            continue;
        }

        std::string type = *beg;
        ++beg;
        if (type.compare("ModelParameter") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 5)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
            std::string name = *beg;
            if (checkDuplicateParameter[name].get<0>()) throw std::runtime_error("\nERROR: ModelParameter " + name + " appears more than once ...!! \n" +
                                                                "1st Occurrence: Line No:" + boost::lexical_cast<std::string>(checkDuplicateParameter[name].get<2>()) + 
                                                                " in file " + checkDuplicateParameter[name].get<1>() + ".\n"
                                                                "2nd Occurrence: Line No:" + boost::lexical_cast<std::string>(lineNo) + " in file " + filename + ".\n");
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
                if (rank == 0) std::cout << "WARNING: unread information in parameter " << name << std::endl;
            checkDuplicateParameter[name] = boost::make_tuple (true, filename, lineNo);
        } else if (type.compare("Observable") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 8)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
            Observable * o = new Observable(ParseObservable(beg));
            ++beg;
            std::string distr = *beg;
            if (distr.compare("file") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 10)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
                std::string fname = filepath + *(++beg);
                std::string histoname = *(++beg);
                o->setLikelihoodFromHisto(fname, histoname);
            } else if (distr.compare("weight") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 11)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
                ++beg;
                o->setAve(atof((*beg).c_str()));
                ++beg;
                o->setErrg(atof((*beg).c_str()));
                ++beg;
                o->setErrf(atof((*beg).c_str()));
                if (o->getErrf() == 0. && o->getErrg() == 0.){
                    throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + o->getName() + " cannot both be 0. in the " + filename + " file.\n");
                }                
            } else if (distr.compare("noweight") == 0) {
            } else
                throw std::runtime_error("ERROR: wrong distribution flag in " + o->getName());
            o->setDistr(distr);
            Observables.push_back(o);
            ++beg;
            if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in observable "
                    << Observables.back().getName() << std::endl;
        } else if (type.compare("BinnedObservable") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 10)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
            Observable * bo = new Observable(ParseObservable(beg));
            bo->setObsType(2);
            ++beg;
            std::string distr = *beg;
            if (distr.compare("file") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 12)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
                std::string fname = filepath + *(++beg);
                std::string histoname = *(++beg);
                bo->setLikelihoodFromHisto(fname, histoname);
            } else if (distr.compare("weight") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 13)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
                ++beg;
                bo->setAve(atof((*beg).c_str()));
                ++beg;
                bo->setErrg(atof((*beg).c_str()));
                ++beg;
                bo->setErrf(atof((*beg).c_str()));
                if (bo->getErrf() == 0. && bo->getErrg() == 0.) {
                    throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + bo->getName() + " cannot both be 0. in the " + filename + " file.");
                }
            } else if (distr.compare("noweight") == 0) {
                ++beg;
                ++beg;
                ++beg;
            } else
                throw std::runtime_error("ERROR: wrong distribution flag in " + bo->getName());
            bo->setDistr(distr);
            ++beg;
            bo->getTho()->setBinMin(atof((*beg).c_str()));
            ++beg;
            bo->getTho()->setBinMax(atof((*beg).c_str()));
            Observables.push_back(bo);
            ++beg;
            if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in observable "
                    << Observables.back().getName() << std::endl;
        } else if (type.compare("Observable2D") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 12)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
            Observable2D o2(ParseObservable(beg));
            ++beg;
            std::string distr = *beg;
            if (distr.compare("file") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 14)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
                std::string fname = filepath + *(++beg);
                std::string histoname = *(++beg);
                o2.setLikelihoodFromHisto(fname, histoname);
            } else if (distr.compare("noweight") == 0) {
            } else
                throw std::runtime_error("ERROR: wrong distribution flag in " + o2.getName());
            o2.setDistr(distr);
            ++beg;
            o2.setThname2(*beg);
            o2.setTho2(myObsFactory.CreateThMethod(*beg, *myModel));
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
                if (rank == 0) std::cout << "WARNING: unread information in observable2D "
                    << Observables2D.back().getName() << std::endl;
        } else if (type.compare("HiggsObservable") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 8)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
            HiggsObservable * ho = new HiggsObservable(ParseObservable(beg));
            ho->setObsType(1);
            ++beg;
            std::string distr = *beg;
            if (distr.compare("parametric") == 0) {
                std::vector<ThObservable*> hthobs;
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
                } else
                    throw std::runtime_error("ERROR: wrong keyword " + distr + " in " + ho->getName());
                ho->setParametricLikelihood(*(++beg), hthobs);
            } else
                throw std::runtime_error("ERROR: wrong distribution flag " + distr + " in " + ho->getName());
            Observables.push_back(ho);
            ++beg;
            if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in HiggsObservable "
                    << Observables.back().getName() << std::endl;
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
                    if (rank == 0) std::cout << "ERROR: no comments or empty lines in CorrelatedGaussianObservables please!"
                            << std::endl;
                    exit(EXIT_FAILURE);
                }
                boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                beg = mytok.begin();
                std::string type = *beg;
                ++beg;
                if (type.compare("Observable") != 0 && type.compare("BinnedObservable") != 0)
                    throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + filename + ", expecting an Observable or BinnedObservable type here...\n");
                Observable * tmp = new Observable(ParseObservable(beg));
                if (type.compare("Observable") == 0) tmp->setObsType(1);
                if (type.compare("BinnedObservable") == 0) tmp->setObsType(2);
                ++beg;
                std::string distr = *beg;
                if (distr.compare("weight") == 0) {
                    ++beg;
                    tmp->setAve(atof((*beg).c_str()));
                    ++beg;
                    tmp->setErrg(atof((*beg).c_str()));
                    ++beg;
                    tmp->setErrf(atof((*beg).c_str()));
                } else if (distr.compare("noweight") == 0) {
                    if (type.compare("BinnedObservable") == 0){
                        ++beg;
                        ++beg;
                        ++beg;
                    }
                } else
                    throw std::runtime_error("ERROR: wrong distribution flag in " + tmp->getName());
                tmp->setDistr(distr);
                if (type.compare("BinnedObservable") == 0){
                    ++beg;
                    tmp->getTho()->setBinMin(atof((*beg).c_str()));
                    ++beg;
                    tmp->getTho()->setBinMax(atof((*beg).c_str()));
                }
                if (tmp->isTMCMC()) {
                    o3.AddObs(*tmp);
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
                    if (rank == 0) std::cout << "ERROR: no comments or empty lines in CorrelatedGaussianObservables please!"
                            << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (lines.at(i)) {
                    boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                    beg = mytok.begin();
                    int nj = 0;
                    for (int j = 0; j < size; j++) {
                        if ((*beg).compare(0, 1, "0") == 0
                                || (*beg).compare(0, 1, "1") == 0
                                || (*beg).compare(0, 1, "-") == 0) {
                            if (lines.at(j)) {
                                myCorr(ni, nj) = atof((*beg).c_str());
                                nj++;
                            }
                            beg++;
                        } else {
                            if (rank == 0) std::cout << "ERROR: invalid correlation matrix for "
                                    << name << ". Check element (" << ni+1 << "," << nj+1 << ")" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                    ni++;
                }
            }
            o3.ComputeCov(myCorr);
            CGO.push_back(o3);
        } else if (type.compare("ModelFlag") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 3)
                throw std::runtime_error("ERROR: lack of information on "
                    + *beg + " in " + filename);
            std::string flagname = *beg;
            ++beg;
            if (boost::iequals(*beg, "true") || boost::iequals(*beg, "false")) {
                /* Boolean flags */
                bool value_bool;
                if (boost::iequals(*beg, "true"))
                    value_bool = 1;
                else
                    value_bool = 0;
                if (!myModel->setFlag(flagname, value_bool))
                    throw std::runtime_error("ERROR: setFlag error for " + flagname);
                else
                    if (rank == 0) std::cout << "set flag " << flagname << "=" << *beg << std::endl;
            } else {
                /* String flags */
                std::string value_str = *beg;
                if (!myModel->setFlagStr(flagname, value_str))
                    throw std::runtime_error("ERROR: setFlag error for " + flagname);
                else
                    if (rank == 0) std::cout << "set flag " << flagname << "=" << value_str << std::endl;
            }
            ++beg;
            if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in Flag " << flagname << std::endl;
        } else if (type.compare("IncludeFile") == 0) {
            std::string IncludeFileName = filepath + *beg;
            if (rank == 0) std::cout << "\nIncluding File: " + IncludeFileName << std::endl;
            ReadParameters(IncludeFileName, rank, ModelPars, Observables, Observables2D, CGO);
            ++beg;
        } else
            throw std::runtime_error("\nERROR: wrong keyword " + type + " in file " + filename + ". Make sure to specify a valid model configuration file.");
    } while (!IsEOF);

    if (modelset == 0)
        throw std::runtime_error("ERROR: Incorrect or missing model name in the model configuration file.\n");
    if (!myModel->CheckFlags())
        throw std::runtime_error("ERROR: incompatible flag(s)\n");

    return (modname);
}


