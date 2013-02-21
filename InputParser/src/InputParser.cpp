/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "InputParser.h"
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <iostream>

InputParser::InputParser() {
    myModel = NULL;
    thf = NULL;
}

InputParser::InputParser(const InputParser& orig) {
    myModel = new StandardModel(*orig.myModel);
    thf = new ThFactory(*orig.thf);
}

InputParser::~InputParser() {
    if (myModel != NULL)
        delete myModel;
    if (thf != NULL)
        delete thf;
}

Observable InputParser::ParseObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg) {
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
    else {
        std::cout << "wrong MCMC flag in " << name << std::endl;
        exit(EXIT_FAILURE);
    }
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
    } else {
        std::cout << "wrong distribution flag in " << name << std::endl;
        exit(EXIT_FAILURE);
    }
    o.setDistr(distr);
    return (o);
}

std::string InputParser::ReadParameters(const std::string filename, std::vector<ModelParameter>&
        ModelPars, std::vector<Observable>& Observables, std::vector<Observable2D>& Observables2D, std::vector<CorrelatedGaussianObservables>& CGO) {
    std::string modname = "";
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open()) {
        std::cout << filename << " does not exist." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    while (!getline(ifile, line).eof()) {
        if (line.empty() || line.at(0) == '#')
            continue;
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
        if (beg->compare("StandardModel") == 0) {
            modname = *beg;
            myModel = new StandardModel();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("NPSTU") == 0) {
            modname = *beg;
            myModel = new NPSTU();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("NPSTUVWXY") == 0) {
            modname = *beg;
            myModel = new NPSTUVWXY();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("NPEpsilons") == 0) {
            modname = *beg;
            myModel = new NPEpsilons();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("NPHiggs") == 0) {
            modname = *beg;
            myModel = new NPHiggs();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("NPHiggsST") == 0) {
            modname = *beg;
            myModel = new NPHiggsST();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("NPZbbbar") == 0) {
            modname = *beg;
            myModel = new NPZbbbar();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("MFV") == 0) {
            modname = *beg;
            myModel = new MFV();
            myModel->InitializeModel();
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("SusyMI") == 0) {
            modname = *beg;
            SUSYMassInsertion* LocalPointer = new SUSYMassInsertion();
            myModel = LocalPointer;
            thf = new ThFactory(*myModel);
            continue;
        } else if (beg->compare("THDM") == 0) {
            modname = *beg;
            myModel = new THDM();
            thf = new ThFactory(*myModel);
            continue;
        }

        std::string type = *beg;
        ++beg;
        if (type.compare("ModelParameter") == 0) {
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
                std::cout << "warning: unread information in parameter " << name << std::endl;
        } else if (type.compare("Observable") == 0) {
            Observables.push_back(ParseObservable(beg));
            ++beg;
            if (beg != tok.end()) std::cout << "warning: unread information in observable "
                    << Observables.back().getName() << std::endl;
        } else if (type.compare("Observable2D") == 0) {
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
            if (beg != tok.end()) std::cout << "warning: unread information in observable "
                    << Observables.back().getName() << std::endl;
        } else if (type.compare("CorrelatedGaussianObservables") == 0) {
            std::string name = *beg;
            ++beg;
            int size = atoi((*beg).c_str());
	    CorrelatedGaussianObservables o3(name);
            int nlines = 0;
            std::vector<bool> lines;
            for (int i = 0; i < size; i++) {
                getline(ifile, line);
                if (line.empty() || line.at(0) == '#')
                    throw std::runtime_error("no comments or empty lines in CorrelatedGaussianObservables please!");
                boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                beg = mytok.begin();
                std::string type = *beg;
                ++beg;
                if (type.compare("Observable") != 0)
                    throw std::runtime_error("Expecting an Observable type here...");
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
	      getline(ifile, line);
	      if (lines.at(i)) {
		boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
		beg = mytok.begin();
		int nj = 0;
		for (int j = 0; j < size; j++) {
		  if (lines.at(j)){
		    myCorr(ni, nj) = atof((*beg).c_str());
		    nj++;
		  }
		  beg++;
		}
		ni++;   
	      }
            }
            o3.ComputeCov(myCorr);
            CGO.push_back(o3);
        } else if (type.compare("ModelFlag") == 0) {

            std::string name = *beg;
                    ++beg;
                    bool value = boost::lexical_cast<bool>((*beg).c_str());
                    ++beg;

            if (!myModel->SetFlag(name, value)) {
                std::stringstream ss;
                        ss << myModel->ModelName() << " SetFlag error for Flag " << name;
                        throw std::runtime_error(ss.str());
            }
            if (beg != tok.end())
                    std::cout << "warning: unread information in Flag " << name << std::endl;
            } else {
            std::cout << "wrong keyword " << *beg << " in config file (first word must be ModelParameter, ModelFlag or Observable)" << std::endl;
                    exit(EXIT_FAILURE);
        }
    }
    return (modname);
}


