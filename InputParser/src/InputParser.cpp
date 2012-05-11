/* 
 * File:   InputParser.cpp
 * Author: silvest
 * 
 * Created on March 15, 2011, 2:36 PM
 */

#include "InputParser.h"

InputParser::InputParser() {
    myModel = NULL;
    myModelMatching = NULL;
    myEWModel = NULL;
    thf = NULL;
}

InputParser::InputParser(const InputParser& orig) {
    myModel = new StandardModel(*orig.myModel);
    myModelMatching = new StandardModelMatching(*orig.myModelMatching);
    myEWModel = new EWSM(*orig.myEWModel);
    thf = new ThFactory(*orig.thf);
}

InputParser::~InputParser() {
    if (myModel != NULL)
        delete myModel;
    if (myModelMatching != NULL)
        delete myModelMatching;
    if (myEWModel != NULL)
        delete myEWModel;
    if (thf != NULL)
        delete thf;
}

std::string InputParser::ReadParameters(const std::string filename, std::vector<ModelParameter>&
        ModelPars, std::vector<Observable>& Observables, std::vector<Observable2D>& Observables2D) {
    std::string modname = "";
    std::ifstream ifile(filename.c_str());
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
            myModelMatching = new StandardModelMatching(*myModel);
            myEWModel = new EWSM(*myModel);
            thf = new ThFactory(*myModel,*myModelMatching,*myEWModel);
            continue;
        }
        else if (beg->compare("MFV") == 0) {
            modname = *beg;
            myModel = new MFV();
            myModelMatching = new StandardModelMatching(*myModel);
            myEWModel = new EWSM(*myModel);            
            thf = new ThFactory(*myModel,*myModelMatching,*myEWModel);
            continue;
        }
        else if (beg->compare("SusyMI") == 0) {
            modname = *beg;
            SUSYMassInsertion* LocalPointer = new SUSYMassInsertion();
            myModel = LocalPointer;
            myModelMatching = new SUSYMassInsertionMatching(*LocalPointer);
            myEWModel = new EWSM(*myModel);
            thf = new ThFactory(*myModel,*myModelMatching,*myEWModel);
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
        } else if (type.compare("Observable") == 0 || type.compare("Observable2D") == 0) {
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
            if (type.compare("Observable") == 0) {
                Observables.push_back(o);
                ++beg;
                if (beg != tok.end()) std::cout << "warning: unread information in observable "
                        << Observables.back().getName() << std::endl;
            } else { // Observable2D
                ++beg;
                Observable2D o2(o);
                o2.setThname2(*beg);
                o2.setTho2(thf->getThMethod(*beg));
                ++beg;
                label = *beg;
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
            }
        } else {
            std::cout << "wrong keyword " << *beg << " in config file (first word must be ModelParameter or Observable)" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    return(modname);
}

