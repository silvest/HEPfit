/* 
 * File:   InputParser.cpp
 * Author: silvest
 * 
 * Created on March 15, 2011, 2:36 PM
 */

#include "InputParser.h"

ThObservable * InputParser::ThFactory(const std::string& name){
    ThObservable * tho;
    if(name.compare("Dmd0")==0)
        tho = new Dmb((StandardModel *) myModel,0);
    else if(name.compare("Dmd1")==0)
        tho = new Dmb((StandardModel *) myModel,1);
    else if(name.compare("Vud")==0)
        tho = new Vud((StandardModel *) myModel);
    else if(name.compare("Vus")==0)
        tho = new Vus((StandardModel *) myModel);
    else if(name.compare("Vub")==0)
        tho = new Vub((StandardModel *) myModel);
    else if(name.compare("Vcb")==0)
        tho = new Vcb((StandardModel *) myModel);
    else if(name.compare("alpha")==0)
        tho = new alpha((StandardModel *) myModel);
    else if(name.compare("gamma")==0)
        tho = new gammac((StandardModel *) myModel);
    else {
        std::cout << "wrong observable " << name <<" in ThFactory" << std::endl;
        exit(EXIT_FAILURE);        
    }
    return tho;
}

InputParser::InputParser() {
    myModel = NULL;
    myFlavour = NULL;
}

InputParser::InputParser(const InputParser& orig) {
    myModel = orig.myModel;
    myFlavour = orig.myFlavour;
}

InputParser::~InputParser() {
}

void InputParser::ReadParameters(const char * filename, std::vector<ModelParameter>&
ModelPars, std::vector<Observable>& Observables){
    std::ifstream ifile(filename);
    std::string line;
    while(!getline(ifile,line).eof()){
        if(line.at(0)=='#') continue;
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > tok(line,sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin();
        if(beg->compare("StandardModel")==0) {
            myModel = new StandardModel();
            myFlavour = new Flavour((StandardModel&) *myModel); 
            continue;
        }
        if(beg->compare("ModelParameter")==0){
            ++beg;
            std::string name = *beg;
            ++beg;
            double mean = atof((*beg).c_str());
            ++beg;
            double errg = atof((*beg).c_str());
            ++beg;
            double errf = atof((*beg).c_str());
            ++beg;
            ModelParameter m(name,mean,errg,errf);
            ModelPars.push_back(m);
            if(beg!=tok.end()) std::cout << "warning: unread information in parameter " << m.name << std::endl;
        }
        else if(beg->compare("Observable")==0){
            ++beg;
            std::string name = *beg;
            ++beg;
            double min = atof((*beg).c_str());
            ++beg;
            double max = atof((*beg).c_str());
             ++beg;
            std::string toMCMC = *beg;
            bool tMCMC;
            if(toMCMC.compare("MCMC")==0)
                tMCMC = true;
            else if(toMCMC.compare("noMCMC")==0)
                tMCMC = false;
            else {
                std::cout << "wrong MCMC flag in " << name << std::endl;
                exit(EXIT_FAILURE);
            }
            ThObservable * tho= ThFactory(name);
            Observable o(name,tMCMC,min,max,tho);
            ++beg;
            std::string distr = *beg;
            if(distr.compare("file")==0) {
                ++beg;
                std::string fname = *beg;
                ++beg;
                std::string histoname = *beg;
                o.Set(distr,fname,histoname);
                Observables.push_back(o);
            }
            else if(distr.compare("weight")==0) {
                ++beg;
                double mean = atof((*beg).c_str());
                ++beg;
                double errg = atof((*beg).c_str());
                ++beg;
                double errf = atof((*beg).c_str());
                o.Set(distr,mean,errg,errf);
                Observables.push_back(o);
            }
            else if(distr.compare("noweight")==0) {
                o.Set(distr);
                Observables.push_back(o);
            }
            else {
                std::cout << "wrong distribution flag in " << name << std::endl;
                exit(EXIT_FAILURE);
            }
            ++beg;
            if(beg!=tok.end()) std::cout << "warning: unread information in observable " << Observables.back().name << std::endl;
        }
        else {
            std::cout << "wrong keyword in config file (first word must be ModelParameter or Observable)" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}
