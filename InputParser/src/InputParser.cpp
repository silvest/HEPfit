/*
 * Copyright (C) 2012 HEPfit Collaboration
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

InputParser::InputParser(ModelFactory& ModF, ThObsFactory& ObsF) : myModelFactory(ModF), myObsFactory(ObsF), filename(""), rank(0)
{
    modelset = 0;
}

InputParser::InputParser(const InputParser& orig) : myModelFactory(orig.myModelFactory), myObsFactory(orig.myObsFactory), filename(""), rank(0)
{
    myModel = new StandardModel(*orig.myModel);
}

InputParser::~InputParser()
{}

std::string InputParser::ReadParameters(const std::string filename_i,
        const int rank_i,
        std::vector<ModelParameter>& ModelPars,
        boost::ptr_vector<Observable>& Observables,
        std::vector<Observable2D>& Observables2D,
        std::vector<CorrelatedGaussianObservables>& CGO,
        std::vector<CorrelatedGaussianParameters>& CGP)
{
    filename = filename_i;
    rank = rank_i;
    modname = "";
    lineNo = 0;
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open()) {
        if(rank == 0) throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid model configuration file.\n");
        else sleep (2);
    }

    if (filename.find("\\/") == std::string::npos) filepath = filename.substr(0, filename.find_last_of("\\/") + 1);
    IsEOF = false;
    do {
        IsEOF = getline(ifile, line).eof();
        lineNo++;
        if (*line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
        if (line.empty() || line.find_first_not_of(' ') == std::string::npos || line.at(0) == '#')
            continue;
        sep = new boost::char_separator<char>(" \t");
        tok = new boost::tokenizer<boost::char_separator<char> >(line, *sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg = tok->begin();

        if (modelset == 0) {
            modname = *beg;
            myModel = myModelFactory.CreateModel(modname);
            myModel->setModelName(modname);
            myModel->InitializeModel();
            if (myModel->IsModelInitialized()) {
                if (rank == 0) std::cout << "\nModel Initialized: " << modname << std::endl;
                modeldefinedinfile = filename;
            } else if (rank == 0)
                throw std::runtime_error("\nERROR: " + modname + " not initialized successfully.\n");
            modelset = 1;
            continue;
        } else if (modelset == 1 && beg->compare(myModel->ModelName()) == 0) {
            continue;
        }

        std::string type = *beg;
        ++beg;
        if (type.compare("ModelParameter") == 0) {
            
            if (std::distance(tok->begin(), tok->end()) < 5) {
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on " + *beg + " in " + filename + ".\n");
                else sleep(2);
            } else {
            ModelParameter tmpMP;
            beg = tmpMP.ParseModelParameter(beg);
            
            if (checkDuplicateParameter[tmpMP.getname()].get<0>()) {
                if(rank == 0) throw std::runtime_error("\nERROR: ModelParameter " + tmpMP.getname() + " appears more than once ...!! \n" +
                "1st Occurrence: Line No:" + boost::lexical_cast<std::string>(checkDuplicateParameter[tmpMP.getname()].get<2>()) +
                " in file " + checkDuplicateParameter[tmpMP.getname()].get<1>() + ".\n"
                "2nd Occurrence: Line No:" + boost::lexical_cast<std::string>(lineNo) + " in file " + filename + ".\n");
                else sleep (2);
            }
            
            if (beg != tok->end())
                if (rank == 0) std::cout << "WARNING: unread information in parameter " << tmpMP.getname() << std::endl;
            checkDuplicateParameter[tmpMP.getname()] = boost::make_tuple(true, filename, lineNo);
               
            ModelPars.push_back(tmpMP);
            }

        } else if (type.compare("CorrelatedGaussianParameters") == 0) {
            
            CorrelatedGaussianParameters tmpCGP;
            lineNo = tmpCGP.ParseCGP(ModelPars, filename, ifile, beg, lineNo, rank);
            IsEOF = tmpCGP.isEOF();
            CGP.push_back(tmpCGP);

        } else if (type.compare("Observable") == 0 || type.compare("BinnedObservable") == 0 || type.compare("FunctionObservable") == 0) {
            
            Observable * tmpObs = new Observable();
            beg = tmpObs->ParseObservable(type, tok, beg, filepath, filename, rank);
            tmpObs->setTho(myObsFactory.CreateThMethod(tmpObs->getThname(), *myModel));
            Observables.push_back(tmpObs);

        } else if (type.compare("Observable2D") == 0) {
            
            Observable2D tmpObs2D;
            lineNo = tmpObs2D.ParseObservable2D(type, tok, beg, filename, ifile, lineNo, rank);
            tmpObs2D.setTho1Tho2(myObsFactory.CreateThMethod(tmpObs2D.getThname(), *myModel), myObsFactory.CreateThMethod(tmpObs2D.getThname2(), *myModel));
            IsEOF = tmpObs2D.isEOF();
            Observables2D.push_back(tmpObs2D);

        } else if (type.compare("HiggsObservable") == 0) {
            
            Observable * tmphObs = new Observable();
            beg = tmphObs->ParseObservable(type, tok, beg, filepath, filename, rank);
            tmphObs->setTho(myObsFactory.CreateThMethod(tmphObs->getThname(), *myModel));
            HiggsObservable * tmpho = new HiggsObservable(*tmphObs);
            beg = tmpho->ParseHiggsObservable(beg, myObsFactory, myModel, rank);
            Observables.push_back(tmpho);
            ++beg;
            if (beg != tok->end() && rank == 0) std::cout << "WARNING: unread information in HiggsObservable " << tmpho->getName() << std::endl;

        } else if (type.compare("CorrelatedGaussianObservables") == 0) {
            CorrelatedGaussianObservables tmpCGO;
            lineNo = tmpCGO.ParseCGO(Observables, ifile, beg, filename, myObsFactory, myModel, lineNo, rank);
            IsEOF = tmpCGO.isEOF();
            if (tmpCGO.getObs().size() > 1) CGO.push_back(tmpCGO);

        } else if (type.compare("CustomObservable") == 0) {
//            if (std::distance(tok->begin(), tok->end()) < 2) {
//                if (rank == 0) throw std::runtime_error("ERROR: lack of information on " + *beg + " in " + filename + ".\n");
//            else sleep(2);
//        }
//            std::string customObsName = *beg;
//            beg++;
//            std::string customParserName = ObservableToParsermap[customObsName];
//            if (rank == 0 && customObservableTypeMap.find(customObsName) == customObservableTypeMap.end()) throw std::runtime_error("\nERROR: No Observable Type defined for " + customObsName + "\n");
//            if (rank == 0 && customParserMap.find(customParserName) == customParserMap.end()) throw std::runtime_error("\nERROR: No parser defined for " + customObsName + "\n");
//            if (rank == 0 && ObservableToParsermap.find(customObsName) == ObservableToParsermap.end()) throw std::runtime_error("\nERROR: No parser linked to the observable " + customObsName + "\n");
//            InputParser* customParser = CreateCustomParser(customObsName);
//            customParser->setModel(myModel);
//            Observable * tmp = customParser->ParseObservable(type, beg);
//            Observable * customObs = CreateObservableType(customObsName, *tmp);
//            customObs->setObsType(customObsName);
//            Observables.push_back(customObs);
//            delete customParser;
        } else if (type.compare("ModelFlag") == 0) {
            if (std::distance(tok->begin(), tok->end()) < 3) {
                if(rank == 0) throw std::runtime_error("ERROR: lack of information on " + *beg + " in " + filename);
                else sleep (2);
            }
            std::string flagname = *beg;
            ++beg;
            if (boost::iequals(*beg, "true") || boost::iequals(*beg, "false")) {
                /* Boolean flags */
                bool value_bool;
                if (boost::iequals(*beg, "true"))
                    value_bool = 1;
                else
                    value_bool = 0;
                if (!myModel->setFlag(flagname, value_bool)) {
                    if(rank == 0) throw std::runtime_error("ERROR: setFlag error for " + flagname);
                    else sleep (2);
                }
                else if (rank == 0) std::cout << "set flag " << flagname << "=" << *beg << std::endl;
            } else {
                /* String flags */
                std::string value_str = *beg;
                if (!myModel->setFlagStr(flagname, value_str)) {
                    if(rank == 0) throw std::runtime_error("ERROR: setFlag error for " + flagname);
                    else sleep (2);
                }
                else if (rank == 0) std::cout << "set flag " << flagname << "=" << value_str << std::endl;
            }
            ++beg;
            if (beg != tok->end() && rank == 0) std::cout << "WARNING: unread information in Flag " << flagname << std::endl;
        } else if (type.compare("IncludeFile") == 0) {
            std::string IncludeFileName = filepath + *beg;
            if (rank == 0) std::cout << "\nIncluding File: " + IncludeFileName << std::endl;
            ReadParameters(IncludeFileName, rank, ModelPars, Observables, Observables2D, CGO, CGP);
            IsEOF = false;
            ++beg;
        } else {
            if (rank == 0) throw std::runtime_error("\nERROR: wrong keyword " + type + " in file " + filename + " line no. " + boost::lexical_cast<std::string>(lineNo) + ". Make sure to specify a valid model configuration file.\n");
            else sleep(2);
        }

    } while (!IsEOF);

    if (modelset == 0 && rank == 0)
        throw std::runtime_error("ERROR: Incorrect or missing model name in the model configuration file.\n");
    if (!myModel->CheckFlags() && rank == 0)
        throw std::runtime_error("ERROR: incompatible flag(s)\n");

    return (modname);
}

void InputParser::addCustomParser(const std::string name, boost::function<InputParser*(ModelFactory& ModF, ThObsFactory& ObsF) > funct)
{
    customParserMap[name] = funct;
}

void InputParser::addCustomObservableType(const std::string name, boost::function<Observable*(Observable& obs_i) > funct)
{
    customObservableTypeMap[name] = funct;
}

void InputParser::linkParserToObservable(std::string name_obs, std::string name_par)
{
    ObservableToParsermap[name_obs] = name_par;
}

InputParser * InputParser::CreateCustomParser(const std::string& name) const
{
    if (customParserMap.find(ObservableToParsermap.at(name)) == customParserMap.end())
        throw std::runtime_error("\nERROR: No parser defined for " + ObservableToParsermap.at(name) + " so it cannot be created.\n");
    return (customParserMap.at(ObservableToParsermap.at(name))(myModelFactory, myObsFactory));
}

Observable * InputParser::CreateObservableType(const std::string& name, Observable& obs_i) const
{
    if (customObservableTypeMap.find(name) == customObservableTypeMap.end())
        throw std::runtime_error("ERROR: No observable defined for " + name + " so it cannot be created");
    return (customObservableTypeMap.at(name)(obs_i));
}