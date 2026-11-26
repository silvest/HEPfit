/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "InputParser.h"
#include "ModelFactory.h"
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <algorithm>
#include <cstdlib>

namespace {
// Canonical name of a file, used to detect recursive IncludeFile directives; the name itself if it cannot be resolved.
std::string canonicalFileName(const std::string& name)
{
    char* resolved = realpath(name.c_str(), NULL);
    if (resolved == NULL) return name;
    std::string canonical(resolved);
    free(resolved);
    return canonical;
}
}

InputParser::InputParser(ModelFactory& ModF, ThObsFactory& ObsF) : myModelFactory(ModF), myObsFactory(ObsF), filename(""), rank(0)
{
    modelset = false;
    myModel = NULL;
}

InputParser::InputParser(const InputParser& orig) : myModelFactory(orig.myModelFactory), myObsFactory(orig.myObsFactory), filename(""), rank(0)
{
    myModel = new StandardModel(*orig.myModel);
}

InputParser::~InputParser()
{
    if (myModel != NULL) {
        delete myModel;
        myModel = NULL;
    }
}

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
    lineNo = 0;
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open()) {
        if(rank == 0) throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid model configuration file.\n");
        else sleep (2);
    }

    // files being read (nested by IncludeFile), to reject a file that includes itself directly or indirectly
    includeStack.push_back(canonicalFileName(filename));
    struct PopIncludeStack {
        std::vector<std::string>& stack;
        ~PopIncludeStack() { stack.pop_back(); }
    } popIncludeStack = {includeStack};

    filepath = filename.substr(0, filename.find_last_of("\\/") + 1);
    const boost::char_separator<char> sep(" \t");
    std::string line;
    while (std::getline(ifile, line)) {
        lineNo++;

        if (!line.empty() && line.back() == '\r') line.pop_back(); // for CR+LF
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        auto beg = tok.begin();
        if (beg == tok.end() || beg->front() == '#')
            continue;

        if (modelset == 0) {
            modname = *beg;
            myModel = myModelFactory.CreateModel(modname);
            myModel->setModelName(modname);
            myModel->InitializeModel();
            if (myModel->IsModelInitialized()) {
                if (rank == 0) std::cout << "\nModel Initialized: " << modname << "\n" << std::endl;
                modeldefinedinfile = filename;
            } else if (rank == 0)
                throw std::runtime_error("\nERROR: " + modname + " not initialized successfully.\n");
            modelset = 1;
            continue;
        } else if (modelset == 1 && beg->compare(myModel->getModelName()) == 0) {
            continue;
        }

        std::string type = *beg;
        ++beg;
        const int keywordLineNo = lineNo;
        const size_t nObservables = Observables.size(), nObservables2D = Observables2D.size(), nCGO = CGO.size();
        if (type.compare("ModelParameter") == 0) {
            
            ModelParameter tmpMP;
            try {
                tmpMP.ParseModelParameter(beg);
            } catch (const std::runtime_error& err) {
                throw std::runtime_error(filename + ":" + std::to_string(lineNo) + ": " + err.what());
            }
            if (checkDuplicateParameter[tmpMP.getname()].get<0>()) {
                if(rank == 0) throw std::runtime_error("\nERROR: ModelParameter " + tmpMP.getname() + " appears more than once ...!! \n" +
                "1st Occurrence: Line No:" + boost::lexical_cast<std::string>(checkDuplicateParameter[tmpMP.getname()].get<2>()) +
                " in file " + checkDuplicateParameter[tmpMP.getname()].get<1>() + ".\n"
                "2nd Occurrence: Line No:" + boost::lexical_cast<std::string>(lineNo) + " in file " + filename + ".\n");
                else sleep (2);
            }
            
            if (tmpMP.getname().compare("lambdaNP") == 0 && !tmpMP.IsFixed())
                if(rank == 0) throw std::runtime_error("\nERROR: ModelParameter " + tmpMP.getname() + " cannot float in the MonteCarlo run ");
            
            if (beg != tok.end() && beg->front() != '#')
                if (rank == 0) std::cout << "WARNING: unread information in parameter " << tmpMP.getname() << std::endl;
            checkDuplicateParameter[tmpMP.getname()] = boost::make_tuple(true, filename, lineNo);
               
            ModelPars.push_back(tmpMP);

        } else if (type.compare("CorrelatedGaussianParameters") == 0) {
            
            CorrelatedGaussianParameters tmpCGP;
            lineNo = tmpCGP.ParseCGP(ModelPars, filename, ifile, beg, lineNo, rank);
            if (!tmpCGP.getDiagPars().empty()) CGP.push_back(tmpCGP);

        } else if (type.compare("Observable") == 0 || type.compare("BinnedObservable") == 0 || type.compare("FunctionObservable") == 0 || type.compare("AsyGausObservable") == 0) {
            
            Observable * tmpObs = new Observable();
            beg = tmpObs->ParseObservable(type, &tok, beg, filepath, filename, rank);
            tmpObs->setTho(myObsFactory.CreateThMethod(tmpObs->getThname(), *myModel));
            Observables.push_back(tmpObs);

        } else if (type.compare("Observable2D") == 0) {
            
            Observable2D tmpObs2D;
            lineNo = tmpObs2D.ParseObservable2D(type, &tok, beg, filename, ifile, lineNo, rank);
            tmpObs2D.setTho1Tho2(myObsFactory.CreateThMethod(tmpObs2D.getThname(), *myModel), myObsFactory.CreateThMethod(tmpObs2D.getThname2(), *myModel));
            Observables2D.push_back(tmpObs2D);

        } else if (type.compare("HiggsObservable") == 0) {
            
            Observable * tmphObs = new Observable();
            beg = tmphObs->ParseObservable(type, &tok, beg, filepath, filename, rank);
            tmphObs->setTho(myObsFactory.CreateThMethod(tmphObs->getThname(), *myModel));
            HiggsObservable * tmpho = new HiggsObservable(*tmphObs);
            beg = tmpho->ParseHiggsObservable(beg, myObsFactory, myModel, rank);
            Observables.push_back(tmpho);
            ++beg;
            if (beg != tok.end() && rank == 0) std::cout << "WARNING: unread information in HiggsObservable " << tmpho->getName() << std::endl;

        } else if (type.compare("CorrelatedHiggsObservables") == 0) {
            
            Observable * tmphObs = new Observable();
            beg = tmphObs->ParseObservable(type, &tok, beg, filepath, filename, rank);
            tmphObs->setTho(myObsFactory.CreateThMethod(tmphObs->getThname(), *myModel));
            HiggsObservable * tmpho = new HiggsObservable(*tmphObs);
            tmpho->setIsCorrelated(true);
            beg = tmpho->ParseHiggsObservable(beg, myObsFactory, myModel, rank);
            Observables.push_back(tmpho);
            ++beg;
            if (beg != tok.end() && rank == 0) std::cout << "WARNING: unread information in HiggsObservable " << tmpho->getName() << std::endl;

        } else if (type.compare("CorrelatedGaussianObservables") == 0) {
            
            CorrelatedGaussianObservables tmpCGO;
            lineNo = tmpCGO.ParseCGO(Observables, ifile, beg, filename, myObsFactory, myModel, lineNo, rank);
            if (tmpCGO.getObs().size() > 1) CGO.push_back(tmpCGO);
            
        } else if (type.compare("ObservablesWithCovarianceInverse") == 0) {
            
            CorrelatedGaussianObservables tmpCGO;
            tmpCGO.setCovarianceFromConfig(true);
            lineNo = tmpCGO.ParseCGO(Observables, ifile, beg, filename, myObsFactory, myModel, lineNo, rank);
            if (tmpCGO.getObs().size() > 1) CGO.push_back(tmpCGO);    

        } else if (type.compare("CorrelatedObservables") == 0) {
            
            CorrelatedGaussianObservables tmpCO;
            tmpCO.setIsPrediction(true);
            lineNo = tmpCO.ParseCGO(Observables, ifile, beg, filename, myObsFactory, myModel, lineNo, rank);
            if (tmpCO.getObs().size() > 1) CGO.push_back(tmpCO); // a single member is booked as a normal Observable by ParseCGO()

        } else if (type.compare("CustomObservable") == 0) {
            
            if (std::distance(tok.begin(), tok.end()) < 2) {
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on " + *beg + " in " + filename + ".\n");
            else sleep(2);
            }
            std::string customObsName = *beg;
            beg++;
            if (customObservableTypeMap.find(customObsName) == customObservableTypeMap.end()) {
                if (rank == 0) throw std::runtime_error("\nERROR: No Observable Type defined for " + customObsName + "\n");
                else sleep(2);
            }
            Observable * tmpcustomObs = CreateObservableType(customObsName);
            tmpcustomObs->setObsType(customObsName);
            beg = tmpcustomObs->ParseObservable(type, &tok, beg, filepath, filename, rank);
            tmpcustomObs->setTho(myObsFactory.CreateThMethod(tmpcustomObs->getThname(), *myModel));
            Observables.push_back(tmpcustomObs);
            ++beg;
            if (beg != tok.end() && rank == 0) std::cout << "WARNING: unread information in CustomObservable " << tmpcustomObs->getName() << std::endl;

        } else if (type.compare("ModelFlag") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 3) {
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
                } else if (rank == 0) std::cout << "set flag " << flagname << "=" << value_str << std::endl;
            }
            ++beg;
            if (beg != tok.end() && rank == 0) std::cout << "WARNING: unread information in Flag " << flagname << std::endl;
        } else if (type.compare("IncludeFile") == 0) {
            if (beg == tok.end()) {
                if (rank == 0) throw std::runtime_error("\nERROR: IncludeFile without a file name in file " + filename + " line no. " + boost::lexical_cast<std::string>(lineNo) + ".\n");
                else sleep(2);
                continue;
            }
            std::string oldfilepath = filepath;
            std::string oldfilename = filename;
            int oldlineNo = lineNo;
            // a relative name is relative to the directory of the including file; an absolute name is used as it is
            std::string IncludeFileName = ((*beg).at(0) == '/' ? *beg : filepath + *beg);
            if (std::find(includeStack.begin(), includeStack.end(), canonicalFileName(IncludeFileName)) != includeStack.end()) {
                if (rank == 0) throw std::runtime_error("\nERROR: recursive IncludeFile of " + IncludeFileName + " in file " + filename + " line no. " + boost::lexical_cast<std::string>(lineNo) + ".\n");
                else sleep(2);
                continue;
            }
            if (rank == 0) std::cout << "Including File: " + IncludeFileName << std::endl;
            ReadParameters(IncludeFileName, rank, ModelPars, Observables, Observables2D, CGO, CGP);
            filepath = oldfilepath;
            filename = oldfilename;
            lineNo = oldlineNo;
            ++beg;
        } else {
            if (rank == 0) throw std::runtime_error("\nERROR: wrong keyword " + type + " in file " + filename + " line no. " + boost::lexical_cast<std::string>(lineNo) + ". Make sure to specify a valid model configuration file.\n");
            else sleep(2);
        }

        // Observable names must be unique (histograms, statistics and chains are keyed by name).
        // The observables of an included file have been checked while reading it.
        if (type.compare("IncludeFile") != 0) {
            for (size_t i = nObservables; i < Observables.size(); i++)
                CheckDuplicateObservable(checkDuplicateObservable, Observables[i].getName(), "Observable", keywordLineNo);
            for (size_t i = nCGO; i < CGO.size(); i++) {
                const std::vector<Observable> members = CGO[i].getObs();
                for (size_t j = 0; j < members.size(); j++)
                    CheckDuplicateObservable(checkDuplicateObservable, members[j].getName(), "Observable", keywordLineNo);
            }
            for (size_t i = nObservables2D; i < Observables2D.size(); i++)
                CheckDuplicateObservable(checkDuplicateObservable2D, Observables2D[i].getName(), "Observable2D", keywordLineNo);
        }
    }

    if (modelset == 0 && rank == 0)
        throw std::runtime_error("ERROR: Incorrect or missing model name in the model configuration file.\n");
    if (rank == 0 && !myModel->CheckFlags())
        throw std::runtime_error("ERROR: incompatible flag(s)\n");
    return (modname);
}

void InputParser::CheckDuplicateObservable(std::map<std::string, std::pair<std::string, int> >& seen, const std::string& name,
                                           const std::string& kind, const int line)
{
    std::map<std::string, std::pair<std::string, int> >::const_iterator it = seen.find(name);
    if (it != seen.end()) {
        // A repeated name is only reported: several shipped configurations reuse a name for different bins.
        if (rank == 0) std::cout << "\nWARNING: " << kind << " " << name << " appears more than once: "
                << "first at line " << it->second.second << " in file " << it->second.first << ", again at line "
                << line << " in file " << filename << ". Histograms, statistics and chain branches of observables"
                << " with the same name are merged; rename one of them to keep them separate.\n" << std::endl;
    } else
        seen[name] = std::make_pair(filename, line);
}

void InputParser::addCustomObservableType(const std::string name, boost::function<Observable*() > funct)
{
    customObservableTypeMap[name] = funct;
}

Observable * InputParser::CreateObservableType(const std::string& name) const
{
    if (customObservableTypeMap.find(name) == customObservableTypeMap.end()) {
        if (rank ==0) throw std::runtime_error("ERROR: No observable defined for " + name + " so it cannot be created");
        else sleep(0);
    }
    return (customObservableTypeMap.at(name)());
}
