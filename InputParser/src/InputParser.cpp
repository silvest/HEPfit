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

void InputParser::ParseHiggsObservable(HiggsObservable * ho, boost::tokenizer<boost::char_separator<char> >::iterator & beg)
{
    std::vector<ThObservable*> hthobs;
    std::string filename = *(++beg);
    std::ifstream ifile(filename.c_str());
    if (!ifile.is_open())
        throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid Higgs parameters configuration file.\n");
    std::string line;
    bool IsEOF;
    int i = 0, nrows = 0;
    do {
        IsEOF = getline(ifile, line).eof();
        if (line.compare(0, 10, "CATEGORIES") == 0) {
            boost::tokenizer<boost::char_separator<char> > tok(line.substr(11));
            boost::tokenizer<boost::char_separator<char> >::iterator beg2 = tok.begin();
            while (beg2 != tok.end()) {
                std::string cat = *beg2;
                hthobs.push_back(myObsFactory.CreateThMethod(cat, *myModel));
                beg2++;
            }
        }
    } while (!IsEOF);
    if (hthobs.size() > 0)
        ho->setParametricLikelihood(filename, hthobs);
    else
        throw std::runtime_error("\nERROR: " + ho->getName() + " does not provide at least one category\n");
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
        if (rank == 0) throw std::runtime_error("\nERROR: " + filename + " does not exist. Make sure to specify a valid model configuration file.\n");
    std::string filepath;
    if (filename.find("\\/") == std::string::npos) filepath = filename.substr(0, filename.find_last_of("\\/") + 1);
    std::string line;
    bool IsEOF = false;
    do {
        IsEOF = getline(ifile, line).eof();
        lineNo += 1;
        if (*line.rbegin() == '\r') line.erase(line.length() - 1); // for CR+LF
        if (line.empty() || line.find_first_not_of(' ') == std::string::npos || line.at(0) == '#')
            continue;
        boost::char_separator<char> sep(" \t");
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
            if (std::distance(tok.begin(), tok.end()) < 5)
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                        + *beg + " in " + filename);
            std::string name = *beg;
            if (checkDuplicateParameter[name].get<0>() && rank == 0) throw std::runtime_error("\nERROR: ModelParameter " + name + " appears more than once ...!! \n" +
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
            checkDuplicateParameter[name] = boost::make_tuple(true, filename, lineNo);
        } else if (type.compare("Observable") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 8)
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                        + *beg + " in " + filename);
            Observable * o = new Observable(ParseObservable(beg));
            ++beg;
            std::string distr = *beg;
            if (distr.compare("file") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 10)
                    if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                            + *beg + " in " + filename);
                std::string fname = filepath + *(++beg);
                std::string histoname = *(++beg);
                o->setLikelihoodFromHisto(fname, histoname);
                if (rank == 0) std::cout << "added input histogram "  << fname << "/" << histoname << std::endl;
            } else if (distr.compare("weight") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 11)
                    if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                            + *beg + " in " + filename);
                ++beg;
                o->setAve(atof((*beg).c_str()));
                ++beg;
                o->setErrg(atof((*beg).c_str()));
                ++beg;
                o->setErrf(atof((*beg).c_str()));
                if (o->getErrf() == 0. && o->getErrg() == 0.) {
                    if (rank == 0) throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + o->getName() + " cannot both be 0. in the " + filename + " file.\n");
                }
            } else if (distr.compare("noweight") == 0) {
            } else if (rank == 0)
                throw std::runtime_error("ERROR: wrong distribution flag in " + o->getName());
            o->setDistr(distr);
            Observables.push_back(o);
            ++beg;
            if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in observable "
                        << Observables.back().getName() << std::endl;
        } else if (type.compare("BinnedObservable") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 10)
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                        + *beg + " in " + filename);
            Observable * bo = new Observable(ParseObservable(beg));
            bo->setObsType(2);
            ++beg;
            std::string distr = *beg;
            if (distr.compare("file") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 12)
                    if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                            + *beg + " in " + filename);
                std::string fname = filepath + *(++beg);
                std::string histoname = *(++beg);
                bo->setLikelihoodFromHisto(fname, histoname);
                if (rank == 0) std::cout << "added input histogram "  << fname << "/" << histoname << std::endl;
            } else if (distr.compare("weight") == 0) {
                if (std::distance(tok.begin(), tok.end()) < 13)
                    if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                            + *beg + " in " + filename);
                ++beg;
                bo->setAve(atof((*beg).c_str()));
                ++beg;
                bo->setErrg(atof((*beg).c_str()));
                ++beg;
                bo->setErrf(atof((*beg).c_str()));
                if (bo->getErrf() == 0. && bo->getErrg() == 0.) {
                    if (rank == 0) throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + bo->getName() + " cannot both be 0. in the " + filename + " file.");
                }
            } else if (distr.compare("noweight") == 0) {
                ++beg;
                ++beg;
                ++beg;
            } else if (rank == 0)
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
            if (std::distance(tok.begin(), tok.end()) < 12) {
                std::string name = *beg;
                ++beg;
                if (std::distance(tok.begin(), tok.end()) < 4 && rank == 0) throw std::runtime_error("ERROR: lack of information on "
                        + name + " in " + filename + " at line number" + boost::lexical_cast<std::string>(lineNo));
                std::string toMCMC = *beg;
                bool tMCMC;
                if (toMCMC.compare("MCMC") == 0)
                    tMCMC = true;
                else if (toMCMC.compare("noMCMC") == 0)
                    tMCMC = false;
                else
                    throw std::runtime_error("ERROR: wrong MCMC flag in Observable2D" + name + " at line number:" + boost::lexical_cast<std::string>(lineNo) + " of file " + filename);
                ++beg;
                std::string distr = *beg;
                std::string fname = filepath + *(++beg);
                std::string histoname = *(++beg);

                std::vector<double> min(2, 0.);
                std::vector<double> max(2, 0.);
                std::vector<double> ave(2, 0.);
                std::vector<double> errg(2, 0.);
                std::vector<double> errf(2, 0.);
                std::vector<double> bin_min(2, 0.);
                std::vector<double> bin_max(2, 0.);
                std::vector<std::string> thname(2, "");
                std::vector<std::string> label(2, "");
                std::vector<std::string> type2D(2, "");
                size_t pos = -1;
                for (int i = 0; i < 2; i++) {
                    IsEOF = getline(ifile, line).eof();
                    if (line.empty() || line.at(0) == '#') {
                        if (rank == 0) std::cout << "ERROR: no comments or empty lines in Observables2D (" + name + ") please!"
                                << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    lineNo++;
                    boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                    beg = mytok.begin();
                    type2D[i] = *beg;
                    if (type2D[i].compare("Observable") != 0 && type2D[i].compare("BinnedObservable") != 0)
                        if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + filename + ", expecting an Observable or BinnedObservable type here...\n");
                    ++beg;
                    thname[i] = *beg;
                    ++beg;
                    label[i] = *beg;
                    while ((pos = label[i].find("~", pos + 1)) != std::string::npos)
                        label[i].replace(pos, 1, " ");
                    ++beg;
                    min[i] = atof((*beg).c_str());
                    ++beg;
                    max[i] = atof((*beg).c_str());
                    if (distr.compare("weight") == 0) {
                        ++beg;
                        ave[i] = atof((*beg).c_str());
                        ++beg;
                        errg[i] = atof((*beg).c_str());
                        ++beg;
                        errf[i] = atof((*beg).c_str());
                        if (errg[i] == 0. && errg[i] == 0.) {
                            if (rank == 0) throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + name + " cannot both be 0. in the " + filename + " file, line number:" + boost::lexical_cast<std::string>(lineNo));
                        }
                    } else if (distr.compare("noweight") == 0) {
                        if (type2D[i].compare("BinnedObservable") == 0) {
                            ++beg;
                            ++beg;
                            ++beg;
                        }
                    }
                    if (type2D[i].compare("BinnedObservable") == 0) {
                        ++beg;
                        bin_min[i] = atof((*beg).c_str());
                        ++beg;
                        bin_max[i] = atof((*beg).c_str());
                    }
                }
                Observable2D o2(name, thname[0], thname[1], label[0], label[1], tMCMC, min[0], max[0], min[1], max[1], myObsFactory.CreateThMethod(thname[0], *myModel), myObsFactory.CreateThMethod(thname[1], *myModel));
                o2.setDistr(distr);
                o2.setObsType(type2D[0]);
                o2.setObsType2(type2D[1]);
                o2.setAve(ave[0]);
                o2.setAve2(ave[1]);
                o2.setErrg(errg[0]);
                o2.setErrg2(errg[1]);
                o2.setErrf(errf[0]);
                o2.setErrf2(errf[1]);
                if (distr.compare("file") == 0) {
                    o2.setLikelihoodFromHisto(fname, histoname);
                    if (rank == 0) std::cout << "added input histogram " << fname << "/" << histoname << std::endl;
                }
                if (o2.getObsType() == 2) {
                    o2.getTho()->setBinMin(bin_min[0]);
                    o2.getTho()->setBinMax(bin_max[0]);
                }
                if (o2.getObsType2() == 2) {
                    o2.getTho2()->setBinMin(bin_min[1]);
                    o2.getTho2()->setBinMax(bin_max[1]);
                }
                Observables2D.push_back(o2);
            } else {
                Observable2D o2(ParseObservable(beg));
                ++beg;
                std::string distr = *beg;
                if (distr.compare("file") == 0) {
                    if (std::distance(tok.begin(), tok.end()) < 14)
                        if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
                                + *beg + " in " + filename);
                    std::string fname = filepath + *(++beg);
                    std::string histoname = *(++beg);
                    o2.setLikelihoodFromHisto(fname, histoname);
                    if (rank == 0) std::cout << "added input histogram "  << fname << "/" << histoname << std::endl;
                } else if (distr.compare("noweight") == 0) {
                } else if (rank == 0)
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
            }
        } else if (type.compare("HiggsObservable") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 8)
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
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
                } else if (rank == 0)
                    throw std::runtime_error("ERROR: wrong keyword " + distr + " in " + ho->getName());
                ho->setParametricLikelihood(*(++beg), hthobs);
            }
            else if (distr.compare("new_parametric") == 0) {
                ParseHiggsObservable(ho, beg);
            } else if (rank == 0)
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
                lineNo++;
                boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                beg = mytok.begin();
                std::string type = *beg;
                ++beg;
                if (type.compare("Observable") != 0 && type.compare("BinnedObservable") != 0)
                    if (rank == 0) throw std::runtime_error("ERROR: in line no." + boost::lexical_cast<std::string>(lineNo) + " of file " + filename + ", expecting an Observable or BinnedObservable type here...\n");
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
                    if (tmp->getErrf() == 0. && tmp->getErrg() == 0.) {
                        if (rank == 0) throw std::runtime_error("ERROR: The Gaussian and flat error in weight for " + tmp->getName() + " cannot both be 0. in the " + filename + " file, line number:" + boost::lexical_cast<std::string>(lineNo));
                    }
                } else if (distr.compare("noweight") == 0) {
                    if (type.compare("BinnedObservable") == 0) {
                        ++beg;
                        ++beg;
                        ++beg;
                    }
                } else if (rank == 0)
                    throw std::runtime_error("ERROR: wrong distribution flag in " + tmp->getName());

                tmp->setDistr(distr);
                if (type.compare("BinnedObservable") == 0) {
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
            if (nlines > 1) {
                gslpp::matrix<double> myCorr(gslpp::matrix<double>::Id(nlines));
                int ni = 0;
                for (int i = 0; i < size; i++) {
                    IsEOF = getline(ifile, line).eof();
                    if (line.empty() || line.at(0) == '#') {
                        if (rank == 0) std::cout << "ERROR: no comments or empty lines in CorrelatedGaussianObservables please!"
                                << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    lineNo++;
                    if (lines.at(i)) {
                        boost::tokenizer<boost::char_separator<char> > mytok(line, sep);
                        beg = mytok.begin();
                        int nj = 0;
                        for (int j = 0; j < size; j++) {
                            if ((*beg).compare(0, 1, "0") == 0
                                    || (*beg).compare(0, 1, "1") == 0
                                    || (*beg).compare(0, 1, "-") == 0) {
                                if (std::distance(mytok.begin(), mytok.end()) < size && rank == 0) throw std::runtime_error(("ERROR: Correlation matrix is of wrong size in Correlated Gaussian Observables: " + name).c_str());
                                if (lines.at(j)) {
                                    myCorr(ni, nj) = atof((*beg).c_str());
                                    nj++;
                                }
                                beg++;
                            } else {
                                if (rank == 0) std::cout << "ERROR: invalid correlation matrix for "
                                        << name << ". Check element (" << ni + 1 << "," << nj + 1 << ") in line number " + boost::lexical_cast<std::string>(lineNo) << std::endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                        ni++;
                    }
                }
                o3.ComputeCov(myCorr);
                CGO.push_back(o3);
            } else {
                if (rank == 0) std::cout << "\nWARNING: Correlated Gaussian Observable " << name.c_str() << " defined with less than two correlated observables. The set is being marked as normal Observables." << std::endl;
                if (o3.getObs().size() == 1) Observables.push_back(new Observable(o3.getObs(0)));
                for (int i = 0; i < size; i++) {
                    IsEOF = getline(ifile, line).eof();
                    lineNo++;
                }
            }
        } else if (type.compare("ModelFlag") == 0) {
            if (std::distance(tok.begin(), tok.end()) < 3)
                if (rank == 0) throw std::runtime_error("ERROR: lack of information on "
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
                if (!myModel->setFlag(flagname, value_bool)) {
                    if (rank == 0) throw std::runtime_error("ERROR: setFlag error for " + flagname);
                } else if (rank == 0) std::cout << "set flag " << flagname << "=" << *beg << std::endl;
            } else {
                /* String flags */
                std::string value_str = *beg;
                if (!myModel->setFlagStr(flagname, value_str)) {
                    if (rank == 0) throw std::runtime_error("ERROR: setFlag error for " + flagname);
                } else if (rank == 0) std::cout << "set flag " << flagname << "=" << value_str << std::endl;
            }
            ++beg;
            if (beg != tok.end())
                if (rank == 0) std::cout << "WARNING: unread information in Flag " << flagname << std::endl;
        } else if (type.compare("IncludeFile") == 0) {
            std::string IncludeFileName = filepath + *beg;
            if (rank == 0) std::cout << "\nIncluding File: " + IncludeFileName << std::endl;
            ReadParameters(IncludeFileName, rank, ModelPars, Observables, Observables2D, CGO);
            ++beg;
        } else if (rank == 0)
            throw std::runtime_error("\nERROR: wrong keyword " + type + " in file " + filename + ". Make sure to specify a valid model configuration file.");

    } while (!IsEOF);

    if (modelset == 0 && rank == 0)
        throw std::runtime_error("ERROR: Incorrect or missing model name in the model configuration file.\n");
    if (!myModel->CheckFlags() && rank == 0)
        throw std::runtime_error("ERROR: incompatible flag(s)\n");

    return (modname);
}


