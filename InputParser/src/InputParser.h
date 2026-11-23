/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef INPUTPARSER_H
#define	INPUTPARSER_H

#include "ThObsFactory.h"
#include "ModelFactory.h"
#include "Observable.h"
#include "Observable2D.h"
#include "CorrelatedGaussianObservables.h"
#include "CorrelatedGaussianParameters.h"
#include "HiggsObservable.h"
#include "ThObservable.h"
#include "ModelParameter.h"
#include "StandardModel.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <string>

/**
 * @addtogroup InputParser 
 * @brief A module for interfaces to model parameters and observables.
 * @details The %InputParser class parses all the input given to the executable
 * through the SomeModel.conf file and the output directives in the same. \n
 * The %ThObsFactory defines and initiates all available observables that can
 * be computed.
 * @{
 */

/**
 * @class InputParser
 * @brief A class for reading input parameters and output directives
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The InputParser parses the SomeModel.conf file for the name of the model, ModelParameters, 
 * %Observables, * %Observables2D, %ModelParaVsObs and %CorrelatedGaussianObservables. The observables 
 * (or set of observables) are parsed for their inclusion in the Monte Carlo run.
 */
class InputParser {
public:
    
    /**
     * @brief Constructor.
     */
    InputParser(ModelFactory& ModF, ThObsFactory& ObsF);
    
    /**
     * @brief The copy constructor.
     * @param[in] orig a reference to an object of type InputParser()
     */
    InputParser(const InputParser& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~InputParser();
    
    /**
     * @brief The member that parses the Observable2D directives from SomeModel.conf file
     * @param[in] Observables a vector of Observables
     * @param[in] ifile the file being read
     * @param[in] beg an iterator over words in a line separated by a specific separator character
     */
//    CorrelatedGaussianObservables ParseCGO(boost::ptr_vector<Observable>& Observables, std::ifstream& ifile, boost::tokenizer<boost::char_separator<char> >::iterator& beg);
    
    /**
     * @brief Responsible for parsing the SomeModel.conf file.
     * @details This method parses the SomeModel.conf file for all input instructions. The algorithm
     * is as follows:
     * \li Search for the predefined names of the available models, initialize the necessary model
     * class and build a ThFactory for that model to define the observables.
     * \li Search and read the ModelParameter list.
     * \li Search and read the CorrelatedGaussianParameters.
     * \li Search and read the Observable list using the ParseObservable() method.
     * \li Search and read the Observable2D list using the ParseObservable() method.
     * \li Search and read the CorrelatedGaussianObservables.
     * \li Return the name of the model initialized.
     *
     * A line "IncludeFile <file>" reads <file> at that point; a relative name is relative to the directory of
     * the file containing the directive, an absolute name is used as it is. A file that includes itself, directly
     * or through other included files, is an error.
     *
     * Observable names must be unique, including the members of correlated sets (histograms, statistics and
     * chains are keyed by name); the same holds for the names of the Observable2D. A repeated name is reported
     * with a warning (unlike repeated ModelParameters, which are an error), since the outputs of observables
     * sharing a name are merged.
     * @param[in] filename name of the SomeModel.conf file passed as an argument to the executable
     * @param[in] rank the rank of the process in an MPI run (set to 0 for a serial run)
     * @param[out] ModelPars the vector of model parameters
     * @param[out] Observables the vector of observables
     * @param[out] Observables2D the vector of observable pairs
     * @param[out] CGO the vector of correlated Gaussian observables
     * @param[out] CGP the vector of correlated Gaussian parameters
     * @return modname the name of the model initialized
     */
    std::string ReadParameters(const std::string filename_i,
                               const int rank,
                               std::vector<ModelParameter>& ModelPars,
                               boost::ptr_vector<Observable>& Observables,
                               std::vector<Observable2D>& Observables2D,
                               std::vector<CorrelatedGaussianObservables>& CGO,
                               std::vector<CorrelatedGaussianParameters>& CGP);

    /**
     * @brief A get method to access the pointer to the object of the StandardModel class.
     * @return myModel a pointer to the object of the StandardModel class
     */
    
    StandardModel* getModel() const
    {
        return myModel;
    }
    
    const ThObsFactory& getObsFactory() const
    {
        return myObsFactory;
    }

    ModelFactory& getModelFactory() const {
        return myModelFactory;
    }
    
    void addCustomObservableType(const std::string name, boost::function<Observable*() > funct);
    
    Observable * CreateObservableType(const std::string& name) const;
    
    void setModel(StandardModel* myModel_i){
        myModel = myModel_i;
    }
    
private:
    /**
     * @brief Warns if an observable name has already been booked, otherwise records where it was booked.
     * @param[in,out] seen the names booked so far, with the file and line number where they were defined
     * @param[in] name the name of the observable
     * @param[in] kind the kind of observable, for the warning
     * @param[in] line the line number of the definition in the current file
     */
    void CheckDuplicateObservable(std::map<std::string, std::pair<std::string, int> >& seen, const std::string& name,
                                  const std::string& kind, const int line);

    StandardModel* myModel;///< Pointer to an object of type StandardModel.
    ModelFactory& myModelFactory;///< Pointer to an object of type ModelFactory.
    ThObsFactory& myObsFactory;///< Reference to an object of type ThObsFactory.
    std::string modname;///< A string to store the model name in.
    bool modelset;
    std::map <std::string, boost::tuple<bool, std::string, int> > checkDuplicateParameter;
    std::map<std::string, std::pair<std::string, int> > checkDuplicateObservable;///< Observable names booked so far (file, line).
    std::map<std::string, std::pair<std::string, int> > checkDuplicateObservable2D;///< Observable2D names booked so far (file, line).
    std::string modeldefinedinfile;
    //std::map<std::string, boost::function<InputParser*(ModelFactory&, ThObsFactory&) > > customParserMap;
    std::map<std::string, boost::function<Observable* ()> > customObservableTypeMap;
    //std::map<std::string, std::string> ObservableToParsermap;
    int lineNo;
    std::string filename;
    int rank;
    std::string filepath;
    std::vector<std::string> includeStack;///< The canonical names of the files being read (nested by IncludeFile).
};

/** 
 * @}
 */

#endif	/* INPUTPARSER_H */
