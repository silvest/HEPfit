/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERATEEVENT_H
#define	GENERATEEVENT_H

#include "InputParser.h"
#include "Observable.h"
#include "Observable2D.h"
#include "CorrelatedGaussianObservables.h"
#include "ModelParameter.h"
#include "Model.h"
#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"
#include <TF1.h>
#ifdef _MPI
#include <mpi.h>
#endif

/**
 * @addtogroup EventGeneration
 * @brief A module for generating events based on the ROOT random number generator.
 * @details This module is for using the implementations of the observables without
 * running a Markov Chain Monte Carlo. It allows generation of events based on the
 * ROOT random number generator.
 * @{
 */

/**
 * @class GenerateEvent
 * @ingroup EventGeneration
 * @brief A class for generating events.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class GenerateEvent {
public:
    /**
     * @brief Constructor.
     * @details The default constructor sets the names of the configuration file
     * and the names of the specific output directory along with a job ID if it is
     * specified as an argument to the executable. Possible arguments to the executable are:
     * \li SomeModel.conf --noMC : For single event (central event) printed to terminal only.
     * \li SomeModel.conf --noMC --it=10 : For the central event and 10 events printed to terminal only.
     * \li SomeModel.conf --noMC --it=10 --output_folder FOLDER: For the central event and 10 
     * events printed to GenerateEvents/FOLDER only.
     * @param[in] ModelConf_i the name of the input configuration file for the model name,
     * the model parameters and observables to be calculated.
     * @param[in] OutDirName_i the name of the root output directory to be given
     * @param[in] JobTag_i optional job tag that might be specified
     * @param[in] noMC_i the noMC specification given to the executable=
     */
    GenerateEvent(ModelFactory& ModelF, ThObsFactory& ThObsF, 
                  const std::string& ModelConf_i,
                  const std::string& OutDirName_i,
                  const std::string& JobTag_i,
                  const bool noMC_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~GenerateEvent();
    
    /**
     * @brief The method used to generate events and format output
     * @details This member creates all out directories that might be needed in case the 
     * output is stored on disk. Output on disk is generated as follows:
     * \li Existence of GenerateEvents folder is checked or else created in the current directory.
     * \li The user specified OutDirName_i and a directory called OLD is created under 
     * GenerateEvents if it does not already exist.
     * \li If OutDirName_i exists, it is moved into OLD.
     * \li If OutDirName_i exits in OLD is deleted before the existing GenerateEvents/OutDirName_i
     * is moved into OLD.
     * \li In the now empty GenerateEvents/OutDirName_i the directories Parameters and Observables
     * are created. The names indicate what they will contain.\n
     * The file Summary.txt contains the list of parameters that are varied and observables that 
     * are generated in a run. This is in the GenerateEvents/OutDirName_i directory.\n
     * @param[in] rank the rank of the process for MPI implementations (0 for serial implementation)
     * @param[in] nIteration number of events generated
     * @param[in] seed seed for the random number generator (default = 0)
     */
    void generate(int unsigned nIteration_i, int seed = 0, bool weight = false);
    
    void addCustomObservableType(const std::string name, boost::function<Observable*() > funct);
    
private:
    /**
     * @brief The parameter distributions as specified in SomeModel.conf is defined here.
     * @details The parameters in the SomeModel.conf either have a flat error or a Gaussian
     * error and might even have both. This member defines the distributions of the parameters
     * using the <a href="http://root.cern.ch/root/htmldoc/TF1.html" target=blank>TF1</a> class
     * in ROOT.
     */
    void defineParameterDistributions();
    
    /**
     * @brief This member generates random numbers for the parameters being varied in the model.
     * The first set is always the central value of the parameter.
     * @param[in] iterationNo the iteration number
     */
    void generateRandomEvent(int iterationNo);
    
    void createDirectories();
    
    void initModel();
    
    void setDParsFromParameters(const std::vector<double>& parameters, std::map<std::string,double>& DPars_i);
    
    InputParser myInputParser; ///< An oject of the InputParser() class.
    std::map<std::string, double> DPars; ///< Map of parameters to be passed to Model().
    std::map<std::string, TF1*> DDist; ///< Map of parameter distributions.
    std::map<std::string, boost::shared_ptr<std::ofstream> > ObsOut;///< Map of output stream for observables.
    std::map<std::string, boost::shared_ptr<std::ofstream> > CGOOut;///< Map of output stream for corellated Gaussian observables.
    std::map<std::string, boost::shared_ptr<std::ofstream> > ParsOut;///< Map of output stream for parameters.
    Model* Mod; ///< Name of the model as defined in SomeModel.conf
    std::vector<ModelParameter> ModPars; ///< Vector for the model parameters defined in SomeModel.conf.
    std::vector<ModelParameter> ModParsVar; ///< Vector for the model parameters varied in SomeModel.conf.
    boost::ptr_vector<Observable> Obs; ///< Vector for the observables defined in SomeModel.conf.
    std::vector<Observable2D> Obs2D; ///< Vector for the Observables2D defined in SomeModel.conf.
    std::vector<CorrelatedGaussianObservables> CGO; ///< vector for the Correlated Gaussian Observables defined in SomeModel.conf.
    std::vector<CorrelatedGaussianParameters> CGP; ///< vector for the Correlated Gaussian Parameters defined in SomeModel.conf.
    std::string ModelConf; ///< String for the name of the SomeModel.conf file.
    std::string OutDirName; ///< String for the name of the output directory.
    std::string OldOutDirName; ///< String for the name of the backup output directory.
    std::string ObsDirName; ///< String for the name of the observables output directory.
    std::string CGODirName; ///< String for the name of the Correlated Gaussian Observables output directory.
    std::string SMDebugDirName; ///< String for the name of the observables output directory.
    std::string SUSYDebugDirName; ///< String for the name of the Correlated Gaussian Observables output directory.
    std::string ParsDirName; ///< String for the name of the parameters output directory.
    std::string JobTag; ///< String for the optional JobTag argument to be passes to the executable.
    bool noMC;///< Flag to initiate noMC mode.
    bool weight;///< Flag to get weight for observables.
    bool outputTerm; ///< Flag to specify output stream storage.
    int rank;
    int procnum;
#ifdef _MPI
    int nameLen;
    char processorName[MPI_MAX_PROCESSOR_NAME];
#endif
    std::ofstream summary;
    int nIteration;
    int positionID;
    int buffersize;
    double *sendbuff;
    double *sendbuff_w;
    double **buff;
    double **buff_w;
    int *sendbuff_int;
    int **buff_int;
    std::string ModelName;
};

#endif	/* GENERATEEVENT_H */

