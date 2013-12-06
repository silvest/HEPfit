/*
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef __MONTECARLOENGINE__H
#define __MONTECARLOENGINE__H

#include <BAT/BCModel.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCModelOutput.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TRandom3.h>
#include <Observable.h>
#include <Observable2D.h>
#include <CorrelatedGaussianObservables.h>
#include <ModelParaVsObs.h>
#include <ModelParameter.h>
#include <Model.h>
#include <map>
#include <string>
#include <sstream>

#define NBINSMODELPARS 100
#define NBINS1D 100
#define NBINS2D 100

/**
 * @class MonteCarloEngine
 * @ingroup MonteCarlo
 * @brief An engine class for Monte Carlo.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class MonteCarloEngine : public BCModel {
public:
    
    /**
     * @brief The default constructor.
     * @param[in] ModPars_i the vector of model parameters as defined in SomeModel.conf
     * @param[in] Obs_i the vector of observables as defined in Somemodel.conf
     * @param[in] Obs2D_i the vector of observables2D as defined in SomeModel.conf
     * @param[in] CGO_i the vector of correlated Gaussian observables as defined in SomeModel.conf
     * @param[in] ParaObs_i the vector of parameter vs. observables as defined in SomeModel.conf
     * @param[in] checkHistRange_i
     */
    MonteCarloEngine(const std::vector<ModelParameter>& ModPars_i,
                     std::vector<Observable>& Obs_i,
                     std::vector<Observable2D>& Obs2D_i,
                     std::vector<CorrelatedGaussianObservables>& CGO_i,
                     std::vector<ModelParaVsObs>& ParaObs_i,
                     const bool checkHistRange_i=false);
    
    /**
     * @brief The default destructor.
     */
    ~MonteCarloEngine();
    
    /**
     * @brief Initialization of the Monte Carlo Engine
     * @param[in] Mod_i the pointer to the model defined in SomeModel.conf
     */
    void Initialize(Model* Mod_i);
    
    /**
     * @brief A member to classify the prior of the model parameters varied in the Monte Carlo
     */
    void DefineParameters();
    
    /**
     * @brief
     * @param[in] parameters the vector of the parameters that are being varied in the Monte Carlo
     * @return the loglikelihood
     */
    double LogLikelihood(const std::vector <double>& parameters);
    
    /**
     * @brief
     * @param[in] hist a reference to a object of type <a href="http://root.cern.ch/root/html/TH1D.html"
     * target=blank>TH1D</a> as defined in the ROOT libraries
     * @param[in] name the name for the histogram
     */
    void CheckHistogram(const TH1D& hist, const std::string name);
    
    /**
     * @brief
     * @param[in] hist a reference to a object of type <a href="http://root.cern.ch/root/html/TH2D.html"
     * target=blank>TH2D</a> as defined in the ROOT libraries
     * @param[in] name the name for the histogram
     */
    void CheckHistogram(const TH2D& hist, const std::string name);
    
    /**
     * @brief
     * @param[in] out a reference to an object of tyoe BCModelOutput as defined in the
     * <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT libraries</a>
     * @param[in] it a iterator of the vector of objects of type Observable()
     * @param[in] OutputDir the name of the output directory
     */
    void PrintHistogram(BCModelOutput & out,
                        std::vector<Observable>::iterator it,
                        const std::string OutputDir);
    
    /**
     * @brief
     * @param[in] out a reference to an object of tyoe BCModelOutput as defined in the
     * <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT libraries</a>
     * @param[in] OutputDir the name of the output directory
     */
    void PrintHistogram(BCModelOutput& out, const std::string OutputDir);
    
    /**
     * @brief
     * @return
     */
    void MCMCIterationInterface();
    
    /**
     * @brief
     * @param[in]
     * @return
     */
    void setNChains(unsigned int i);
    
    /**
     * @brief
     * @return
     */
    void AddChains();
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
     */
    double Weight(const Observable& obs, const double& th);
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @param[in]
     * @return
     */
    double Weight(const Observable2D& obs, const double& th1, const double& th2);
    
    /**
     * @brief
     * @param[in]
     * @return
     */
    double Weight(const CorrelatedGaussianObservables& obs);
    
    /**
     * @brief
     * @return
     */
    std::string getHistoLog() const
    {
        return HistoLog.str().c_str();
    }
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
     */
    void PrintCorrelationMatrix(const std::string filename);
    
    /**
     * @brief
     * @return
     */
    int getNumOfDiscardedEvents() const
    {
        return NumOfDiscardedEvents;
    }
    
    /**
     * @brief
     * @return
     */
    int getNumOfUsedEvents() const
    {
        return NumOfUsedEvents;
    }
    
private:
    const std::vector<ModelParameter>& ModPars; /**< */
    std::vector<Observable> Obs_MCMC; /**< */
    std::vector<Observable>& Obs_ALL;  /**< */
    std::vector<Observable2D> Obs2D_MCMC; /**< */
    std::vector<Observable2D>& Obs2D_ALL; /**< */
    std::vector<CorrelatedGaussianObservables>& CGO; /**< */
    std::vector<ModelParaVsObs>& ParaObs; /**< */
    Model* Mod; /**< */
    std::map<std::string, double> DPars; /**< */
    std::map<std::string, BCH1D * > Histo1D; /**< */
    std::map<std::string, TH1D * > InHisto1D; /**< */
    std::map<std::string, BCH2D * > Histo2D; /**< */
    std::map<std::string, TH2D * > InHisto2D; /**< */
    std::map<std::string, double> thMin; /**< */
    std::map<std::string, double> thMax; /**< */
    double *obval; /**< */
    double *obweight; /**< */
    unsigned int kwmax; /**< */
    unsigned int kmax; /**< */
    std::ostringstream HistoLog; /**< */
    int NumOfUsedEvents; /**< */
    int NumOfDiscardedEvents; /**< */
    bool checkTheoryRange; /**< */
    
};

#endif

