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
     * @param[in]
     * @param[in]
     * @param[in]
     * @param[in]
     * @param[in]
     * @param[in]
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
     * @brief
     * @param[in]
     * @return
     */
    void Initialize(Model* Mod_i);
    
    /**
     * @brief
     * @return
     */
    void DefineParameters();
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
     */
    double LogLikelihood(const std::vector <double>& parameters);
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
     */
    void CheckHistogram(const TH1D& hist, const std::string name);
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
     */
    void CheckHistogram(const TH2D& hist, const std::string name);
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @param[in]
     * @return
     */
    void PrintHistogram(BCModelOutput & out, std::vector<Observable>::iterator it,
                        const std::string OutputDir);
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
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

