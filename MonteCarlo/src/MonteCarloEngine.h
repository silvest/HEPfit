/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef __MONTECARLOENGINE__H
#define __MONTECARLOENGINE__H

#include "Observable.h"
#include "Observable2D.h"
#include "CorrelatedGaussianObservables.h"
#include "CorrelatedGaussianParameters.h"
#include "ModelParameter.h"
#include "Model.h"
#include <BAT/BCModel.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <TFile.h>
#include <TPrincipal.h>
#include <TColor.h>
#include <map>
#include <string>
#include <sstream>
#include <boost/ptr_container/ptr_vector.hpp>

#define NBINS1D 100
#define NBINS2D 100
#define NSTEPS 1e5

/**
 * @class MonteCarloEngine
 * @ingroup MonteCarlo
 * @brief An engine class for Monte Carlo.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The MonteCarlo engine is has overloaded from BCEngineMCMC class in the
 * <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a> libraries 
 * and user defined fumctions to facilitate the Monte Carlo run.
 */
class MonteCarloEngine : public BCModel {
public:

    /**
     * @brief Constructor.
     * @param[in] ModPars_i the vector of model parameters as defined in SomeModel.conf
     * @param[in] Obs_i the vector of observables as defined in SomeModel.conf
     * @param[in] Obs2D_i the vector of observables2D as defined in SomeModel.conf
     * @param[in] CGO_i the vector of correlated Gaussian observables as defined in SomeModel.conf
     */
    MonteCarloEngine(const std::vector<ModelParameter>& ModPars_i,
            boost::ptr_vector<Observable>& Obs_i,
            std::vector<Observable2D>& Obs2D_i,
            std::vector<CorrelatedGaussianObservables>& CGO_i,
            std::vector<CorrelatedGaussianParameters>& CGP_i);

    /**
     * @brief The default destructor. Some pointers defined in this class are explicitly freed.
     */
    ~MonteCarloEngine();

    /**
     * @brief Initialization of the Monte Carlo Engine.
     * @details The initialization of the Monte Carlo Engine performs the following tasks
     * \li The observables and observables2D lists are checked and sorted into the ones which are included
     * in the MCMC and the ones which are not.
     * \li The experimental likelihood for observables and observalbes2D are read if provided through a .root
     * file.
     * \li The correlated Gaussian observables are sorted for inclusion in the MCMC.
     * @param[in] Mod_i the pointer to the model defined in SomeModel.conf
     */
    void Initialize(StandardModel* Mod_i);
    
    /**
     * @brief Creation of the histogram maps for Observables, Observable2D and Correlated Gaussian Observable.
     */
    void CreateHistogramMaps();

    /**
     * @brief A member to classify the prior of the model parameters varied in the Monte Carlo
     * @details The model parameters being varied are first sorted out checking for the existence of 
     * of Gaussian error (\f$\delta_{g}\f$) or flat error (\f$\delta_{f}\f$). The SetPriorGaus method 
     * and SetPriorConstant method in BCEngineMCMC 
     * (<a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a>) is used to set a Gaussian
     * or flat prior for the model parameters respectively. In case a model parameter has both a Gaussian
     * and a flat error then the combined function is built using <a href="http://root.cern.ch/root/html/TF1.html"
     * target=blank>TF1</a>. The definition of the combined error is
     * \f[
     * \delta_{combined}(x) = \frac{1}{4\delta_f}\left(\frac{Erf(x - x_{ave}+\delta_f)}{\sqrt{2}\delta_g}
       - \frac{Erf(x - x_{ave}-\delta_f)}{\sqrt{2}\delta_g}\right),
     * \f]
     * with
     * \f[
     * Erf(x)=\frac{1}{\sqrt{2}\pi}\int_0^x \exp^{-\frac{t^2}{2}}dt.
     * \f]
     */
    void DefineParameters();

    /**
     * @brief This member calculates the loglikelihood for the observables included in the MCMC run.
     * @details The model is updated with the new set of parameters through the Model::Update() method.
     * If this is successful this even is counted as a used event. Otherwise it is counted as a 
     * discarded event. The log probability is calculated using Weight(const Observable&, const double& )
     * function for the observables, Weight(const Observable2D&, const double&, const double& ) for 
     * the observable2D and Weight(const CorrelatedGaussianObservables& ) for the correlated Gaussian
     * observables. Overloaded from BCEngineMCMC in <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a>.
     * @param[in] parameters the vector of the parameters that are being varied in the Monte Carlo
     * @return the loglikelihood
     */
    double LogLikelihood(const std::vector <double>& parameters);

    /**
     * @brief This member checks if there is overflow of the 1D histogram.
     * @param[in] hist a reference to an object of type <a href="http://root.cern.ch/root/html/TH1D.html"
     * target=blank>TH1D</a> as defined in the ROOT libraries
     * @param[in] name the name for the histogram
     */
    void CheckHistogram(TH1& hist, const std::string name);

    /**
     * @brief This member checks if there is overflow of the 2D histogram.
     * @param[in] hist a reference to an object of type <a href="http://root.cern.ch/root/html/TH2D.html"
     * target=blank>TH2D</a> as defined in the ROOT libraries
     * @param[in] name the name for the histogram
     */
    void CheckHistogram(TH2& hist, const std::string name);
    
    void Print1D(BCH1D hist1D, const char * filename, int ww=0, int wh=0);
    
    void Print2D(BCH2D hist2D, const char * filename, int ww=0, int wh=0);

    /**
     * @brief Overloaded from PrintHistogram(BCModelOutput&, const std::string) to print
     * histogram for observables.
     * @param[in] out a reference to an object of type BCModelOutput as defined in the
     * <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT libraries</a>
     * @param[in] it a reference to an object of type Observable
     * @param[in] OutputDir the name of the output directory
     */
    void PrintHistogram(std::string& OutFile, Observable & it, const std::string OutputDir);

    /**
     * @brief Member used for printing histograms for observables, observable2D, correlated Gaussian observables
     * and model parameters vs. observables.
     * @param[in] out a reference to an object of type BCModelOutput as defined in the
     * <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT libraries</a>
     * @param[in] OutputDir the name of the output directory
     */
    void PrintHistogram(std::string& OutFile, const std::string OutputDir);

    /**
     * @brief Overloaded from BCEngineMCMC in <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a>
     * @details The interface is used to update the model parameters using the Model::Update() method. Then the
     * values of the observables are computed and the respective histograms are filled. If the checkTheoryRange
     * flag is true then the minimum and maximum of the theory value is checked and reset to include the current
     * theory values.
     */
    void MCMCUserIterationInterface();

    /**
     * @brief A set method to fix the number of chains.
     * @details The number of chains are set using the MCMCSetNChains() from
     * BCEngineMCMC class in <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a>.
     * This also creates pointers to the vector of observable values and weights for all the chains.
     * @param[in] i the number of chains
     */
    void setNChains(unsigned int i);

    /**
     * @brief A method to add the observable values and weights to the chain information.
     * @details If the WriteChain flag is set to true in the MonteCarlo.conf file then chain information
     * is written. This method adds the observable values and weights (for the observables included in the MCMC)
     * to the chain information.
     */
    void AddChains();
    
    void InChainFillObservablesTree();

    /**
     * @brief A get method to access the stream that stores the log messages coming from histogram printing and checking.
     * @return a string containing the log messages
     */
    std::string getHistoLog() const
    {
        return HistoLog.str().c_str();
    }
    
    /**
     * @brief A get method to access the number of Observable chains requested.
     * @return an integer for the number of Observable chains
     */
    int getchainedObsSize() const
    {
        return kchainedObs;
    }

    /**
     * @brief A get method to compute the mean and rms of the computed observables.
     * @return a string containing the statistics
     */
    std::string computeStatistics();

    /**
     * @brief A method to write in a text file the best fit parameters and the prerun scale factors.
     * @return a string containing the data
     */
    std::string writePreRunData();

    /**
     * @brief This member generates the correlation matrix using BCH2D from the
     * <a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT libraries</a>.
     * @param[in] filename the name of the file where the correlation matrix is printed
     */
    void PrintCorrelationMatrixToLaTeX(const std::string filename);

    /**
     * @brief A get method to access the number of events discarded due to failure to update model.
     * These events are not used for the MCMC run.
     * @return the number of discarded events
     */
    int getNumOfDiscardedEvents() const
    {
        return NumOfDiscardedEvents;
    }

    /**
     * @brief A get method to access the number of events used in the MCMC run.
     * @return the number of events used in the MCMC run
     */
    int getNumOfUsedEvents() const
    {
        return NumOfUsedEvents;
    }
    
    /**
     * @brief A get method to access the map of 1D histograms.
     * @return the map of 1D histograms
     */
    std::map<std::string, BCH1D> getHistograms1D() const
    {
        return Histo1D;
    }

    /**
     * @brief A get method to access the map of 2D histograms.
     * @return the map of 2D histograms
     */
    std::map<std::string, BCH2D> getHistograms2D() const
    {
        return Histo2D;
    }
    
    /**
     * @brief A method to calculate the normalization based on the Laplace-Metropolis Estimator
     * @return the normalization
     */
    double computeNormalizationLME();
    
    /**
     * @brief A method to calculate the normalization based on the Monte Carlo Simulation
     * @return the normalization
     */
    std::vector<double> computeNormalizationMC(int NIterationNormalizationMC);
    
    /**
     * @brief A method to calculate the second derivative.
     * @return the second derivative
     */
    double SecondDerivative(BCParameter par1, BCParameter par2, std::vector<double> point);
    
    /**
     * @brief A method to calculate the first derivative.
     * @return the first derivative
     */
    double FirstDerivative(BCParameter par, std::vector<double> point);
    
    /**
     * @brief A method to calculate the LogLikelihood + LogAprioriProbability.
     * @param[in] point the set of points in the parameter space
     * @return LogLikelihood + LogAprioriProbability
     */
    double Function_h(std::vector<double> point);
    
    /**
     * @brief A method to rotate the diagonalized parameters to the original basis for correlated parameters.
     * @param[in] parameters the parameters in the diagonal basis
     * @param[in] DPars_i the parameters in the in the input basis
     */
    void setDParsFromParameters(const std::vector<double>& parameters, std::map<std::string,double>& DPars_i);
    
    /**
     * @brief A set method to toggle the printing of logo on the histogram pdf
     * @param[in] print a boolean to toggle the printing
     */
    void setPrintLogo(bool print)
    {
        printLogo = print;
    };
    
    /**
     * @brief A set method to toggle the printing of legends in 1D and 2D histograms.
     * @param[in] legend a boolean to toggle the the printing of legends
     */
    void setNoLegend(bool legend)
    {
        noLegend = legend;
    };
    
    /**
     * @brief A set method to toggle the printing of Parameters vs. Loglikelihood
     * @param[in] legend a boolean to toggle the the printing of the 2D plots
     */
    void setPrintLoglikelihoodPlots(bool plot)
    {
        PrintLoglikelihoodPlots = plot;
    };
    
    /**
     * @brief A set method to toggle the printing of legends in 1D and 2D histograms.
     * @param[in] legend a boolean to toggle the the printing of legends
     */
    void setAlpha2D(double alpha)
    {
        alpha2D = alpha;
    };
    
    /**
     * @brief A set method to set the number of smoothing passes for ROOT in 1D histograms
     * @param[in] int_N number of smoothing passes for ROOT
     */
    void setSmooth(int int_N)
    {
        nSmooth = int_N;
    };
    
    /**
     * @brief A set method to set the band fill type for 2D histograms
     * @param[in] type 2D Histogram type, 1001 -> Lego, 101 -> Filled, 0 -> Contour
     */
    void setHistogram2DType(int type)
    {
        histogram2Dtype = type;
    };
    
    /**
     * @brief A set method to set the number of bins for a 1D histograms
     * @param[in] nbins the number of bins in the 1D histogram
     */
    void setNBins1D(unsigned int nbins)
    {
        nBins1D = nbins;
    };
    
    /**
     * @brief A set method to set the number of bins for a 2D histograms
     * @param[in] nbins the number of bins in the 2D histogram
     */
    void setNBins2D(unsigned int nbins)
    {
        nBins2D = nbins;
    };
    
    /**
     * @brief A set method to set the number of significant digits in the output
     * @param[in] significants_i the number of significant digits in the output
     */
    void setSignificants(unsigned int significants_i)
    {
        significants = significants_i;
    };
    
    /**
     * @brief A set method to set the size of the buffer used by the histograms
     * @param[in] histogramBufferSize_i the size of the buffer used by the histograms
     */
    void setHistogramBufferSize(unsigned int histogramBufferSize_i)
    {
        histogramBufferSize = histogramBufferSize_i;
    };
    
    /**
     * @brief A get method to get the number of MPI world size.
     * @return the size of the MPI world
     */    
    int getMPIWorldSize() 
    {
#ifdef _MPI
        return procnum;
#else
        return 1;
#endif
    }
    
private:
    
    /**
     * @brief A get method to calculate the precision of the numbers in the Statistics file based on the rms.
     * @param[in] value the mean value
     * @param[in] rms the rms value
     * @return the number of digits of precision
     */
    int getPrecision(double value, double rms);
    
    const std::vector<ModelParameter>& ModPars; ///< A vector of model parameters.
    const std::vector<CorrelatedGaussianParameters>& CGP; ///< A vector of correlated Gaussian parameters.
    boost::ptr_vector<Observable>& Obs_ALL; ///< A vector of all observables.
    std::vector<Observable2D>& Obs2D_ALL; ///< A vector of all pairs of observable for Observable2D.
    std::vector<CorrelatedGaussianObservables>& CGO; ///< A vector of correlated Gaussian observables.
    StandardModel* Mod; ///< A pointer to an abject of type Model.
    std::map<std::string, double> DPars; ///< A map between parameter names and their values.
    std::vector<std::map<std::string, double> > DPars_allChains; ///< A vector of maps between parameter names and their values for all chains.
    std::map<std::string, BCH1D> Histo1D; ///< A map between pointers to objects of type BCH1D (<a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a>) and their given names.
    std::map<std::string, BCH2D> Histo2D; ///< A map between pointers to objects of type BCH2D (<a href="https://www.mppmu.mpg.de/bat/?page=home" target=blank>BAT</a>) and their given names.
    std::map<std::string, double> thMin; ///< A map between the name of a theory observable and its minimum value.
    std::map<std::string, double> thMax; ///< A map between the name of a theory observable and its maximum value.
    std::map<std::string, TPrincipal *> CorrelationMap; ///< A map between the name of a theory observable and its maximum value.
    double *obval; ///< A pointer to the vector of observable values.
    double *obweight; ///< A pointer to the vector of observable weights.
    unsigned int kwmax; ///< The number of observables whose weights are used for the MCMC.
    unsigned int kmax; ///< The number of observables.
    unsigned int kchainedObs; ///< The number of observables for which the chains should be written.
    std::ostringstream HistoLog; ///< A stream to store the output messages from printing and checking histograms.
    int NumOfUsedEvents; ///< The number of events for which the model is successfully updated and hence used for the MCMC run.
    int NumOfDiscardedEvents; ///< The number of events for which the update of the model fails and these events are not used for the MCMC run.
    int rank; ///< Rank of the process for a MPI run. Value is 0 for a serial run.
    TTree * hMCMCObservableTree; ///< A ROOT tree that contains the observables values and weight when the chains are written.
    std::vector<std::vector<double> > hMCMCObservables; ///< A vector of vectors containing the observables values of all the chains to be put into the ROOT tree.
    std::vector<std::vector<double> > hMCMCObservables_weight; ///< A vector of vectors containing the observables weight of all the chains to be put into the ROOT tree.
    std::vector<double> hMCMCTree_Observables; ///< A vector containing the observables values to be put into the ROOT tree.
    std::vector<double> hMCMCTree_Observables_weight; ///< A vector containing the observables weight to be put into the ROOT tree.
    unsigned int cindex;///< An index to distinguish between succesive canvases used to draw histograms.
    std::ofstream ofi;
    std::vector<std::string> unknownParameters; ///< A vector to contain the unkenown parameters passed in the configuration file.
    bool printLogo; ///< A flag that is set to true for printing the logo on the histogram pdf.
    int nSmooth; ///< The number of times a 1D  histogram should be smoothed.
    int histogram2Dtype; ///< Type of 2D Histogram 1001 -> box pixel, 101 -> filled, 1 -> contour.
    bool noLegend; ///< A flag to toggle the histogram legends.
    bool PrintLoglikelihoodPlots; ///< A flag to toggle the printing of Parameters vs. Loglikelihood.
    double alpha2D;///< A number between 0. and 1. that sets the opacity level of 2D Histograms, 1. being fully opaque.
    int gIdx;
    int rIdx;
    unsigned int nBins1D; ///< The number of bins in a 1D histogram.
    unsigned int nBins2D; ///< The number of bins in a 2D histogram.
    unsigned int significants; /// < The number of significant digits in the Statistics File.
    TColor * HEPfit_green; /// < The colour green for HEPfit.
    TColor * HEPfit_red; /// < The colour red for HEPfit.
    unsigned int histogramBufferSize; /// < The size of the buffer used for the histograms.
    std::vector<double> par_at_LL_max; /// < vector to store the value of the parameters at maximum LogLikelihood;
    double LogLikelihood_max; /// < Maximum value of LogLikelihood
};

#endif

