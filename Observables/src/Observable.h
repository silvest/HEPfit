/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBSERVABLE_H
#define	OBSERVABLE_H

#include "ThObservable.h"
#include <string>
#include <iostream>
#include <TH1D.h>
#include <boost/tokenizer.hpp>

/**
 * @class Observable
 * @ingroup Observable
 * @brief A class for observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The class for building an observable and storing its different 
 * parameters read from the SomeModel.conf file or specified by the user. The name (thname) of the observable has
 * to correspond to the allowed name of observables listed in the ThFactory class.
 */
class Observable {
public:
    /**
     * @brief Constructor.
     * @param[in] name_i a given name for the observable
     * @param[in] thname_i the thname for the observable fixed in ThFactory()
     * @param[in] label_i the label assigned to the observable
     * @param[in] tMCMC_i boolean flag to indicate inclusion in MCMC
     * @param[in] min_i minimum value for the observable
     * @param[in] max_i maximum value for the observable
     * @param[in] tho_i a pointer to an object of type ThObservable
     */
    Observable(const std::string name_i,
               const std::string thname_i,
               const std::string label_i,
               const bool tMCMC_i,
               const double min_i,
               const double max_i,
               ThObservable * tho_i);
    
    /**
     * @brief The copy constructor.
     */
    Observable(const Observable& orig);
    
    /**
     * @brief The default constructor.
     */
    Observable();
    
    /**
     * @brief The parser for Observables.
     * @param[in] type the string specifying the type of the observable
     * @param[in] tok the tokenizer containing the line being parsed
     * @param[in] beg the iterator that parses a line in the config file
     * @param[in] filepath the path to the config file being parsed
     * @param[in] filename the name of the config file being parsed
     * @return the line number (integer) after the parsing is done
     */
    boost::tokenizer<boost::char_separator<char> >::iterator &  ParseObservable(std::string& type,
                                                                                boost::tokenizer<boost::char_separator<char> >* tok, 
                                                                                boost::tokenizer<boost::char_separator<char> >::iterator & beg, 
                                                                                std::string& filepath, 
                                                                                std::string& infilename,
                                                                                int rank);
    
    /**
     * @brief The default destructor.
     */
    virtual ~Observable();
    
    /**
     * @brief A method to access the computed theory value of the observable.
     */
    double computeTheoryValue();

    /**
     * @brief A method to compute the weight associated with the observable.
     * @param[in] th the theoretical value of the observable
     */
    virtual double computeWeight(double th);
    
    /**
     * @brief A method to compute the weight associated with the observable.
     * @param[in] th the theoretical value of the observable
     * @param[in] ave_i the average value of the observable
     * @param[in] errg_i the Gaussian error of the observable
     * @param[in] errf_i the flat error of the observable
     */
    virtual double computeWeight(double th, double ave_i, double errg_i, double errf_i);
    
    /**
     * @brief A method to compute the weight associated with the observable.
     * @param[in] th1 the theoretical value of the first observable
     * @param[in] th2 the theoretical value of the second observable
     */
    virtual double computeWeight(double th1, double th2)
    {
        return 0.0;
    };

    /**
     * @brief A method to compute the weight associated with the observable.
     */
    virtual double computeWeight()
    {
        return computeWeight(computeTheoryValue());
    }

    /**
     * @brief A get method to access the average value of the observable.
     * @return the average value of the observable
     */
    double getAve() const
    {
        return ave;
    }

    /**
     * @brief A set method to fix the average value of the observable.
     * @param[in] ave the average value of the observable
     */
    void setAve(double ave)
    {
        this->ave = ave;
    }

    /**
     * @brief A get method to access the name of the distribution of the observable.
     * @return the name of the distribution of the observable
     */
    std::string getDistr() const
    {
        return distr;
    }

    /**
     * @brief A set method to fix the name of the distribution of the observable.
     * @param[in] distr the name of the distribution of the observable
     */
    void setDistr(std::string distr)
    {
        this->distr = distr;
    }

    /**
     * @brief A get method to access the flat error of the observable.
     * @return the flat error of the observable
     */
    double getErrf() const
    {
        return errf;
    }

    /**
     * @brief A set method to fix the flat error of the observable.
     * @param[in] errf the flat error of the observable
     */
    void setErrf(double errf)
    {
        this->errf = errf;
    }

    /**
     * @brief A get method to access the Gaussian error of the observble.
     * @return the Gauissian error of the observable
     */
    double getErrg() const
    {
        return errg;
    }

    /**
     * @brief A set method to fix the gaussian error of the observable.
     * @param[in] errg the Gaussian error of the observable
     */
    void setErrg(double errg)
    {
        this->errg = errg;
    }

    /**
     * @brief A get method to access the filename of the observables experimental likelihood file.
     * @return the name of the file
     */
    std::string getFilename() const
    {
        return filename;
    }
    
    void setFilename(std::string filename_i)
    {
        filename = filename_i;
    }
    
    /**
     * @brief A set method to set the likelihood from which the experimental likelihood of the observable will
     * be read.
     * @param filename the name of the file
     * @param histoname the name of the histogram
     */
    virtual void setLikelihoodFromHisto(std::string filename, std::string histoname);

    /**
     * @brief A set method to set a parametric likelihood reading parameters from a file
     * @param filename the name of the file
     */
    virtual void setParametricLikelihood(std::string filename)
    {
        this->filename=filename; //real implementation will be done in extension
    }
    
    /**
     * Set the parametric likelihood to be overloaded by HiggsObservable.
     * @param filename the name of the config file
     * @param thObsV a vector of ThObservables
     */
    virtual void setParametricLikelihood(std::string filename, std::vector<ThObservable*> thObsV)
    {};

    /**
     * @brief A get method to access the name for the histogram of the observable.
     * @return the name of the histogram for the observable
     */
    std::string getHistoname() const
    {
        return histoname;
    }
    
    /**
     * @brief A set method to set the name of the histogram containing the likelihood
     * @param[in] histoname_i a string that contains the name of the histogram
     */
    void setHistoname(std::string histoname_i)
    {
        histoname = histoname_i;
    }

    /**
     * @brief A get method to access the label for the observable.
     * @return the label for the observable
     */
    std::string getLabel() const
    {
        return label;
    }

    /**
     * @brief A set method to fix the label for the observable.
     * @param[in] label the label for the observable
     */
    void setLabel(std::string label)
    {
        this->label = label;
    }

    /**
     * @brief A get method to access the maximum value of the observable
     * @return the maximum value of the observable
     */
    double getMax() const
    {
        return max;
    }

    /**
     * @brief A set method to fix the maximum value for the observable.
     * @param[in] max the maximum value for the observable
     */
    void setMax(double max)
    {
        this->max = max;
    }
    
    /**
     * @brief A get method to access the minimum value of the observable.
     * @return the minimum value of the observable
     */
    double getMin() const 
    {
        return min;
    }
    
    /**
     * @brief A set method to fix the minimum value for the observable.
     * @param[in] min the minimum value for the observable
     */
    void setMin(double min) 
    {
        this->min = min;
    }

    /**
     * @brief A get method to access the name of the observable.
     * @return the name of the observable
     */
    std::string getName() const
    {
        return name;
    }

    /**
     * @brief A set method to fix the name for the observable.
     * @param name for the observable
     */
    void setName(std::string name)
    {
        this->name = name;
    }

    /**
     * @brief A method to check if the observable is listed for MCMC.
     * @return true or false
     */
    bool isTMCMC() const
    {
        return tMCMC;
    }

    /**
     * @brief A set method to fix the observable's inclusion in the MCMC listing.
     * @param[in] tMCMC true or false
     */
    void setTMCMC(bool tMCMC)
    {
        this->tMCMC = tMCMC;
    }

    /**
     * @brief A get method to access the thname of the observable as defined in ThFactory class.
     * @return thname the name of the observable as listed in ThFactory class
     */
    std::string getThname() const
    {
        return thname;
    }

    /**
     * @brief A set method to fix the name of the observable as listed in ThFactory class.
     * @param[in] thname the name of the observable as listed in ThFactory class
     */
    void setThname(std::string thname)
    {
        this->thname = thname;
    }

    /**
     * @brief A get method to access the pointer to the object of the ThObservable class.
     * @return pointer to the object of type ThObservable
     */
    ThObservable* getTho() const
    {
        return tho;
    }
    
    /**
     * @brief A set method to set the Observable type
     * @param[in] obsType_s a string that contains the parameter name
     */
    void setObsType(std::string& obsType_s)
    {
        obsType = obsType_s;
    }
    
    /**
     * @brief A get method to get the Observable type.
     * @return a string containing the observable type
     */
    std::string getObsType() const
    {
        return obsType;
    }

    /**
     * @brief A set method to fix the pointer to object of type ThObservable.
     * @param[in] tho pointer to the object of type ThObservable
     */
    void setTho(ThObservable* tho_i)
    {
        tho = tho_i;
        tho->setBinMin(bin_min);
        tho->setBinMax(bin_max);
    }
    
    /**
     * @brief A set method to fix the pointer to object of type ThObservable.
     * @param[in] tho pointer to the object of type ThObservable
     */
    void setTho(ThObservable* tho_i, double bmin, double bmax)
    {
        tho = tho_i;
        tho->setBinMin(bmin);
        tho->setBinMax(bmax);
    }
    
    /**
     * A method to compute the log of a split Gaussian likelihood
     * @param x the value of the split-Gaussian distributed variable
     * @param ave the average value
     * @param errl the left-side error
     * @param errr the right-side error
     * @return the log likelihood at point x
     */
    double LogSplitGaussian(double x, double ave, double errl, double errr);
    
    /**
     * A method to compute the log of a Gaussian likelihood
     * @param x the value of the Gaussian distributed variable
     * @param ave the average value
     * @param sigma the error
     * @return the log likelihood at point x
     */
    double LogGaussian(double x, double ave, double sigma);
    
    

    /**
     * @brief Befriending of the std::ostream operator << to generate an
     * output stream for printing the observables details.
     * @param[out] output the formatted output stream to print the model parameters
     * @param[in] o a reference to an object of type Observable()
     */
    friend std::ostream& operator<<(std::ostream& output, const Observable& o);

protected:
    ThObservable * tho; ///< A pointer of to the object of the ThObservables class.
    std::string name; ///< A name for the observable.
    std::string thname; ///< The name for the oservable as fixed in the ThObservable class.
    std::string label; ///< A label for the observable.
    std::string distr; ///< The name of the distribution of the the observable.
    std::string filename; ///< The name of the file containing the experimental likelihood for the observable.
    std::string histoname; ///< The name of the histogram for the observable.
    double ave; ///< The average value of the observable.
    double errg; ///< The gaussian error of the observable.
    double errf; ///< the flat error of the observable.
    double min; ///< The minimum value of the observable.
    double max; ///< The maximum valus of the observable.
    bool tMCMC; ///< The flag to include or exclude the observable from the MCMC run.
    TH1D * inhisto; ///< 1D Histogram containing the experimental likelihood for the observable
    std::string obsType; ///< Type of the Observable. 0: Observable, 1: HiggsObservable, 2: BinnedObservable, 3: FunctionObservable
    double bin_min; ///< The minimum value of the observable bin.
    double bin_max; ///< The maximum valus of the observable bin.
    int iterationNo; ///< A counter for the interation that helps with the observable caching.
    double thValue; ///< The theory value of the first observable.
};


#endif	/* OBSERVABLE_H */

