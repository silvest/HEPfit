/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBSERVABLE2D_H
#define	OBSERVABLE2D_H

#include "Observable.h"
#include <TH2D.h>

/**
 * @class Observable2D
 * @ingroup Observable
 * @brief A class for analyzing observables pairwise
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The class for building a pair of observables and storing their different
 * parameters read from the SomeModel.conf file or specified by the user. The names (thname) of the observables have
 * to correspond to the allowed names of observables listed in the ThFactory class.
 */
class Observable2D : public Observable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] name_i a given name for the observable pair
     * @param[in] thname_i the thname for the first observable fixed in ThFactory
     * @param[in] thname2_i the thname for the second observable fixed in ThFactory
     * @param[in] label_i the label assigned to the first observable
     * @param[in] label2_i the label assigned to the second observable
     * @param[in] tMCMC_i boolean flag to indicate inclusion in MCMC
     * @param[in] min_i minimum value for the first observable
     * @param[in] max_i maximum value for the first observable
     * @param[in] min2_i minimum value for the second observable
     * @param[in] max2_i maximum value for the second observable
     * @param[in] tho_i a pointer to an object of type ThObservable for the first observable
     * @param[in] tho2_i a pointer to an object of type ThObservable for the second observable
     */
    Observable2D(const std::string name_i,
                 const std::string thname_i,
                 const std::string thname2_i,
                 const std::string label_i,
                 const std::string label2_i,
                 const bool tMCMC_i,
                 const double min_i,
                 const double max_i,
                 const double min2_i,
                 const double max2_i,
                 ThObservable * tho_i,
                 ThObservable * tho2_i);
    
    /**
     * @brief The default constructor.
     */
    Observable2D();
    
    /**
     * @brief A conversion constructor. Constructs
     * Observable2D with just one observable.
     */
    Observable2D(const Observable& o1d);
    
    /**
     * @brief The copy constructor.
     */
    Observable2D(const Observable2D& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~Observable2D();
    
    int ParseObservable2D(std::string& type, 
                           boost::tokenizer<boost::char_separator<char> >* tok, 
                           boost::tokenizer<boost::char_separator<char> >::iterator& beg,
                           std::string& infilename,
                           std::ifstream& ifile, 
                           int lineNo,
                           int rank);

    /**
     * @brief A set method to set the likelihood from which the experimental likelihood of the observable will
     * be read.
     * @param filename the name of the file
     * @param histoname the name of the histogram
     */
    virtual void setLikelihoodFromHisto(std::string filename, std::string histoname);

    /**
     * @brief A method to access the computed theory value of the second observable.
     */
    double computeTheoryValue2();

    /**
     * @brief A get method to access the label for the second observable.
     * @return the label for the second observable
     */
    std::string getLabel2() const
    {
        return label2;
    }

    /**
     * @brief A set method to fix the label for the second observable.
     * @param[in] label2 the label for the second observable
     */
    void setLabel2(std::string label2)
    {
        this->label2 = label2;
    }

    /**
     * @brief A get method to access the maximum value of the second observable.
     * @return the maximum value of the second observable
     */
    double getMax2() const
    {
        return max2;
    }

    /**
     * @brief A set method to fix the maximum value for the second observable.
     * @param[in] max2 the maximum value for the second observable
     */
    void setMax2(double max2)
    {
        this->max2 = max2;
    }

    /**
     * @brief A get method to access the minimum value of the second observable.
     * @return the minimum value of the second observable
     */
    double getMin2() const
    {
        return min2;
    }

    /**
     * @brief A set method to fix the minimum value for the second observable.
     * @param[in] min2 the minimum value for the second observable
     */
    void setMin2(double min2)
    {
        this->min2 = min2;
    }

    /**
     * @brief A get method to access the thname of the second observable as defined in ThFactory class.
     * @return thname the name of the second observable as listed in ThFactory() class
     */
    std::string getThname2() const
    {
        return thname2;
    }

    /**
     * @brief A set method to fix the name of the second observable as listed in ThFactory class.
     * @param[in] thname2 the name of the second observable as listed in ThFactory class
     */
    void setThname2(std::string thname2)
    {
        this->thname2 = thname2;
    }

    /**
     * @brief A get method to access the pointer to the object of the ThObservable class for
     * the second observable.
     * @return pointer to the object of type ThObservable() for the second observable
     */
    ThObservable* getTho2() const
    {
        return tho2;
    }

    /**
     * @brief A set method to fix the pointer to object of type ThObservable class for
     * the second observable.
     * @param[in] tho2 pointer to the object of type ThObservable() for the second observable
     */
    void setTho2(ThObservable* tho2_i)
    {
        tho2 = tho2_i;
        tho2->setBinMin(bin_min[1]);
        tho2->setBinMax(bin_max[1]);
    }
    
    /**
     * @brief A set method to fix the pointer to object of type ThObservable class for
     * the second observable.
     * @param[in] tho2 pointer to the object of type ThObservable() for the second observable
     */
    void setTho1Tho2(ThObservable* tho1_i, ThObservable* tho2_i)
    {
        setTho(tho1_i, bin_min[0], bin_max[0]);
        setTho2(tho2_i);
    }

    /**
     * @brief A method to compute the weight associated with the observable.
     * @param[in] th1 the theoretical value of the first observable
     * @param[in] th2 the theoretical value of the second observable
     */
    virtual double computeWeight(double th1, double th2);

    /**
     * @brief A method to compute the weight associated with the observable.
     */
    virtual double computeWeight()
    {
        return computeWeight(computeTheoryValue(),computeTheoryValue2());
    }
    
    void setObsType2(std::string& obsType2_i)
    {
        obsType2 = obsType2_i;
    }
    
    std::string getObsType2() const
    {
        return obsType2;
    }
    
    /**
     * @brief A get method to access the average value of the second observable.
     * @return the average value of the observable
     */
    double getAve2() const
    {
        return ave2;
    }

    /**
     * @brief A set method to fix the average value of the second observable.
     * @param[in] ave the average value of the observable
     */
    void setAve2(double ave2)
    {
        this->ave2 = ave2;
    }

    /**
     * @brief A get method to access the flat error of the second observable.
     * @return the flat error of the observable
     */
    double getErrf2() const
    {
        return errf2;
    }

    /**
     * @brief A set method to fix the flat error of the second observable.
     * @param[in] errf the flat error of the observable
     */
    void setErrf2(double errf2)
    {
        this->errf2 = errf2;
    }

    /**
     * @brief A get method to access the Gaussian error of the second observble.
     * @return the Gauissian error of the observable
     */
    double getErrg2() const
    {
        return errg2;
    }

    /**
     * @brief A set method to fix the gaussian error of the observable.
     * @param[in] errg the Gaussian error of the observable
     */
    void setErrg2(double errg2)
    {
        this->errg2 = errg2;
    }
    
    /**
     * @brief A method to check if the end of file has been reached
     * @return a boolean which is true if the end of file has been reached
     */
    bool isEOF()
    {
        return IsEOF;
    }
    
private:
    std::string thname2; ///< The name for the second oservable as fixed in the ThObservable() class.
    std::string label2; ///< A label for the second observable.
    double min2; ///< The minimum value of the second observable.
    double max2; ///< The maximum valus of the second observable.
    double ave2; ///< The average value of the second observable.
    double errg2; ///< The gaussian error of the second observable.
    double errf2; ///< the flat error of the second observable.
    std::string obsType2; ///< Type of the second Observable. 0: Observable, 1: HiggsObservable, 2: BinnedObservable, 3: FunctionObservable.
    ThObservable * tho2; ///< A pointer to an object of the ThObservable class.
    TH2D * inhisto2d;  ///< 2D Histogram containing the experimental likelihood for the observable.
    std::vector<double> bin_min; ///< The minimum value of the bin.
    std::vector<double> bin_max; ///< The maximum value of the bin.
    std::string filepath; ///< The path to the file being parsed.
    int iterationNo2; ///< Counts the iteration to help with caching.
    double thValue2; ///< The theory value of the second observable.
    bool IsEOF; ///< A bolean that is true if the end of file is reached.
};

#endif	/* OBSERVABLE2D_H */

