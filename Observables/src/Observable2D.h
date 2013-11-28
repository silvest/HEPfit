/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OBSERVABLE2D_H
#define	OBSERVABLE2D_H

#include "Observable.h"

/**
 * @class Observable2D
 * @ingroup Observable
 * @brief A class for analysing observables pairwise
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The class for building a pair of observables and storing their different
 * parameters read from the SomeModel.conf file. The names (thname) of the observables have
 * to correspond to the allowed names of observables listed in the ThFactory() class.
 */
class Observable2D : public Observable {
public:
    
    /**
     * @brief The default constructor.
     * @param[in] name_i a given name for the observable pair
     * @param[in] thname_i the thname for the first observable fixed in ThFactory()
     * @param[in] thname2_i the thname for the second observable fixed in ThFactory()
     * @param[in] label_i the label assigned to the first observable
     * @param[in] label2_i the label assigned to the second observable
     * @param[in] tMCMC_i boolean flag to indicate inclusion in MCMC
     * @param[in] min_i minimum value for the first observable
     * @param[in] max_i maximum value for the first observable
     * @param[in] min2_i minimum value for the second observable
     * @param[in] max2_i maximum value for the second observable
     * @param[in] tho_i a pointer to an object of type ThObservable() for the first observable
     * @param[in] tho2_i a pointer to an object of type ThObservable() for the second observable
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

    /**
     * @brief A method to access the computed theory value of the second observable
     */
    double computeTheoryValue2();

    /**
     * @brief A get method to access the label for the second observable
     * @return the label for the second observable
     */
    std::string getLabel2() const
    {
        return label2;
    }

    /**
     * @brief A set method to fix the label for the second observable
     * @param[in] label the label for the second observable
     */
    void setLabel2(std::string label2)
    {
        this->label2 = label2;
    }

    /**
     * @brief A get method to access the maximum value of the second observable
     * @return the maximum value of the second observable
     */
    double getMax2() const
    {
        return max2;
    }

    /**
     * @brief A set method to fix the maximum value for the second observable
     * @param[in] the maximum value for the second observable
     */
    void setMax2(double max2)
    {
        this->max2 = max2;
    }

    /**
     * @brief A get method to access the minimum value of the second observable
     * @return the minimum value of the second observable
     */
    double getMin2() const
    {
        return min2;
    }

    /**
     * @brief A set method to fix the minimum value for the second observable
     * @param[in] the minimum value for the second observable
     */
    void setMin2(double min2)
    {
        this->min2 = min2;
    }

    /**
     * @brief A get method to access the thname of the second observable as defined in ThFactory() class
     * @return thname the name of the second observable as listed in ThFactory() class
     */
    std::string getThname2() const
    {
        return thname2;
    }

    /**
     * @brief A set method to fix the name of the second observable as listed in ThFactory() class
     * @param[in] thname the name of the second observable as listed in ThFactory() class
     */
    void setThname2(std::string thname2)
    {
        this->thname2 = thname2;
    }

    /**
     * @brief A get method to access the pointer to the object of the ThObservable() class for 
     * the second observable.
     * @return pointer to the object of type ThObservable() for the second observable
     */
    ThObservable* getTho2() const
    {
        return tho2;
    }

    /**
     * @brief A set method to fix the pointer to object of type ThObservable()class for
     * the second observable.
     * @param[in] pointer to the object of type ThObservable() for the second observable
     */
    void setTho2(ThObservable* tho2)
    {
        this->tho2 = tho2;
    }

private:
    std::string thname2; /**< The name for the second oservable as fixed in the ThObservable() class.*/
    std::string label2; /**< A label for the second observable. */
    double min2; /**< The minimum value of the second observable. */
    double max2; /**< The maximum valus of the second observable. */
    ThObservable * tho2;
};

#endif	/* OBSERVABLE2D_H */

