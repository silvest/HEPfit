/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELPARAVSOBS_H
#define	MODELPARAVSOBS_H

#include "Observable.h"

/**
 * @class ModelParaVsObs
 * @ingroup Observable
 * @brief A class for the correlation between a model parameter and an observable. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class builds the model parameter vs. observable comparisons 
 * that are specified as output in the SomeModel.conf file. The name (ParaName) of the 
 * parameter has to correspond to the list of names in the SomeModel.cpp file in the
 * Model src directories.\n
 * e.g. GF, which can be found in StandardModel.cpp file in the
 * StandardModel project source folder %StandardModel/src.
 * The name (ObsName) of the observable has
 * to correspond to the allowed name of observables listed in the ThFactory class.
 * A list of both can also be found on the main page of the documentation website.
 */
class ModelParaVsObs : public Observable {
public:

    /**
     * @brief Constructor.
     * @param[in] name_i a given name for model parameter vs. observable comparison
     * @param[in] ParaName_i the name of the parameter to be compared
     * @param[in] ParaLabel_i the label for the parameter to be compared
     * @param[in] ParaMin_i the minimum value of the parameter
     * @param[in] ParaMax_i the maximum value for the parameter
     * @param[in] ObsName_i the name of the observable
     * @param[in] ObsLabel_i the label for the observable
     * @param[in] ObsMin_i the minimum value for the observable
     * @param[in] ObsMax_i the maxinmum value for the observable
     * @param[in] Obs_i a pointer of type ThObservable for the observable
     */
    ModelParaVsObs(const std::string name_i,
                   const std::string ParaName_i,
                   const std::string ParaLabel_i,
                   const double ParaMin_i,
                   const double ParaMax_i,
                   const std::string ObsName_i,
                   const std::string ObsLabel_i,
                   const double ObsMin_i,
                   const double ObsMax_i,
                   ThObservable * Obs_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ModelParaVsObs();

    /**
     * @brief A get method to access the parameter name.
     * @return the parameter name
     */
    std::string getParaName() const
    {
        return ParaName;
    }

    /**
     * @brief A set method to fix the parameter name.
     * @param[in] ParaName the parameter name
     */
    void setParaName(std::string ParaName)
    {
        this->ParaName = ParaName;
    }

    /**
     * @brief A get method to access the parameter label.
     * @return the parameter label
     */
    std::string getParaLabel() const
    {
        return ParaLabel;
    }

    /**
     * @brief A set method to fix the parameter label.
     * @param[in] ParaLabel the parameter label
     */
    void setParaLabel(std::string ParaLabel)
    {
        this->ParaLabel = ParaLabel;
    }

    /**
     * @brief A get methos to access the minimum value of the parameter.
     * @return the minimum value of the parameter
     */
    double getParaMin() const
    {
        return ParaMin;
    }

    /**
     * @brief A set method to fix the minimum value of the parameter.
     * @param[in] ParaMin the minimum value of the parameter
     */
    void setParaMin(double ParaMin)
    {
        this->ParaMin = ParaMin;
    }

    /**
     * @brief A get methos to access the maximum value of the parameter.
     * @return the maximum value of the parameter
     */
    double getParaMax() const
    {
        return ParaMax;
    }

    /**
     * @brief A set methos to fix the maximum value of the parameter.
     * @param[in] ParaMax the maximum value of the parameter
     */
    void setParaMax(double ParaMax)
    {
        this->ParaMax = ParaMax;
    }

private:
    std::string ParaName; ///< The name of the parameter.
    std::string ParaLabel; ///< The label for the parameter.
    double ParaMin; ///< The minimum value for the parameter.
    double ParaMax; ///< The maximum value for the parameter.
};

#endif	/* MODELPARAVSOBS_H */

