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
 * @details
 */
class ModelParaVsObs : public Observable {
public:

    ModelParaVsObs(const std::string name_i,
                   const std::string ParaName_i, const std::string ParaLabel_i,
                   const double ParaMin_i, const double ParaMax_i,
                   const std::string ObsName_i, const std::string ObsLabel_i,
                   const double ObsMin_i, const double ObsMax_i,
                   ThObservable * Obs_i);

    std::string getParaName() const
    {
        return ParaName;
    }

    void setParaName(std::string ParaName)
    {
        this->ParaName = ParaName;
    }

    std::string getParaLabel() const
    {
        return ParaLabel;
    }

    void setParaLabel(std::string ParaLabel)
    {
        this->ParaLabel = ParaLabel;
    }

    double getParaMin() const
    {
        return ParaMin;
    }

    void setParaMin(double ParaMin)
    {
        this->ParaMin = ParaMin;
    }

    double getParaMax() const
    {
        return ParaMax;
    }

    void setParaMax(double ParaMax)
    {
        this->ParaMax = ParaMax;
    }

private:
    std::string ParaName, ParaLabel;
    double ParaMin, ParaMax;    
};

#endif	/* MODELPARAVSOBS_H */

