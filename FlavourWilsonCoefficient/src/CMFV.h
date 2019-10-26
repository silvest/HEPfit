/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CMFV_H
#define CMFV_H

#include "StandardModel.h"
#include "CMFVMatching.h"

/**
 * @addtogroup FlavourWilsonCoefficent
 * @brief A module for EFT models for flavour observaqbles.
 * @details This module has implementations of several EFTs for flavour
 * analyses like the Weak Effective Theory at low energy and generalized
 * EFT models for @f$\Delta F = 1@f$ and @f$\Delta F = 2@f$.
 * @{
 */

/**
 * @class CMFV
 * @brief Model for CMFV contributions to flavour. 
 */

class CMFV : public StandardModel {
public:
    
    static const int NCMFVvars = 1;

    static const std::string CMFVvars[NCMFVvars];
    
    /**
     * @brief CMFV constructor
     */
    CMFV();
    
    /**
     * @brief A method to check if all the mandatory parameters for %CMFV
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A get method to access the member reference of type %CMFVMatching
     * @return a reference to a %CMFVMatching object
     */
    virtual CMFVMatching& getMatching() const
    {
        return CMFVM.getObj();
    }
    
    /**
     *
     * @return Ftt
     */
    double getFtt() const
    {
        return Ftt;
    }


protected: 
    
    virtual void setParameter(const std::string, const double&);
    mutable Matching<CMFVMatching,CMFV> CMFVM;

private:
    double Ftt;

};

#endif /* CMFV_H */

