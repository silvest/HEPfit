/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE1_H
#define	NPEFFECTIVE1_H

#include "NPEffective.h"

/**
 * @class NPEffective1
 * @brief A base class for new physics in the form of contributions to the 
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details Base class for an effective Lagrangian extension in the basis of 
 * \cite Barbieri:1999tm, assuming diagonal and family-universal \f${\cal O}_{HQ,HL,HE,HU,HD}\f$
 * and \f${\cal O}_{HQ,HL}^\prime\f$ interactions.
 */
class NPEffective1 : public NPEffective {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NNPEffectiveVars = 11;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string NPEffectiveVars[NNPEffectiveVars];

   /**
     * @brief Constructor.
     */
    NPEffective1()
    : NPEffective()
    {
    }

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPEffective1";
    }

    /**
     * @brief A method to check if all the mandatory parameters for the model have been
     * provided in the model configuration file.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


protected:
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
      
};

#endif	/* NPEFFECTIVE1_H */

