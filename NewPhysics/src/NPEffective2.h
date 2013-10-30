/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE2_H
#define	NPEFFECTIVE2_H

#include "NPEffective.h"

/**
 * @class NPEffective2
 * @brief A class for new physics with the effective Lagrangian approach.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NPEffective2 : public NPEffective {
public:
    static const int NNPEffectiveVars = 19;
    static const std::string NPEffectiveVars[NNPEffectiveVars];

    /**
     * @brief NPEffective2 constructor.
     */
    NPEffective2()
    : NPEffective()
    {
    }

    virtual std::string ModelName() const
    {
        return "NPEffective2";
    }

    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


protected:
    virtual void parseParameters(const std::string name, const double& value);

};

#endif	/* NPEFFECTIVE2_H */

