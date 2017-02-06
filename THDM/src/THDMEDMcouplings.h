/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMEDMCOUPLINGS_H
#define THDMEDMCOUPLINGS_H

#include "gslpp.h"
#include "EffectiveEDMcouplings.h"

class THDM;

/**
 * @class THDMEDMcouplings
 * @ingroup THDM
 * @brief A class the effective %THDM couplings for EDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class THDMEDMcouplings : public EffectiveEDMcouplings {
public:
    /**
     * @brief The default constructor.
     */
    THDMEDMcouplings(const THDM & THDM_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~THDMEDMcouplings();
    
    /**
     * 
     * @brief \f$ \kappa_t \f$ 
     * @return the value of \f$ \kappa_t \f$
     */
    virtual double kappa_t() ;
    
    /**
     * 
     * @brief \f$ \kappa_b \f$ 
     * @return the value of \f$ \kappa_b \f$
     */
    virtual double kappa_b() ;
    
private:
    const THDM & myTHDM;

};

#endif /* THDMEDMCOUPLINGS_H */

