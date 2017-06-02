/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYEDMCOUPLINGS_H
#define SUSYEDMCOUPLINGS_H

#include "gslpp.h"
#include "EffectiveEDMcouplings.h"

class SUSY;

/**
 * @class SUSYEDMcouplings
 * @ingroup SUSY
 * @brief A class the effective %SUSY couplings for EDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class SUSYEDMcouplings : public EffectiveEDMcouplings {
public:
    /**
     * @brief The default constructor.
     */
    SUSYEDMcouplings(const SUSY & SUSY_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~SUSYEDMcouplings();
    
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
//    const SUSY & mySUSY;

};

#endif /* SUSYEDMCOUPLINGS_H */

