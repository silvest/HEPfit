/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMEDMCOUPLINGS_H
#define GENERALTHDMEDMCOUPLINGS_H

#include "gslpp.h"
#include "EffectiveEDMcouplings.h"

class GeneralTHDM;

/**
 * @class GeneralTHDMEDMcouplings
 * @ingroup THDM
 * @brief A class the effective %GeneralTHDM couplings for EDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GeneralTHDMEDMcouplings : public EffectiveEDMcouplings {
public:
    /**
     * @brief The default constructor.
     */
    GeneralTHDMEDMcouplings(const GeneralTHDM & GeneralTHDM_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~GeneralTHDMEDMcouplings();
    
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
    const GeneralTHDM & myGeneralTHDM;

};

#endif /* GENERALTHDMEDMCOUPLINGS_H */

