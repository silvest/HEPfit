/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EFFECTIVEEDMCOUPLINGS_H
#define	EFFECTIVEEDMCOUPLINGS_H

#include <gslpp.h>

#define LEPS 1.e-5 // tolerance in the limit of S(x,y) to S(x)

class StandardModel;

/**
 * @class EffectiveEDMcouplings
 * @ingroup StandardModel
 * @brief A class for effective EDM couplings. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */

class EffectiveEDMcouplings {
public:
    /**
     * @brief The default constructor.
     */
    EffectiveEDMcouplings(const StandardModel & SM_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~EffectiveEDMcouplings();
    
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
  
protected:
    
private:
    
    const StandardModel & SM;
    
    
};

#endif	/* EFFECTIVEEDMCOUPLINGS_H */

