/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBSTRANGE_H
#define	AFBSTRANGE_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class AFBstrange
 * @ingroup EW 
 * @brief An observable class for the forward-backward asymmetry in
 * @f$e^+ e^-\to Z\to s\bar{s}@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry for
 * the @f$e^+ e^-\to Z\to s\bar{s}@f$ process at the @f$Z@f$ pole in terms of
 * the asymmetry parameters @f$\mathcal{A}_f@f$, computed in the EW class:
 * @f[
 * A_{\mathrm{FB}}^{0,c}=\frac 34 \mathcal{A}_e \mathcal{A}_s\,.
 * @f]
 * 
 * @sa EW_NPZff::AFBstrange() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class AFBstrange : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AFBstrange(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };
    
    /**
     * @brief The forward-backward asymmetry for @f$e^+ e^-\to Z \to s\bar{s}@f$
     * at the @f$Z@f$ pole, @f$A^{0,s}_{\mathrm{FB}}@f$
     * @return @f$A^{0,s}_{\mathrm{FB}}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* AFBSTRANGE_H */

