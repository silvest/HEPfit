/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMINUS2_H
#define	GMINUS2_H

#include "ThObservable.h"




//class gminus2 : public ThObservable {
//public:
//    /**
//     * constructor
//     * @param LeptonFlavour
//     */
//    gminus2(const StandardModel& SM_i);
//
//    /**
//     *
//     * @brief 
//     * @return
//     */
//    double computeThValue();
//
//protected:
//
//private:
//    
//};


/**
 * @class gminus2_mu
 * @ingroup LeptonFlavour
 * @brief A class for calculating the \f$ (g-2)_{\mu} \f$ at one-loop. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The gminus2_mu class calculates the contribution to \f$ (g-2)_{\mu} \f$ at one-loop generated in the model.
 */
class gminus2_mu : public ThObservable {
public:
    
    /**
     * @brief Constructor of the class gminus2_mu
     */
    gminus2_mu(const StandardModel& SM_i);
    
    /**
     * @brief Calculates the value of \f$ (g-2)_{\mu} \f$ at one-loop.
     * @return value of \f$ (g-2)_{\mu} \f$
     */
    double computeThValue();
    
private:


};

#endif	/* GMINUS2_H */
