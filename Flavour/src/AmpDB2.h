/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDB2_H
#define	AMPDB2_H

class StandardModel;
#include "gslpp_complex.h"
#include "OrderScheme.h"
#include "gslpp.h"

/**
 * @addtogroup Flavour
 * @brief A module for all the flavour observables implemented in HEPfit.
 * @details This module has several flavour physics observables which include
 * quark flavour violation in the beauty, charm and strange sectors. This includes 
 * the evolutors, Hamiltonians and low energy observables.
 * @{
 */


/**
 * @class AmpDB2
 * @ingroup Flavour
 * @brief \f$ | \Delta B = 2 | \f$ Amplitude Class
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is related to the calculation of the \f$ B_{d,s}-\bar{B}_{d,s}\f$
 * mixing.
 *
 */

class AmpDB2 {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    AmpDB2(const StandardModel& SM_i);

    /**
    * @brief The value of @f$M_{12}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{12}^{bd}@f$
    */
    gslpp::complex getAmpBd(orders order){
        return AmpBd(order);
    }

    /**
    * @brief The value of @f$M_{12}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{12}^{bs}@f$
    */
    gslpp::complex getAmpBs(orders order){
        return AmpBs(order);
    }

    gslpp::complex getPBd(){
        return PBd();
    }

    gslpp::complex getPBs(){
        return PBs();
    }

protected:
    /**
    * @brief A method to compute @f$M_{12}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{12}^{bd}@f$
    */
    gslpp::complex AmpBd(orders order);
    
    /**
    * @brief A method to compute @f$M_{12}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{12}^{bs}@f$
    */
    gslpp::complex AmpBs(orders order);
    
    /**
    * @brief A method to compute the ratio of the absolute value of the $B_s$ mixing amplitude over the Standard Model value.
    * @param[in] order the %QCD order of the computation
    * @return @f$\vert (M_{12}^{bs})_\mathrm{full}/(M_{12}^{bs})_\mathrm{SM}\vert@f$
    */
    gslpp::complex RBs(orders order);
    gslpp::complex PBd();
    gslpp::complex PBs();

private:

    const StandardModel& mySM;/**< Model type */
    
    gslpp::complex C_1_SM;/**<Wilson coeffients @f$C_1@f$*/

};

/**
 * @}
 */

#endif	/* AMPDB2_H */

