/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLLOBSERVABLES_H
#define	MPLLOBSERVABLES_H

#include "MPll.h"
#include <StandardModel.h>
#include <ThObservable.h>




/**
 * @class Branching Fraction for electron channel
 * @ingroup flavour
 * @brief A class for the clean observable BR_e. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BR_MPll : public ThObservable{
public:
    
    /**
    * @brief \f$ BR_e \f$ 
    */
    BR_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable BR_e
    */
    double computeBR_MPll(double qmin, double qmax, StandardModel::lepton lep_i);
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
};


/**
 * @class ratio between BR for electron and muon channels
 * @ingroup flavour
 * @brief A class for the Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class R_MPll : public BR_MPll{
public:
    
    /**
    * @brief \f$ BR_e/BR_mu \f$ 
    */
    R_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @return the ratio between branching fractions of \f$ B\to K e^+ e^- \f$ and \f$ B\to K \mu^+ \mu^- \f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep1;
    StandardModel::lepton lep2;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
};


/**
 * @class ACP for electron channel
 * @ingroup flavour
 * @brief A class for the clean observable ACP_e. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ACP_MPll : public BR_MPll{
public:
    
    /**
    * @brief \f$ A_{CP} \f$ 
    */
    ACP_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable ACP_e
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
};


#endif	/* MPLLOBSERVABLES_H */

    