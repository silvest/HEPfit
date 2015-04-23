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
 * @class Branching Fraction
 * @ingroup flavour
 * @brief A class for the BR. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BR_MPll : public ThObservable{
public:
    
    /**
    * @brief \f$ BR \f$ 
    */
    BR_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the BR
    */
    double computeBR_MPll(double qmin, double qmax, StandardModel::lepton lep_i);
    double computeThValue ();
    
private:
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
};


/**
 * @class R_K
 * @ingroup flavour
 * @brief A class for the Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class R_MPll : public BR_MPll{
public:
    
    /**
    * @brief \f$ BR_lep1/BR_lep2 \f$ 
    */
    R_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @return the ratio between branching fractions of \f$ B\to K \l_1^+ \l_1^- \f$ and \f$ B\to K \l_2^+ \l_2^- \f$ 
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep1;
    StandardModel::lepton lep2;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
};


/**
 * @class ACP
 * @ingroup flavour
 * @brief A class for the ACP. 
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
    * @return return the ACP
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
};


#endif	/* MPLLOBSERVABLES_H */

    