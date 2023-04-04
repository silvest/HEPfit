/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GeneralTHDMSTU_H
#define	GeneralTHDMSTU_H

#include "ThObservable.h"

class GeneralTHDM;
class GeneralTHDMcache;


/**
 * @class GeneralTHDMSTU
 * @ingroup GeneralTHDM 
 * @brief An observable class for the electroweak Peskin-Takeuchi pseudo-observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observables @f$S@f$, @f$T@f$ and @f$U@f$. Formulae from equations (21), (22) and (23) in @cite Baak:2011ze.
 */
class GeneralTHDMSTU : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   GeneralTHDMSTU(const StandardModel& SM_i);
     
    /**
     * @brief Empty constructor.
     */
   double computeThValue();

    const GeneralTHDM * myGTHDM;



protected:
// We do not need the cache anymore
//    GeneralTHDMcache * mycache;
          
};

/**
 * @class GTHDMDeltaS
 * @ingroup GeneralTHDM
 * @brief An observable class for the %THDM contribution to the electroweak Peskin-Takeuchi pseudo-observable @f$S@f$.
 */
class GTHDMDeltaS: public GeneralTHDMSTU {
public:

    /**
     * @brief Constructor for DeltaS.
     */
    GTHDMDeltaS(const StandardModel& SM_i);

    /**
     * @brief GeneralTHDM contribution to @f$S@f$.
     * @return @f$\Delta S@f$
     */
    double computeThValue ();
};

/**
 * @class DeltaT
 * @ingroup GeneralTHDM
 * @brief An observable class for the %THDM contribution to the electroweak Peskin-Takeuchi pseudo-observable @f$T@f$.
 */
class GTHDMDeltaT: public GeneralTHDMSTU {
public:

    /**
     * @brief Constructor for DeltaT.
     */
    GTHDMDeltaT(const StandardModel& SM_i);

    /**
     * @brief GeneralTHDM contribution to @f$T@f$.
     * @return @f$\Delta T@f$
     */
    double computeThValue ();
};
/**
 * @class DeltaU
 * @ingroup GeneralTHDM
 * @brief An observable class for the GeneralTHDM contribution to the electroweak Peskin-Takeuchi pseudo-observable @f$U@f$.
 */
class GTHDMDeltaU: public GeneralTHDMSTU {
public:

    /**
     * @brief THDM contribution to @f$U@f$.
     * @brief Constructor for DeltaU.
     */
    GTHDMDeltaU(const StandardModel& SM_i);

    /**
     * @return @f$\Delta U@f$
     */
    double computeThValue ();

protected:
    GTHDMDeltaS * myDeltaS;
};



    #endif	/* GeneralTHDMSTU_H */
