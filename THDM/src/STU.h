/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STU_H
#define	STU_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class STU
 * @ingroup THDM 
 * @brief An observable class for the electroweak Peskin-Takeuchi pseudo-observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observables S, T and U.
 */
class STU : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   STU(const StandardModel& SM_i, int obsFlag);
     
    /**
     * @brief S, T and U
     * @return
     */
    double computeThValue();

    
    ///////////////////////////////////////////////////////////////////////////
    /* Functions for EW precision observables */

//    double obliqueS() const;
//    double obliqueT() const;
//    double obliqueU() const;
    
    
    ///////////////////////////////////////////////////////////////////////////
protected: 
    
    THDMcache * mycache;

    private:
        const THDM * myTHDM;
        int obs;
    double DeltaS, DeltaT, DeltaU;
    double Mz2;
    double s_W2;//\sin^2(\theta_W)   
    double Mw_i, Mw2;
    double cos2_ba, sin2_ba;   
    
    bool requireCKM, requireYe, requireYn;
          
    ////////////////////////////////////////////////////////////////////////////
    /*One-loop functions*/
    
    /**
     * @brief function F(m0,m1) used for THDM. Remember that this function is
     * defined for THDM while for SUSY we have a multiplicative factor 2.
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the function F for THDM 
     */
    double F(const double m0, const double m1) const;

};

#endif	/* STU_H */
