/* 
 * Copyright (C) 2015 HEPfit Collaboration
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
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observables S, T and U. Formulae from http://arxiv.org/pdf/1107.0975v2.pdf.
 */
class STU : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   STU(const StandardModel& SM_i);
     
    /**
     * @brief S, T and U
     * @return
     */
    double computeThValue();

    const THDM * myTHDM;

    /**
     * @brief function F(m02,m12) used for THDM. Remember that this function is
     * defined for THDM while for SUSY we have a multiplicative factor 2.
     * @param[in] m02 mass square m_0^2
     * @param[in] m12 mass square m_1^2
     * @return the function F for THDM 
     */
    double F(const double m02, const double m12) const;

protected:
    THDMcache * mycache;

private:
//    double DeltaS, DeltaT, DeltaU;
//    double Mz2;
//    double s_W2;//\sin^2(\theta_W)   
//    double Mw_i, Mw2;
    
//    bool requireCKM, requireYe, requireYn;
          
};

class DeltaS: public STU {
public:

    /**
     * @brief Constructor.
     */
    DeltaS(const StandardModel& SM_i);

    /**
     * @return DeltaS
     */
    double computeThValue ();
};

class DeltaT: public STU {
public:

    /**
     * @brief Constructor.
     */
    DeltaT(const StandardModel& SM_i);

    /**
     * @return DeltaT
     */
    double computeThValue ();
};

class DeltaU: public STU {
public:

    /**
     * @brief Constructor.
     */
    DeltaU(const StandardModel& SM_i);

    /**
     * @return DeltaU
     */
    double computeThValue ();

protected:
    DeltaS * myDeltaS;
};

#endif	/* STU_H */
