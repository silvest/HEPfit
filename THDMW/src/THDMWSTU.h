/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef THDMWSTU_H
#define THDMWSTU_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDMW.h"
#include "THDMWcache.h"
#include "gslpp.h"

/**
 * @class THDMWSTU
 * @ingroup THDMW 
 * @brief An observable class for the electroweak Peskin-Takeuchi pseudo-observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observables @f$S@f$, @f$T@f$ and @f$U@f$. Formulae from the appendix A of arxiv:0907.2696.
 */
class THDMWSTU : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   THDMWSTU(const StandardModel& SM_i);
     
    /**
     * @brief Empty constructor.
     */
    double computeThValue();

    const THDMW * myTHDMW;

    /**
     * @brief Function @f$F(m02,m12)@f$ used for %THDMW. Remember that this function is
     * defined for %THDMW and SUSY while for %THDM we have a multiplicative factor 0.5.
     * @param[in] m02 mass square @f$m_0^2@f$
     * @param[in] m12 mass square @f$m_1^2@f$
     * @return @f$F(m02,m12)@f$
     */
    double F(const double m02, const double m12) const;


protected:
    THDMWcache * mycache;
          
};






class THDMWDeltaS: public THDMWSTU {
public:

    /**
     * @brief Constructor for DeltaS.
     */
    THDMWDeltaS(const StandardModel& SM_i);

    /**
     * @brief THDMW contribution to @f$S@f$.
     * @return @f$\Delta S@f$
     */
    double computeThValue ();
};


class THDMWDeltaT: public THDMWSTU {
public:

    /**
     * @brief Constructor for DeltaT.
     */
    THDMWDeltaT(const StandardModel& SM_i);

    /**
     * @brief THDMW contribution to @f$T@f$.
     * @return @f$\Delta T@f$
     */
    double computeThValue ();
};

class THDMWDeltaU: public THDMWSTU {
public:

    /**
     * @brief Constructor for DeltaU.
     */
    THDMWDeltaU(const StandardModel& SM_i);

    /**
     * @brief THDMW contribution to @f$U@f$.
     * @return @f$\Delta U@f$
     */
    double computeThValue ();
};



#endif /* THDMWSTU_H */
