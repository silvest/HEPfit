/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWPO_H
#define	EWPO_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class EWPO
 * @ingroup THDM 
 * @brief An observable class to calculate the electroweak precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details An observable class to calculate ??? in context of the 2HDM.
 */
class EWPO : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    EWPO(const StandardModel& SM_i);

    /**
     * @brief EWPO.
     * @return 
     */
    double computeThValue();
    double dDelta_r();
    void computeTHDMcouplings();

    const THDM * myTHDM;
    THDMcache * mycache;

    private:
};

class  AlTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    AlTHDM(const StandardModel& SM_i);

    /**
     * @return AlTHDM
     */
    double computeThValue ();
private:
};

class  PpoltauTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    PpoltauTHDM(const StandardModel& SM_i);

    /**
     * @return PpoltauTHDM
     */
    double computeThValue ();
private:
};

class  AcTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    AcTHDM(const StandardModel& SM_i);

    /**
     * @return AcTHDM
     */
    double computeThValue ();
private:
};

class  AbTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    AbTHDM(const StandardModel& SM_i);

    /**
     * @return AbTHDM
     */
    double computeThValue ();
private:
};

class  AFBl0THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    AFBl0THDM(const StandardModel& SM_i);

    /**
     * @return AFBl0THDM
     */
    double computeThValue ();
private:
};

class  AFBc0THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    AFBc0THDM(const StandardModel& SM_i);

    /**
     * @return AFBc0THDM
     */
    double computeThValue ();
private:
};

class  AFBb0THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    AFBb0THDM(const StandardModel& SM_i);

    /**
     * @return AFBb0THDM
     */
    double computeThValue ();
private:
};

class  GammaZTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    GammaZTHDM(const StandardModel& SM_i);

    /**
     * @return GammaZTHDM
     */
    double computeThValue ();
private:
};

class  Rl0THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    Rl0THDM(const StandardModel& SM_i);

    /**
     * @return Rl0THDM
     */
    double computeThValue ();
private:
};

class  Rc0THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    Rc0THDM(const StandardModel& SM_i);

    /**
     * @return Rc0THDM
     */
    double computeThValue ();
private:
};

class  Rb0THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    Rb0THDM(const StandardModel& SM_i);

    /**
     * @return Rb0THDM
     */
    double computeThValue ();
private:
};

class  SigmahadTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    SigmahadTHDM(const StandardModel& SM_i);

    /**
     * @return SigmahadTHDM
     */
    double computeThValue ();
private:
};

class  GammaWTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    GammaWTHDM(const StandardModel& SM_i);

    /**
     * @return GammaWTHDM
     */
    double computeThValue ();
private:
};

class  sinthetaeffl_2THDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    sinthetaeffl_2THDM(const StandardModel& SM_i);

    /**
     * @return sinthetaeffl_2THDM
     */
    double computeThValue ();
private:
};

class  MWTHDM: public EWPO {
public:

    /**
     * @brief Constructor.
     */
    MWTHDM(const StandardModel& SM_i);

    /**
     * @return MWTHDM
     */
    double computeThValue ();
private:
};

#endif	/* EWPO_H */
