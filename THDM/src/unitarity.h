/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef UNITARITY_H
#define	UNITARITY_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class unitarity
 * @ingroup THDM 
 * @brief An observable class for the requirement of tree level perturbative 
 * unitarity.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require the unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 */
class unitarity : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   unitarity(const StandardModel& SM_i);
     
    /**
     * @brief The unitarity conditions for all the tree level scalar-scalar 
     * scattering amplitudes.
     * @return
     */
    double computeThValue();

    const THDM * myTHDM;
    
};

class unitarity1: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity1(const StandardModel& SM_i);

    /**
     * @return unitarity1
     */
    double computeThValue();
};

class unitarity2: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity2(const StandardModel& SM_i);

    /**
     * @return unitarity2
     */
    double computeThValue();
};

class unitarity3: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity3(const StandardModel& SM_i);

    /**
     * @return unitarity3
     */
    double computeThValue();
};

class unitarity4: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity4(const StandardModel& SM_i);

    /**
     * @return unitarity4
     */
    double computeThValue();
};

class unitarity5: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity5(const StandardModel& SM_i);

    /**
     * @return unitarity5
     */
    double computeThValue();
};

class unitarity6: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity6(const StandardModel& SM_i);

    /**
     * @return unitarity6
     */
    double computeThValue();
};

class unitarity7: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity7(const StandardModel& SM_i);

    /**
     * @return unitarity7
     */
    double computeThValue();
};

class unitarity8: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity8(const StandardModel& SM_i);

    /**
     * @return unitarity8
     */
    double computeThValue();
};

class unitarity9: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity9(const StandardModel& SM_i);

    /**
     * @return unitarity9
     */
    double computeThValue();
};

class unitarity10: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity10(const StandardModel& SM_i);

    /**
     * @return unitarity10
     */
    double computeThValue();
};

class unitarity11: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity11(const StandardModel& SM_i);

    /**
     * @return unitarity11
     */
    double computeThValue();
};

class unitarity12: public unitarity {
public:

    /**
     * @brief Constructor.
     */
    unitarity12(const StandardModel& SM_i);

    /**
     * @return unitarity12
     */
    double computeThValue();
};

#endif	/* UNITARITY_H */

