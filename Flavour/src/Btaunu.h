/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BTAUNU_H
#define	BTAUNU_H

#include "ThObservable.h"
#include "QCD.h"

class StandardModel;

/**
 * @class Btaunu
 * @ingroup Flavour
 * @brief A class for @f$S_{J/\psi\phi}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * branching ratio of \f$B\to \tau\nu\f.
 */
class Btaunu : public ThObservable {
public:   
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Btaunu(const StandardModel& SM_i, QCD::meson meson_i);
    
    /**
     * 
     * @brief arXiv:1206.2634v2
     * @return theoretical value of |\f$ BR(B \rightarrow \tau \nu) \f$|
     */
    double computeThValue();
    
protected:
    
private:
    
    QCD::meson meson;
    
};

#endif	/* BTAUNU_H */