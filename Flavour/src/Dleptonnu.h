/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#ifndef DLEPTONNU_H
#define DLEPTONNU_H



#include "ThObservable.h"
#include "QCD.h"

class StandardModel;

/**
 * @class Dleptonnu
 * @ingroup Flavour
 * @brief A class for @f$D\to \lepton \nu @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * branching ratio of \f$D\to \lepton \nu\f.
 * The SM amplitude includes the leading-log short-distance electroweak (Sirlin) factor
 * \f$1+\frac{\alpha}{\pi}\ln\frac{M_Z}{M_D}\f$ through HeffDF1_diujlknu.
 */
class Dleptonnu : public ThObservable {
public:   
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Dleptonnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::lepton lepton_i);
    
    /**
     * 
     * @brief arXiv:1206.2634v2
     * @return theoretical value of |\f$ BR(D \rightarrow \lepton \nu) \f$|
     */
    double computeThValue();
    
protected:
    
private:
    
    QCD::meson meson;
    QCD::lepton lepton;
    
};


#endif /* DLEPTONNU_H */

