/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef TAUKNU_O_TAUPNU_H
#define TAUKNU_O_TAUPNU_H



#include "ThObservable.h"
#include "QCD.h"


class StandardModel;

/**
 * @class Kmunu_o_Pmunu
 * @ingroup Flavour
 * @brief A class for @f$ (K\to \lepton \nu) / (\pi \to \lepton \nu) @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * branching ratio of @f$ (K\to \lepton \nu) / (\pi \to \lepton \nu) @f$.
 */
class tauKnu_o_tauPnu : public ThObservable {
public:   
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    tauKnu_o_tauPnu(const StandardModel& SM_i);
    
    /**
     * 
     * @brief arXiv:1206.2634v2
     * @return theoretical value of |\f$ BR(D \rightarrow \lepton \nu) \f$|
     */
    double computeThValue();
    
protected:
    
    
private:
    

    
};


#endif /* TAUKNU_O_TAUPNU_H */

