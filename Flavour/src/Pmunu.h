/* 
 * Copyright (C) 2024 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#ifndef PMUNU_H
#define PMUNU_H



#include "ThObservable.h"
#include "QCD.h"


class StandardModel;

/**
 * @class Pmunu
 * @ingroup Flavour
 * @brief A class for @f$ \pi \to \lepton \nu @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * the width of @f$ \pi \to \lepton \nu @f$.
 */
class Pmunu : public ThObservable {
public:   
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Pmunu(const StandardModel& SM_i);
    
    /**
     * 
     * @brief arXiv:1605.07114
     * @return theoretical value of |\f$ \Gamma(\pi \rightarrow \lepton \nu) \f$|
     */
    double computeThValue();
    
protected:
    
    
private:
    

    
};



#endif /* KMUNU_O_PMUNU_H */

