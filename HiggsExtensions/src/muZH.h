/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MUZH_H
#define	MUZH_H
#include <ThObservable.h>

/**
 * @class muZH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{ZH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ZH}@f] between the Z Higgs associated production cross-section
 * in the current model and in the Standard Model
 */

class muZH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExensionModel object or to any extension of it
     */
    muZH(const HiggsExtensionModel& HESM_i)
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")==0) 
        HESM=static_cast<HiggsExtensionModel&>(HESM_i);
        else 
            throw std::runtime_error("ERROR: the muZH constructor can only be used with a HiggsExtensionModel reference, "
                   + "while I got " + HESM_i.ModelName() + " as argument");
    };
    muZH(const muZH& orig);
    virtual ~muZH();
    
    /**
     * method to compute the value of  @f[\mu_{ZH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKZ()*HESM.computeKZ());
    }
private:
    HiggsExtensionModel& HESM;
};

#endif	/* MUZH_H */
