/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MUttH_H
#define	MUttH_H
#include <ThObservable.h>

/**
 * @class muttH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{ttH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ttH}@f] between the t-tbar-Higgs associated production cross-section
 * in the current model and in the Standard Model
 */

class muttH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    muttH(const HiggsExtensionModel& HESM_i)
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")==0) 
        HESM=static_cast<HiggsExtensionModel&>(HESM_i);
        else 
            throw std::runtime_error("ERROR: the muttH constructor can only be used with a HiggsExtensionModel reference, "
                   + "while I got " + HESM_i.ModelName() + " as argument");
    };
    muttH(const muttH& orig);
    virtual ~muttH();
    
    /**
     * method to compute the value of  @f[\mu_{ttH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKt()*HESM.computeKt());
    }
private:
    HiggsExtensionModel& HESM;
};

#endif	/* MUttH_H */
