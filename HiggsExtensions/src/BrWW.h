/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BrWW_H
#define	BrWW_H
#include <ThObservable.h>

/**
 * @class BrWW
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to WW@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to WW@f]
 * in the current model and in the Standard Model
 */

class BrWW : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrWW(const StandardModel& HESM_i)
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")==0) 
        HESM=static_cast<HiggsExtensionModel&>(HESM_i);
        else 
            throw std::runtime_error("ERROR: the BrWW constructor can only be used with a HiggsExtensionModel reference, "
                   + "while I got " + HESM_i.ModelName() + " as argument");
    };
    
    BrWW(const BrWW& orig);
    virtual ~BrWW();
    
    /**
     * method to compute the the ratio of the @f[BR(H\to WW@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKW()*HESM.computeKW()/HESM.computeGTotalRatio();
    }
    
private:
    HiggsExtensionModel& HESM;
};

#endif	/* BrWW_H */
