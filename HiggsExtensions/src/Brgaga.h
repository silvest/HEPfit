/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef Brgaga_H
#define	Brgaga_H
#include <ThObservable.h>

/**
 * @class Brgaga
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to\gamma\gamma@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to\gamma\gamma@f]
 * in the current model and in the Standard Model
 */

class Brgaga : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    Brgaga(const HiggsExtensionModel& HESM_i)
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")==0) 
        HESM=static_cast<HiggsExtensionModel&>(HESM_i);
        else 
            throw std::runtime_error("ERROR: the Brgaga constructor can only be used with a HiggsExtensionModel reference, "
                   + "while I got " + HESM_i.ModelName() + " as argument");
    };
    Brgaga(const Brgaga& orig);
    virtual ~Brgaga();
    
    /**
     * method to compute the the ratio of the @f[BR(H\to\gamma\gamma@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKgaga()*HESM.computeKgaga()/HESM.computeGTotalRatio();
    }
private:
    HiggsExtensionModel& HESM;
};

#endif	/* Brgaga_H */
