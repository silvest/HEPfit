/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef Brtautau_H
#define	Brtautau_H
#include <ThObservable.h>
#include <HiggsExtensionModel.h>

/**
 * @class Brtautau
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to\tau\tau@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to\tau\tau@f]
 * in the current model and in the Standard Model
 */

class Brtautau : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    Brtautau(const StandardModel& HESM_i) : ThObservable(HESM_i), HESM(static_cast<const HiggsExtensionModel&> (HESM_i))
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")!=0)
            throw std::runtime_error("ERROR: the Brtautau constructor can only be used with a HiggsExtensionModel reference, while I got " +
                    HESM_i.ModelName() + " as argument");
    };
    Brtautau(const Brtautau& orig);
    virtual ~Brtautau();
    
    /**
     * method to compute the the ratio of the @f[BR(H\to\tau\tau@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKtau()*HESM.computeKtau()/HESM.computeGTotalRatio();
    }
private:
    const HiggsExtensionModel& HESM;
};

#endif	/* Brtautau_H */
