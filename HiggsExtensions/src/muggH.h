/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MUGGH_H
#define	MUGGH_H
#include <ThObservable.h>
#include <HiggsExtensionModel.h>

/**
 * @class muggH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{ggH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ggH}@f] between the gluon-gluon fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */

class muggH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    muggH(const StandardModel& HESM_i) : ThObservable(HESM_i), HESM(static_cast<const HiggsExtensionModel&> (HESM_i))
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")!=0)
            throw std::runtime_error("ERROR: the muggH constructor can only be used with a HiggsExtensionModel reference, while I got " +
                    HESM_i.ModelName() + " as argument");

    };
    muggH(const muggH& orig);
    virtual ~muggH();
    
    /**
     * method to compute the value of  @f[\mu_{ggH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return HESM.computeKglgl()*HESM.computeKglgl();
    }
    
private:
    const HiggsExtensionModel& HESM;
};

#endif	/* MUGGH_H */

