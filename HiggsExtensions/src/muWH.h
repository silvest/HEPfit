/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MUWH_H
#define	MUWH_H
#include <ThObservable.h>
#include <HiggsExtensionModel.h>

/**
 * @class muWH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{WH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{WH}@f] between the W Higgs associated production cross-section
 * in the current model and in the Standard Model
 */

class muWH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    muWH(const StandardModel& HESM_i) : ThObservable(HESM_i), HESM(static_cast<const HiggsExtensionModel&> (HESM_i))
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")!=0)
            throw std::runtime_error("ERROR: the muWH constructor can only be used with a HiggsExtensionModel reference, while I got " +
                    HESM_i.ModelName() + " as argument");
    };
    muWH(const muWH& orig);
    virtual ~muWH();
    
    /**
     * method to compute the value of  @f[\mu_{WH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKW()*HESM.computeKW());
    }
private:
    const HiggsExtensionModel& HESM;
};

#endif	/* MUWH_H */
