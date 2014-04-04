/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MUVBF_H
#define	MUVBF_H
#include <ThObservable.h>

/**
 * @class muVBF
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{VBF}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{VBF}@f] between the vector-boson fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */

class muVBF : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    muVBF(const HiggsExtensionModel& HESM_i)
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")==0) 
        HESM=static_cast<HiggsExtensionModel&>(HESM_i);
        else 
            throw std::runtime_error("ERROR: the muVBF constructor can only be used with a HiggsExtensionModel reference, "
                   + "while I got " + HESM_i.ModelName() + " as argument");
    };
    muVBF(const muVBF& orig);
    virtual ~muVBF();
    
    /**
     * method to compute the value of  @f[\mu_{VBF}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKW()*HESM.computeKW()*HESM.computeSigmaWF()+HESM.computeKZ()*HESM.computeKZ()*HESM.computeSigmaZF()+
                HESM.computeKW()*HESM.computeKZ()*HESM.computeSigmaZWF())/
                (HESM.computeSigmaWF()+HESM.computeSigmaZF()+HESM.computeSigmaZF());
    }
private:
    HiggsExtensionModel& HESM;
};

#endif	/* MUVBF_H */
