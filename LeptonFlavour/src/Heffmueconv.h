/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFMUECONV_H
#define	HEFFMUECONV_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>

using namespace gslpp;

class Heffmueconv {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    Heffmueconv(const StandardModel & SM_i);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~Heffmueconv();
    
    /**
     * 
     * @param 
     * @param 
     * @return
     */
    vector<complex>** ComputeCoeffmueconv();

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;

    WilsonCoefficient coeffmueconv;

};

#endif	/* HEFFMUECONV_H */
