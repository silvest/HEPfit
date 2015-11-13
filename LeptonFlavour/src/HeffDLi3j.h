/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDLI3J_H
#define	HEFFDLI3J_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>

using namespace gslpp;

/**
 * @addtogroup LeptonFlavour
 * @brief A module for lepton flavour observables.
 * @{
 */

class HeffDLi3j {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDLi3j(const StandardModel & SM_i);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDLi3j();
    
    /**
     * 
     * @param 
     * @param 
     * @return
     */
    vector<complex>** ComputeCoeffDLi3j(int li_lj);

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;

    WilsonCoefficient coeffDLi3j_1;
    WilsonCoefficient coeffDLi3j_2;
    WilsonCoefficient coeffDLi3j_3;
    WilsonCoefficient coeffDLi3j_4;

};

/**
 * @}
 */

#endif	/* HEFFDLI3J_H */

