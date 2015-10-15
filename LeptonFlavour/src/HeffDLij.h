/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDLIJ_H
#define	HEFFDLIJ_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>

using namespace gslpp;

/**
 * @addtogroup LeptonFlavour
 * @brief A module for lepton flavour observables.
 * @{
 */

class HeffDLij {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDLij(const StandardModel & SM_i);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDLij();
    
    /**
     * 
     * @param 
     * @param 
     * @return
     */
    vector<complex>** ComputeCoeffDLij(int li_lj);
//    vector<complex>** ComputeCoeffDL1_2();
//    vector<complex>** ComputeCoeffDL1_3();

//    WilsonCoefficient getCoeffDL1_1() const {
//        return coeffDL1_1;
//    }
//
//    WilsonCoefficient getCoeffDL1_2() const {
//        return coeffDL1_2;
//    }
//
//    WilsonCoefficient getCoeffDL1_3() const {
//        return coeffDL1_3;
//    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;

    WilsonCoefficient coeffDLij_1;
    WilsonCoefficient coeffDLij_2;
    WilsonCoefficient coeffDLij_3;

    //gslpp::vector<complex> nlep, nlep2, nlepCC;
};

/**
 * @}
 */

#endif	/* HEFFDLIJ_H */

