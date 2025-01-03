/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2THOBSERVABLE_H
#define	LEP2THOBSERVABLE_H

#include <stdexcept>
#include <cstring>
#include <ThObservable.h>
//#include <NPSTUVWXY.h>
#include <StandardModel.h>
#include "LEP2test.h"

// TEST: use the Zfitter outputs defined in LEP2test class for SM predictions
//#define LEP2TEST

/**
 * @class LEP2ThObservable
 * @ingroup EW 
 * @brief A class for the LEP2 inclusive observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class LEP2ThObservable : public ThObservable  {
public:

    // Radiative Corrections for the LEP2 observables
//    enum LEP2RCs {Weak=0, WeakBox, ISR, QEDFSR, QCDFSR, NUMofLEP2RCs};
        
    /**
     * @brief LEP2ThObservable constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     * @param[in] bSigmaForR_i true for the denominator of R_b or R_c
     */
    LEP2ThObservable(const StandardModel& SM_i, const double sqrt_s_i) 
            : ThObservable(SM_i),sqrt_s(sqrt_s_i), s(sqrt_s_i*sqrt_s_i)
    {
    }

    const LEP2test myTEST;
    const double sqrt_s, s;
    
    
private:   
      
};


/**
 * @class LEP2ThDiffObservable
 * @ingroup EW 
 * @brief A class for the LEP2 differential observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class LEP2ThDiffObservable : public ThObservable  {
public:
        
    /**
     * @brief LEP2ThDiffObservable constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] cos_i the polar angle of the final state particle wrt e^-
     */
    LEP2ThDiffObservable(const StandardModel& SM_i, const double sqrt_s_i, const double cos_i) 
            : ThObservable(SM_i),sqrt_s(sqrt_s_i), s(sqrt_s_i*sqrt_s_i), cos(cos_i)
    {
    }

    const LEP2test myTEST;
    const double sqrt_s, s, cos;
    
    
private:   
      
};


/**
 * @class eeffThObservable
 * @ingroup EW 
 * @brief A class for the two to two fermion inclusive observables in electron positron colliders. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class eeffThObservable : public ThObservable  {
public:
        
    /**
     * @brief eeffThObservable constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffThObservable(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
            : ThObservable(SM_i), pol_e(pol_e_i), pol_p(pol_p_i), sqrt_s(sqrt_s_i), s(sqrt_s_i*sqrt_s_i)
    {
    }

    const double pol_e, pol_p, sqrt_s, s;
     
private:   
      
};


#endif	/* LEP2THOBSERVABLE_H */

