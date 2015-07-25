/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2TEST_H
#define	LEP2TEST_H

#include <iostream>
#include <stdexcept>

/**
 * @class LEP2test
 * @ingroup EW 
 * @brief A test class for the LEP-II observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2test {
public:
    
    LEP2test() 
    {};
    
    int checkSqrtS(const double sqrt_s_i) const 
    {
        for (int i=0; i<12; i++)
            if (sqrt_s_i == sqrt_s[i]) 
                return i; 
        throw std::runtime_error("Error in LEP2test::checkSqrtS()");   
    }
    
    int checkSqrtS_HF(const double sqrt_s_i) const 
    {
        for (int i=0; i<10; i++)
            if (sqrt_s_i == sqrt_s_HF[i]) 
                return i; 
        throw std::runtime_error("Error in LEP2test::checkSqrtS()");   
    }
    
    double sigmaHadronTEST(const double sqrt_s_i) const 
    {
        return sigmaHadron[checkSqrtS(sqrt_s_i)];
    } 

    double sigmaMuTEST(const double sqrt_s_i) const 
    {
        return sigmaMu[checkSqrtS(sqrt_s_i)];
    } 

    double sigmaTauTEST(const double sqrt_s_i) const
    {
        return sigmaTau[checkSqrtS(sqrt_s_i)];
    } 

    double AFBmuTEST(const double sqrt_s_i) const 
    {
        return AFBmu[checkSqrtS(sqrt_s_i)];
    } 

    double AFBtauTEST(const double sqrt_s_i) const 
    {
        return AFBtau[checkSqrtS(sqrt_s_i)];
    } 

    double AFBbottomTEST(const double sqrt_s_i) const
    {
        return AFBbottom[checkSqrtS_HF(sqrt_s_i)];
    } 
    
    double AFBcharmTEST(const double sqrt_s_i) const
    {
        return AFBcharm[checkSqrtS_HF(sqrt_s_i)];
    } 

    double RbottomTEST(const double sqrt_s_i) const 
    {
        return Rbottom[checkSqrtS_HF(sqrt_s_i)];
    } 
    
    double RcharmTEST(const double sqrt_s_i) const 
    {
        return Rcharm[checkSqrtS_HF(sqrt_s_i)];
    } 

    
private:    
    static const double sqrt_s[12];
    static const double sqrt_s_HF[10];
    
    static const double sigmaHadron[12];
    static const double sigmaMu[12];
    static const double sigmaTau[12];
    static const double AFBmu[12];
    static const double AFBtau[12];
    static const double AFBbottom[10];
    static const double AFBcharm[10];
    static const double Rbottom[10];
    static const double Rcharm[10];
    
};

#endif	/* LEP2TEST_H */

