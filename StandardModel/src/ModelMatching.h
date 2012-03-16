/* 
 * File:   ModelMatching.h
 * Author: Mauro_87
 *
 * Created on 13 dicembre 2011, 13.15
 */

#ifndef MODELMATCHING_H
#define	MODELMATCHING_H

#include"WilsonCoefficient.h"
#include"StandardModel.h"

class ModelMatching {
public:
    
    virtual const std::vector<WilsonCoefficient>& CMdbd2() = 0;
    virtual const std::vector<WilsonCoefficient>& CMdbs2() = 0;
    virtual const std::vector<WilsonCoefficient>& CMdd2() = 0;

    virtual const std::vector<WilsonCoefficient>& CMbnlep( const int& a) = 0;
    virtual const std::vector<WilsonCoefficient>& CMbnlepCC( const int& a) = 0;
    virtual const std::vector<WilsonCoefficient>& CMbsg()= 0;
    virtual const std::vector<WilsonCoefficient>& CMd1()= 0;
    virtual const std::vector<WilsonCoefficient>& CMd1Buras()= 0;

    
};

#endif	/* MODELMATCHING_H */

