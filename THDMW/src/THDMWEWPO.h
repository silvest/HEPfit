/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   THDMWEWPO.h
 * Author: victormirallesaznar
 *
 * Created on September 25, 2018, 1:02 PM
 */

#ifndef THDMWEWPO_H
#define THDMWEWPO_H

#include "ThObservable.h"
#include "THDMW.h"
#include "THDMWcache.h"
#include "StandardModel.h"

/**
 * @class Rb0
 * @ingroup THDMW
 * @brief 
 */
class  Rb0THDMW: public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Rb0THDMW(const StandardModel& SM_i);

    /**
     * @return Rb0GTHDM
     */
    double computeThValue ();
private:
    const THDMW& myTHDMW;
};




#endif /* THDMWEWPO_H */
