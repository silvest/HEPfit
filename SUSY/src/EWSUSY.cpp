/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSUSY.h"

EWSUSY::EWSUSY(const SUSY& SUSY_in)
: mySUSY(SUSY_in),
        Yu_JR(3,3,0.0), Yd_JR(3,3,0.0), Yl_JR(3,3,0.0),
        Au_JR(3,3,0.0), Ad_JR(3,3,0.0), Al_JR(3,3,0.0)
{
}

void EWSUSY::SetRosiekParameters()
{
    Yu_JR = mySUSY.getYu();
    Yd_JR = - mySUSY.getYd();
    Yl_JR = - mySUSY.getYe();

    Au_JR = - mySUSY.getTU();
    Ad_JR = mySUSY.getTD();
    Al_JR = mySUSY.getTE();
}