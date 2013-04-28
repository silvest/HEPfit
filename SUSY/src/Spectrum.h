/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SPECTRUM_H
#define	SPECTRUM_H

#include "SUSY.h"

class Spectrum {
public:
    Spectrum();
    
    void CalcSpectrum(SUSY & SUSY);
//    void Sflav_init(SUSY & SUSY);
//    void Sflav_Spectrum(const std::string);
    
private:

};

#endif	/* SPECTRUM_H */

