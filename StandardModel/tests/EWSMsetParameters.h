/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMSETPARAMETERS_H
#define	EWSMSETPARAMETERS_H

#include "StandardModel.h"

void setSMparameters(StandardModel& SM_i) {
    std::map<std::string, double> Parameters;
    // parameters defined in StandardModel
    Parameters["GF"] = 1.16637E-5;
    Parameters["mneutrino_1"] = 0.0;
    Parameters["mneutrino_2"] = 0.0;
    Parameters["mneutrino_3"] = 0.0;
    Parameters["melectron"] = 0.54857990943e-3;
    Parameters["mmu"] = 0.1134289256;
    Parameters["mtau"] = 1.77682;
    Parameters["lambda"] = 0.2253;
    Parameters["A"] = 0.808;
    Parameters["rhob"] = 0.132;
    Parameters["etab"] = 0.341;
    Parameters["ale"] = 1.0/137.035999679;
    Parameters["dAle5Mz"] = 0.02758;
    Parameters["mHl"] = 130.0;
    Parameters["muw"] = 80.0;
    
    // parameters defined in QCD    
    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    Parameters["mtop"] = 174.0;
    Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 175.0;
    Parameters["mub"] = 5.0;
    Parameters["muc"] = 1.5;
    
    /** Test for alpha_lep **/
    //Parameters["melectron"] = 0.00051099907;
    //Parameters["mmu"] = 0.105658389;
    //Parameters["mtau"] = 1.777;    
    //Parameters["AlsMz"] = 1.0/137.0359895;
    //Parameters["Mz"] = 91.187;    
    
    /** To make comparisons with ZFITTER codes **/
    Parameters["GF"] = 1.16637E-5;  // for GFER=2
    Parameters["melectron"] = 0.51099907e-3;
    Parameters["mmu"] = 0.105658389;
    Parameters["mtau"] = 1.77705;
    Parameters["mup"] = 0.062;
    Parameters["mdown"] = 0.083;
    Parameters["mcharm"] = 1.50; // In ZFitter, 1.5 is the pole mass, but we input it to loop functions by hand for comparisons. 
    Parameters["mstrange"] = 0.215;
    Parameters["mbottom"] = 4.70; // In ZFitter, 4.7 is the pole mass, but we input it to loop functions by hand for comparisons. 
    Parameters["ale"] = 1.0/137.0359895;

    /* TEST for Table 6.1 in hep-ph/0507146*/
    /* flags: AMT4=6, ALEM=2 */
    Parameters["Mz"] = 91.1875;
    Parameters["mtop"] = 175.0;    
    Parameters["mHl"] = 150.0;
    Parameters["AlsMz"] = 0.118;
    Parameters["dAle5Mz"] = 0.02758;    

    /* to make mb(Mz) and mc(Mz) similar to ZFitter ones */
    /* mcMz = 0.56381685, mbMz = 2.8194352 */
    /* muMz = 0.062, mdMz = 0.083, msMz = 0.215 */
    //Parameters["mbottom"] = 4.122;
    //Parameters["mub"] = 4.122;    
    //Parameters["mcharm"] = 1.171;
    //Parameters["muc"] = 1.171;    
    
    //SM_i.Init(Parameters);
    SM_i.Update(Parameters);
}

#endif	/* EWSMSETPARAMETERS_H */

