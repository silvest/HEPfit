/* 
 * File:   ZFitterTest-2.cpp
 * Author: mishima
 */

/*
 *  Test for Table 6.1 in hep-ph/0507146
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <StandardModel.h>
#include "ZFitter.h"
#include "ZFitterObservables.h"


void setSMparameters(StandardModel& SM_i) {
    std::map<std::string, double> Parameters;
    // 17 parameters defined in StandardModel
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
    //Parameters["dAle5Mz"] = 0.02758;
    //Parameters["mHl"] = 130.0;
    Parameters["muw"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
    // 26 parameters defined in QCD    
    //Parameters["AlsMz"] = 0.1184;
    //Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    //Parameters["mtop"] = 174.0;
    Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
    Parameters["MBd"] = 0.0;
    Parameters["MBs"] = 0.0;
    Parameters["MBp"] = 0.0;
    Parameters["MK0"] = 0.0;
    Parameters["MKp"] = 0.0;
    Parameters["FBs"] = 0.0;
    Parameters["FBsoFBd"] = 0.0;
    Parameters["BBsoBBd"] = 0.0;
    Parameters["BBs1"] = 0.0;
    Parameters["BBs2"] = 0.0;
    Parameters["BBs3"] = 0.0;
    Parameters["BBs4"] = 0.0;
    Parameters["BBs5"] = 0.0;
    Parameters["BBsscale"] = 0.0;
    Parameters["BBsscheme"] = 0.0;

    /* TEST for Table 6.1 in hep-ph/0507146*/
    Parameters["Mz"] = 91.1875;
    Parameters["mtop"] = 175.0;    
    Parameters["mHl"] = 150.0;
    Parameters["AlsMz"] = 0.118;
    Parameters["dAle5Mz"] = 0.02758;
    
    SM_i.Init(Parameters);
}


int main(int argc, char** argv) {
 
    StandardModel* myModel;
    myModel = new StandardModel();
    setSMparameters(*myModel);
    
    /* call a constructor of ZFitter with the above parameters (essential)
     * set flags and cuts to be their defalut values */
    ZFitterObservables ZFO(*myModel);
    
    int indexFermion;

    /* sqrt(s) = Mz */
    //double sqrt_s = ZFO.getZMASS();
    //double s = sqrt_s*sqrt_s;

    /* print input parameters (not essential) */
    ZFO.printInputs();

    /* set flags (not necessary if using the default flags) */
    // default: AMT4=4, ALEM=3
    int AFBC = 1, SCAL = 0, SCRE = 0, AMT4 = 6, BORN = 0,
        BOXD = 1, CONV = 1, FINR = 1, FOT2 = 3, GAMS = 1,
        DIAG = 1, INTF = 1, BARB = 2, PART = 0, POWR = 1,
        PRNT = 0, ALEM = 2, QCDC = 3, VPOL = 1, WEAK = 1,
        FTJR = 1, EXPR = 0, EXPF = 0, HIGS = 0, AFMT = 3,
        CZAK = 1, PREC = 10,HIG2 = 0, ALE2 = 3, GFER = 2,
        ISPP = 2, FSRS = 1, MISC = 0, MISD = 1, IPFC = 5,
        IPSC = 0, IPTO = 3, FBHO = 0, FSPP = 0, FUNA = 0,
        ASCR = 1, SFSR = 1, ENUE = 1, TUPV = 1, DMWW = 0,
        DSWW = 0;
    int flags[46] = {AFBC, SCAL, SCRE, AMT4, BORN, BOXD, CONV, FINR, FOT2, GAMS,
                     DIAG, INTF, BARB, PART, POWR, PRNT, ALEM, QCDC, VPOL, WEAK,
                     FTJR, EXPR, EXPF, HIGS, AFMT, CZAK, PREC, HIG2, ALE2, GFER,
                     ISPP, FSRS, MISC, MISD, IPFC, IPSC, IPTO, FBHO, FSPP, FUNA,
                     ASCR, SFSR, ENUE, TUPV, DMWW, DSWW};
    ZFO.setAllFlags(flags, 1);

    /* set cuts (not necessary if using the default cuts) */
    /*
    int ICUT[12];
    double ACOL[12], EMIN[12], S_PR[12], ANG0[12], ANG1[12], SPP[12];
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        ICUT[indexFermion] = -1;
        ACOL[indexFermion] = 0.0;
        EMIN[indexFermion] = 0.0;
        S_PR[indexFermion] = 0.01*s;
        ANG0[indexFermion] = 0.0;
        ANG1[indexFermion] = 180.0;
        SPP[indexFermion]  = 0.25*s;
    }
    ICUT[11] = 0;
    ACOL[11] = 15.0;
    EMIN[11] = 10.0;
    S_PR[11] = 0.0;
    ANG0[11] = 35.0;
    ANG1[11] = 145.0;
    SPP[11]  = 0.0;
    ZFO.setAllCuts(ICUT, ACOL, EMIN, S_PR, ANG0, ANG1, SPP, 1);
    std::cout << std::endl;
     */

    /* calculate EW common blocks (necessary before computing observables) */
    ZFO.calcCommonBlocks();
    std::cout << std::endl;

    /* ptint constants defiend in ZFITTER */
    ZFO.printConstants();

    /* calcualte and print EW precision obserbables */
    std::cout << "###  see Table 6.1 in hep-ph/0507146  ###" << std::endl;
    ZFO.printPO();

    return (EXIT_SUCCESS);
}

