/* 
 * File:   ZFitter.cpp
 * Author: mishima
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ZFitter.h"


ZFitter::ZFitter(const StandardModel& mySM) : SM(mySM) {
    ZMASS = mySM.getMz();
    TMASS = mySM.getQuarks(mySM.TOP).getMass();
    HMASS = mySM.getMHl();
    ALFAS = mySM.getAlsMz();
    DAL5H = mySM.getDAle5Mz();
    V_TB = 1.0; // V_{tb} = 1
    UMASS = 0.1;// constituent u-quark mass
    DMASS = 0.1;// constituent d-quark mass
    
    init(0);
}

//ZFitter::ZFitter(const ZFitter& orig) {
//}

ZFitter::~ZFitter() {
}


////////////////////////////////////////////////////////////////////////

double ZFitter::getCommonALPHST() { return zupars_.ALPHST; }

double ZFitter::getCommonSIN2TW() { return zupars_.SIN2TW; }

double ZFitter::getCommonS2TEFF(const int INDF) {
    if (INDF < 12) {
        return zupars_.S2TEFF[INDF];
    } else {
        std::cout << "S2TEFF[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonALLCH(const int INDF) {
    if (INDF < 12) {
        return zfchms_.ALLCH[INDF];
    } else {
        std::cout << "ALLCH[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonALLMS(const int INDF) {
    if (INDF < 12) {
        return zfchms_.ALLMS[INDF];
    } else {
        std::cout << "ALLMS[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


////////////////////////////////////////////////////////////////////////

double ZFitter::getCommonARROFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARROFZ[INDF];
    } else {
        std::cout << "ARROFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonARKAFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARKAFZ[INDF];
    } else {
        std::cout << "ARKAFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonARVEFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARVEFZ[INDF];
    } else {
        std::cout << "ARVEFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonARSEFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARSEFZ[INDF];
    } else {
        std::cout << "ARSEFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAROTFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AROTFZ[INDF];
    } else {
        std::cout << "AROTFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAIROFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AIROFZ[INDF];
    } else {
        std::cout << "AIROFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAIKAFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AIKAFZ[INDF];
    } else {
        std::cout << "AIKAFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAIVEFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AIVEFZ[INDF];
    } else {
        std::cout << "AIVEFZ[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


////////////////////////////////////////////////////////////////////////

void ZFitter::init(const int IPRINT) { zuinit_(&IPRINT); }

void ZFitter::flag(const std::string CHFLAG, const int IVALUE) {
    zuflag_(CHFLAG.c_str(), &IVALUE, CHFLAG.size());
}

void ZFitter::setAllFlags(const int flags[46], const int flagPrint) {

    std::string flagNames[46]
            = {"AFBC", "SCAL", "SCRE", "AMT4", "BORN",
               "BOXD", "CONV", "FINR", "FOT2", "GAMS",
               "DIAG", "INTF", "BARB", "PART", "POWR",
               "PRNT", "ALEM", "QCDC", "VPOL", "WEAK",
               "FTJR", "EXPR", "EXPF", "HIGS", "AFMT",
               "CZAK", "PREC", "HIG2", "ALE2", "GFER",
               "ISPP", "FSRS", "MISC", "MISD", "IPFC",
               "IPSC", "IPTO", "FBHO", "FSPP", "FUNA",
               "ASCR", "SFSR", "ENUE", "TUPV", "DMWW",
               "DSWW" };

    /* set flags */
    int i;
    for (i=0; i<46; i++) {
        flag(flagNames[i].c_str(), flags[i]);
    }

    /* print flags if flagPrint=1 */
    if (flagPrint==1) {
        for (i=0; i<46; i++) {
            std::cout << "  " << flagNames[i].c_str() << ": " << flags[i];
            if ((i+1)%5==0) std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}

void ZFitter::FlagInfo() { 
    const int mode = 0;
    zuinfo_(&mode); 
}

void ZFitter::calcCommonBlocks() {
    zvweak_(&ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS);

    alphaMZ = calqed_.ALQEDZ;
    Mw = ZMASS*sqrt(1.0 - getCommonSIN2TW());
  
    rhoZ_l[SM.NEUTRINO_1].real() = getCommonAROTFZ(0);
    rhoZ_l[SM.ELECTRON].real() = getCommonAROTFZ(1);
    rhoZ_l[SM.NEUTRINO_2].real() = getCommonAROTFZ(0);
    rhoZ_l[SM.MU].real() = getCommonAROTFZ(2);
    rhoZ_l[SM.NEUTRINO_3].real() = getCommonAROTFZ(0);
    rhoZ_l[SM.TAU].real() = getCommonAROTFZ(3);            

    rhoZ_l[SM.NEUTRINO_1].imag() = getCommonAIROFZ(0);    
    rhoZ_l[SM.ELECTRON].imag() = getCommonAIROFZ(1);
    rhoZ_l[SM.NEUTRINO_2].imag() = getCommonAIROFZ(0);
    rhoZ_l[SM.MU].imag() = getCommonAIROFZ(2);
    rhoZ_l[SM.NEUTRINO_3].imag() = getCommonAIROFZ(0);
    rhoZ_l[SM.TAU].imag() = getCommonAIROFZ(3);            
    
    rhoZ_q[SM.UP].real() = getCommonAROTFZ(4);
    rhoZ_q[SM.DOWN].real() = getCommonAROTFZ(5);
    rhoZ_q[SM.CHARM].real() = getCommonAROTFZ(6);
    rhoZ_q[SM.STRANGE].real() = getCommonAROTFZ(7);
    rhoZ_q[SM.TOP].real() = 0.0;
    rhoZ_q[SM.BOTTOM].real() = getCommonAROTFZ(9);
    
    rhoZ_q[SM.UP].imag() = getCommonAIROFZ(4);
    rhoZ_q[SM.DOWN].imag() = getCommonAIROFZ(5);
    rhoZ_q[SM.CHARM].imag() = getCommonAIROFZ(6);
    rhoZ_q[SM.STRANGE].imag() = getCommonAIROFZ(7);
    rhoZ_q[SM.TOP].imag() = 0.0;
    rhoZ_q[SM.BOTTOM].imag() = getCommonAIROFZ(9);    
    
    kappaZ_l[SM.NEUTRINO_1].real() = getCommonARKAFZ(0);
    kappaZ_l[SM.ELECTRON].real() = getCommonARKAFZ(1);
    kappaZ_l[SM.NEUTRINO_2].real() = getCommonARKAFZ(0);
    kappaZ_l[SM.MU].real() = getCommonARKAFZ(2);
    kappaZ_l[SM.NEUTRINO_3].real() = getCommonARKAFZ(0);
    kappaZ_l[SM.TAU].real() = getCommonARKAFZ(3);            

    kappaZ_l[SM.NEUTRINO_1].imag() = getCommonAIKAFZ(0);    
    kappaZ_l[SM.ELECTRON].imag() = getCommonAIKAFZ(1);
    kappaZ_l[SM.NEUTRINO_2].imag() = getCommonAIKAFZ(0);
    kappaZ_l[SM.MU].imag() = getCommonAIKAFZ(2);
    kappaZ_l[SM.NEUTRINO_3].imag() = getCommonAIKAFZ(0);
    kappaZ_l[SM.TAU].imag() = getCommonAIKAFZ(3);            
    
    kappaZ_q[SM.UP].real() = getCommonARKAFZ(4);
    kappaZ_q[SM.DOWN].real() = getCommonARKAFZ(5);
    kappaZ_q[SM.CHARM].real() = getCommonARKAFZ(6);
    kappaZ_q[SM.STRANGE].real() = getCommonARKAFZ(7);
    kappaZ_q[SM.TOP].real() = 0.0;
    kappaZ_q[SM.BOTTOM].real() = getCommonARKAFZ(9);
    
    kappaZ_q[SM.UP].imag() = getCommonAIKAFZ(4);
    kappaZ_q[SM.DOWN].imag() = getCommonAIKAFZ(5);
    kappaZ_q[SM.CHARM].imag() = getCommonAIKAFZ(6);
    kappaZ_q[SM.STRANGE].imag() = getCommonAIKAFZ(7);
    kappaZ_q[SM.TOP].imag() = 0.0;
    kappaZ_q[SM.BOTTOM].imag() = getCommonAIKAFZ(9);    
}






