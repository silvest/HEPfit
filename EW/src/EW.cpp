/* 
 * File:   EW.cpp
 * Author: mishima
 */

#include "EW.h"


EW::EW(const StandardModel& SM_i) : ThObsType(SM_i), myEWSM(SM_i), myZFitter(SM_i) {
}

//EW::EW(const EW& orig) : ThObsType(orig.SM), EWSM(orig.SM) {
//}

EW::~EW() {
}


////////////////////////////////////////////////////////////////////////

void EW::ComputeEWSM() {
    myEWSM.Compute();
    
    alphaMZ = myEWSM.getAlphaMZ();
    Mw = myEWSM.getMw();    
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        rhoZ_l[i] = myEWSM.getRhoZ_l(i_l);
        rhoZ_q[i] = myEWSM.getRhoZ_q(i_q);
        kappaZ_l[i] = myEWSM.getKappaZ_l(i_l);
        kappaZ_q[i] = myEWSM.getKappaZ_q(i_q);
    }
}

void EW::ComputeZFitter() {
    /* set flags */
    // default: AMT4=4, ISPP=2
    int AFBC = 1, SCAL = 0, SCRE = 0, AMT4 = 6, BORN = 0,
        BOXD = 1, CONV = 1, FINR = 1, FOT2 = 3, GAMS = 1,
        DIAG = 1, INTF = 1, BARB = 2, PART = 0, POWR = 1,
        PRNT = 0, ALEM = 3, QCDC = 3, VPOL = 1, WEAK = 1,
        FTJR = 1, EXPR = 0, EXPF = 0, HIGS = 0, AFMT = 3,
        CZAK = 1, PREC = 10,HIG2 = 0, ALE2 = 3, GFER = 2,
        ISPP = 1, FSRS = 1, MISC = 0, MISD = 1, IPFC = 5,
        IPSC = 0, IPTO = 3, FBHO = 0, FSPP = 0, FUNA = 0,
        ASCR = 1, SFSR = 1, ENUE = 1, TUPV = 1, DMWW = 0,
        DSWW = 0;
    int flags[46] = {AFBC, SCAL, SCRE, AMT4, BORN, BOXD, CONV, FINR, FOT2, GAMS,
                     DIAG, INTF, BARB, PART, POWR, PRNT, ALEM, QCDC, VPOL, WEAK,
                     FTJR, EXPR, EXPF, HIGS, AFMT, CZAK, PREC, HIG2, ALE2, GFER,
                     ISPP, FSRS, MISC, MISD, IPFC, IPSC, IPTO, FBHO, FSPP, FUNA,
                     ASCR, SFSR, ENUE, TUPV, DMWW, DSWW};
    myZFitter.setAllFlags(flags, 1);    

    myZFitter.calcCommonBlocks();
    
    alphaMZ = myZFitter.getAlphaMZ();
    Mw = myZFitter.getMw();
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        rhoZ_l[i] = myZFitter.getRhoZ_l(i_l);
        rhoZ_q[i] = myZFitter.getRhoZ_q(i_q);
        kappaZ_l[i] = myZFitter.getKappaZ_l(i_l);
        kappaZ_q[i] = myZFitter.getKappaZ_q(i_q);
    }
}



