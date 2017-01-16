/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMCACHE_H
#define GENERALTHDMCACHE_H

#include "GeneralTHDM.h"
#include "GeneralTHDMquantities.h"

class GeneralTHDMcache {

public:

    /**
     * @brief THDMcache constructor.
     * @details Reads all the tables values and stores them in the memory.
     */
    GeneralTHDMcache(const StandardModel& SM_i);

    /**
     * @brief THDMcache destructor.
     * @details Reads all the tables values and stores them in the memory.
     */
    ~GeneralTHDMcache();

    void updateCache();

    double mH1_2;
    double mH2_2;
    double mH3_2;
    double M11_2;
    double M12_2;
    double M13_2;
    double M22_2;
    double M23_2;
    double M33_2;
    double M2_GTHDM;
    double m11_2_GTHDM;
    double m22_2_GTHDM;
    double Imm12_2_GTHDM;
    double lambda1_GTHDM;
    double lambda2_GTHDM;
    double lambda3_GTHDM;
    double lambda4_GTHDM;
    double Relambda5_GTHDM;

private:

    const GeneralTHDM * myGTHDM;

    double mHl;
    double vev;
    double tanb;
    double cosb;
    double sinb;
    double cosalpha1;
    double sinalpha1;
    double cosalpha2;
    double sinalpha2;
    double cosalpha3;
    double sinalpha3;
    double Imlambda5;
    double Imlambda6;
    double Imlambda7;
    double mHp2;
    double Relambda6;
    double Relambda7;
//    double Q_THDM;
//    double bma;
//    double m12_2;
//    double mHh2;
//    double mA2;
//    double MW;
//    double cW2;
//    double Ale;
//    double Als;
//    double Mt;
//    double Mb;
//    double Mtau;
//    double Mc;
//    double Ms;
//    double Mmu;
//    double Mu;
//    double Md;
//    double Me;
//    double MZ;
};

#endif /* GENERALTHDMCACHE_H */

