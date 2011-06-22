/* 
 * File:   StandardModel.cpp
 * Author: silvest
 * 
 * Created on November 30, 2010, 1:27 PM
 */


//#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include "StandardModel.h"
#include "CKM.h"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>

const std::string StandardModel::SMvars[NSMvars] = {"GF", "mneutrino_1", "mneutrino_2",
    "mneutrino_3", "melectron", "mmu", "mtau", "lambda", "A", "rhob", "etab", "ale",
    "dAle5Mz", "mHl", "muw", "mub", "muc"};

void StandardModel::Update(const std::map<std::string, double>& DPars) {
    computeCKM = false;
    computeYe = false;
    computeYn = false;
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetSMParameter(it->first, it->second);
    QCD::Update(DPars);
    if (computeCKM) {
        myCKM.setWolfenstein(lambda, A, rhob, etab);
        myCKM.getCKM(VCKM);
    }
    if (computeYu || computeCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v() * sqrt(2.));
        Yu = Yu*VCKM;
    }
    if (computeYd) {
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (computeYe) {
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (computeYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }
}

void StandardModel::SetSMParameter(std::string name, double value) {
    if (name.compare("GF") == 0)
        GF = value;
    else if (name.compare("ale") == 0)
        ale = value;
    else if (name.compare("dAle5Mz") == 0)
        dAle5Mz = value;
    else if (name.compare("mHl") == 0)
        mHl = value;
    else if (name.compare("muw") == 0)
        muw = value;
    else if (name.compare("mub") == 0)
        mub = value;
    else if (name.compare("muc") == 0)
        muc = value;
    else if (name.compare("mneutrino_1") == 0) {
        leptons[NEUTRINO_1].setMass(value);
        computeYn = true;
    } else if (name.compare("mneutrino_2") == 0) {
        leptons[NEUTRINO_2].setMass(value);
        computeYn = true;
    } else if (name.compare("mneutrino_3") == 0) {
        leptons[NEUTRINO_3].setMass(value);
        computeYn = true;
    } else if (name.compare("melectron") == 0) {
        leptons[ELECTRON].setMass(value);
        computeYe = true;
    } else if (name.compare("mmu") == 0) {
        leptons[MU].setMass(value);
        computeYe = true;
    } else if (name.compare("mtau") == 0) {
        leptons[TAU].setMass(value);
        computeYe = true;
    } else if (name.compare("lambda") == 0) {
        lambda = value;
        computeCKM = true;
    } else if (name.compare("A") == 0) {
        A = value;
        computeCKM = true;
    } else if (name.compare("rhob") == 0) {
        rhob = value;
        computeCKM = true;
    } else if (name.compare("etab") == 0) {
        etab = value;
        computeCKM = true;
    } else
        SetQCDParameter(name, value);
}

bool StandardModel::Init(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NSMvars; i++) {
        if (DPars.find(SMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory SM parameter " << SMvars[i] << std::endl;
            return false;
        }
    }
    for (int i = 0; i < NQCDvars; i++) {
        if (DPars.find(QCDvars[i]) == DPars.end()) {
            std::cout << "missing mandatory QCD parameter " << QCDvars[i] << std::endl;
            return false;
        }
    }
    Update(DPars);
    return true;
}

StandardModel::StandardModel() : QCD(), VCKM(3, 3, 0.), UPMNS(3, 3, 0.), Yu(3, 3, 0.),
Yd(3, 3, 0.), Yn(3, 3, 0.), Ye(3, 3, 0.) {
}

///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const {
    return 1. / sqrt(sqrt(2.) * GF);
}


////////////////////////////////////////////////////////////////////////

double StandardModel::dAleLepMz() const {
    /*
     * Eqs.(14-15) in J.H.Kuhn, M.Steinhauser, PLB437,425(1998) [hep-ph/9802241]
     * Eqs.(5-10) in M.Steinhauser, PLB429,158(1998) [hep-ph/9803313]
     *
     * oneLoop=314.19007, twoLoop=0.77617, threeLoop=0.01063
     * sum=314.97686
     *
     * Notes: oneLoop and twoLoop are OK for me=0.00051099907, mmu=0.105658389,
     *        mtau=1.777, ale=1.0/137.0359895 and mZ=91.187, but only threeLoop
     *        differs from the above value by 5% (Why?).
     */
    double xl[3] = {mZ * mZ / leptons[ELECTRON].getMass() / leptons[ELECTRON].getMass(),
        mZ * mZ / leptons[MU].getMass() / leptons[MU].getMass(),
        mZ * mZ / leptons[TAU].getMass() / leptons[TAU].getMass()};
    double log_l[3] = {log(xl[0]), log(xl[1]), log(xl[2])};
    double log2 = log(2.0);

    /* zeta functions */
    double zeta2 = gsl_sf_zeta_int(2);
    double zeta3 = gsl_sf_zeta_int(3);
    double zeta5 = gsl_sf_zeta_int(5);

    double oneLoop[3], twoLoop[3], threeLoop[3];
    for (int i = 0; i < 3; i++) {
        oneLoop[i] = -5.0 / 9.0 + log_l[i] / 3.0 - 2.0 / xl[i];
        twoLoop[i] = -5.0 / 24.0 + zeta3 + log_l[i] / 4.0 + 3.0 / xl[i] * log_l[i];
        threeLoop[i] = -121.0 / 48.0 + (-5.0 + 8.0 * log2) * zeta2 - 99.0 / 16.0 * zeta3
                + 10.0 * zeta5 + log_l[i] / 8.0;
        for (int j = 0; j < 3; j++) {
            if (i > j) { /* Pi^{(2)}_l */
                threeLoop[i] += -116.0 / 27.0 + 4.0 / 3.0 * zeta2 + 38.0 / 9.0 * zeta3
                        + 14.0 / 9.0 * log_l[i]
                        + (5.0 / 18.0 - 4.0 / 3.0 * zeta3) * log_l[j]
                        + log_l[i] * log_l[i] / 6.0
                        - log_l[i] * log_l[j] / 3.0;
            } else if (i == j) { /* Pi^{(2)}_F */
                threeLoop[i] += -307.0 / 216.0 - 8.0 / 3.0 * zeta2 + 545.0 / 144.0 * zeta3
                        + (11.0 / 6.0 - 4.0 / 3.0 * zeta3) * log_l[i]
                        - log_l[i] * log_l[i] / 6.0;
            } else { /* Pi^{(2)}_h */
                threeLoop[i] += -37.0 / 6.0 + 38.0 / 9.0 * zeta3
                        + (11.0 / 6.0 - 4.0 / 3.0 * zeta3) * log_l[j]
                        - log_l[j] * log_l[j] / 6.0;
            }
        }
        threeLoop[i] /= -4.0;
    }

    /* TESTS */
    //for (int i=0; i<3; i++) {
    //    std::cout << ale/M_PI*oneLoop[i] << "  "
    //              << ale/M_PI*ale/M_PI*twoLoop[i] << "  "
    //              << ale/M_PI*ale/M_PI*ale/M_PI*threeLoop[i] << std::endl;
    //}

    return ( ale / M_PI * (oneLoop[0] + oneLoop[1] + oneLoop[2])
            + ale / M_PI * ale / M_PI * (twoLoop[0] + twoLoop[1] + twoLoop[2])
            + ale / M_PI * ale / M_PI * ale / M_PI
            * (threeLoop[0] + threeLoop[1] + threeLoop[2]));
}

double StandardModel::dAleTopMz() const {
    /*
     * Eq.(12) in J.H.Kuhn, M.Steinhauser, PLB437,425(1998) [hep-ph/9802241]
     *
     * (-0.70+-0.05)*10^{-4} for mt=175.6+-5.5 and alpha_s(mZ)=0.118
     */
    double xt = mZ * mZ / quarks[TOP].getMass() / quarks[TOP].getMass();
    double log_t = log(xt);

    return ( -4.0 / 45.0 * ale / M_PI * xt
            * (1.0 + 5.062 * alsMz / M_PI
            + (28.220 + 9.702 * log_t) * alsMz / M_PI * alsMz / M_PI
            + xt * (0.1071 + 0.8315 * alsMz / M_PI
            + (6.924 + 1.594 * log_t) * alsMz / M_PI * alsMz / M_PI)));
}

double StandardModel::dAleTotalMz() const {
    return ( dAle5Mz + dAleLepMz() + dAleTopMz());
}

double StandardModel::aleMz() const {
    return ( ale / (1.0 - dAleTotalMz()));
}

double StandardModel::mcMz() const {
    double mc_at_mb = Mrun(quarks[BOTTOM].getMass(), quarks[CHARM].getMass(), 4.0);
    double mc_at_mZ = Mrun(mZ, mc_at_mb, 5.0);

    /* TEST */
    //std::cout << "mc(mc)= " << quarks[CHARM].getMass() << std::endl;
    //std::cout << "mc(mb)_LO+NLO= " << mc_at_mb << std::endl;
    //std::cout << "mc(mZ)_LO+NLO= " << mc_at_mZ << std::endl;

    return ( mc_at_mZ);
    //return ( 0.563817 );// <--- used in ZFITTER with the effective mass mc=1.5
}

double StandardModel::mbMz() const {
    /* TEST */
    //std::cout << "mb(mb)= " << quarks[BOTTOM].getMass() << std::endl;
    //std::cout << "mb(mZ)_LO= " << mrun(mZ, quarks[BOTTOM].getMass(), 5.0, 0) << std::endl;
    //std::cout << "mb(mZ)_LO+NLO= " << mrun(mZ, quarks[BOTTOM].getMass(), 5.0, 1) << std::endl;

    return ( Mrun(mZ, quarks[BOTTOM].getMass(), 5.0));
    //return (2.819440);// <--- used in ZFITTER with the effective mass mb=4.7
}

double StandardModel::mW() const {

    //std::cout << "Write codes for StandardModel::mW() " << std::endl;
    return (80.3613);
}

complex StandardModel::gZf(const int INDF) const {

    std::cout << "Write codes for StandardModel::gZf() " << std::endl;
    complex tmp(0.0738065, -0.0120949, false);
    return (tmp);
}

complex StandardModel::rhoZf(const int INDF) const {

    std::cout << "Write codes for StandardModel::rhoZf() " << std::endl;
    complex tmp(1.00516, -0.00473674, false);
    return (tmp);
}

double StandardModel::Delta_r() const {

    std::cout << "Write codes for StandardModel::Delta_r() " << std::endl;
    return (0.0378211);
}


////////////////////////////////////////////////////////////////////////

// Angles

double StandardModel::getBeta() const {
    return (-VCKM(1, 0) * VCKM(1, 2).conjugate() / (VCKM(2, 0) * VCKM(2, 2).conjugate())).arg();
}

double StandardModel::getGamma() const {
    return (-VCKM(0, 0) * VCKM(0, 2).conjugate() / (VCKM(1, 0) * VCKM(1, 2).conjugate())).arg();
}

double StandardModel::getAlpha() const {
    return (-VCKM(2, 0) * VCKM(2, 2).conjugate() / (VCKM(0, 0) * VCKM(0, 2).conjugate())).arg();
}

double StandardModel::getBetas() const {
    return (-VCKM(2, 1) * VCKM(2, 2).conjugate() / (VCKM(1, 1) * VCKM(1, 2).conjugate())).arg();
}

// Lambda_q

complex StandardModel::getlamt() const {
    return VCKM(2, 0) * VCKM(2, 1).conjugate();
}

complex StandardModel::getlamc() const {
    return VCKM(1, 0) * VCKM(1, 1).conjugate();
}

complex StandardModel::getlamu() const {
    return VCKM(0, 0) * VCKM(0, 1).conjugate();
}

complex StandardModel::getlamt_d() const {
    return VCKM(2, 0) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_d() const {
    return VCKM(1, 0) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_d() const {
    return VCKM(0, 0) * VCKM(0, 2).conjugate();
}

complex StandardModel::getlamt_s() const {
    return VCKM(2, 1) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_s() const {
    return VCKM(1, 1) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_s() const {
    return VCKM(0, 1) * VCKM(0, 2).conjugate();
}

