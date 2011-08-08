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

double StandardModel::mW() const {

    //std::cout << "Write codes for EWSM::mW() " << std::endl;
    return (80.3613);
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

