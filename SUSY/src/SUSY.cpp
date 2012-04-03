/* 
 * File:   SUSY.cpp
 * Author: marco
 * 
 * Created on December 2, 2010, 3:32 PM
 */

#include "SUSY.h"
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <FHCouplings.h>

const std::string SUSY::SUSYvars[NSUSYvars] = {"m1r", "m1i", "m2r", "m2i", "m3r", 
      "m3i", "muHr", "muHi", "mHp", "tanb", "Q"};

SUSY::SUSY() :
        StandardModel(), Ru(6,6,0.), Rd(6,6,0.), Rl(6,6,0.), Rn(6,6,0.),
        U(2,2,0.), V(2,2,0.), N(4,4,0.), Msu2(6,0.), Msd2(6,0.), Msl2(6,0.),
        Msn2(6,0.), Mch(2,0.), Mneu(4,0.), MsQ2(3,3,0.), MsU2(3,3,0.), 
        MsD2(3,3,0.), MsL2(3,3,0.), MsE2(3,3,0.), MsN2(3,3,0.), TU(3,3,0.), 
        TD(3,3,0.), TE(3,3,0.), TN(3,3,0.), CMsQ2(3,3,0.), CMsU2(3,3,0.), CMsD2(3,3,0.), 
        CMsL2(3,3,0.), CMsE2(3,3,0.), CMsN2(3,3,0.), CTU(3,3,0.), CTD(3,3,0.), 
        CTE(3,3,0.), CTN(3,3,0.), UH(3,3,0.), ZH(3,3,0.) {
    int err;
    FHSetFlags(&err,
            4, // Full MSSM 
            0, //DRbar
            0, //DRbar
            3, //cMSSM
            4, //UHiggs at q^2=0
            2, //two loops where available
            1, //running top mass, 
            1, //resum tan beta contributions
            3 //interpolation in phases for missing 2-loop corrections
            );
    if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHSetFlags error " << err;
        throw std::runtime_error(ss.str());
    }
}

void SUSY::SetParameter(const std::string name, const double& value){
    if(name.compare("m1r") == 0)
        m1.real() = value;
    else if(name.compare("m1i") == 0)
        m1.imag() = value;
    else if(name.compare("m2r") == 0)
        m2.real() = value;
    else if(name.compare("m2i") == 0)
        m2.imag() = value;
    else if(name.compare("m3") == 0)
        m3 = value;
    else if(name.compare("muHr") == 0)
        muH.real() = value;
    else if(name.compare("muHi") == 0)
        muH.imag() = value;
    else if(name.compare("mHp") == 0)
        mHp = value;
    else if(name.compare("tanb") == 0)
        tanb = value;
    else if(name.compare("Q") == 0)
        Q = value;
    else
        StandardModel::SetParameter(name,value);
}
///////////////////////////////////////////////////////////////////////////

bool SUSY::CheckParameters(const std::map<std::string, double>& DPars) {
    for(int i=0;i<NSUSYvars;i++)
        if(DPars.find(SUSYvars[i])==DPars.end()) {
            std::cout << "missing mandatory SUSY parameter " << SUSYvars[i] << std::endl;
            return false;
        }
    return(StandardModel::CheckParameters(DPars));
}

void SUSY::SetFeynHiggsPars(void) {
    int err;
    FHSetSMPara(&err,
            (1.-314.98e-4-dAle5Mz)/ale,
            AlsMz, GF, 
            leptons[ELECTRON].getMass(), quarks[UP].getMass(), quarks[DOWN].getMass(),
            leptons[MU].getMass(), quarks[CHARM].getMass(), quarks[STRANGE].getMass(),
            leptons[TAU].getMass(),                           quarks[BOTTOM].getMass(),
            Mw_tree(), Mz,
            lambda, A, rhob, etab);
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHSetSMPara error " << err;
        throw std::runtime_error(ss.str());
     }
    double x1 = v1()/sqrt(2.);
    double x2 = v2()/sqrt(2.);
    FHSetPara(&err, mut/quarks[TOP].getMass(),
            mtpole, tanb, -1, //using MHp
            mHp,
            sqrt(MsL2(2,2).real()), sqrt(MsE2(2,2).real()), sqrt(MsQ2(2,2).real()),
            sqrt(MsU2(2,2).real()), sqrt(MsD2(2,2).real()),
            sqrt(MsL2(1,1).real()), sqrt(MsE2(1,1).real()), sqrt(MsQ2(1,1).real()),
            sqrt(MsU2(1,1).real()), sqrt(MsD2(1,1).real()),
            sqrt(MsL2(0,0).real()), sqrt(MsE2(0,0).real()), sqrt(MsQ2(0,0).real()),
            sqrt(MsU2(0,0).real()), sqrt(MsD2(0,0).real()),
            ToComplex2(muH.real(),muH.imag()),
            ToComplex2(TE(2,2).real(),TE(2,2).imag())*x1/leptons[TAU].getMass(),
            ToComplex2(TU(2,2).real(),TU(2,2).imag())*x2/MS2DRqmass(Mrun(Q,quarks[TOP].getMass())),
            ToComplex2(TD(2,2).real(),TD(2,2).imag())*x1/MS2DRqmass(Mrun(Q,quarks[BOTTOM].getMass())),
            ToComplex2(TE(1,1).real(),TE(1,1).imag())*x1/leptons[MU].getMass(),
            ToComplex2(TU(1,1).real(),TU(1,1).imag())*x2/MS2DRqmass(Mrun(Q,quarks[CHARM].getMass())),
            ToComplex2(TD(1,1).real(),TD(1,1).imag())*x1/MS2DRqmass(Mrun(Q,quarks[STRANGE].getMass())),
            ToComplex2(TE(0,0).real(),TE(0,0).imag())*x1/leptons[ELECTRON].getMass(),
            ToComplex2(TU(0,0).real(),TU(0,0).imag())*x2/MS2DRqmass(Mrun(Q,quarks[UP].getMass())),
            ToComplex2(TD(0,0).real(),TD(0,0).imag())*x1/MS2DRqmass(Mrun(Q,quarks[DOWN].getMass())),
            ToComplex2(m1.real(),m1.imag()),ToComplex2(m2.real(),m2.imag()),ToComplex2(m3,0.),
            Q,Q,Q);
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHSetPara error " << err;
        throw std::runtime_error(ss.str());
     }
    FHSetNMFV(&err, ToComplex2(MsL2(0,1).real(),MsL2(0,1).imag())/sqrt(MsL2(0,0).real()*MsL2(1,1).real()),
            ToComplex2(MsL2(1,2).real(),MsL2(1,2).imag())/sqrt(MsL2(1,1).real()*MsL2(2,2).real()),
            ToComplex2(MsL2(0,2).real(),MsL2(0,2).imag())/sqrt(MsL2(0,0).real()*MsL2(2,2).real()),
            ToComplex2(TU(0,1).real(),-TU(0,1).imag())*x2/sqrt(MsL2(0,0).real()*MsU2(1,1).real()),
            ToComplex2(TU(1,2).real(),-TU(1,2).imag())*x2/sqrt(MsL2(1,1).real()*MsU2(2,2).real()),
            ToComplex2(TU(0,2).real(),-TU(0,2).imag())*x2/sqrt(MsL2(0,0).real()*MsU2(2,2).real()),
            ToComplex2(TU(1,0).real(),TU(1,0).imag())*x2/sqrt(MsU2(0,0).real()*MsL2(1,1).real()),
            ToComplex2(TU(2,1).real(),TU(2,1).imag())*x2/sqrt(MsU2(1,1).real()*MsL2(2,2).real()),
            ToComplex2(TU(2,0).real(),TU(2,0).imag())*x2/sqrt(MsU2(0,0).real()*MsL2(2,2).real()),
            ToComplex2(MsU2(0,1).real(),MsU2(0,1).imag())/sqrt(MsU2(0,0).real()*MsU2(1,1).real()),
            ToComplex2(MsU2(1,2).real(),MsU2(1,2).imag())/sqrt(MsU2(1,1).real()*MsU2(2,2).real()),
            ToComplex2(MsU2(0,2).real(),MsU2(0,2).imag())/sqrt(MsU2(0,0).real()*MsU2(2,2).real()),
            ToComplex2(TD(0,1).real(),-TD(0,1).imag())*x1/sqrt(MsL2(0,0).real()*MsD2(1,1).real()),
            ToComplex2(TD(1,2).real(),-TD(1,2).imag())*x1/sqrt(MsL2(1,1).real()*MsD2(2,2).real()),
            ToComplex2(TD(0,2).real(),-TD(0,2).imag())*x1/sqrt(MsL2(0,0).real()*MsD2(2,2).real()),
            ToComplex2(TD(1,0).real(),TD(1,0).imag())*x1/sqrt(MsD2(0,0).real()*MsL2(1,1).real()),
            ToComplex2(TD(2,1).real(),TD(2,1).imag())*x1/sqrt(MsD2(1,1).real()*MsL2(2,2).real()),
            ToComplex2(TD(2,0).real(),TD(2,0).imag())*x1/sqrt(MsD2(0,0).real()*MsL2(2,2).real()),
            ToComplex2(MsD2(0,1).real(),MsD2(0,1).imag())/sqrt(MsD2(0,0).real()*MsD2(1,1).real()),
            ToComplex2(MsD2(1,2).real(),MsD2(1,2).imag())/sqrt(MsD2(1,1).real()*MsD2(2,2).real()),
            ToComplex2(MsD2(0,2).real(),MsD2(0,2).imag())/sqrt(MsD2(0,0).real()*MsD2(2,2).real()));
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHSetNMFV error " << err;
        throw std::runtime_error(ss.str());
     }
}

void SUSY::CalcHiggsSpectrum(void){
    int err;
    double_complex SAeff;
    double_complex UHiggs[3][3];
    double_complex ZHiggs[3][3];
    FHHiggsCorr(&err, mh, &SAeff, UHiggs, ZHiggs);
    saeff = complex(SAeff.real(),SAeff.imag()); 
    mHp = mh[3];
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++){
            UH.assign(i,j,complex(UHiggs[i][j].real(),UHiggs[i][j].imag()));
            ZH.assign(i,j,complex(ZHiggs[i][j].real(),ZHiggs[i][j].imag()));            
        }
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHHiggsCorr error " << err;
        throw std::runtime_error(ss.str());
     }
}

void SUSY::CalcHiggsCouplings(void){
    int err;
    double_complex couplings[ncouplings];
    double_complex couplingsms[ncouplingsms];
    double gammas[ngammas];
    double gammasms[ngammasms];
    FHCouplings(&err, couplings, couplingsms, gammas, gammasms, 0);
    if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHCouplings error " << err;
        throw std::runtime_error(ss.str());
     }
}

void SUSY::CalcHiggsProd(const double& sqrts){
    int err;
    double prodxs[nprodxs];
    FHHiggsProd(&err, sqrts, prodxs);
    if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHHiggsProd error " << err;
        throw std::runtime_error(ss.str());
     }
}

void SUSY::CalcConstraints(){
    int err;
    FHConstraints(&err, &FHgm2, &FHdeltarho, &FHMWMSSM, &FHMWSM, &FHSW2MSSM, 
            &FHSW2SM, &FHedmeTh, &FHedmn, &FHedmHg);
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHConstraints error " << err;
        throw std::runtime_error(ss.str());
     }   
}

void SUSY::CalcFlavour(){
    int err;
    FHFlavour(&err, &FHbsgMSSM, &FHbsgSM, &FHdeltaMsMSSM, &FHdeltaMsSM,
            &FHbsmumuMSSM, &FHbsmumuSM);
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHFlavour error " << err;
        throw std::runtime_error(ss.str());
     }   
}

void SUSY::CalcSpectrum(){
    int err, nmfv;
    double MASf[4][6], MCha[2], MNeu[4];
    double_complex UASf[4][6][6], UCha[2][2], VCha[2][2], ZNeu[4][4], Deltab;
    FHGetPara(&err, &nmfv, MASf, UASf, MCha, UCha, VCha, MNeu, ZNeu, &Deltab, 
            &FHMGl, FHMHtree, &FHSAtree);
    for(int i = 0; i < 6; i++){
        Msn2(i) = MASf[0][i]*MASf[0][i];
        Msl2(i) = MASf[1][i]*MASf[1][i];
        Msu2(i) = MASf[2][i]*MASf[2][i];
        Msd2(i) = MASf[3][i]*MASf[3][i];
        for(int j = 0; j < 6; j++){
            Rn.assign(i,j,complex(UASf[0][i][j].real(),UASf[0][i][j].imag()));
            Rl.assign(i,j,complex(UASf[1][i][j].real(),UASf[1][i][j].imag()));
            Ru.assign(i,j,complex(UASf[2][i][j].real(),UASf[2][i][j].imag()));
            Rd.assign(i,j,complex(UASf[3][i][j].real(),UASf[3][i][j].imag()));
        }
    }
    for(int i = 0; i < 4; i++){
        Mneu(i) = MNeu[i];
        for(int j = 0; j < 4; j++)
            N.assign(i,j,complex(ZNeu[i][j].real(),ZNeu[i][j].imag()));
    }
    for(int i = 0; i < 2; i++){
        Mch(i) = MCha[i];
        for(int j = 0; j < 2; j++){
            U.assign(i,j,complex(UCha[i][j].real(),UCha[i][j].imag()));
            V.assign(i,j,complex(VCha[i][j].real(),VCha[i][j].imag()));
        }
    }
    FHDeltab = complex(Deltab.real(),Deltab.imag());
    if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHGetPara error " << err;
        throw std::runtime_error(ss.str());
     }   
   
}

double SUSY::v1() {
    return v()*cosb;
}

double SUSY::v2() {
    return v()*sinb;
}


///////////////////////////////////////////////////////////////////////////

void SUSY::setY(double tanb_i) {

    setTanb(tanb_i);

//    Yd.assign(0,0,md/v1()*sqrt(2.));
//    Yd.assign(1,1,ms/v1()*sqrt(2.));
//    Yd.assign(2,2,mb/v1()*sqrt(2.));
//    Yu.assign(0,0,mu/v2()*sqrt(2.));
//    Yu.assign(1,1,mc/v2()*sqrt(2.));
//    Yu.assign(2,2,mt/v2()*sqrt(2.));
//    Yu = Yu*VCKM;
}

void SUSY::setTanb(double tanb) {
    sinb = tanb * sqrt(1. / (1. + tanb * tanb));
    cosb = sqrt(1. / (1. + tanb * tanb));
    this->tanb = tanb;
}

void SUSY::setSinb(double sinb) {
    cosb = sqrt(1. - sinb * sinb);
    tanb = sinb / cosb;
    this->sinb = sinb;
}

void SUSY::setCosb(double cosb) {
    sinb = sqrt(1. - cosb * cosb);
    tanb = sinb / cosb;
    this->cosb = cosb;
}


///////////////////////////////////////////////////////////////////////////

double SUSY::Mw() const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::Mw() " << std::endl;
    return (80.3613);
}

gslpp::complex SUSY::gZf(const int INDF) const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::gZf() " << std::endl;
    gslpp::complex tmp(0.0738065, -0.0120949, false);
    return (tmp);
}

gslpp::complex SUSY::rhoZf(const int INDF) const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::rhoZf() " << std::endl;
    gslpp::complex tmp(1.00516, -0.00473674, false);
    return (tmp);
}

double SUSY::Delta_r() const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::Delta_r() " << std::endl;
    return (0.0378211);
}





