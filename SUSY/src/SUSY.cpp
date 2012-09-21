/* 
 * File:   SUSY.cpp
 * Author: marco
 * 
 * Created on December 2, 2010, 3:32 PM
 */

#include "SUSY.h"
#include "SUSYMatching.h"
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <FHCouplings.h>


const std::string SUSY::SUSYvars[NSUSYvars] = {"m1r", "m1i", "m2r", "m2i", "m3" , "muHr", "muHi", "mHptree", "tanb", "Q"};
const std::string SUSY::SUSYFlags[NSUSYFlags] = {"Flag_H","Flag_g","Flag_Chi","Flag_Chi0"};
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


bool SUSY::SetFlag(const std::string name, const bool& value){

    bool res = false;
    if(name.compare("Flag_H") == 0){
        Fh = value;
        res = true;
    }
    else if(name.compare("Flag_g") == 0){
        Fg = value;
        res = true;
    }
    else if(name.compare("Flag_Chi") == 0){
        FChi = value;
        res = true;
    }
    else if(name.compare("Flag_Chi0") == 0){
        FChi0 = value;
        res = true;
    } 
    else {
        res = StandardModel::SetFlag(name,value);
    }
        return(res);

}


void SUSY::SetParameter(const std::string name, const double& value) {
    if (name.compare("m1r") == 0)
        m1.real() = value;
    else if (name.compare("m1i") == 0)
        m1.imag() = value;
    else if (name.compare("m2r") == 0)
        m2.real() = value;
    else if (name.compare("m2i") == 0)
        m2.imag() = value;
    else if (name.compare("m3") == 0)
        m3 = value;
    else if (name.compare("muHr") == 0)
        muH.real() = value;
    else if (name.compare("muHi") == 0)
        muH.imag() = value;
    else if (name.compare("mHptree") == 0)
        mHptree = value;
    else if (name.compare("tanb") == 0)
        setTanb(value);
    else if (name.compare("Q") == 0)
        Q = value;
    else
        StandardModel::SetParameter(name, value);
        
}

bool SUSY::InitializeMatching(){
    
    mySUSYMatching = new SUSYMatching(*this);
    SetMatchingInitialized(true);
    return(true);
}

bool SUSY::PreUpdate(){
    
     if(!StandardModel::PreUpdate())  return (false);
    
     return (true);
}


bool SUSY::PostUpdate(){
    
    mySUSYMatching->Comp_mySUSYMQ();
    
    mySUSYMatching->Comp_VdDNL(0);
    mySUSYMatching->Comp_VdDNR(0);
    mySUSYMatching->Comp_VdUCL();
    mySUSYMatching->Comp_VdUCR(0);
    
    mySUSYMatching->Comp_DeltaMd();
    mySUSYMatching->Comp_mySUSY_CKM();

    if (IsFh()) mySUSYMatching->Comp_VUDHH();
    if (IsFChi0()) {
        
        mySUSYMatching->Comp_VdDNL(1);
        mySUSYMatching->Comp_VdDNR(1);
        mySUSYMatching->Comp_VuUN();
    }
    
    
    if (IsFChi()) {
        
        mySUSYMatching->Comp_VdUCR(1);
        mySUSYMatching->Comp_VuDCL();
        mySUSYMatching->Comp_VuDCR();
    }
    
    
    //mySUSYMatching->Test();
    
     return (true);
}

bool SUSY::Update(const std::map<std::string, double>& DPars) {

    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
    SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
     return (true);
     
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

bool SUSY::SetFeynHiggsPars(void) {
    
    int err;
    //FHSetDebug(3);
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
        return (false);
        //throw std::runtime_error(ss.str());
     }
    double x1 = v1()/sqrt(2.);
    double x2 = v2()/sqrt(2.);
    //FHSetDebug(2);  
        
    
    /* FeynHiggs input: the FeynHiggs's notation is related to SLHA notation by hconjugate 
       operation for T matrix and conjugate operation for mu parameter */
    
    matrix<complex> TUFH(3.,3.,0.), TDFH(3.,3.,0.), TEFH(3.,3.,0.);
    complex muHFH(0.,0.,false);
    
    TUFH = TU.hconjugate();
    TDFH = TD.hconjugate();
    TEFH = TE.hconjugate();
    muHFH = muH.conjugate(); 
    
    /** test - lines **/
//    std::cout << "TUFH = " << TUFH << std::endl;
//    std::cout << "TU = " << TU << std::endl;
//    std::cout << "TDFH = " << TDFH << std::endl;
//    std::cout << "TD = " << TD << std::endl;
//    std::cout << "TEFH = " << TEFH << std::endl;
//    std::cout << "TE = " << TE << std::endl;
//    std::cout << "muHFH = " << muHFH << std::endl;
//    std::cout << "mu = " << muH << std::endl;
    /** end - test **/
    
    FHSetPara(&err, mut/quarks[TOP].getMass(),
            mtpole, tanb, -1, //using MHptree
            mHptree,
            sqrt(MsL2(2,2).real()), sqrt(MsE2(2,2).real()), sqrt(MsQ2(2,2).real()),
            sqrt(MsU2(2,2).real()), sqrt(MsD2(2,2).real()),
            sqrt(MsL2(1,1).real()), sqrt(MsE2(1,1).real()), sqrt(MsQ2(1,1).real()),
            sqrt(MsU2(1,1).real()), sqrt(MsD2(1,1).real()),
            sqrt(MsL2(0,0).real()), sqrt(MsE2(0,0).real()), sqrt(MsQ2(0,0).real()),
            sqrt(MsU2(0,0).real()), sqrt(MsD2(0,0).real()),
            ToComplex2(muHFH.real(),muHFH.imag()),
            ToComplex2(TEFH(2,2).real(),TEFH(2,2).imag())*x1/leptons[TAU].getMass(),
            ToComplex2(TUFH(2,2).real(),TUFH(2,2).imag())*x2/MS2DRqmass(Q,Mrun(Q,quarks[TOP].getMass())),
            ToComplex2(TDFH(2,2).real(),TDFH(2,2).imag())*x1/MS2DRqmass(Q,Mrun(Q,quarks[BOTTOM].getMass())),
            ToComplex2(TEFH(1,1).real(),TEFH(1,1).imag())*x1/leptons[MU].getMass(),
            ToComplex2(TUFH(1,1).real(),TUFH(1,1).imag())*x2/MS2DRqmass(Q,Mrun(Q,quarks[CHARM].getMass())),
            ToComplex2(TDFH(1,1).real(),TDFH(1,1).imag())*x1/MS2DRqmass(Q,Mrun(Q,quarks[STRANGE].getMass_scale(),quarks[STRANGE].getMass())),  
            ToComplex2(TEFH(0,0).real(),TEFH(0,0).imag())*x1/leptons[ELECTRON].getMass(),
            ToComplex2(TUFH(0,0).real(),TUFH(0,0).imag())*x2/MS2DRqmass(Q,Mrun(Q,quarks[UP].getMass_scale(),quarks[UP].getMass())),    
            ToComplex2(TDFH(0,0).real(),TDFH(0,0).imag())*x1/MS2DRqmass(Q,Mrun(Q,quarks[DOWN].getMass_scale(),quarks[DOWN].getMass())),  
            ToComplex2(m1.real(),m1.imag()),ToComplex2(m2.real(),m2.imag()),ToComplex2(m3,0.),
            Q,Q,Q);
    
    
     if (err != 0) {
        std::stringstream ss;
         std::cout << "FeynHiggs FHSetPara error " << std::endl;
        ss << "FeynHiggs FHSetPara error " << err;
        return (false);
        //throw std::runtime_error(ss.str());
     }
    FHSetNMFV(&err, ToComplex2(MsL2(0,1).real(),MsL2(0,1).imag())/sqrt(MsL2(0,0).real()*MsL2(1,1).real()),
            ToComplex2(MsL2(1,2).real(),MsL2(1,2).imag())/sqrt(MsL2(1,1).real()*MsL2(2,2).real()),
            ToComplex2(MsL2(0,2).real(),MsL2(0,2).imag())/sqrt(MsL2(0,0).real()*MsL2(2,2).real()),
            ToComplex2(TUFH(0,1).real(),-TUFH(0,1).imag())*x2/sqrt(MsL2(0,0).real()*MsU2(1,1).real()),
            ToComplex2(TUFH(1,2).real(),-TUFH(1,2).imag())*x2/sqrt(MsL2(1,1).real()*MsU2(2,2).real()),
            ToComplex2(TUFH(0,2).real(),-TUFH(0,2).imag())*x2/sqrt(MsL2(0,0).real()*MsU2(2,2).real()),
            ToComplex2(TUFH(1,0).real(),TUFH(1,0).imag())*x2/sqrt(MsU2(0,0).real()*MsL2(1,1).real()),
            ToComplex2(TUFH(2,1).real(),TUFH(2,1).imag())*x2/sqrt(MsU2(1,1).real()*MsL2(2,2).real()),
            ToComplex2(TUFH(2,0).real(),TUFH(2,0).imag())*x2/sqrt(MsU2(0,0).real()*MsL2(2,2).real()),
            ToComplex2(MsU2(0,1).real(),MsU2(0,1).imag())/sqrt(MsU2(0,0).real()*MsU2(1,1).real()),
            ToComplex2(MsU2(1,2).real(),MsU2(1,2).imag())/sqrt(MsU2(1,1).real()*MsU2(2,2).real()),
            ToComplex2(MsU2(0,2).real(),MsU2(0,2).imag())/sqrt(MsU2(0,0).real()*MsU2(2,2).real()),
            ToComplex2(TDFH(0,1).real(),-TDFH(0,1).imag())*x1/sqrt(MsL2(0,0).real()*MsD2(1,1).real()),
            ToComplex2(TDFH(1,2).real(),-TDFH(1,2).imag())*x1/sqrt(MsL2(1,1).real()*MsD2(2,2).real()),
            ToComplex2(TDFH(0,2).real(),-TDFH(0,2).imag())*x1/sqrt(MsL2(0,0).real()*MsD2(2,2).real()),
            ToComplex2(TDFH(1,0).real(),TDFH(1,0).imag())*x1/sqrt(MsD2(0,0).real()*MsL2(1,1).real()),
            ToComplex2(TDFH(2,1).real(),TDFH(2,1).imag())*x1/sqrt(MsD2(1,1).real()*MsL2(2,2).real()),
            ToComplex2(TDFH(2,0).real(),TDFH(2,0).imag())*x1/sqrt(MsD2(0,0).real()*MsL2(2,2).real()),
            ToComplex2(MsD2(0,1).real(),MsD2(0,1).imag())/sqrt(MsD2(0,0).real()*MsD2(1,1).real()),
            ToComplex2(MsD2(1,2).real(),MsD2(1,2).imag())/sqrt(MsD2(1,1).real()*MsD2(2,2).real()),
            ToComplex2(MsD2(0,2).real(),MsD2(0,2).imag())/sqrt(MsD2(0,0).real()*MsD2(2,2).real()));
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHSetNMFV error " << err;
        return (false);
        //throw std::runtime_error(ss.str());
     }
    
    return (true);
}

bool SUSY::CalcHiggsSpectrum(void){
    int err;
    Complex SAeff;
    Complex UHiggs[3][3];
    Complex ZHiggs[3][3];
    //FHSetDebug(2);
    FHHiggsCorr(&err, mh, &SAeff, UHiggs, ZHiggs);  
    saeff = complex(SAeff.real(),SAeff.imag()); 
    
    // test - lines
//    int i,j;
//    for( i= 0; i< 3; i++){
//        for(j =0 ; j < 3 ;j++){
//            
//            std::cout << "UHiggs = " << UHiggs[i][j] << "i,j = " << i << " " << j << std::endl;
//        }
//    }
    
    // end - test
    
    
//    std::cout << "mh[0] = mh = " << mh[0] << std::endl;
//    std::cout << "mh[1] = mH = " << mh[1] << std::endl;
//    std::cout << "mh[2] = mA = " << mh[2] << std::endl;
//    std::cout << "mh[3] = mH+ =" << mh[3] << std::endl;
    
    for(int i = 0; i < 4; i++){
        if(std::isnan(mh[i])){
            std::cout << "FeynHiggs FHCorr error " << std::endl;
            return (false);
        }
    }

    mHp = mh[3];

    if ((mh[0] < 10.) || (mh[0] > 1000.)) {

        std::cout << "Evento scartato Mh valore sballato = " << mh[0] << std::endl;

        return (false);
    }
    
            
//    for(int i=0; i<3; i++)
//        for(int j=0; j<3; j++){
//            UH.assign(i,j,complex(UHiggs[i][j].real(),UHiggs[i][j].imag()));
//            ZH.assign(i,j,complex(ZHiggs[i][j].real(),ZHiggs[i][j].imag()));            
//        }
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHHiggsCorr error " << err;
        std::cout << "FeynHiggs FHCorr error " << std::endl;
        return (false);
        //throw std::runtime_error(ss.str());
     }
    
    return (true);
}

void SUSY::CalcHiggsCouplings(void){
    int err;
    Complex couplings[ncouplings];
    Complex couplingsms[ncouplingsms];    
    double gammas[ngammas];
    double gammasms[ngammasms];
    //FHSetDebug(2);
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

bool SUSY::CalcConstraints(){
    int err, ccb;
    //FHSetDebug(2);
    FHConstraints(&err, &FHgm2, &FHdeltarho, &FHMWMSSM, &FHMWSM, &FHSW2MSSM, 
            &FHSW2SM, &FHedmeTh, &FHedmn, &FHedmHg, &ccb);
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHConstraints error " << err;
        return (false);
        //throw std::runtime_error(ss.str());
     }   
    if (ccb != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHConstraints: colour breaking minimum " << err;
        //std::cout << "Error in Constraints.f in line: " << ccb << std::endl;
        return (false);       
    }
    return (true);
}

bool SUSY::CalcFlavour(){
    int err;
    FHFlavour(&err, &FHbsgMSSM, &FHbsgSM, &FHdeltaMsMSSM, &FHdeltaMsSM,
            &FHbsmumuMSSM, &FHbsmumuSM);
     if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs FHFlavour error " << err;
        return (false);
        //throw std::runtime_error(ss.str());
     }   
    /**** test -lines ****/
    
//      std::cout << " Delta MB_s MSSM FeynHiggs = " << " " << FHdeltaMsMSSM << std::endl; 
//      std::cout << " Delta MB_s SM FeynHiggs = " << " " << FHdeltaMsSM << std::endl;
//      std::cout << " " << std::endl;
    
    /**** end - test ***/
    return (true);
}

bool SUSY::CalcSpectrum(){
    int err, nmfv;
    double MASf[4][6], MCha[2], MNeu[4];
    Complex UASf[4][6][6], UCha[2][2], VCha[2][2], ZNeu[4][4], Deltab;
    FHGetPara(&err, &nmfv, MASf, UASf, MCha, UCha, VCha, MNeu, ZNeu, &Deltab, 
            &FHMGl, FHMHtree, &FHSAtree);
//            FHMHtree, &FHSAtree);
    for(int i = 0; i < 6; i++){
        
        Msn2(i) = MASf[0][i]*MASf[0][i];
        Msl2(i) = MASf[1][i]*MASf[1][i];
        Msu2(i) = MASf[2][i]*MASf[2][i];
        Msd2(i) = MASf[3][i]*MASf[3][i];
        
//        std::cout << "i : " << i << "Msu = " << MASf[2][i] << std::endl;
        
        for(int j = 0; j < 6; j++){
            Rn.assign(i,j,complex(UASf[0][i][j].real(),UASf[0][i][j].imag()));
            Rl.assign(i,j,complex(UASf[1][i][j].real(),UASf[1][i][j].imag()));
            Ru.assign(i,j,complex(UASf[2][i][j].real(),UASf[2][i][j].imag()));
            Rd.assign(i,j,complex(UASf[3][i][j].real(),UASf[3][i][j].imag()));
            
//          std::cout << "i, j " << i << j << "Msu = " << UASf[2][i][j] << std::endl;
//            
            
//            
//            if(std::isnan(complex(UASf[0][i][j].real()))) return (false);
//            if(std::isnan(complex(UASf[0][i][j].imm()))) return (false);
            
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
        std::cout << "FeynHiggs FHGetPara error " << std::endl;
        return (false);
        //throw std::runtime_error(ss.str());
     }   
       
    return (true);
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
    return (80.385);
}

double SUSY::cW2() const {

    std::cout << "Write codes for SUSY::cW2() " << std::endl;
    return (Mw() / 91.1876);
}

double SUSY::sW2() const {


    std::cout << "Write codes for SUSY::sW2() " << std::endl;
    return (1 - cW2());

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





