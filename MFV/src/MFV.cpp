/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MFV.h"
#include <math.h>

    /**
     * @brief MFV constructor
     * @param mQtilde_i squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     * @param mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     * @param mDtilde_i right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     * @param Au_i universal trilinear up-type squark coupling
     * @param Ad_i universal trilinear down-type squark coupling
     * @param mLtilde_i slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     * @param mEtilde_i right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     * @param mNtilde_i right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     * @param Ae_i universal trilinear charged slepton coupling
     * @param An_i universal trilinear sneutrino coupling
     * @param m1_i bino soft DRbar mass @f$m_{1}(m_{1})@f$
     * @param m2_i wino soft DRbar mass @f$m_{2}(m_{2})@f$
     * @param m3_i gluino soft DRbar mass @f$m_{3}(m_{3})@f$
     * @param muH_i superpotential @f$\mu@f$ term
     * @param tanb_i @f$\tan \beta @f$
     * @param mHp_i charged Higgs mass @f$m_{H^+}@f$
     * @param VCKM_i The CKM matrix
     * @param mu_i up quark mass at 2 GeV
     * @param md_i down quark mass at 2 GeV
     * @param mc_i charm quark mass mc(mc)
     * @param ms_i strange quark mass at 2 GeV
     * @param mt_i top quark mass mt(mt)
     * @param mb_i bottom quark mass mb(mb)
     * @param UPMNS_i The PMNS matrix
     * @param me_i electron mass
     * @param mmu_i muon mass
     * @param mtau_i tau mass
     * @param mnu1_i lightest neutrino mass
     * @param mnu2_i middle neutrino mass
     * @param mnu3_i hevier neutrino mass
     */

const std::string MFV::MFVvars[NMFVvars] = {"a1", "a2", "a3", "a4r", "a4i", 
        "a5r", "a5i", "a6", "a7", "a8r", "a8i", "x1", "x2", "y1", "y2r", "y2i", 
        "y3", "y4r", "y4i", "y5r", "y5i", "y6", "y7", "w1r", "w1i", "w2r", 
        "w2i", "w3r", "w3i", "w4r", "w4i", "w5r", "w5i"};

MFV::MFV() :
        SUSY(), X(), FHPost_cache(51,0.), FHPre_cache(51,0.)
{
    NCache = 0;
    Npassi = 0;
}


bool MFV::InitializeModel(){
    
    SetModelInitialized(SUSY::InitializeModel());
    return (IsModelInitialized());
    
}


bool MFV::PreUpdate(){
    
    if(!SUSY::PreUpdate())  return (false);
    
     return (true);
    
}

bool MFV::PostUpdate(){
   
    /** Cache **/
    
//    for (int i = 0; i < 6; i++) FHPost_cache(i) = quarks[i].getMass();
//    FHPost_cache(6) = lambda;
//    FHPost_cache(7) = A;
//    FHPost_cache(8) = rhob;
//    FHPost_cache(9) = etab;
//    FHPost_cache(10) = leptons[TAU].getMass();
//    FHPost_cache(11) = mut;
//    FHPost_cache(12) = mtpole;
//    FHPost_cache(13) = tanb;
//    FHPost_cache(14) = mHptree;
//    FHPost_cache(15) = muH.real();
//    FHPost_cache(16) = muH.imag();
//    FHPost_cache(17) = a1;
//    FHPost_cache(18) = a2;
//    FHPost_cache(19) = a3;
//    FHPost_cache(20) = a4.real();
//    FHPost_cache(21) = a4.imag();
//    FHPost_cache(22) = a5.real();
//    FHPost_cache(23) = a5.imag();
//    FHPost_cache(24) = a6;
//    FHPost_cache(25) = a7;
//    FHPost_cache(26) = a8.real();
//    FHPost_cache(27) = a8.imag();
//    FHPost_cache(28) = x1;
//    FHPost_cache(29) = x2;
//    FHPost_cache(30) = y1;
//    FHPost_cache(31) = y2.real();
//    FHPost_cache(32) = y2.imag();
//    FHPost_cache(33) = y3;
//    FHPost_cache(34) = y4.real();
//    FHPost_cache(35) = y4.imag();
//    FHPost_cache(36) = y5.real();
//    FHPost_cache(37) = y5.imag();
//    FHPost_cache(38) = y6;
//    FHPost_cache(39) = y7;
//    FHPost_cache(40) = w1.real();
//    FHPost_cache(41) = w1.imag();
//    FHPost_cache(42) = w2.real();
//    FHPost_cache(43) = w2.imag();
//    FHPost_cache(44) = w3.real();
//    FHPost_cache(45) = w3.imag();
//    FHPost_cache(46) = w4.real();
//    FHPost_cache(47) = w4.imag();
//    FHPost_cache(48) = w5.real();
//    FHPost_cache(49) = w5.imag();
//    FHPost_cache(50) = Q;
    /** end - cache **/

    
    
//    std::cout << "Post cache = " << FHPost_cache  << std::endl;
//    std::cout << "Pre cache = " << FHPre_cache  << std::endl;
//    std::cout << "Delta cache = " << FHPost_cache - FHPre_cache << std::endl;
//    if((FHPost_cache - FHPre_cache).mod() == 0){
//        
//        NCache++;
//        std::cout << "NusingCache = " << NCache << std::endl;
//        
//    }
//    
//    Npassi++;
//    std::cout << "NPassi = " << Npassi << std::endl;
    
    
    if (!StandardModel::PostUpdate()) return (false);
    
    //if((FHPost_cache - FHPre_cache).mod() != 0){
    SetSoftTerms();
    if(!SetFeynHiggsPars())  return (false);   //USO SPECTRUM FATTO DA ME
    if(!CalcHiggsSpectrum()) return (false);   //USO SPECTRUM FATTO DA ME
    
    
    //CalcHiggsCouplings();  
    //CalcHiggsProd(7.);     //cross sections at 7 TeV
    //if(!CalcConstraints()) return (false);
    //if(!CalcFlavour()) return (false);
    if(!CalcSpectrum()) return (false);
    
    
    
    if(!SUSY::PostUpdate())  return (false);
   // }
    
    /** Cache **/
    
//    for (int i = 0; i < 6; i++) FHPre_cache(i) = quarks[i].getMass();
//    FHPre_cache(6) = lambda;
//    FHPre_cache(7) = A;
//    FHPre_cache(8) = rhob;
//    FHPre_cache(9) = etab;
//    FHPre_cache(10) = leptons[TAU].getMass();
//    FHPre_cache(11) = mut;
//    FHPre_cache(12) = mtpole;
//    FHPre_cache(13) = tanb;
//    FHPre_cache(14) = mHptree;
//    FHPre_cache(15) = muH.real();
//    FHPre_cache(16) = muH.imag();
//    FHPre_cache(17) = a1;
//    FHPre_cache(18) = a2;
//    FHPre_cache(19) = a3;
//    FHPre_cache(20) = a4.real();
//    FHPre_cache(21) = a4.imag();
//    FHPre_cache(22) = a5.real();
//    FHPre_cache(23) = a5.imag();
//    FHPre_cache(24) = a6;
//    FHPre_cache(25) = a7;
//    FHPre_cache(26) = a8.real();
//    FHPre_cache(27) = a8.imag();
//    FHPre_cache(28) = x1;
//    FHPre_cache(29) = x2;
//    FHPre_cache(30) = y1;
//    FHPre_cache(31) = y2.real();
//    FHPre_cache(32) = y2.imag();
//    FHPre_cache(33) = y3;
//    FHPre_cache(34) = y4.real();
//    FHPre_cache(35) = y4.imag();
//    FHPre_cache(36) = y5.real();
//    FHPre_cache(37) = y5.imag();
//    FHPre_cache(38) = y6;
//    FHPre_cache(39) = y7;
//    FHPre_cache(40) = w1.real();
//    FHPre_cache(41) = w1.imag();
//    FHPre_cache(42) = w2.real();
//    FHPre_cache(43) = w2.imag();
//    FHPre_cache(44) = w3.real();
//    FHPre_cache(45) = w3.imag();
//    FHPre_cache(46) = w4.real();
//    FHPre_cache(47) = w4.imag();
//    FHPre_cache(48) = w5.real();
//    FHPre_cache(49) = w5.imag();
//    FHPre_cache(50) = Q;
  
    /** end - cache **/
    
    
    return (true);
}




bool MFV::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate())  return (false);
    
    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate())  return (false);
    
     return (true);
    
}

void MFV::SetParameter(const std::string name, const double& value) {
    if(name.compare("a1") == 0)
        a1 = value;
    else if(name.compare("a2") == 0)
        a2 = value;
    else if(name.compare("a3") == 0)
        a3 = value;
    else if(name.compare("a4r") == 0)
        a4.real() = value;
    else if(name.compare("a4i") == 0)
        a4.imag() = value;
    else if(name.compare("a5r") == 0)
        a5.real() = value;
    else if(name.compare("a5i") == 0)
        a5.imag() = value;
    else if(name.compare("a6") == 0)
        a6 = value;
    else if(name.compare("a7") == 0)
        a7 = value;
    else if(name.compare("a8r") == 0)
        a8.real() = value;
    else if(name.compare("a8i") == 0)
        a8.imag() = value;
    else if(name.compare("x1") == 0)
        x1 = value;
    else if(name.compare("x2") == 0)
        x2 = value;
    else if(name.compare("y1") == 0)
        y1 = value;
    else if(name.compare("y2r") == 0)
        y2.real() = value;
    else if(name.compare("y2i") == 0)
        y2.imag() = value;
    else if(name.compare("y3") == 0)
        y3 = value;
    else if(name.compare("y4r") == 0)
        y4.real() = value;
    else if(name.compare("y4i") == 0)
        y4.imag() = value;
    else if(name.compare("y5r") == 0)
        y5.real() = value;
    else if(name.compare("y5i") == 0)
        y5.imag() = value;
    else if(name.compare("y6") == 0)
        y6 = value;
    else if(name.compare("y7") == 0)
        y7 = value;
    else if(name.compare("w1r") == 0)
        w1.real() = value;
    else if(name.compare("w1i") == 0)
        w1.imag() = value;
    else if(name.compare("w2r") == 0)
        w2.real() = value;
    else if(name.compare("w2i") == 0)
        w2.imag() = value;
    else if(name.compare("w3r") == 0)
        w3.real() = value;
    else if(name.compare("w3i") == 0)
        w3.imag() = value;
    else if(name.compare("w4r") == 0)
        w4.real() = value;
    else if(name.compare("w4i") == 0)
        w4.imag() = value;
    else if(name.compare("w5r") == 0)
        w5.real() = value;
    else if(name.compare("w5i") == 0)
        w5.imag() = value;
    else
        SUSY::SetParameter(name, value);
}

bool MFV::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NMFVvars; i++) {
        if (DPars.find(MFVvars[i]) == DPars.end()) {
            std::cout << "missing mandatory MFV parameter " << MFVvars[i] << std::endl;
            return false;
        }
    }
    return(SUSY::CheckParameters(DPars));
}

bool MFV::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return (CheckParameters(DPars));
}

void MFV::SetSoftTerms(void){
    
    /// test lines
    
//    std::cout << "a1 = " << a1 << std::endl;
//    std::cout << "a2 = " << a2 << std::endl;
//    std::cout << "a3 = " << a3 << std::endl;
//    std::cout << "a6 = " << a6 << std::endl;
//    std::cout << "a7 = " << a7 << std::endl;
//    std::cout << "m1 = " << m1 << std::endl;
//    std::cout << "m2 = " << m2 << std::endl;
//    std::cout << "m3 = " << m3 << std::endl;
//    std::cout << "muH = " << muH << std::endl;
//    std::cout << "mHp = " << mHp << std::endl;
//    std::cout << "tanb = " << tanb << std::endl;
//    std::cout << "Q = " << Q << std::endl;
//    std::cout << "mup = " << quarks[UP].getMass() << std::endl;
//    std::cout << "mdown = " << quarks[DOWN].getMass() << std::endl;
//    std::cout << "mstrange = " << quarks[STRANGE].getMass() << std::endl;
//    std::cout << "mcharm = " << quarks[CHARM].getMass() << std::endl;
//    std::cout << "mtop = " << quarks[TOP].getMass() << std::endl;
//    std::cout << "mbottom = " << quarks[BOTTOM].getMass() << std::endl;
//    std::cout << "AlsMz = " << AlsMz << std::endl;
    
    
    
    //// end-test
    
    
    
    // Colangelo's expression in Colangelo's basis
    X.Update(myCKM);
    MsQ2 = matrix<complex>::Id(3) * a1 + X.GetX13() * x1 + X.GetX1() * y1 + 
            X.GetX5() * y2 + X.GetX9() * y2.conjugate();
    MsU2 = matrix<complex>::Id(3) * a2 + X.GetX1() * x2;
    MsD2 = matrix<complex>::Id(3) * a3 + X.GetX1() * y3 + X.GetX3() * w1 + 
            X.GetX4() * w1.conjugate();
    TU = (X.GetX5() * a4 + X.GetX1() * y4 + 
            X.GetX6() * w2).transpose();
    TD = (X.GetX1() * a5 + X.GetX5() * y5 + 
            X.GetX2() * w3 + X.GetX4() * w4).transpose();
    MsL2 = matrix<complex>::Id(3) * a6 + X.GetX1() * y6;
    MsE2 = matrix<complex>::Id(3) * a7 + X.GetX1() * y7;
    TE = (X.GetX1() * a8 + X.GetX2() * w5).transpose();
    
    //rotation to the SCKM basis according to SLHA notation
    
    matrix<complex> ckm(3,3,0.);
    myCKM.getCKM(ckm);
  
    TU = sqrt(2.) * TU * ckm.hconjugate();  
    TD = sqrt(2.) * TD;
    TE = sqrt(2.) * TE;
    
//    std::cout << "TU = " << TU << std::endl;
//    std::cout << "TD = " << TD << std::endl;
    
    
    // test - lines
    
    //std::cout << "X1 = " << X.GetX1() << std::endl;
    
//    std::cout << "TU = " << TU << std::endl;
//    std::cout << "TD = " << TD << std::endl;
//    std::cout << "TE = " << TE << std::endl;
//    std::cout << "MsU2 = " << MsU2 << std::endl;
//    std::cout << "MsD2 = " << MsU2 << std::endl;
//    std::cout << "MsL2 = " << MsU2 << std::endl;
//    std::cout << "MsE2 = " << MsU2 << std::endl;
}

