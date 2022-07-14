//AG: The following functions were defined:
// virtual double deltaGV_f_2(const Particle f) const;
// virtual double deltaGA_f_2(const Particle f) const;
// virtual double deltaGamma_Zf_Test(const Particle f) const;
// virtual double deltaGamma_Zf_2(const Particle f) const;
// virtual double deltaGamma_Z_Test() const;
// virtual double deltaGamma_Z_2() const;
// virtual double deltaGamma_Zhad_Test() const;
// virtual double deltaGamma_Zhad_2() const;
// virtual double deltaSigmaHadron_Test() const;
// virtual double deltaSigmaHadron_2() const;
// virtual double deltaSin2thetaEff_e_2() const;
// virtual double deltaSin2thetaEff_mu_2() const;
// virtual double deltaA_f_2(const Particle f) const;
// virtual double deltaAFB_Test(const Particle f) const;
// virtual double deltaAFB_2(const Particle f) const;
// virtual double deltaR0_f_Test(const Particle f) const;
// virtual double deltaR0_f_2(const Particle f) const;

/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPbase.h"

NPbase::NPbase()
: StandardModel(), trueSM()
{
    trueSM.InitializeModel();
    trueSM.setSliced(true);
}

bool NPbase::Update(const std::map<std::string, double>& DPars)
{
    if (!trueSM.Update(DPars)) return (false);
    return StandardModel::Update(DPars);
}

bool NPbase::PostUpdate()
{
    return StandardModel::PostUpdate();
}

//AG:begin
bool NPbase::NumericCheck() const
{
    return false;
}
int NPbase::Output() const {
    // 0 Linear
    // 1 Linear_Test
    // 2 Linear + Quadratic
    // 3 Linear_Test + Quadratic
    // 4 Quadratic
    // 5 SM
    return 3;
}
//AG:end

double NPbase::alphaMz() const
{
    /*AG:begin
    // FlagMWinput flag to be implemented in NPbase, especially for alphaMz and Mw!
    double dalphaMz = 0.0;
    if(FlagMWinput){
        double alpha = trueSM.alphaMz();
        double c2 = trueSM.cW2();
        double s2 = trueSM.sW2();
        
        dalphaMz = sqrt(2.0)*GF*Mz*Mz/4.0/M_PI*c2*( alpha/2.0*(obliqueS() - 2.0*c2*obliqueT() - (c2-s2)/2.0/s2*obliqueU()) - s2/2.0/(c2-s2)*DeltaGF() );
    }
    return (trueSM.alphaMz() + dalphaMz);
    //AG:end */
    
    double myAlphaMz = trueSM.alphaMz();
    
    return myAlphaMz;
}

double NPbase::Mw() const
{
    double myMw = trueSM.Mw();

    double alpha = trueSM.alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    myMw *= 1.0 - alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT() - (c2 - s2) * obliqueU() / 2.0 / s2)
            - s2 / 2.0 / (c2 - s2) * DeltaGF();

    //std::cout << "Mw: c_S=" << - alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Mw: c_T=" << - alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Mw: c_U=" << - alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return myMw;
}

double NPbase::GammaW(const Particle fi, const Particle fj) const
{
    double Gamma_Wij = trueSM.GammaW(fi, fj);

    double alpha = trueSM.alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    Gamma_Wij *= 1.0 - 3.0 * alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT() - (c2 - s2) * obliqueU() / 2.0 / s2)
            - (1.0 + c2) / 2.0 / (c2 - s2) * DeltaGF();

    //std::cout << "Gw: c_S=" << - 3.0*alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Gw: c_T=" << - 3.0*alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Gw: c_U=" << - 3.0*alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return Gamma_Wij;
}

double NPbase::GammaW() const
{
    double Gamma_W = trueSM.GammaW();

    double alpha = trueSM.alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    Gamma_W *= 1.0 - 3.0 * alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT() - (c2 - s2) * obliqueU() / 2.0 / s2)
            - (1.0 + c2) / 2.0 / (c2 - s2) * DeltaGF();

    //std::cout << "Gw: c_S=" << - 3.0*alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Gw: c_T=" << - 3.0*alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Gw: c_U=" << - 3.0*alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return Gamma_W;
}

double NPbase::BrW(const Particle fi, const Particle fj) const
{
    double GammW = GammaW();
    double GammWij = GammaW(fi, fj);

    return GammWij/GammW;
}


double NPbase::RWlilj(const Particle li, const Particle lj) const
{
    double GammWli, GammWlj;
    
    if (li.is("ELECTRON"))
        GammWli = GammaW(leptons[NEUTRINO_1],li);
    else if (li.is("MU"))
        GammWli = GammaW(leptons[NEUTRINO_2],li);        
    else if (li.is("TAU"))
        GammWli = GammaW(leptons[NEUTRINO_3],li);        
    else
        throw std::runtime_error("Error in NPbase::RWlilj. li must be a charged lepton");
    
    if (lj.is("ELECTRON"))
        GammWlj = GammaW(leptons[NEUTRINO_1],lj);
    else if (lj.is("MU"))
        GammWlj = GammaW(leptons[NEUTRINO_2],lj);        
    else if (lj.is("TAU"))
        GammWlj = GammaW(leptons[NEUTRINO_3],lj);        
    else
        throw std::runtime_error("Error in NPbase::RWlilj. lj must be a charged lepton");
    
    return GammWli/GammWlj;
}

double NPbase::RWc() const
{  
    double GammWcX, GammWhad;
    
//  Basic NPBase implementation is only NP universal and respects CKM unitarity.
//  As it directly uses SM implementation of W width, use same definition here
//  (with the modified NPBase implementation of GammaW)

//  Add all the  W-> cX decays
//  In SM GammaW fermion masses are ignored and CKM=1 but uses that SM CKM is unitary => I only need W->cs
    GammWcX = GammaW(quarks[CHARM], quarks[STRANGE]);
    
//  For the same reasons, I only need to add the W-> ud decays into the hadronic part
    GammWhad = GammWcX
            + GammaW(quarks[UP], quarks[DOWN]);

    return GammWcX/GammWhad;
}

double NPbase::deltaGV_f(const Particle f) const
{
    if (f.is("TOP")) return 0.;

    /* SM values */
    double alpha = trueSM.alphaMz();
    double sW2SM = trueSM.sW2();
    double cW2SM = trueSM.cW2();
    double gVSM = trueSM.gV_f(f).real();
    double gASM = trueSM.gA_f(f).real();

    return ( gVSM * (alpha * obliqueT() - DeltaGF()) / 2.0
            + (gVSM - gASM) / 4.0 / sW2SM / (cW2SM - sW2SM)
            *(alpha * (obliqueS() - 4.0 * cW2SM * sW2SM * obliqueT())
            + 4.0 * cW2SM * sW2SM * DeltaGF()));
}

gslpp::complex NPbase::gV_f(const Particle f) const
{
    //AG: deltaGV_f_2(f) added below.
    return ( trueSM.gV_f(f) + deltaGV_f(f) + deltaGV_f_2(f) );
}

double NPbase::deltaGA_f(const Particle f) const
{
    if (f.is("TOP")) return 0.;
    /* SM values */
    double alpha = trueSM.alphaMz();
    double gASM = trueSM.gA_f(f).real();

    return ( gASM * (alpha * obliqueT() - DeltaGF()) / 2.0);
}

gslpp::complex NPbase::gA_f(const Particle f) const
{
    //AG: deltaGA_f_2(f) added below
    return ( trueSM.gA_f(f) + deltaGA_f(f) + deltaGA_f_2(f) );
}

gslpp::complex NPbase::rhoZ_f(const Particle f) const
{
    return ( gA_f(f) * gA_f(f) / f.getIsospin() / f.getIsospin());

}

gslpp::complex NPbase::kappaZ_f(const Particle f) const
{
    return ( (1.0 - gV_f(f) / gA_f(f)) / (4.0 * fabs(f.getCharge()) * sW2()));
}

////////////////////////////////////////////////////////////////////////
//AG:begin

double NPbase::deltaGamma_Zf_Test(const Particle f) const   //c_alpha
{    
    double DeltaGamma_Zf_Test = 0.0;
    
    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f); 
    
    if (delGVf != 0.0 || delGAf != 0.0) {
        double Nf;
        if (f.is("LEPTON")) {
            Nf = 1.0;
        } else {
            Nf = 3.0;
        }

        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real(); 
        double DelGammaZfTest;
        
        DelGammaZfTest = 2.0 * Nf * (gVf * delGVf + gAf * delGAf);
        //DeltaGamma_Zf_Test = 4.0 * GF * pow(Mz,3.0) / 24.0 / sqrt(2.0) / M_PI * DelGammaZfTest;
        DeltaGamma_Zf_Test = alphaMz()*Mz / 12.0 / trueSM.sW2()/trueSM.cW2() * DelGammaZfTest;
        
        // As in deltaGamma_Zf(f):
        //delGammaZf = 2.0 * Nf * (gVf * delGVf + gAf * delGAf);
        //deltaGamma_Zf = alphaMz() * Mz / 12.0 / sW2_SM / cW2_SM * delGammaZf;
    }

    return DeltaGamma_Zf_Test;
}
double NPbase::deltaGamma_Zf_2(const Particle f) const        //c_alpha
{    
    double DeltaGamma_Zf_2=0.0;
    
    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f);  
    double delGVf2 = deltaGV_f_2(f);
    double delGAf2 = deltaGA_f_2(f);
    
    bool nonZeroNP = false;
    if (delGVf2!=0.0 || delGAf2!=0.0) {nonZeroNP = true;}
    
    if (nonZeroNP) {
        double Nf;
        if (f.is("LEPTON")) {
            Nf = 1.0;
        } else {
            Nf = 3.0;
        }

        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real(); 
        double DelGammaZf2;
        
        DelGammaZf2 = Nf * ( 2.0*(gVf*delGVf2 + gAf*delGAf2) + delGVf*delGVf + delGAf*delGAf );
        //DeltaGamma_Zf_2 = 4.0 * GF * pow(Mz,3.0) / 24.0 / pow(2,0.5) / M_PI * DelGammaZf2;
        DeltaGamma_Zf_2 = alphaMz()*Mz / 12.0 / trueSM.sW2()/trueSM.cW2() * DelGammaZf2;
    }

    return DeltaGamma_Zf_2;
}

double NPbase::deltaGamma_Z_Test() const
{
    double deltaGamma_Z_Test = 0.0;
    
    bool nonZeroNP = false;
    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p] != 0.0 || delGAl[p] != 0.0 || delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        for (int p=0; p<6; p++){
            deltaGamma_Z_Test += deltaGamma_Zf_Test(leptons[p]) + deltaGamma_Zf_Test(quarks[p]);
        }
    }

    return deltaGamma_Z_Test;
}

double NPbase::deltaGamma_Z_2() const
{
    double deltaGamma_Z_2 = 0.;
    
    bool nonZeroNP = false;
    double delGVl2[6], delGAl2[6], delGVq2[6], delGAq2[6];
    for (int p = 0; p < 6; ++p) {
        delGVl2[p] = deltaGV_f_2(leptons[p]);
        delGAl2[p] = deltaGA_f_2(leptons[p]);
        delGVq2[p] = deltaGV_f_2(quarks[p]);
        delGAq2[p] = deltaGA_f_2(quarks[p]);
        if (delGVq2[p]!=0.0 or delGAq2[p]!=0.0 or delGVl2[p]!=0.0 or delGAl2[p]!=0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        for(int p=0; p<6; p++){
            deltaGamma_Z_2 += deltaGamma_Zf_2(leptons[p]) + deltaGamma_Zf_2(quarks[p]);
        }
    }

    return deltaGamma_Z_2;
}

double NPbase::deltaGamma_Zhad_Test() const
{
    double DeltaGamma_Zhad_Test = 0.;
    
    bool nonZeroNP = false;
    double delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        for(int p=0; p<6; p++){
            DeltaGamma_Zhad_Test += deltaGamma_Zf_Test(quarks[p]);
        }
    }
    
    return DeltaGamma_Zhad_Test;
}

double NPbase::deltaGamma_Zhad_2() const
{
    double DeltaGamma_Zhad_2 = 0.;
    bool nonZeroNP = false;
    double delGVq2[6], delGAq2[6];
    for (int p = 0; p < 6; ++p) {
        delGVq2[p] = deltaGV_f_2(quarks[p]);
        delGAq2[p] = deltaGA_f_2(quarks[p]);
        if (delGVq2[p] != 0.0 || delGAq2[p] != 0.0) {nonZeroNP = true;}
    }
    
    if (nonZeroNP) {
        for(int p=0; p<6; p++){
            DeltaGamma_Zhad_2 += deltaGamma_Zf_2(quarks[p]);
        }
    }
    
    return DeltaGamma_Zhad_2;
}
/*
double NPbase::deltaSigmaHadron_Test() const
{
    double sigma_had_Test = 0.0;
    
    bool nonZeroNP = false;
    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p]!=0.0 || delGAl[p]!=0.0 || delGVq[p]!=0.0 || delGAq[p]!=0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double prefactor = 4.0 * GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI; // It will cancel in the ratio either way
        
        // SM contributions (without higher order SM corrections)
        double Gamma_e_SMtree = 1.0 * prefactor * ( pow(trueSM.gV_f(leptons[ELECTRON]).real(),2.0) + pow(trueSM.gA_f(leptons[ELECTRON]).real(),2.0) );
        
        double Gamma_had_SMtree = 0.0;
        for (int q = 0; q < 6; ++q) {
            if (quarks[q].getName()!="TOP") {
                Gamma_had_SMtree += 3.0 * prefactor * ( pow(trueSM.gV_f(quarks[q]).real(),2.0) + pow(trueSM.gA_f(quarks[q]).real(),2.0) );
            }
        };
        
        double Gamma_lep_SMtree = 0.0;
        for (int l = 0; l < 6; ++l) {
            //std::cout << "l=" << l << ", leptons[l] = " << leptons[l].getName() << std::endl;
            Gamma_lep_SMtree += 1.0 * prefactor * ( pow(trueSM.gV_f(leptons[l]).real(),2.0) + pow(trueSM.gA_f(leptons[l]).real(),2.0) );
        };
        
        double Gamma_Z_SMtree = Gamma_had_SMtree + Gamma_lep_SMtree;
        
        // Linear contributions: (TO BE CHANGED TO NON-TEST AFTER CHECKING PREFACTOR) 
        double dGamma_e = deltaGamma_Zf_Test(leptons[ELECTRON]);        
        double dGamma_had = deltaGamma_Zhad_Test();
        double dGamma_Z = deltaGamma_Z_Test();
       
        // Then,
        sigma_had_Test = 12.0*M_PI/pow(Mz,2.0) * Gamma_e_SMtree*Gamma_had_SMtree/pow(Gamma_Z_SMtree,2.0) * (
                dGamma_e/Gamma_e_SMtree 
                + dGamma_had/Gamma_had_SMtree 
                - 2.0*dGamma_Z/Gamma_Z_SMtree
                );
        //std::cout << dGamma_e << '\t' << Gamma_e_SMtree << '\t' << dGamma_e/Gamma_e_SMtree << std::endl;
        //std::cout << dGamma_had << '\t' << Gamma_had_SMtree << '\t' << dGamma_had/Gamma_had_SMtree << std::endl;
        //std::cout << dGamma_Z << '\t' << Gamma_Z_SMtree << '\t' << dGamma_Z/Gamma_Z_SMtree << std::endl;
        //std::cout << " " << std::endl;
    }

    return sigma_had_Test;
}
*/

double NPbase::deltaSigmaHadron_Test() const        //c_alpha
{
    double sigma_had_Test = 0.0;
    
    bool nonZeroNP = false;
    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p]!=0.0 || delGAl[p]!=0.0 || delGVq[p]!=0.0 || delGAq[p]!=0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        //double prefactor = 4.0 * GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI; // It will cancel in the ratio either way
        double prefactor = alphaMz()*Mz / 12.0 / trueSM.sW2()/trueSM.cW2();
        
        // SM contributions (without higher order SM corrections)
        double Gamma_e_SMtree = 1.0 * prefactor * ( pow(trueSM.gV_f(leptons[ELECTRON]).real(),2.0) + pow(trueSM.gA_f(leptons[ELECTRON]).real(),2.0) );
        
        double Gamma_lep_SMtree = 0.0;
        double Gamma_had_SMtree = 0.0;
        for (int p = 0; p < 6; ++p) {
            Gamma_lep_SMtree += 1.0 * prefactor * ( pow(trueSM.gV_f(leptons[p]).real(),2.0) + pow(trueSM.gA_f(leptons[p]).real(),2.0) );
            if (quarks[p].getName()!="TOP") {
                Gamma_had_SMtree += 3.0 * prefactor * ( pow(trueSM.gV_f(quarks[p]).real(),2.0) + pow(trueSM.gA_f(quarks[p]).real(),2.0) );
            }
        }
        double Gamma_Z_SMtree = Gamma_had_SMtree + Gamma_lep_SMtree;
        
        // Linear contributions: (TO BE CHANGED TO NON-TEST AFTER CHECKING PREFACTOR) 
        double dGamma_e = deltaGamma_Zf_Test(leptons[ELECTRON]);        
        double dGamma_had = deltaGamma_Zhad_Test();
        double dGamma_Z = deltaGamma_Z_Test();
       
        // Then,
        sigma_had_Test = 12.0*M_PI/pow(Mz,2.0) * Gamma_e_SMtree*Gamma_had_SMtree/pow(Gamma_Z_SMtree,2.0) * (
                dGamma_e/Gamma_e_SMtree 
                + dGamma_had/Gamma_had_SMtree 
                - 2.0*dGamma_Z/Gamma_Z_SMtree
                );   
        
    }
    return sigma_had_Test;
}

double NPbase::deltaSigmaHadron_2() const             //c_alpha
{
    double sigma_had_2 = 0.;

    bool nonZeroNP = false;
    double delGVl2[6], delGAl2[6], delGVq2[6], delGAq2[6];
    for (int p = 0; p < 6; ++p) {
        delGVl2[p] = deltaGV_f_2(leptons[p]);
        delGAl2[p] = deltaGA_f_2(leptons[p]);
        delGVq2[p] = deltaGV_f_2(quarks[p]);
        delGAq2[p] = deltaGA_f_2(quarks[p]);
        if (delGVl2[p]!=0.0 || delGAl2[p]!=0.0 || delGVq2[p]!=0.0 || delGAq2[p]!=0.0) 
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        //double prefactor = 4.0 * GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI; // It will cancel in the ratio either way
        double prefactor = alphaMz()*Mz / 12.0 / trueSM.sW2()/trueSM.cW2();
        
        // tree-level SM:
        double Gamma_e_SM = 1.0 * prefactor * ( pow(trueSM.gV_f(leptons[ELECTRON]).real(),2.0) + pow(trueSM.gA_f(leptons[ELECTRON]).real(),2.0) );
        
        double Gamma_lep_SM = 0.0;
        double Gamma_had_SM = 0.0;
        for (int p = 0; p < 6; ++p) {
            Gamma_lep_SM += 1.0 * prefactor * ( pow(trueSM.gV_f(leptons[p]).real(),2.0) + pow(trueSM.gA_f(leptons[p]).real(),2.0) );
            if (quarks[p].getName()!="TOP") {
                Gamma_had_SM += 3.0 * prefactor * ( pow(trueSM.gV_f(quarks[p]).real(),2.0) + pow(trueSM.gA_f(quarks[p]).real(),2.0) );
            }
        };      
        double Gamma_Z_SM = Gamma_had_SM + Gamma_lep_SM;
        
        double dGamma_e = deltaGamma_Zf(leptons[ELECTRON]);        
        double dGamma_had = deltaGamma_Zhad();
        double dGamma_Z = deltaGamma_Z();
        
        double dGamma_e_2 = deltaGamma_Zf_2(leptons[ELECTRON]);        
        double dGamma_had_2 = deltaGamma_Zhad_2();
        double dGamma_Z_2 = deltaGamma_Z_2();
       
        // Then,
        sigma_had_2 = 12.0*M_PI/pow(Mz,2.0) * Gamma_e_SM*Gamma_had_SM/pow(Gamma_Z_SM,2.0) * (
                  dGamma_e_2/Gamma_e_SM 
                + dGamma_had_2/Gamma_had_SM 
                - 2.0*dGamma_Z_2/Gamma_Z_SM
                + dGamma_e*dGamma_had/Gamma_e_SM/Gamma_had_SM
                - 2.0*dGamma_e*dGamma_Z/Gamma_e_SM/Gamma_Z_SM
                - 2.0*dGamma_had*dGamma_Z/Gamma_had_SM/Gamma_Z_SM
                + 3.0*pow(dGamma_Z,2.0)/pow(Gamma_Z_SM,2.0)
                );
    }

    return sigma_had_2;
}

double NPbase::deltaSin2thetaEff_e_2() const
{
    double sin2_theta_eff_2=0.0;
    double delGVf = deltaGV_f(leptons[ELECTRON]);
    double delGAf = deltaGA_f(leptons[ELECTRON]);
    double delGVf2 = deltaGV_f_2(leptons[ELECTRON]);
    double delGAf2 = deltaGA_f_2(leptons[ELECTRON]);
    
    bool nonZeroNP = false;
    if (delGVf2!=0.0 || delGAf2!=0.0) {nonZeroNP = true;}
    
    if (nonZeroNP) {
        double gVf = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAf = trueSM.gA_f(leptons[ELECTRON]).real();
        sin2_theta_eff_2 = 1.0/4.0 * delGVf * delGAf / pow(gAf,2.0)
                           - 1.0/4.0 * gVf * pow(delGAf,2.0) / pow(gAf,3.0)
                           - 1.0/4.0 * ( gAf*delGVf2 - gVf*delGAf2) / pow(gAf,2.0) ; 
    }
    
    return sin2_theta_eff_2;
}

double NPbase::deltaSin2thetaEff_mu_2() const
{
    double sin2_theta_eff_2=0.0;
    double delGVf = deltaGV_f(leptons[MU]);
    double delGAf = deltaGA_f(leptons[MU]);
    double delGVf2 = deltaGV_f_2(leptons[MU]);
    double delGAf2 = deltaGA_f_2(leptons[MU]);
    
    bool nonZeroNP = false;
    if (delGVf2!=0.0 || delGAf2!=0.0) {nonZeroNP = true;}
    
    if (nonZeroNP) {
        double gVf = trueSM.gV_f(leptons[MU]).real();
        double gAf = trueSM.gA_f(leptons[MU]).real();
        sin2_theta_eff_2 = 1.0/4.0 * delGVf * delGAf / pow(gAf,2.0)
                           - 1.0/4.0 * gVf * pow(delGAf,2.0) / pow(gAf,3.0)
                           - 1.0/4.0 * ( gAf*delGVf2 - gVf*delGAf2) / pow(gAf,2.0) ; 
    }
    
    return sin2_theta_eff_2;
}

double NPbase::deltaA_f_2(const Particle f) const
{     
    double dA_2 = 0.0;
   
    bool nonZeroNP = false;
    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f);
    double delGVf2 = deltaGV_f_2(f);
    double delGAf2 = deltaGA_f_2(f);
    if (delGVf2!=0.0 || delGAf2!=0.0) {nonZeroNP = true;}
    
    if (nonZeroNP) {
        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real();
        double Gf = gVf*gVf + gAf*gAf;

        double f2 = -2.0 * ( gVf*gVf - gAf*gAf ) * ( gAf*delGVf2 - gVf*delGAf2) / Gf / Gf;
        double f1 =  2.0 * ( gVf*gAf*( gAf*gAf - 3.0*gVf*gVf )*delGAf*delGAf 
                    + gVf*gAf*( gVf*gVf - 3.0*gAf*gAf )*delGVf*delGVf 
                    - ( pow(gAf,4.0) - 6.0*pow(gAf,2.0)*pow(gVf,2.0) + pow(gVf,4.0) )*delGVf*delGAf 
                    ) / pow(Gf,3.0);

        dA_2 = f1+f2;
    }
    return dA_2;
}

double NPbase::deltaAFB_Test(const Particle f) const
{
    double dAFB_Test=0.0;
    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f);
    double delGVe = deltaGV_f(leptons[ELECTRON]);
    double delGAe = deltaGA_f(leptons[ELECTRON]);
    
    if (delGVe != 0.0 || delGAe != 0.0 || delGVf != 0.0 || delGAf != 0.0) {
        double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real();
        
        double Ge = gVe*gVe + gAe*gAe;
        double delGVeOverGAe = (gVe*delGAe-gAe*delGVe) ;
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gVf*delGAf-gAf*delGVf) ;
        
        double Ae = 2.0*gVe*gAe/(gVe*gVe+gAe*gAe);
        double deltaAe = deltaA_f(leptons[ELECTRON]);
        double Af = 2.0*gVf*gAf/(gVf*gVf+gAf*gAf);
        double deltaAf = deltaA_f(f);
        
        dAFB_Test = 3.0/4.0 * ( Af*deltaAe + Ae*deltaAf ) ;
    }
    
    return dAFB_Test;
}

double NPbase::deltaAFB_2(const Particle f) const
{
    double dAFB_2=0.0;
    
    bool nonZeroNP = false;
    double delGVf2 = deltaGV_f_2(f);
    double delGAf2 = deltaGA_f_2(f);
    if (delGVf2!=0.0 || delGAf2!=0.0) {nonZeroNP = true;}
    
    if (nonZeroNP) {
        /*
        double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real();
        
        double Ge = gVe*gVe + gAe*gAe;
        double delGVeOverGAe = (gVe*delGAe-gAe*delGVe) ;
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gVf*delGAf-gAf*delGVf) ;
        
        double Ae = 2.0*gVe*gAe/(gVe*gVe+gAe*gAe);
        double deltaAe = deltaA_f(leptons[ELECTRON]);
        double Af = 2.0*gVf*gAf/(gVf*gVf+gAf*gAf);
        double deltaAf = deltaA_f(f);
        
        dAFB_Test = 3.0 * Af/2.0 * deltaAe/2.0 + 3.0 * Ae/2.0 * deltaAf/2.0 ;*/
        
        double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real();
        
        double Ae = 2.0*gVe*gAe/(gVe*gVe+gAe*gAe);
        double deltaAe = deltaA_f(leptons[ELECTRON]);
        double deltaAe2 = deltaA_f_2(leptons[ELECTRON]);
        double Af = 2.0*gVf*gAf/(gVf*gVf+gAf*gAf);
        double deltaAf = deltaA_f(f);
        double deltaAf2 = deltaA_f_2(f);
        
        if (f.is("ELECTRON")) 
            dAFB_2 = 3.0/4.0 * ( deltaAe*deltaAe + 2.0*Ae*deltaAe2 );
        else 
            dAFB_2 = 3.0/4.0 * ( deltaAe*deltaAf + Ae*deltaAf2 + Af*deltaAe2 );
    }
    return dAFB_2;
}

double NPbase::deltaR0_f_Test(const Particle f) const        //c_alpha
{
    double dR0_f_Test = 0.;
    double delGVl = 0.0, delGAl = 0.0;
    double delGVq[6], delGAq[6];
    bool nonZeroNP = false;
    if (f.is("LEPTON")) {
        delGVl = deltaGV_f(f);
        delGAl = deltaGA_f(f);
        if (delGVl != 0.0 || delGAl != 0.0) {nonZeroNP = true;}
    }
    for (int q = 0; q < 6; ++q) {
        delGVq[q] = deltaGV_f(quarks[q]);
        delGAq[q] = deltaGA_f(quarks[q]);
        if (delGVq[q] != 0.0 || delGAq[q] != 0.0) {nonZeroNP = true;}
    }

    if (nonZeroNP) {
        //double prefactor = 4.0 * GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI; // It will cancel in the ratio either way
        double prefactor = alphaMz()*Mz / 12.0 / trueSM.sW2()/trueSM.cW2();
        
        double Gamma_l_SMtree = 1.0 * prefactor * ( pow(trueSM.gV_f(f).real(),2.0) + pow(trueSM.gA_f(f).real(),2.0) );
        double Gamma_q_SMtree = 3.0 * prefactor * ( pow(trueSM.gV_f(f).real(),2.0) + pow(trueSM.gA_f(f).real(),2.0) );
        
        double Gamma_had_SMtree = 0.0;
        for (int q = 0; q < 6; ++q) {
            if (quarks[q].getName()!="TOP") {
                Gamma_had_SMtree += 3.0 * prefactor * ( pow(trueSM.gV_f(quarks[q]).real(),2.0) + pow(trueSM.gA_f(quarks[q]).real(),2.0) );
            }
        };
        
        double deltaGamma_f = deltaGamma_Zf_Test(f);        
        double deltaGamma_had = deltaGamma_Zhad_Test();
        
        if(f.is("ELECTRON") || f.is("MU") || f.is("TAU")){
            dR0_f_Test = ( Gamma_l_SMtree*deltaGamma_had - Gamma_had_SMtree*deltaGamma_f ) / pow(Gamma_l_SMtree,2.0);
        }
        if(f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3")){
            dR0_f_Test = ( Gamma_had_SMtree*deltaGamma_f - Gamma_l_SMtree*deltaGamma_had ) / pow(Gamma_had_SMtree, 2.0);
        }
        if(f.is("QUARK")){
            dR0_f_Test = ( Gamma_had_SMtree*deltaGamma_f - Gamma_q_SMtree*deltaGamma_had ) / pow(Gamma_had_SMtree, 2.0);
        }
    }
    
    return dR0_f_Test;
   
}
double NPbase::deltaR0_f_2(const Particle f) const           //c_alpha
{
    double dR0_f_2 = 0.;
    double delGVl2=0.0, delGAl2=0.0, delGVq2[6], delGAq2[6];
    bool nonZeroNP = false;
    if (f.is("LEPTON")) {
        delGAl2 = deltaGA_f_2(f);
        delGVl2 = deltaGV_f_2(f);
        if (delGVl2!=0.0 || delGAl2!=0.0) {nonZeroNP = true;}
    }
    for (int q = 0; q < 6; ++q) {
        delGVq2[q] = deltaGV_f_2(quarks[q]);
        delGAq2[q] = deltaGA_f_2(quarks[q]);
        if (delGVq2[q]!=0.0 || delGAq2[q]!=0.0) {nonZeroNP = true;}
    }
   
    if (nonZeroNP) {
        //double prefactor = 4.0 * GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI; // It will cancel in the ratio either way
        double prefactor = alphaMz()*Mz / 12.0 / trueSM.sW2()/trueSM.cW2();
        
        double Gamma_l_SMtree = 1.0 * prefactor * ( pow(trueSM.gV_f(f).real(),2.0) + pow(trueSM.gA_f(f).real(),2.0) );
        double Gamma_q_SMtree = 3.0 * prefactor * ( pow(trueSM.gV_f(f).real(),2.0) + pow(trueSM.gA_f(f).real(),2.0) );
        
        double Gamma_had_SMtree = 0.0;
        for (int q = 0; q < 6; ++q) {
            Gamma_had_SMtree += 3.0 * prefactor * ( pow(trueSM.gV_f(quarks[q]).real(),2.0) + pow(trueSM.gA_f(quarks[q]).real(),2.0) );
        };
        
        double deltaGamma_f = deltaGamma_Zf(f);        
        double deltaGamma_had = deltaGamma_Zhad();
        
        double deltaGamma_f_2 = deltaGamma_Zf_2(f);        
        double deltaGamma_had_2 = deltaGamma_Zhad_2();
        
        // There may be s mistake in the charged leptons case:
        if(f.is("ELECTRON") || f.is("MU") || f.is("TAU")){
            dR0_f_2 =  Gamma_had_SMtree*pow(deltaGamma_f,2.0) / pow(Gamma_l_SMtree,3.0)
                     - deltaGamma_had*deltaGamma_f / pow(Gamma_l_SMtree,2.0)
                     + (Gamma_l_SMtree*deltaGamma_had_2-Gamma_had_SMtree*deltaGamma_f_2) / pow(Gamma_l_SMtree,2.0);
        }
        if(f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3")){
            dR0_f_2 = Gamma_l_SMtree*pow(deltaGamma_had,2.0) / pow(Gamma_had_SMtree,3.0)
                     - deltaGamma_f*deltaGamma_had / pow(Gamma_had_SMtree,2.0)
                     + (Gamma_had_SMtree*deltaGamma_f_2-Gamma_l_SMtree*deltaGamma_had_2) / pow(Gamma_had_SMtree,2.0);
        }
        if(f.is("QUARK")){
            dR0_f_2 = Gamma_q_SMtree*pow(deltaGamma_had,2.0) / pow(Gamma_had_SMtree,3.0)
                     - deltaGamma_f*deltaGamma_had / pow(Gamma_had_SMtree,2.0)
                     + (Gamma_had_SMtree*deltaGamma_f_2-Gamma_q_SMtree*deltaGamma_had_2) / pow(Gamma_had_SMtree,2.0);
        }
    }
    
    return dR0_f_2;
}

//AG:end
////////////////////////////////////////////////////////////////////////

double NPbase::deltaGamma_Zf(const Particle f) const
{
    double deltaGamma_Zf = 0.;
    bool nonZeroNP = false;

    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f);
    
    double gVf = trueSM.gV_f(f).real();
    double gAf = trueSM.gA_f(f).real(); 
    
    double Nf;
    
    if (f.is("LEPTON")) {
        Nf = 1.0;
    } else {
        Nf = 3.0;
    }
    
    double alpha = trueSM.alphaMz();
    double sW2_SM = trueSM.sW2();
    double cW2_SM = trueSM.cW2();
    
    if (delGVf != 0.0 || delGAf != 0.0)
            nonZeroNP = true;
 
    if (nonZeroNP) {
        double delGammaZf = 0.0;
        delGammaZf = 2.0 * Nf * (gVf * delGVf + gAf * delGAf);
        
        deltaGamma_Zf = alpha * Mz / 12.0 / sW2_SM / cW2_SM * delGammaZf;
    }

    return deltaGamma_Zf;
}

double NPbase::Gamma_Zf(const Particle f) const                   //AG: Output()
{
    //AG:begin
    if(Output()==0){ return (trueSM.GammaZ(f) + deltaGamma_Zf(f)); }
    if(Output()==1){ return (trueSM.GammaZ(f) + deltaGamma_Zf_Test(f) ); }
    if(Output()==2){ return (trueSM.GammaZ(f) + deltaGamma_Zf(f) + deltaGamma_Zf_2(f) ); }
    if(Output()==3){ return (trueSM.GammaZ(f) + deltaGamma_Zf_Test(f) + deltaGamma_Zf_2(f) ); }
    if(Output()==4){ return (trueSM.GammaZ(f) + deltaGamma_Zf_2(f) ); }
    if(Output()==5){ return (trueSM.GammaZ(f) ); }
    //AG:end
    
    //return (trueSM.GammaZ(f) + deltaGamma_Zf(f));
}


double NPbase::deltaGamma_Z() const
{
    double deltaGamma_Z = 0.;
    bool nonZeroNP = false;

    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p] != 0.0 || delGAl[p] != 0.0
                || delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double deltaGl[6], deltaGq[6];
        double delGammaZ = 0.0;
        for (int p = 0; p < 6; ++p) {
            gVf = trueSM.gV_f(leptons[p]).real();
            gAf = trueSM.gA_f(leptons[p]).real();
            deltaGl[p] = 2.0 * (gVf * delGVl[p] + gAf * delGAl[p]);

            gVf = trueSM.gV_f(quarks[p]).real();
            gAf = trueSM.gA_f(quarks[p]).real();
            deltaGq[p] = 2.0 * (gVf * delGVq[p] + gAf * delGAq[p]);

            delGammaZ += deltaGl[p] + 3.0 * deltaGq[p];
        }

        double alpha = trueSM.alphaMz();       
        double sW2_SM = trueSM.sW2();
        double cW2_SM = trueSM.cW2();
        deltaGamma_Z = alpha * Mz / 12.0 / sW2_SM / cW2_SM
                * delGammaZ;
    }

    return deltaGamma_Z;
}

double NPbase::Gamma_Z() const                   //AG: NumericCheck() & Output()
{
    //AG:begin
    if(NumericCheck()) { 
        //// Testing the vector-axial Z-couplings:
       std::cout << " " << std::endl;
       std::cout << "Numeric Values of V-A couplings (from GammaZ)" << std::endl;
       double trueSMGV[6], trueSMGA[6], deltaGV[6], deltaGA[6], deltaGV2[6], deltaGA2[6];
       for (int p = 0; p < 6; ++p) {
           trueSMGV[p] = trueSM.gV_f(leptons[p]).real();
           trueSMGA[p] = trueSM.gA_f(leptons[p]).real();
           deltaGV[p] = deltaGV_f(leptons[p]);
           deltaGA[p] = deltaGA_f(leptons[p]);
           deltaGV2[p] = deltaGV_f_2(leptons[p]);
           deltaGA2[p] = deltaGA_f_2(leptons[p]);
       }
       std::cout << "  " << leptons[0].getName() << " " << leptons[2].getName() << " " << leptons[4].getName() << std::endl;
       std::cout << "    trueSM.gV_f(p) = " << trueSMGV[0] << " " << trueSMGV[2] << " " << trueSMGV[4] << std::endl;
       std::cout << "    deltaGV_f(p) = " << deltaGV[0] << " " << deltaGV[2] << " " << deltaGV[4] << std::endl;
       std::cout << "    deltaGV_f_2(p) = " << deltaGV2[0] << " " << deltaGV2[2] << " " << deltaGV2[4] << std::endl;
       std::cout << "    trueSM.gA_f(p) = " << trueSMGA[0] << " " << trueSMGA[2] << " " << trueSMGA[4] << std::endl;
       std::cout << "    deltaGA_f(p) = " << deltaGA[0] << " " << deltaGA[2] << " " << deltaGA[4] << std::endl;
       std::cout << "    deltaGA_f_2(p) = " << deltaGA2[0] << " " << deltaGA2[2] << " " << deltaGA2[4] << std::endl;
       std::cout << "  " << leptons[1].getName() << " " << leptons[3].getName() << " " << leptons[5].getName() << std::endl;
       std::cout << "    trueSM.gV_f(p) = " << trueSMGV[1] << " " << trueSMGV[3] << " " << trueSMGV[5] << std::endl;
       std::cout << "    deltaGV_f(p) = " << deltaGV[1] << " " << deltaGV[3] << " " << deltaGV[5] << std::endl;
       std::cout << "    deltaGV_f_2(p) = " << deltaGV2[1] << " " << deltaGV2[3] << " " << deltaGV2[5] << std::endl;
       std::cout << "    trueSM.gA_f(p) = " << trueSMGA[1] << " " << trueSMGA[3] << " " << trueSMGA[5] << std::endl;
       std::cout << "    deltaGA_f(p) = " << deltaGA[1] << " " << deltaGA[3] << " " << deltaGA[5] << std::endl;
       std::cout << "    deltaGA_f_2(p) = " << deltaGA2[1] << " " << deltaGA2[3] << " " << deltaGA2[5] << std::endl;
       for (int p = 0; p < 6; ++p) {
           trueSMGV[p] = trueSM.gV_f(quarks[p]).real();
           trueSMGA[p] = trueSM.gA_f(quarks[p]).real();
           deltaGV[p] = deltaGV_f(quarks[p]);
           deltaGA[p] = deltaGA_f(quarks[p]);
           deltaGV2[p] = deltaGV_f_2(quarks[p]);
           deltaGA2[p] = deltaGA_f_2(quarks[p]);
       }
       std::cout << "  " << quarks[0].getName() << " " << quarks[2].getName() << " " << quarks[4].getName() << std::endl;
       std::cout << "    trueSM.gV_f(p) = " << trueSMGV[0] << " " << trueSMGV[2] << " " << trueSMGV[4] << std::endl;
       std::cout << "    deltaGV_f(p) = " << deltaGV[0] << " " << deltaGV[2] << " " << deltaGV[4] << std::endl;
       std::cout << "    deltaGV_f_2(p) = " << deltaGV2[0] << " " << deltaGV2[2] << " " << deltaGV2[4] << std::endl;
       std::cout << "    trueSM.gA_f(p) = " << trueSMGA[0] << " " << trueSMGA[2] << " " << trueSMGA[4] << std::endl;
       std::cout << "    deltaGA_f(p) = " << deltaGA[0] << " " << deltaGA[2] << " " << deltaGA[4] << std::endl;
       std::cout << "    deltaGA_f_2(p) = " << deltaGA2[0] << " " << deltaGA2[2] << " " << deltaGA2[4] << std::endl;
       std::cout << "  " << quarks[1].getName() << " " << quarks[3].getName() << " " << quarks[5].getName() << std::endl;
       std::cout << "    trueSM.gV_f(p) = " << trueSMGV[1] << " " << trueSMGV[3] << " " << trueSMGV[5] << std::endl;
       std::cout << "    deltaGV_f(p) = " << deltaGV[1] << " " << deltaGV[3] << " " << deltaGV[5] << std::endl;
       std::cout << "    deltaGV_f_2(p) = " << deltaGV2[1] << " " << deltaGV2[3] << " " << deltaGV2[5] << std::endl;
       std::cout << "    trueSM.gA_f(p) = " << trueSMGA[1] << " " << trueSMGA[3] << " " << trueSMGA[5] << std::endl;
       std::cout << "    deltaGA_f(p) = " << deltaGA[1] << " " << deltaGA[3] << " " << deltaGA[5] << std::endl;
       std::cout << "    deltaGA_f_2(p) = " << deltaGA2[1] << " " << deltaGA2[3] << " " << deltaGA2[5] << std::endl;
       
       std::cout << " " << std::endl;
       std::cout << "ale = " << ale << std::endl;
       std::cout << "DeltaAlpha() = " << DeltaAlpha() << std::endl;
       std::cout << "alphaMz = " << alphaMz() << std::endl;
       std::cout << "trueSM.alphaMz() = " << trueSM.alphaMz() << std::endl;
       std::cout << "trueSM.sW2() = " << trueSM.sW2() << std::endl;
       std::cout << "trueSM.cW2() = " << trueSM.cW2() << std::endl;
       std::cout << "trueSM.alphaMz()*Mz/48/trueSM.sW2()/trueSM.cW2() = " << trueSM.alphaMz()*Mz/48/trueSM.sW2()/trueSM.cW2() << std::endl;

       //// Testing Gamma_Zf:
       //double prefactorAlpha = alphaMz() * Mz / 48.0 / trueSM.sW2() / trueSM.cW2();
       //double prefactorUni = GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI;
       //std::cout << prefactorAlpha << '\t' << prefactorUni << std::endl;
       std::cout << " " << std::endl;
       std::cout << "Numeric Values of the partial decay width Gamma_Zf:" << std::endl;
       std::cout << "  trueSM e: " << trueSM.GammaZ(leptons[ELECTRON]) << " " << trueSM.GammaZ(leptons[MU]) << " " << trueSM.GammaZ(leptons[TAU]) << std::endl;
       std::cout << "  trueSM nu: " << trueSM.GammaZ(leptons[NEUTRINO_1]) << " " << trueSM.GammaZ(leptons[NEUTRINO_2]) << " " << trueSM.GammaZ(leptons[NEUTRINO_3]) << std::endl;
       std::cout << "  trueSM d: " << trueSM.GammaZ(quarks[DOWN])<< " " << trueSM.GammaZ(quarks[STRANGE])<< " " << trueSM.GammaZ(quarks[BOTTOM])<<std::endl;
       std::cout << "  trueSM u: " << trueSM.GammaZ(quarks[UP]) << " " << trueSM.GammaZ(quarks[CHARM]) << " " << trueSM.GammaZ(quarks[TOP]) << std::endl;
       std::cout << "   deltaGamma_Zf e: " << deltaGamma_Zf(leptons[ELECTRON]) << " " << deltaGamma_Zf(leptons[MU]) << " " << deltaGamma_Zf(leptons[TAU]) << std::endl;
       std::cout << "   deltaGamma_Zf nu: " << deltaGamma_Zf(leptons[NEUTRINO_1]) << " " << deltaGamma_Zf(leptons[NEUTRINO_2]) << " " << deltaGamma_Zf(leptons[NEUTRINO_3]) << std::endl;
       std::cout << "   deltaGamma_Zf d: " << deltaGamma_Zf(quarks[DOWN]) << " " << deltaGamma_Zf(quarks[STRANGE]) << " " << deltaGamma_Zf(quarks[BOTTOM]) <<std::endl;
       std::cout << "   deltaGamma_Zf u: " << deltaGamma_Zf(quarks[UP]) << " " << deltaGamma_Zf(quarks[CHARM]) << " " << deltaGamma_Zf(quarks[TOP]) << std::endl;   
       std::cout << "  deltaGamma_Zf_Test e: " << deltaGamma_Zf_Test(leptons[ELECTRON]) << " " << deltaGamma_Zf_Test(leptons[MU]) << " " << deltaGamma_Zf_Test(leptons[TAU]) << std::endl;
       std::cout << "  deltaGamma_Zf_Test nu: " << deltaGamma_Zf_Test(leptons[NEUTRINO_1]) << " " << deltaGamma_Zf_Test(leptons[NEUTRINO_2]) << " " << deltaGamma_Zf_Test(leptons[NEUTRINO_3]) << std::endl;
       std::cout << "  deltaGamma_Zf_Test d: " << deltaGamma_Zf_Test(quarks[DOWN]) << " " << deltaGamma_Zf_Test(quarks[STRANGE]) << " " << deltaGamma_Zf_Test(quarks[BOTTOM])<<std::endl;
       std::cout << "  deltaGamma_Zf_Test u: " << deltaGamma_Zf_Test(quarks[UP]) << " " << deltaGamma_Zf_Test(quarks[CHARM]) << " " << deltaGamma_Zf_Test(quarks[TOP]) << std::endl; 
       std::cout << "   deltaGamma_Zf_2 e: " << deltaGamma_Zf_2(leptons[ELECTRON]) << " " << deltaGamma_Zf_2(leptons[MU]) << " " << deltaGamma_Zf_2(leptons[TAU]) << std::endl;
       std::cout << "   deltaGamma_Zf_2 nu: " << deltaGamma_Zf_2(leptons[NEUTRINO_1]) << " " << deltaGamma_Zf_2(leptons[NEUTRINO_2]) << " " << deltaGamma_Zf_2(leptons[NEUTRINO_3]) << std::endl;
       std::cout << "   deltaGamma_Zf_2 d: " << deltaGamma_Zf_2(quarks[DOWN])<< " " << deltaGamma_Zf_2(quarks[STRANGE])<< " " << deltaGamma_Zf_2(quarks[BOTTOM])<<std::endl;
       std::cout << "   deltaGamma_Zf_2 u: " << deltaGamma_Zf_2(quarks[UP]) << " " << deltaGamma_Zf_2(quarks[CHARM]) << " " << deltaGamma_Zf_2(quarks[TOP]) << std::endl; 

       ////Testing Gamma_Zhad:
       std::cout << " " << std::endl;
       std::cout << "Numeric values of GammaZhad:" << std::endl;
       std::cout << "  trueSM = " << trueSM.Gamma_had() << std::endl;
       std::cout << "  deltaGamma_Zhad = " << deltaGamma_Zhad() << std::endl;
       std::cout << "  deltaGamma_Zhad_Test = " << deltaGamma_Zhad_Test() << std::endl;
       std::cout << "  deltaGamma_Zhad_2 = " << deltaGamma_Zhad_2() << std::endl;

       ////Testing total GammaZ:
       std::cout << " " << std::endl;
       std::cout << "Numeric values of GammaZ:" << std::endl;
       std::cout << "  trueSM = " << trueSM.Gamma_Z() << std::endl;
       std::cout << "  deltaGamma_Z = " << deltaGamma_Z() << std::endl;
       std::cout << "  deltaGamma_Z_Test() = " << deltaGamma_Z_Test() << std::endl;
       std::cout << "  deltaGamma_Z_2() = " << deltaGamma_Z_2() << std::endl;

       //// Testing prefactor:
       /*std::cout << " " << std::endl;
       std::cout << "(NPbase) Gamma_Z Prefactor: " << std::endl;
       std::cout << "  Mw_inp = " << Mw_inp  << std::endl;
       std::cout << "  trueSM.Mw() = " << trueSM.Mw() << std::endl;
       std::cout << "  Mw() = " << Mw() << std::endl;
       std::cout << "  ale = " << ale << std::endl;
       std::cout << "  trueSM.alphaMz() = " << trueSM.alphaMz() << std::endl;
       std::cout << "  alphaMz() = " << alphaMz() << std::endl;
       std::cout << "    (trueMw/Mz)^2 = " << pow(trueSM.Mw()/Mz,2.0) << std::endl;
       double sW2_SM = trueSM.sW2();
       double cW2_SM = trueSM.cW2();
       std::cout << "    sW2_SM = " << sW2_SM << std::endl;
       std::cout << "    cW2_SM = " << cW2_SM << std::endl;
       std::cout << "    trueSM.sW2() * trueSM.cW2() = " << trueSM.sW2() * trueSM.cW2() << std::endl;
       std::cout << "    alphaMz() * Mz / 48.0 / sW2_SM / cW2_SM = " << alphaMz() * Mz / 48.0 / sW2_SM / cW2_SM << std::endl;
       std::cout << "    GF * pow(Mz,3.0) / 24.0 / sqrt(2.0) / M_PI = " << GF * pow(Mz,3.0) / 24.0 / sqrt(2.0) / M_PI << std::endl;
       */
    }

    if(Output()==0){ return (trueSM.Gamma_Z() + deltaGamma_Z()); }
    if(Output()==1){ return (trueSM.Gamma_Z() + deltaGamma_Z_Test() ); }
    if(Output()==2){ return (trueSM.Gamma_Z() + deltaGamma_Z() + deltaGamma_Z_2() ); }
    if(Output()==3){ return (trueSM.Gamma_Z() + deltaGamma_Z_Test() + deltaGamma_Z_2() ); }
    if(Output()==4){ return (trueSM.Gamma_Z() + deltaGamma_Z_2() ); }
    if(Output()==5){ return (trueSM.Gamma_Z() ); }
    //AG:end
    
    //return (trueSM.Gamma_Z() + deltaGamma_Z());
}

//AG:begin
double NPbase::deltaRuc_2() const
{
    double DeltaRuc_2 = 0.0;
    
    bool nonZeroNP=false;
    if(deltaR0_f_2(quarks[UP])!=0.0 || deltaR0_f_2(quarks[CHARM])!=0.0) { nonZeroNP = true;}
    
    if(nonZeroNP){
        // This keeps the same structure in SM. Should it be modified with a more general CKM assumption?
        DeltaRuc_2 = 0.5 * ( deltaR0_f_2(quarks[UP]) + deltaR0_f_2(quarks[CHARM]) );
    }
    
    return DeltaRuc_2;
}

double NPbase::deltaRuc() const
{
    double DeltaRuc = 0.0;
    
    bool nonZeroNP=false;
    if(deltaR0_f(quarks[UP])!=0.0 || deltaR0_f(quarks[CHARM])!=0.0) { nonZeroNP = true;}
    
    if(nonZeroNP){
        // This keeps the same structure in SM. Should it be modified with a more general CKM assumption?
        DeltaRuc = 0.5 * ( deltaR0_f(quarks[UP]) + deltaR0_f(quarks[CHARM]) );
    }
    
    return DeltaRuc;
}

double NPbase::Ruc() const                       //AG: NumericCheck() & Output()
{
    if(NumericCheck()){
        std::cout << " " << std::endl;
        std::cout << "Numeric Values of Ruc():" << std::endl;
        std::cout << "  trueSM.Ruc() = " << trueSM.Ruc() << std::endl;
        std::cout << "  deltaRuc() = " << deltaRuc() << std::endl;
        std::cout << "  deltaRuc_2() = " << deltaRuc_2() << std::endl;
    }
    
    if(Output()==0 || Output()==1){ return (trueSM.Ruc() + deltaRuc()); }
    if(Output()==2 || Output()==3){ return (trueSM.Ruc() + deltaRuc() + deltaRuc_2() ); }
    if(Output()==4){ return (trueSM.Ruc() + deltaRuc_2() ); }
    if(Output()==5){ return (trueSM.Ruc() ); }
    
    //return ( trueSM.Ruc() + deltaRuc() );
}
//AG:end

double NPbase::RZlilj(const Particle li, const Particle lj) const
{
    double GammZli, GammZlj;
    
    if ( li.is("ELECTRON") || li.is("MU") || li.is("TAU") )
        GammZli = Gamma_Zf(li);        
    else
        throw std::runtime_error("Error in NPbase::RZlilj. li must be a charged lepton");
    
    if ( lj.is("ELECTRON") || lj.is("MU") || lj.is("TAU") )
        GammZlj = Gamma_Zf(lj);        
    else
        throw std::runtime_error("Error in NPbase::RZlilj. lj must be a charged lepton");
    
    return GammZli/GammZlj;
}


double NPbase::deltaGamma_Zhad() const
{
    double deltaGamma_Zhad = 0.;
    bool nonZeroNP = false;

    double delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double deltaGq[6];
        double delGammaZhad = 0.0;
        for (int p = 0; p < 6; ++p) {

            gVf = trueSM.gV_f(quarks[p]).real();
            gAf = trueSM.gA_f(quarks[p]).real();
            deltaGq[p] = 2.0 * (gVf * delGVq[p] + gAf * delGAq[p]);

            delGammaZhad += 3.0 * deltaGq[p];
        }

        double alpha = trueSM.alphaMz();
        double sW2_SM = trueSM.sW2();
        double cW2_SM = trueSM.cW2();
        deltaGamma_Zhad = alpha * Mz / 12.0 / sW2_SM / cW2_SM
                * delGammaZhad;
    }

    return deltaGamma_Zhad;
}

double NPbase::Gamma_had() const
{
    return (trueSM.Gamma_had() + deltaGamma_Zhad());
}

double NPbase::BR_Zf(const Particle f) const
{
    double delGammaZTot = deltaGamma_Z();
    double delGammaZf = deltaGamma_Zf(f);
    
    double GammaZTotSM = trueSM.Gamma_Z();
    double GammaZfSM = trueSM.GammaZ(f);
    
    return (GammaZfSM/GammaZTotSM + delGammaZf/GammaZTotSM - GammaZfSM * delGammaZTot /GammaZTotSM/GammaZTotSM);
}

double NPbase::deltaSigmaHadron() const
{
    double sigma_had = 0.;
    bool nonZeroNP = false;

    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p] != 0.0 || delGAl[p] != 0.0
                || delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double Gl[6], deltaGl[6], Gq[6], deltaGq[6];
        double Gq_sum = 0.0, delGq_sum = 0.0;
        double Gf_sum = 0.0, delGf_sum = 0.0;
        for (int p = 0; p < 6; ++p) {
            gVf = trueSM.gV_f(leptons[p]).real();
            gAf = trueSM.gA_f(leptons[p]).real();
            Gl[p] = gVf * gVf + gAf*gAf;
            deltaGl[p] = 2.0 * (gVf * delGVl[p] + gAf * delGAl[p]);

            gVf = trueSM.gV_f(quarks[p]).real();
            gAf = trueSM.gA_f(quarks[p]).real();
            Gq[p] = gVf * gVf + gAf*gAf;
            deltaGq[p] = 2.0 * (gVf * delGVq[p] + gAf * delGAq[p]);

            Gq_sum += 3.0 * Gq[p];
            Gf_sum += Gl[p] + 3.0 * Gq[p];
            delGq_sum += 3.0 * deltaGq[p];
            delGf_sum += deltaGl[p] + 3.0 * deltaGq[p];
        }

        sigma_had = 12.0 * M_PI / Mz / Mz
                * Gl[ELECTRON] * Gq_sum / Gf_sum / Gf_sum
                * (deltaGl[ELECTRON] / Gl[ELECTRON]
                + delGq_sum / Gq_sum - 2.0 * delGf_sum / Gf_sum);
    }

    return sigma_had;
}

double NPbase::sigma0_had() const                //AG: NumericCheck() & Output()
{
    //AG:begin
    if(NumericCheck()){
        std::cout << "" << std::endl;
        std::cout << "Numeric values of sigma0_had : " << std::endl;
        std::cout << "  trueSM.sigma0_had = " << trueSM.sigma0_had() << std::endl;
        std::cout << "  deltaSigmaHadron = " << deltaSigmaHadron() << std::endl;
        std::cout << "  deltaSigmaHadron_Test = " << deltaSigmaHadron_Test() << std::endl;
        std::cout << "  deltaSigmaHadron_2 = " << deltaSigmaHadron_2()  << std::endl;
    }
    
    if(Output()==0){ return (trueSM.sigma0_had() + deltaSigmaHadron()); }
    if(Output()==1){ return (trueSM.sigma0_had() + deltaSigmaHadron_Test() ); }
    if(Output()==2){ return (trueSM.sigma0_had() + deltaSigmaHadron() + deltaSigmaHadron_2() ); }
    if(Output()==3){ return (trueSM.sigma0_had() + deltaSigmaHadron_Test() + deltaSigmaHadron_2() ); }
    if(Output()==4){ return (trueSM.sigma0_had() + deltaSigmaHadron_2() ); }
    if(Output()==5){ return (trueSM.sigma0_had() ); }
    //AG:end
    
    //return (trueSM.sigma0_had() + deltaSigmaHadron());
}

double NPbase::deltaSin2thetaEff_e() const
{
    double sin2_theta_eff = 0.;
    double delGVf = deltaGV_f(leptons[ELECTRON]);
    double delGAf = deltaGA_f(leptons[ELECTRON]);
    if (delGVf != 0.0 || delGAf != 0.0) {
        double gVf = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAf = trueSM.gA_f(leptons[ELECTRON]).real();
        double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

        sin2_theta_eff = -delGVfOverGAf / 4.0;
    }
    return sin2_theta_eff;
}

double NPbase::deltaSin2thetaEff_mu() const
{
    double sin2_theta_eff = 0.;
    double delGVf = deltaGV_f(leptons[MU]);
    double delGAf = deltaGA_f(leptons[MU]);
    if (delGVf != 0.0 || delGAf != 0.0) {
        double gVf = trueSM.gV_f(leptons[MU]).real();
        double gAf = trueSM.gA_f(leptons[MU]).real();
        double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

        sin2_theta_eff = -delGVfOverGAf / 4.0;
    }
    return sin2_theta_eff;
}

double NPbase::sin2thetaEff(const Particle f) const     //AG: NumericCheck() & Output()
{
    //AG:begin
    if(NumericCheck()){
        std::cout << " "  << std::endl;
        std::cout << "Numeric Values of sin2thetaEff: " << std::endl;
        std::cout << "  trueSM.sin2thetaEff_e = " << trueSM.sin2thetaEff(leptons[ELECTRON]) << std::endl;
        std::cout << "  trueSM.sin2thetaEff_mu = " << trueSM.sin2thetaEff(leptons[MU]) << std::endl;
        std::cout << "  trueSM.sin2thetaEff_tau = " << trueSM.sin2thetaEff(leptons[TAU]) << std::endl;
        std::cout << "  deltaSin2thetaEff_e = " << deltaSin2thetaEff_e() << std::endl;
        std::cout << "  deltaSin2thetaEff_mu = " << deltaSin2thetaEff_mu() << std::endl;
        std::cout << "  deltaSin2thetaEff_e_2 = " << deltaSin2thetaEff_e_2() << std::endl;
        std::cout << "  deltaSin2thetaEff_mu_2 = " << deltaSin2thetaEff_mu_2() << std::endl;
    }
    //AG:end    
    
    if (f.is("ELECTRON")){
        //AG:begin
        if(Output()==0 || Output()==1){ return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_e()); }
        if(Output()==2 || Output()==3){ return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_e() + deltaSin2thetaEff_e_2() ); }
        if(Output()==4){ return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_e_2() ); }
        if(Output()==5){ return (trueSM.sin2thetaEff(f)); }
        //AG:end
        //return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_e() );
    }
    else if (f.is("MU")){
        //AG:begin
        if(Output()==0 || Output()==1){ return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_mu()); }
        if(Output()==2 || Output()==3){ return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_mu() + deltaSin2thetaEff_mu_2() ); }
        if(Output()==4){ return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_mu_2() ); }
        if(Output()==5){ return (trueSM.sin2thetaEff(f) ); }
        //AG:end
        //return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_mu() );
    }
    else
        return (trueSM.sin2thetaEff(f));
}

double NPbase::deltaA_f(const Particle f) const
{
    double dAf = 0.;
    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f);
    if (delGVf != 0.0 || delGAf != 0.0) {
        double gVf = trueSM.gV_f(f).real();
        double gAf = trueSM.gA_f(f).real();
        double Gf = gVf * gVf + gAf*gAf;
        double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

        dAf = -2.0 * (gVf * gVf - gAf * gAf) * gAf * gAf / Gf / Gf*delGVfOverGAf;
    }

    return dAf;
}

double NPbase::A_f(const Particle f) const       //AG: NumericCheck() & Output()
{
    //AG: begin
    if(NumericCheck()){
        std::cout << " " << std::endl;
        std::cout << "Numeric Values of Af:" << std::endl;
        std::cout << "  trueSM electron-type :" << trueSM.A_f(leptons[ELECTRON]) << " " << trueSM.A_f(leptons[MU]) << " " << trueSM.A_f(leptons[TAU]) << std::endl;
        std::cout << "  trueSM neutrino-type :" << trueSM.A_f(leptons[NEUTRINO_1]) << " " << trueSM.A_f(leptons[NEUTRINO_2]) << " " << trueSM.A_f(leptons[NEUTRINO_3]) << std::endl;       
        std::cout << "  trueSM down-type :" << trueSM.A_f(quarks[DOWN]) << " " << trueSM.A_f(quarks[STRANGE]) << " " << trueSM.A_f(quarks[BOTTOM]) << std::endl;
        std::cout << "  trueSM up-type :" << trueSM.A_f(quarks[UP]) << " " << trueSM.A_f(quarks[CHARM]) << " " << trueSM.A_f(quarks[TOP]) << std::endl;
        std::cout << "    deltaAf electron-type :" << deltaA_f(leptons[ELECTRON]) << " " << deltaA_f(leptons[MU]) << " " << deltaA_f(leptons[TAU]) << std::endl;
        std::cout << "    deltaAf neutrino-type :" << deltaA_f(leptons[NEUTRINO_1]) << " " << deltaA_f(leptons[NEUTRINO_2]) << " " << deltaA_f(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "    deltaAf up-type :" << deltaA_f(quarks[UP]) << " " << deltaA_f(quarks[CHARM]) << " " << deltaA_f(quarks[TOP]) << std::endl;
        std::cout << "    deltaAf down-type :" << deltaA_f(quarks[DOWN]) << " " << deltaA_f(quarks[STRANGE]) << " " << deltaA_f(quarks[BOTTOM]) << std::endl;
        std::cout << "  deltaAf2 electron-type :" << deltaA_f_2(leptons[ELECTRON]) << " " << deltaA_f_2(leptons[MU]) << " " << deltaA_f_2(leptons[TAU]) << std::endl;
        std::cout << "  deltaAf2 neutrino-type :" << deltaA_f_2(leptons[NEUTRINO_1]) << " " << deltaA_f_2(leptons[NEUTRINO_2]) << " " << deltaA_f_2(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "  deltaAf2 up-type :" << deltaA_f_2(quarks[UP]) << " " << deltaA_f_2(quarks[CHARM]) << " " << deltaA_f_2(quarks[TOP]) << std::endl;
        std::cout << "  deltaAf2 down-type :" << deltaA_f_2(quarks[DOWN]) << " " << deltaA_f_2(quarks[STRANGE]) << " " << deltaA_f_2(quarks[BOTTOM]) << std::endl;
    }
    
    if(Output()==0 || Output()==1){ return (trueSM.A_f(f) + deltaA_f(f)); }
    if(Output()==2 || Output()==3){ return (trueSM.A_f(f) + deltaA_f(f) + deltaA_f_2(f) ); }
    if(Output()==4){ return (trueSM.A_f(f) + deltaA_f_2(f) ); }
    if(Output()==5){ return (trueSM.A_f(f) ); }
    //AG:end
    
    //return (trueSM.A_f(f) + deltaA_f(f));
}

double NPbase::deltaAFB(const Particle f) const
{
    double dAFB = 0.;
    double delGVf = deltaGV_f(f);
    double delGAf = deltaGA_f(f);
    if (f.is("ELECTRON")) {
        if (delGVf != 0.0 || delGAf != 0.0) {
            double gVe = trueSM.gV_f(f).real();
            double gAe = trueSM.gA_f(f).real();
            double Ge = gVe * gVe + gAe*gAe;
            double delGVeOverGAe = (gAe * delGVf - gVe * delGAf) / gAe / gAe;
            dAFB = -6.0 * gVe * gAe * (gVe * gVe - gAe * gAe) * gAe * gAe / Ge / Ge / Ge*delGVeOverGAe;
        }
    } else {
        double delGVe = deltaGV_f(leptons[ELECTRON]);
        double delGAe = deltaGA_f(leptons[ELECTRON]);
        if (delGVe != 0.0 || delGAe != 0.0 || delGVf != 0.0 || delGAf != 0.0) {
            double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
            double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
            double Ge = gVe * gVe + gAe*gAe;
            double delGVeOverGAe = (gAe * delGVe - gVe * delGAe) / gAe / gAe;
            //
            double gVf = trueSM.gV_f(f).real();
            double gAf = trueSM.gA_f(f).real();
            double Gf = gVf * gVf + gAf*gAf;
            double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

            dAFB = -(3.0 * gVf * gAf * (gVe * gVe - gAe * gAe) * gAe * gAe / Gf / Ge / Ge * delGVeOverGAe
                    + 3.0 * gVe * gAe * (gVf * gVf - gAf * gAf) * gAf * gAf / Ge / Gf / Gf * delGVfOverGAf);
        }
    }

    return dAFB;
}

double NPbase::AFB(const Particle f) const       //AG: NumericCheck() & Output()
{
    //AG: begin
    if(NumericCheck()){
        std::cout << " " << std::endl;
        std::cout << "Numeric Values of AFB:" << std::endl;
        std::cout << "  trueSM electron-type :" << trueSM.AFB(leptons[ELECTRON]) << " " << trueSM.AFB(leptons[MU]) << " " << trueSM.AFB(leptons[TAU]) << std::endl;
        std::cout << "  trueSM neutrino-type :" << trueSM.AFB(leptons[NEUTRINO_1]) << " " << trueSM.AFB(leptons[NEUTRINO_2]) << " " << trueSM.AFB(leptons[NEUTRINO_3]) << std::endl;       
        std::cout << "  trueSM down-type :" << trueSM.AFB(quarks[DOWN]) << " " << trueSM.AFB(quarks[STRANGE]) << " " << trueSM.AFB(quarks[BOTTOM]) << std::endl;
        std::cout << "  trueSM up-type :" << trueSM.AFB(quarks[UP]) << " " << trueSM.AFB(quarks[CHARM]) << " " << trueSM.AFB(quarks[TOP]) << std::endl;
        std::cout << "    deltaAFB electron-type :" << deltaAFB(leptons[ELECTRON]) << " " << deltaAFB(leptons[MU]) << " " << deltaAFB(leptons[TAU]) << std::endl;
        std::cout << "    deltaAFB neutrino-type :" << deltaAFB(leptons[NEUTRINO_1]) << " " << deltaAFB(leptons[NEUTRINO_2]) << " " << deltaAFB(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "    deltaAFB up-type :" << deltaAFB(quarks[UP]) << " " << deltaAFB(quarks[CHARM]) << " " << deltaAFB(quarks[TOP]) << std::endl;
        std::cout << "    deltaAFB down-type :" << deltaAFB(quarks[DOWN]) << " " << deltaAFB(quarks[STRANGE]) << " " << deltaAFB(quarks[BOTTOM]) << std::endl;
        std::cout << "  deltaAFBTest electron-type :" << deltaAFB_Test(leptons[ELECTRON]) << " " << deltaAFB_Test(leptons[MU]) << " " << deltaAFB_Test(leptons[TAU]) << std::endl;
        std::cout << "  deltaAFBTest neutrino-type :" << deltaAFB_Test(leptons[NEUTRINO_1]) << " " << deltaAFB_Test(leptons[NEUTRINO_2]) << " " << deltaAFB_Test(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "  deltaAFBTest up-type :" << deltaAFB_Test(quarks[UP]) << " " << deltaAFB_Test(quarks[CHARM]) << " " << deltaAFB_Test(quarks[TOP]) << std::endl;
        std::cout << "  deltaAFBTest down-type :" << deltaAFB_Test(quarks[DOWN]) << " " << deltaAFB_Test(quarks[STRANGE]) << " " << deltaAFB_Test(quarks[BOTTOM]) << std::endl;
        std::cout << "    deltaAFB2 electron-type :" << deltaAFB_2(leptons[ELECTRON]) << " " << deltaAFB_2(leptons[MU]) << " " << deltaAFB_2(leptons[TAU]) << std::endl;
        std::cout << "    deltaAFB2 neutrino-type :" << deltaAFB_2(leptons[NEUTRINO_1]) << " " << deltaAFB_2(leptons[NEUTRINO_2]) << " " << deltaAFB_2(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "    deltaAFB2 up-type :" << deltaAFB_2(quarks[UP]) << " " << deltaAFB_2(quarks[CHARM]) << " " << deltaAFB_2(quarks[TOP]) << std::endl;
        std::cout << "    deltaAFB2 down-type :" << deltaAFB_2(quarks[DOWN]) << " " << deltaAFB_2(quarks[STRANGE]) << " " << deltaAFB_2(quarks[BOTTOM]) << std::endl;   
    }

    if(Output()==0){ return (trueSM.AFB(f) + deltaAFB(f)); }
    if(Output()==1){ return (trueSM.AFB(f) + deltaAFB_Test(f) ); }
    if(Output()==2){ return (trueSM.AFB(f) + deltaAFB(f) + deltaAFB_2(f) ); }
    if(Output()==3){ return (trueSM.AFB(f) + deltaAFB_Test(f) + deltaAFB_2(f) ); }
    if(Output()==4){ return (trueSM.AFB(f) + deltaAFB_2(f) ); }
    if(Output()==5){ return (trueSM.AFB(f) ); }
    //AG:end
    
    //return (trueSM.AFB(f) + deltaAFB(f));
}

double NPbase::deltaR0_f(const Particle f) const
{
    double dR0_f = 0., delGVl = 0., delGAl = 0., deltaGl = 0., Gl = 0.;
    bool nonZeroNP = false;
    if (f.is("LEPTON")) {
        delGVl = deltaGV_f(f);
        delGAl = deltaGA_f(f);
        if (delGVl != 0.0 || delGAl != 0.0) nonZeroNP = true;
    }

    double delGVq[6], delGAq[6];
    for (int q = 0; q < 6; ++q) {
        delGVq[q] = deltaGV_f(quarks[q]);
        delGAq[q] = deltaGA_f(quarks[q]);
        if (delGVq[q] != 0.0 || delGAq[q] != 0.0) nonZeroNP = true;
    }

    if (nonZeroNP) {
        double CF = 1.;
        if (f.is("LEPTON")) {
            double gVl = trueSM.gV_f(f).real();
            double gAl = trueSM.gA_f(f).real();
            Gl = gVl * gVl + gAl*gAl;
            deltaGl = 2.0 * (gVl * delGVl + gAl * delGAl);
            CF = 3.;
        }
        double Gq[6], deltaGq[6];
        double gVq, gAq;
        double Gq_sum = 0.0, delGq_sum = 0.0;
        for (int q = 0; q < 6; ++q) {
            gVq = trueSM.gV_f(quarks[q]).real();
            gAq = trueSM.gA_f(quarks[q]).real();
            Gq[q] = gVq * gVq + gAq*gAq;
            deltaGq[q] = 2.0 * (gVq * delGVq[q] + gAq * delGAq[q]);

            Gq_sum += CF * Gq[q];
            delGq_sum += CF * deltaGq[q];
        }
        if (f.is("LEPTON"))
            if ( f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3")  ) {
                dR0_f = deltaGl / Gq_sum - Gl * delGq_sum / Gq_sum / Gq_sum;                
            } else {
                dR0_f = delGq_sum / Gl - Gq_sum * deltaGl / Gl / Gl;
            }
        else
            dR0_f = deltaGq[f.getIndex() - 6] / Gq_sum
                - Gq[f.getIndex() - 6] * delGq_sum / Gq_sum / Gq_sum;
    }
    return dR0_f;
}

double NPbase::R0_f(const Particle f) const      //AG: NumericCheck() & Output()
{
    //AG:begin
    if(NumericCheck()){
        std::cout << "" << std::endl;
        std::cout << "Numeric values of R0_l : " << std::endl;
        std::cout << "  trueSM.R0_e,mu,tau = " << trueSM.R0_f(leptons[ELECTRON]) <<" "
                                               << trueSM.R0_f(leptons[MU]) <<" "
                                               << trueSM.R0_f(leptons[TAU]) << std::endl;
        std::cout << "  trueSM.R0_nu1,nu2,nu3 = " << trueSM.R0_f(leptons[NEUTRINO_1]) <<" "
                                                  << trueSM.R0_f(leptons[NEUTRINO_2]) <<" "
                                                  << trueSM.R0_f(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "  trueSM.R0_u,c = " << trueSM.R0_f(quarks[UP]) << " " 
                                          << trueSM.R0_f(quarks[CHARM])<< std::endl;
        std::cout << "  trueSM.R0_b = " << trueSM.R0_f(quarks[BOTTOM]) << std::endl;
        std::cout << "    deltaR0_e,mu,tau = " << deltaR0_f(leptons[ELECTRON]) <<" "
                                               << deltaR0_f(leptons[MU]) <<" "
                                               << deltaR0_f(leptons[TAU]) << std::endl;
        std::cout << "    deltaR0_nu1,nu2,nu3 = " << deltaR0_f(leptons[NEUTRINO_1]) <<" "
                                                  << deltaR0_f(leptons[NEUTRINO_2]) <<" "
                                                  << deltaR0_f(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "    deltaR0_u,c = " << deltaR0_f(quarks[UP]) << " " 
                                          << deltaR0_f(quarks[CHARM]) << std::endl;
        std::cout << "    deltaR0_d,s,b = " << deltaR0_f(quarks[DOWN]) << " " 
                                            << deltaR0_f(quarks[STRANGE]) << " " 
                                            << deltaR0_f(quarks[BOTTOM]) << std::endl;
        std::cout << "  deltaR0_Test_e,mu,tau = " << deltaR0_f_Test(leptons[ELECTRON]) <<" "
                                                  << deltaR0_f_Test(leptons[MU]) <<" "
                                                  << deltaR0_f_Test(leptons[TAU]) << std::endl;
        std::cout << "  deltaR0_Test_nu1,nu2,nu3 = " << deltaR0_f_Test(leptons[NEUTRINO_1]) <<" "
                                                     << deltaR0_f_Test(leptons[NEUTRINO_2]) <<" "
                                                     << deltaR0_f_Test(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "  deltaR0_Test_u,c = " << deltaR0_f_Test(quarks[UP]) << " " 
                                             << deltaR0_f_Test(quarks[CHARM]) << std::endl;
        std::cout << "  deltaR0_Test_d,s,b = " << deltaR0_f_Test(quarks[DOWN]) << " " 
                                               << deltaR0_f_Test(quarks[STRANGE]) << " " 
                                               << deltaR0_f_Test(quarks[BOTTOM]) << std::endl;
        std::cout << "    deltaR0_2_e,mu,tau = " << deltaR0_f_2(leptons[ELECTRON]) <<" "
                                                 << deltaR0_f_2(leptons[MU]) <<" "
                                                 << deltaR0_f_2(leptons[TAU]) << std::endl;
        std::cout << "    deltaR0_2_nu1,nu2,nu3 = " << deltaR0_f_2(leptons[NEUTRINO_1]) <<" "
                                                    << deltaR0_f_2(leptons[NEUTRINO_2]) <<" "
                                                    << deltaR0_f_2(leptons[NEUTRINO_3]) << std::endl;
        std::cout << "    deltaR0_2_u,c = " << deltaR0_f_2(quarks[UP]) << " " 
                                            << deltaR0_f_2(quarks[CHARM]) << std::endl;
        std::cout << "    deltaR0_2_d,s,b = " << deltaR0_f_2(quarks[DOWN]) << " " 
                                              << deltaR0_f_2(quarks[STRANGE]) << " " 
                                              << deltaR0_f_2(quarks[BOTTOM]) << std::endl;
    }

    if(Output()==0){ return (trueSM.R0_f(f) + deltaR0_f(f)); }
    if(Output()==1){ return (trueSM.R0_f(f) + deltaR0_f_Test(f) ); }
    if(Output()==2){ return (trueSM.R0_f(f) + deltaR0_f(f) + deltaR0_f_2(f) ); }
    if(Output()==3){ return (trueSM.R0_f(f) + deltaR0_f_Test(f)  + deltaR0_f_2(f) ); }
    if(Output()==4){ return (trueSM.R0_f(f) + deltaR0_f_2(f) ); }
    if(Output()==5){ return (trueSM.R0_f(f) ); }
    //AG:end
    
    //return (trueSM.R0_f(f) + deltaR0_f(f));
}

double NPbase::deltaR_inv() const
{
    double dR_inv = 0., delGVe = 0., delGAe = 0., deltaGe = 0., Ge = 0.;
    bool nonZeroNP = false;

    delGVe = deltaGV_f(leptons[ELECTRON]);
    delGAe = deltaGA_f(leptons[ELECTRON]);
    if (delGVe != 0.0 || delGAe != 0.0) nonZeroNP = true;

    double delGVnu[3], delGAnu[3];
    for (int p = 0; p < 3; ++p) {
        delGVnu[p] = deltaGV_f(leptons[2*p]);
        delGAnu[p] = deltaGA_f(leptons[2*p]);
        if (delGVnu[p] != 0.0 || delGAnu[p] != 0.0 ) nonZeroNP = true;
    }

    if (nonZeroNP) {

        double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
        Ge = gVe * gVe + gAe * gAe;
        deltaGe = 2.0 * (gVe * delGVe + gAe * delGAe);

        double Gnu[3], deltaGnu[3];
        double gVnu, gAnu;
        double Gnu_sum = 0.0, delGnu_sum = 0.0;
        for (int p = 0; p < 3; ++p) {
            gVnu = trueSM.gV_f(leptons[2*p]).real();
            gAnu = trueSM.gA_f(leptons[2*p]).real();
              

            Gnu[p] = gVnu * gVnu + gAnu * gAnu;

            deltaGnu[p] = 2.0 * (gVnu * delGVnu[p] + gAnu * delGAnu[p]);

            Gnu_sum += Gnu[p];
            delGnu_sum += deltaGnu[p];
        }

        dR_inv = delGnu_sum / Ge - Gnu_sum * deltaGe / Ge / Ge;
    }
    return dR_inv;
}

double NPbase::R_inv() const
{
    return ( trueSM.R_inv() + deltaR_inv() );
}


double NPbase::deltaN_nu() const
{   
    double dNnu = 0.0;
    double dGl1, dGl2, dGl3, dGl, dGinv;
    double Gl1, Gl2, Gl3, Gl, Ginv;
    double dRl1, dRl2, dRl3, dRl;
    double Rl1, Rl2, Rl3, Rl;
    double shad0;
    
    dGl1 = deltaGamma_Zf(leptons[ELECTRON]);
    dGl2 = deltaGamma_Zf(leptons[MU]);
    dGl3 = deltaGamma_Zf(leptons[TAU]);
    
    dGl = (1.0/3.0) * (dGl1 + dGl2 + dGl3);
    
    Gl1 = trueSM.GammaZ(leptons[ELECTRON]);
    Gl2 = trueSM.GammaZ(leptons[MU]);
    Gl3 = trueSM.GammaZ(leptons[TAU]);
    
    Gl = (1.0/3.0) * (Gl1 + Gl2 + Gl3);
    
    dGinv = deltaGamma_Zf(leptons[NEUTRINO_1]) + 
            deltaGamma_Zf(leptons[NEUTRINO_2]) +
            deltaGamma_Zf(leptons[NEUTRINO_3]);
    
    Ginv = trueSM.GammaZ(leptons[NEUTRINO_1]) + 
            trueSM.GammaZ(leptons[NEUTRINO_2]) +
            trueSM.GammaZ(leptons[NEUTRINO_3]);
    
    dRl1 = deltaR0_f(leptons[ELECTRON]);
    dRl2 = deltaR0_f(leptons[MU]);
    dRl3 = deltaR0_f(leptons[TAU]);
    
    dRl = (1.0/3.0) * (dRl1 + dRl2 + dRl3);
    
    Rl1 = trueSM.R0_f(leptons[ELECTRON]);
    Rl2 = trueSM.R0_f(leptons[MU]);
    Rl3 = trueSM.R0_f(leptons[TAU]);
    
    Rl = (1.0/3.0) * (Rl1 + Rl2 + Rl3);
    
    shad0 = trueSM.sigma0_had();
    
    dNnu = (trueSM.N_nu())*( dGl/Gl - dGinv/Ginv ) -
            3.0*(Gl/Ginv)*dRl +
            (Gl/Ginv)*sqrt(3.0*M_PI*Rl/Mz/Mz/shad0)*(-3.0*deltaSigmaHadron()/shad0 + 3.0* dRl/Rl);

    return dNnu;
}

double NPbase::N_nu() const
{
    return ( trueSM.N_nu() + deltaN_nu() );
}