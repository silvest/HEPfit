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

double NPbase::alphaMz() const
{
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
    return ( trueSM.gV_f(f) + deltaGV_f(f));
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
    return ( trueSM.gA_f(f) + deltaGA_f(f));
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
    
    double sW2_SM = trueSM.sW2();
    double cW2_SM = trueSM.cW2();
    
    if (delGVf != 0.0 || delGAf != 0.0)
            nonZeroNP = true;
 
    if (nonZeroNP) {
        double delGammaZf = 0.0;
        delGammaZf = 2.0 * Nf * (gVf * delGVf + gAf * delGAf);
        
        deltaGamma_Zf = alphaMz() * Mz / 12.0 / sW2_SM / cW2_SM * delGammaZf;
    }

    return deltaGamma_Zf;
}

double NPbase::Gamma_Zf(const Particle f) const
{
    return (trueSM.GammaZ(f) + deltaGamma_Zf(f));
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

        double sW2_SM = trueSM.sW2();
        double cW2_SM = trueSM.cW2();
        deltaGamma_Z = alphaMz() * Mz / 12.0 / sW2_SM / cW2_SM
                * delGammaZ;
    }

    return deltaGamma_Z;
}

double NPbase::Gamma_Z() const
{
    return (trueSM.Gamma_Z() + deltaGamma_Z());
}


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

        double sW2_SM = trueSM.sW2();
        double cW2_SM = trueSM.cW2();
        deltaGamma_Zhad = alphaMz() * Mz / 12.0 / sW2_SM / cW2_SM
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

double NPbase::sigma0_had() const
{
    return (trueSM.sigma0_had() + deltaSigmaHadron());
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

double NPbase::sin2thetaEff(const Particle f) const
{
    if (f.is("ELECTRON"))
        return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_e());
    else if (f.is("MU"))
        return (trueSM.sin2thetaEff(f) + deltaSin2thetaEff_mu());
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

double NPbase::A_f(const Particle f) const
{
    return (trueSM.A_f(f) + deltaA_f(f));
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

double NPbase::AFB(const Particle f) const
{
    return (trueSM.AFB(f) + deltaAFB(f));
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

double NPbase::R0_f(const Particle f) const
{
    return (trueSM.R0_f(f) + deltaR0_f(f));
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
