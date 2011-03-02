/* 
 * File:   EWphysics.h
 * Author: mishima
 *
 * Created on February 28, 2011, 2:43 PM
 */

#ifndef EWPHYSICS_H
#define	EWPHYSICS_H

/*
 *  To do:
 *    - Gamma_W()
 *    - checks and tests -- errors in Gamma_q()
 *    - constructors with StandardModel, SUSY and MFV objects
 *
 */

#include <cstring>
#include <gslpp_complex.h>
#include <QCD.h>
//#include <StandardModel.h>
//#include <SUSY.h>
//#include <MFV.h>


class EWphysics {
public:

    /**
     * @brief EWphysics constructor
     * @param[in] gZf_i the ratio of the effective vector coupling constants g_Z^f=g_V^f/g_A^f [0-9]
     * @param[in] rho_i the weak form factors [0-9]
     * @param[in] Delta_r_i the radiative corrections Delta r
     * @param[in] mu_i the up quark mass at 2 GeV
     * @param[in] md_i the down quark mass at 2 GeV
     * @param[in] mc_i the charm quark mass mc(mc)
     * @param[in] ms_i the strange quark mass at 2 GeV
     * @param[in] mt_i the top quark mass mt(mt)
     * @param[in] mb_i the bottom quark mass mb(mb)
     * @param[in] me_i the electron mass
     * @param[in] mmu_i the muon mass
     * @param[in] mtau_i the tau mass
     * @param[in] mZ_i the Z boson mass
     * @param[in] mHl_i the Higgs mass
     * @param[in] alsMz_i the strong coupling constant alpha_s(M_Z^2)
     * @param[in] GF_i the Fermi constant
     * @param[in] ale_i the electromagnetic coupling at alpha(0)
     * @param[in] aleMz_i the electromagnetic coupling at alpha(mZ^2)
     */
    EWphysics(gslpp::complex gZf_i[10], gslpp::complex rhoZf_i[10],
              double Delta_r_i,
              double mu_i, double md_i, double mc_i,
              double ms_i, double mt_i, double mb_i,
              double me_i, double mmu_i, double mtau_i,
              double mZ_i, double mHl_i, double alsMz_i, double GF_i,
              double ale_i, double aleMz_i);

    /**
     * @brief EWphysics constructor
     * @param[in] StandardModel_i reference to a StandardModel object
     */
    //EWphysics(StandardModel& StandardModel_i);

    /**
     * @brief EWphysics constructor
     * @param[in] SUSY_i reference to a SUSY object
     */
    //EWphysics(SUSY& SUSY_i);

    /**
     * @brief EWphysics constructor
     * @param[in] MFV_i reference to an MFV object
     */
    //EWphysics(MFV& MFV_i);

    /**
     * @brief EWphysics copy constructor
     * @param[in] orig reference to an EWphysics object
     */
    EWphysics(const EWphysics& orig);

    /**
     * @brief EWphysics destructor
     */
    virtual ~EWphysics();

    
    ///////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the ratio of the effective vector coupling constants @f$g_Z^f=g_V^f/g_A^f@f$ for INDF
     */
    gslpp::complex getGZf(int INDF) const {
        return gZf[INDF];
    }

    /**
     * @brief set the ratio of the effective coupling constants @f$g_Z^f=g_V^f/g_A^f@f$ for INDF
     * @param[in] INDF fermion index [0-9]
     * @param[in] gZf_INDF the ratio of the effective coupling constants for INDF
     */
    void setGZf(int INDF, gslpp::complex gZf_INDF) {
        this->gZf[INDF] = gZf_INDF;
    }

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the weak form factor for INDF
     */
    gslpp::complex getRhoZf(int INDF) const {
        return rhoZf[INDF];
    }

    /**
     * @brief set the weak form factor for INDF
     * @param[in] INDF fermion index [0-9]
     * @param[in] rhoZf_INDF the weak form factor for INDF
     */
    void setRhoZf(int INDF, gslpp::complex rhoZf_INDF) {
        this->rhoZf[INDF] = rhoZf_INDF;
    }

    /**
     * @return @f$\Delta r@f$
     */
    double getDelta_r() const {
        return Delta_r;
    }

    /**
     * @brief set @f$\Delta r@f$
     * @param[in] Delta_r the radiative-correction factor @f$\Delta r@f$
     */
    void setDelta_r(double Delta_r) {
        this->Delta_r = Delta_r;
    }

    /**
     * @return the charm quak masses at mZ, mc(mZ)
     */
    double getMcMz() const {
        return mcMz;
    }

    /**
     * @brief set the charm quak masses at mZ, mc(mZ)
     * @param[in] mcMz the charm quak masses at mZ, mc(mZ)
     */
    void setMcMz(double mcMz) {
        this->mcMz = mcMz;
    }

    /**
     * @return the bottom quak masses at mZ, mb(mZ)
     */
    double getMbMz() const {
        return mbMz;
    }

    /**
     * @brief set the bottom quak masses at mZ, mb(mZ)
     * @param[in] mbMz the bottom quak masses at mZ, mb(mZ)
     */
    void setMbMz(double mbMz) {
        this->mbMz = mbMz;
    }

    /**
     * @return the top quark mass mt(mt)
     */
    double getMt() const {
        return mt;
    }

    /**
     * @brief set the top quark mass mt(mt)
     * @param[in] mt the top quark mass mt(mt)
     */
    void setMt(double mt) {
        this->mt = mt;
    }

    /**
     * @return the electron mass
     */
    double getMe() const {
        return me;
    }

    /**
     * @brief set the electron mass
     * @return the electron mass
     */
    void setMe(double me) {
        this->me = me;
    }

    /**
     * @return the muon mass
     */
    double getMmu() const {
        return mmu;
    }

    /**
     * @brief set the muon mass
     * @param[in] mmu the muon mass
     */
    void setMmu(double mmu) {
        this->mmu = mmu;
    }

    /**
     * @return the tau mass
     */
    double getMtau() const {
        return mtau;
    }

    /**
     * @brief set the tau mass
     * @param[in] mtau the tau mass
     */
    void setMtau(double mtau) {
        this->mtau = mtau;
    }

    /**
     * @return the Z boson mass
     */
    double getMZ() const {
        return mZ;
    }

    /**
     * @brief set the Z boson mass
     * @param[in] mZ the Z boson mass
     */
    void setMZ(double mZ) {
        this->mZ = mZ;
    }

    /**
     * @return the Higgs mass
     */
    double getMHl() const {
        return mHl;
    }

    /**
     * @brief set the Higgs mass
     * @param[in] mHl the Higgs mass
     */
    void setMHl(double mHl) {
        this->mHl = mHl;
    }

    /**
     * @return the strong coupling constant alpha_s(M_Z^2)
     */
    double getAleMz() const {
        return aleMz;
    }

    /**
     * @brief set the strong coupling constant alpha_s(M_Z^2)
     * @param[in] aleMz the strong coupling constant alpha_s(M_Z^2)
     */
    void setAleMz(double aleMz) {
        this->aleMz = aleMz;
    }

    /**
     * @return the Fermi constant
     */
    double getGF() const {
        return GF;
    }

    /**
     * @brief set the Fermi constant
     * @param[in] GF the Fermi constant
     */
    void setGF(double GF) {
        this->GF = GF;
    }

    /**
     * @return the electromagnetic coupling at alpha(0)
     */
    double getAle() const {
        return ale;
    }

    /**
     * @brief set the electromagnetic coupling at alpha(0)
     * @param[in] ale the electromagnetic coupling at alpha(0)
     */
    void setAle(double ale) {
        this->ale = ale;
    }

    /**
     * @return the electromagnetic coupling at alpha(mZ^2)
     */
    double getAlsMz() const {
        return alsMz;
    }

    /**
     * @brief set the electromagnetic coupling at alpha(mZ^2)
     * @param[in] alsMz the electromagnetic coupling at alpha(mZ^2)
     */
    void setAlsMz(double alsMz) {
        this->alsMz = alsMz;
    }


    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * @param[in] flavour the flavour of the final states [nu, e, mu, tau, u, d, c, s, t, b]
     * @return the integer value associated with the flavour
     */
    int flavour_st_to_int(const std::string flavour);

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the string associated with the integer INDF
     */
    std::string flavour_int_to_st(const int INDF);

    /**
     * @param[in] INDF fermion index [0-9]
     * @return electric charge for INDF
     */
    double Qf(const int INDF);

    
    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * @return the mass of the W boson in GeV
     */
    double mW();

    /**
     * @return the total width of the W boson in GeV
     */
    double Gamma_W();

    /**
     * @return the weak mixing angle
     */
    double sw2();

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the effective weak mixing angle
     */
    double s2teff_f(const int INDF);

    /**
     * @param[in] INDF_l fermion index [0-3]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_l(const int INDF_l);

    /**
     * @param[in] INDF_q fermion index [4-9]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_q(const int INDF_q);

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_f(const int INDF);

    /**
     * @return the invisible width of the Z boson in GeV
     */
    double Gamma_inv();

    /**
     * @return the hadronic width of the Z boson in GeV
     */
    double Gamma_had();

    /**
     * @return the total width of the Z boson in GeV
     */
    double Gamma_Z();

    /**
     * @param[in] INDF_l fermion index [0-3]
     * @return pole cross section of the Z boson in GeV^-2
     */
    double sigma0_l(const int INDF_l);

    /**
     * @return hadronic pole cross section of the Z boson in GeV^-2
     */
    double sigma0_had();

    /**
     * @param[in] INDF_l fermion index [0-3]
     * @return Gamma_had/Gamma_l
     */
    double R0_l(const int INDF_l);

    /**
     * @param[in] INDF_q fermion index [4-9]
     * @return Gamma_q/Gamma_had
     */
    double R0_q(const int INDF_q);

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the asymmetry parameter
     */
    double A_f(const int INDF);

    /**
     * @param[in] INDF fermion index [0-9]
     * @return the forward-backward asymmetry
     */
    double AFB0_f(const int INDF);

    /**
     * @return the oblique parameter epsilon_1
     */
    double obliqueEpsilon1();

    /**
     * @return the oblique parameter epsilon_2
     */
    double obliqueEpsilon2();

    /**
     * @return the oblique parameter epsilon_3
     */
    double obliqueEpsilon3();

    /**
     * @return the oblique parameter S
     */
    double obliqueS();

    /**
     * @return the oblique parameter T
     */
    double obliqueT();

    /**
     * @return the oblique parameter U
     */
    double obliqueU();


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] flavour the flavour of the final states [nu, e, mu, tau, u, d, c, s, t, b]
     * @return the effective weak mixing angle
     */
    double s2teff_f(const std::string flavour);

    /**
     * @param[in] flavour_l the flavour of the final states [nu, e, mu, tau]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_l(const std::string flavour_l);

    /**
     * @param[in] flavour_q the flavour of the final states [u, d, c, s, t, b]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_q(const std::string flavour_q);

    /**
     * @param[in] flavour the flavour of the final states [nu, e, mu, tau, u, d, c, s, t, b]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_f(const std::string flavour);

    /**
     * @param[in] flavour_l the flavour of the final states [nu, e, mu, tau]
     * @return pole cross section of the Z boson in GeV^-2
     */
    double sigma0_l(const std::string flavour_l);

    /**
     * @param[in] flavour_l the flavour of the final states [nu, e, mu, tau]
     * @return Gamma_had/Gamma_l
     */
    double R0_l(const std::string flavour_l);

    /**
     * @param[in] flavour_q the flavour of the final states [u, d, c, s, t, b]
     * @return Gamma_q/Gamma_had
     */
    double R0_q(const std::string flavour_q);

    /**
     * @param[in] flavour the flavour of the final states [nu, e, mu, tau, u, d, c, s, t, b]
     * @return the asymmetry parameter
     */
    double A_f(const std::string flavour);

    /**
     * @param[in] flavour the flavour of the final states [nu, e, mu, tau, u, d, c, s, t, b]
     * @return the forward-backward asymmetry
     */
    double AFB0_f(const std::string flavour);


    ///////////////////////////////////////////////////////////////////////////

private:
    gslpp::complex gZf[10];   // gZf = gVf/gAf
    gslpp::complex rhoZf[10];
    double Delta_r; 

    double mcMz, mbMz; /* charm and bottom quak masses at mZ */
    double mt;         /* top quark mass mt(mt) */
    double me, mmu, mtau;
    double mZ, mHl, alsMz, GF, ale;
    double aleMz;      /* alpha at mZ */
   
};

#endif	/* EWPHYSICS_H */

