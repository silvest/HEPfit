/* 
 * File:   EWphysics.h
 * Author: mishima
 *
 * Created on February 28, 2011, 2:43 PM
 */

#ifndef EWPHYSICS_H
#define	EWPHYSICS_H

#include <gslpp_complex.h>
#include <iostream>
#include <cstring>
//#include <QCD.h>
//#include <StandardModel.h>
#include <EWPOSM.h>


class EWphysics {
public:

    /**
     * @brief EWphysics constructor
     * @param[in] gZf_i the ratio of the effective vector coupling constants @f$g_Z^f=g_V^f/g_A^f@f$ [0-9]
     * @param[in] rho_i the weak form factors [0-9]
     * @param[in] Delta_r_i the radiative corrections @f$\Delta r@f$
     * @param[in] mu
     * @param[in] mc
     * @param[in] mt
     * @param[in] md
     * @param[in] ms
     * @param[in] mb
     * @param[in] mnu1
     * @param[in] mnu2
     * @param[in] mnu3
     * @param[in] me
     * @param[in] mmu
     * @param[in] mtau
     * @param[in] mZ
     * @param[in] mHl
     * @param[in] alsMz
     * @param[in] GF
     * @param[in] ale
     * @param[in] aleMz
     * @return
     */
    EWphysics(gslpp::complex gZf_i[10], gslpp::complex rhoZf_i[10],
              double Delta_r_i,
              double mu, double mc, double mt,
              double md, double ms, double mb,
              double mnu1, double mnu2, double mnu3,
              double me, double mmu, double mtau,
              double mZ, double mHl, double alsMz, double GF, 
              double ale, double aleMz);


    /**
     * @brief EWphysics constructor
     * @param StandardModel_i reference to a StandardModel object
     */
//    EWphysics(StandardModel& StandardModel_i);

    /**
     * @brief EWphysics constructor
     * @param SUSY_i reference to a SUSY object
     */
//    EWphysics(SUSY& SUSY_i);

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
     * @brief set @f$\Delta r@f$ for INDF
     * @param[in] INDF fermion index [0-9]
     */
    void setDelta_r(double Delta_r_INDF) {
        this->Delta_r = Delta_r_INDF;
    }

    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * @param[in] flavour the flavour of the final states [nu, e, mu, tau, u, d, c, s, t, b]
     * @return the integer value associated with the flavour
     */
    int flavour_st_to_int(const std::string flavour);


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
     * @param flavour the flavour of the final states [e, mu, tau, b, c, s]
     * @return the effective weak mixing angle
     */
    double s2teff_f(const std::string flavour);

    /**
     * @param[in] flavour the flavour of the final states [e, mu, tau, b, c, s]
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_f(const std::string flavour);

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
     * @param[in] flavour_l the flavour of the final states [e, mu, tau]
     * @return pole cross section of the Z boson in GeV^-2
     */
    double sigma0_l(const std::string flavour_l);

    /**
     * @return hadronic pole cross section of the Z boson in GeV^-2
     */
    double sigma0_had();

    /**
     * @param[in] flavour_l the flavour of the final states [e, mu, tau]
     * @return Gamma_had/Gamma_l
     */
    double R0_l(const std::string flavour_l);

    /**
     * @param[in] flavour_q the flavour of the final states [b, c, s]
     * @return Gamma_q/Gamma_had
     */
    double R0_q(const std::string flavour_q);

    /**
     * @param[in] flavour the flavour of the final states [e, mu, tau, b, c, s]
     * @return the asymmetry parameter
     */
    double A_f(const std::string flavour);

    /**
     * @param flavour the flavour of the final states [e, mu, tau, b, c, s]
     * @return the forward-backward asymmetry
     */
    double AFB0_f(const std::string flavour);

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

private:
    double mu, mc, mt, md, ms, mb;
    double mnu1, mnu2, mnu3, me, mmu, mtau;
    double mZ, mHl, alsMz, GF, ale;
    double aleMz;

    gslpp::complex gZf[10], rhoZf[10]; // gZf = gVf/gAf
    double Delta_r;


    
};

#endif	/* EWPHYSICS_H */

