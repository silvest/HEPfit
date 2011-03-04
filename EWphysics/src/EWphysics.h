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
 *    - Gamma_W(), if necessary
 *
 */

#include <cstring>
#include <gslpp_complex.h>
#include <QCD.h>
#include <StandardModel.h>
//#include <SUSY.h>
//#include <MFV.h>


class EWphysics {
public:

    /** 
     * @brief EWphysics constructor
     * @param[in] StandardModel_i pointer to an object of StandardModel or its extension
     */
    EWphysics(StandardModel* StandardModel_i);

    /**
     * @brief EWphysics constructor
     * @param[in] gZf_i the ratio of the effective vector coupling constants g_Z^f=g_V^f/g_A^f [0-9]
     * @param[in] rho_i the weak form factors [0-9]
     * @param[in] Delta_r_i the radiative corrections Delta r
     * @param[in] mcMz_i the charm quark mass at mZ, mc(mZ)
     * @param[in] mbMz_i the bottom quark mass at mZ, mb(mZ)
     * @param[in] mt_i the top quark mass mt(mt)
     * @param[in] me_i the electron mass
     * @param[in] mmu_i the muon mass
     * @param[in] mtau_i the tau mass
     * @param[in] mZ_i the Z boson mass
     * @param[in] alsMz_i the strong coupling constant at mZ, alpha_s(M_Z^2)
     * @param[in] GF_i the Fermi constant
     * @param[in] ale_i the electromagnetic coupling alpha(0)
     * @param[in] aleMz_i the electromagnetic coupling at mZ, alpha(mZ^2)
     */
    EWphysics(gslpp::complex gZf_i[10], gslpp::complex rhoZf_i[10],
              double Delta_r_i,
              double mcMz_i, double mbMz_i, double mt, 
              double me_i, double mmu_i, double mtau_i,
              double mZ_i, double alsMz_i, double GF_i,
              double ale_i, double aleMz_i);

    /**
     * @brief EWphysics copy constructor
     * @param[in] orig reference to an EWphysics object
     */
    //EWphysics(const EWphysics& orig);

    /**
     * @brief EWphysics destructor
     */
    virtual ~EWphysics();

    
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
    StandardModel *MyModel;

    gslpp::complex gZf[10];   // gZf = gVf/gAf
    gslpp::complex rhoZf[10];
    double Delta_r;

    double mcMz, mbMz; /* charm and bottom quak masses at mZ */
    double mt;         /* top quark mass mt(mt) */
    double me, mmu, mtau; /* charged lepton masses */
    double mZ, alsMz, GF, ale;
    double aleMz; /* alpha at mZ */
   
};

#endif	/* EWPHYSICS_H */

