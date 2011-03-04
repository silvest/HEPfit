/*
 * File:   StandardModel.h
 * Author: silvest
 *
 * Created on November 30, 2010, 1:27 PM
 */

#ifndef STANDARDMODEL_H
#define	STANDARDMODEL_H

#include <gslpp_complex.h>
#include <gslpp_vector_double.h>
#include <gslpp_vector_complex.h>
#include <gslpp_matrix_double.h>
#include <gslpp_matrix_complex.h>
#include <vector>
#include <Parameters.h>
#include "QCD.h"

using namespace gslpp;
/**
 * @class StandardModel
 * @brief Standard Model Class
 */
class StandardModel: public QCD {
public:
    enum lepton {NEUTRINO_1=quark(BOTTOM)+1,ELECTRON,NEUTRINO_2,MU,
    NEUTRINO_3,TAU};
    static const int NSMvars = 13;
    static const std::string SMvars[NSMvars];
///**
//     * @brief StandardModel constructor
//     * @param VCKM_i The CKM matrix
//     * @param mu_i up quark mass at 2 GeV
//     * @param md_i down quark mass at 2 GeV
//     * @param mc_i charm quark mass mc(mc)
//     * @param ms_i strange quark mass at 2 GeV
//     * @param mt_i top quark mass mt(mt)
//     * @param mb_i bottom quark mass mb(mb)
//     * @param UPMNS_i The PMNS matrix
//     * @param me_i electron mass
//     * @param mmu_i muon mass
//     * @param mtau_i tau mass
//     * @param mnu1_i lightest neutrino mass
//     * @param mnu2_i middle neutrino mass
//     * @param mnu3_i hevier neutrino mass
//     * @param GF_i the Fermi constant
//     * @param alsMz_i @f$\alpha_s(M_Z)@f$
//     * @param ale_i the electromagnetic coupling
//     * @param mZ_i the Z boson mass
//     * @param dAle5Mz_i @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$
//     * @param mHl_i the Higgs mass
//     */
//    StandardModel(const matrix<complex>& VCKM_i, double mu_i,
//            double md_i, double ms_i, double mc_i, double mb_i,
//            double mt_i, const matrix<complex>& UPMNS_i,
//            double me_i, double mmu_i, double mtau_i,
//            double mnu1_i, double mnu2_i, double mnu3_i, double GF_i,
//            double alsMz_i, double ale_i, double mZ_i, double dAle5Mz_i,
//            double mHl_i, double mu1_i, double mu2_i, double mu3_i);
//

    /**
     * updates the SM parameters found in the argument
     * @param a Parameters object containing the parameters to be updated
     */
    void update(Parameters&);

    /**
     * StandardModel constructor taking as input a Parameters object
     * @param Par a Parameters object containing all SM parameters listed in the explicit SM constructor
     */
    StandardModel(Parameters& Par);

//    /**
//     * @brief copy constructor
//     * @param orig reference to a StandardModel object
//     */
//    StandardModel(const StandardModel& orig);
//
    /**
     * @brief StandardModel destructor
     */
    virtual ~StandardModel();


    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * @return the VEV
     */
    double v() const;


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @return the CKM matrix
     */
    matrix<complex> getVCKM() const { return VCKM; }

    /**
     * @brief set the CKM matrix
     * @param VCKM the CKM matrix
     */
    void setVCKM(matrix<complex> VCKM) { this->VCKM = VCKM; }

    /**
     * @return the PMNS matrix
     */
    matrix<complex> getUPMNS() const { return UPMNS; }

    /**
     * @brief set the PMNS matrix
     * @param UPMNS the PMNS matrix
     */
    void setUPMNS(matrix<complex> UPMNS) { this->UPMNS = UPMNS; }

    /**
     *
     * @return up Yukawa matrix
     */
    matrix<complex> getYu() const {
        return Yu;
    }

    /**
     *
     * @return down Yukawa matrix
     */
    matrix<complex> getYd() const {
        return Yd;
    }

    /**
     *
     * @return neutrino Yukawa matrix
     */
    matrix<complex> getYn() const {
        return Yn;
    }

    /**
     *
     * @return charged lepton Yukawa matrix
     */
    matrix<complex> getYe() const {
        return Ye;
    }

    /**
     *
     * @return the Fermi constant
     */
    double getGF() const {
        return GF;
    }

    /**
     * @brief set the Fermi constant
     * @param GF the Fermi constant
     */
    void setGF(double GF) {
        this->GF = GF;
    }

   /**
     *
     * @return @f$\alpha_s(M_Z)@f$
     */
    double getAlsMz() const {
        return alsMz;
    }

    /**
     * set @f$\alpha_s(M_Z)@f$
     * @param alsMz @f$\alpha_s(M_Z)@f$
     */
    void setAlsMz(double alsMz) {
        this->alsMz = alsMz;
    }

    /**
     *
     * @return the electromagnetic coupling
     */
    double getAle() const {
        return ale;
    }

    /**
     * @brief set the electromagnetic coupling
     * @param ale the electromagnetic coupling
     */
    void setale(double ale) {
        this->ale = ale;
    }

    /**
     *
     * @return @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$
     */
    double getDAle5Mz() const {
        return dAle5Mz;
    }

    /**
     * set @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$
     * @param dAle5Mz @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$
     */
    void setDAle5Mz(double dAle5Mz) {
        this->dAle5Mz = dAle5Mz;
    }

    /**
     *
     * @return the Z boson mass
     */
    double getMZ() const {
        return mZ;
    }

    /**
     * @brief set the Z boson mass
     * @param mZ the Z boson mass
     */
    void setMZ(double mZ) {
        this->mZ = mZ;
    }

    /**
     *
     * @return the Higgs mass
     */
    double getMHl() const {
        return mHl;
    }

    /**
     * @brief set the Higgs mass
     * @param mHl the Higgs mass
     */
    void setMHl(double mHl) {
        this->mHl = mHl;
    }


    ///////////////////////////////////////////////////////////////////////////
    /* Functions for EW precision observables */
    
    /**
     * @return leptonic correction to alpha at mZ at the three-loop level
     */
    double dAleLepMz() const;

    /**
     * @return top-quark corection to alpha at mZ up to second order in alpha_s
     */
    double dAleTopMz() const;

    /**
     * @return total correction to alpha at mZ
     */
    double dAleTotalMz() const;

    /**
     * @return the electromagnetic coupling alpha at mZ
     */
    double aleMz() const;

    /**
     * @return the charm-quak mass at mZ, mc(mZ)
     */
    double mcMz() const;

    /**
     * @return the bottom-quak mass at mZ, mb(mZ)
     */
    double mbMz() const;

    /**
     * @return the W boson mass, including radiative corrections
     */
    virtual double mW() const;

    /**
     * @param[in] INDF fermion index [0-9] (see EWphysics::flavour_st_to_int())
     * @return the ratio of the effective vector coupling constants @f$g_Z^f=g_V^f/g_A^f@f$ for INDF
     */
    virtual gslpp::complex gZf(const int INDF) const; // gZf = gVf/gAf

    /**
     * @param[in] INDF fermion index [0-9] (see EWphysics::flavour_st_to_int())
     * @return the weak form factor for INDF
     */
    virtual gslpp::complex rhoZf(const int INDF) const;

    /**
     * @return the radiative-correction factor @f$\Delta r@f$
     */
    virtual double Delta_r() const;


    ////////////////////////////////////////////////////////////////////////
    
    // Angles
    double getBeta() const;
    double getGamma() const;
    double getAlpha() const;
    double getBetas() const;

    // Lambda_q
    gslpp::complex getlamt() const;
    gslpp::complex getlamc() const;
    gslpp::complex getlamu() const;

    gslpp::complex getlamt_d() const;
    gslpp::complex getlamc_d() const;
    gslpp::complex getlamu_d() const;

    gslpp::complex getlamt_s() const;
    gslpp::complex getlamc_s() const;
    gslpp::complex getlamu_s() const;

    // Sides
    double getRt() const;
    double getRts() const;
    double getRb() const;

    /**
     * get the @f$\Delta B=\Delta D=2@f$ amplitude
     * @return @f$\langle \bar B_d \vert \mathcal{H}_\mathrm{eff}\vert B_d\rangle@f$ //CHECK!!
     */    
    gslpp::complex getDBD2Amplitude(const int LE) const;


    ///////////////////////////////////////////////////////////////////////////

protected:
    matrix<complex> VCKM, UPMNS, Yu, Yd, Yn, Ye;
    double GF, alsMz, ale, dAle5Mz, mZ, mHl;
    double muw;
    Particle particles[lepton(TAU)+1];
//    static const std::vector<std::string> pino;
//    mutable std::map<std::string,double> Hashes;
//    mutable std::map<std::string,double> DValues;

private:
    double eta2bbar(const int) const;
    double kt_sing_a(const double x) const;
    double kt_sing(const double x) const;
    double kt_oct(const double x) const;
    double kt(const double x, const double mu) const;
    double S(double, double) const;
//    void init(double mu_i,
//            double md_i, double ms_i, double mc_i, double mb_i,
//            double mt_i,
//            double me_i, double mmu_i, double mtau_i,
//            double mnu1_i, double mnu2_i, double mnu3_i, double GF_i,
//            double alsMz_i, double ale_i, double mZ_i, double dAle5Mz_i,
//            double mHl_i, double mu1_i, double mu2_i, double mu3_i);
};

#endif	/* STANDARDMODEL_H */
