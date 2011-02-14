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
#include <Parameters.h>
#include <vector>

/**
 * @class StandardModel
 * @brief Standard Model Class
 */
class StandardModel {
public:
    /**
     * @brief StandardModel constructor
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
     * @param GF_i the Fermi constant
     * @param alsMz_i @f$\alpha_s(M_Z)@f$
     * @param ale_i the electromagnetic coupling
     * @param mZ_i the Z boson mass
     * @param dAle5Mz_i @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$
     */
    StandardModel(const gslpp::matrix<gslpp::complex>& VCKM_i, double mu_i,
            double md_i, double mc_i, double ms_i, double mt_i,
            double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i,
            double me_i, double mmu_i, double mtau_i,
            double mnu1_i, double mnu2_i, double mnu3_i, double GF_i,
            double alsMz_i, double ale_i, double mZ_i, double dAle5Mz_i);

    /**
     * StandardModel constructor taking as input a Parameters object
     * @param Par a Parameters object containing all SM parameters listed in the explicit SM constructor
     */
    StandardModel(Parameters& Par);

    /**
     * @brief copy constructor
     * @param orig reference to a StandardModel object
     */
    StandardModel(const StandardModel& orig);

    /**
     * @brief StandardModel destructor
     */
    virtual ~StandardModel();

    /**
     * @return the VEV
     */
    double v();

    /**
     *
     * @return sin theta_weak^2
     */
    double sin2tw() {return 0.;};

    /**
     * @return the PMNS matrix
     */
    gslpp::matrix<gslpp::complex> getUPMNS() const { return UPMNS; }

    /**
     * @brief set the PMNS matrix
     * @param UPMNS the PMNS matrix
     */
    void setUPMNS(gslpp::matrix<gslpp::complex> UPMNS) { this->UPMNS = UPMNS; }

    /**
     * @return the CKM matrix
     */
    gslpp::matrix<gslpp::complex> getVCKM() const { return VCKM; }

    /**
     * @brief set the CKM matrix
     * @param VCKM the CKM matrix
     */
    void setVCKM(gslpp::matrix<gslpp::complex> VCKM) { this->VCKM = VCKM; }

    /**
     * @return the bottom mass
     */
    double getMb() const { return mb; }

    /**
     * @brief set the bottom mass
     * @param mb the bottom mass mb(mb)
     */
    void setMb(double mb) { this->mb = mb; }

    /**
     * @return the charm mass
     */
    double getMc() const { return mc; }

    /**
     * @brief set the charm mass
     * @param mc the charm mass mc(mc)
     */
    void setMc(double mc) { this->mc = mc; }

    /**
     * @return the down mass
     */
    double getMd() const { return md; }

    /**
     * @brief set the down mass
     * @param md the down mass at 2 GeV
     */
     void setMd(double md) { this->md = md; }

    /**
     * @return the electron mass
     */
    double getMe() const { return me; }

    /**
     * @brief set the electron mass
     * @param me the electron mass
     */
     void setMe(double me) { this->me = me; }

    /**
     * @return the muon mass
     */
    double getMmu() const { return mmu; }

    /**
     * @brief set the muon mass
     * @param mmu the muon mass
     */
     void setMmu(double mmu) { this->mmu = mmu; }

    /**
     * @return the lightest neutrino mass
     */
    double getMnu1() const { return mnu1; }

    /**
     * @brief set the lightest neutrino mass
     * @param mnu1 the lightest neutrino mass
     */
     void setMnu1(double mnu1) { this->mnu1 = mnu1; }

    /**
     * @return the middle neutrino mass
     */
    double getMnu2() const { return mnu2; }

    /**
     * @brief set the middle neutrino mass
     * @param mnu2 the middle neutrino mass
     */
     void setMnu2(double mnu2) { this->mnu2 = mnu2; }

    /**
     * @return the heaviest neutrino mass
     */
    double getMnu3() const { return mnu3; }

    /**
     * @brief set the heaviest neutrino mass
     * @param mnu3 the heaviest neutrino mass
     */
     void setMnu3(double mnu3) { this->mnu3 = mnu3; }

    /**
     * @return the strange mass
     */
    double getMs() const { return ms; }

    /**
     * @brief set the strange mass
     * @param ms the strange mass at 2 GeV
     */
     void setMs(double ms) { this->ms = ms; }

    /**
     * @return the top mass
     */
    double getMt() const { return mt; }

    /**
     * @brief set the top mass
     * @param mt the top mass mt(mt)
     */
     void setMt(double mt) { this->mt = mt; }

    /**
     * @return the tau mass
     */
    double getMtau() const { return mtau; }

    /**
     * @brief set the tau mass
     * @param mtau the tau mass
     */
     void setMtau(double mtau) { this->mtau = mtau; }

    /**
     * @return the up mass
     */
    double getMu() const { return mu; }

     /**
     * @brief set the up mass
     * @param mu the up mass at 2 GeV
     */
    void setMu(double mu) { this->mu = mu; }

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
     * @return the W boson mass
     */
    double mW() const;

    /**
     *
     * @return down Yukawa matrix
     */
    gslpp::matrix<gslpp::complex> getYd() const {
        return Yd;
    }

    /**
     *
     * @return charged lepton Yukawa matrix
     */
    gslpp::matrix<gslpp::complex> getYe() const {
        return Ye;
    }

    /**
     *
     * @return neutrino Yukawa matrix
     */
    gslpp::matrix<gslpp::complex> getYn() const {
        return Yn;
    }

    /**
     *
     * @return up Yukawa matrix
     */
    gslpp::matrix<gslpp::complex> getYu() const {
        return Yu;
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


protected:
    gslpp::matrix<gslpp::complex> VCKM, UPMNS, Yd, Yu, Ye, Yn;
    double mu, md, mc, ms, mt, mb, me, mmu, mtau, mnu1, mnu2, mnu3;
    double mHl, alsMz, ale, mZ, GF, dAle5Mz;
    static const std::map<std::string,std::string[] > Deps;

};

#endif	/* STANDARDMODEL_H */
