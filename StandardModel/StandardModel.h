/* 
 * File:   StandardModel.h
 * Author: silvest
 *
 * Created on November 30, 2010, 1:27 PM
 */

#ifndef STANDARDMODEL_H
#define	STANDARDMODEL_H

#include <gslpp_matrix_complex.h>

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
     */
    StandardModel(const gslpp::matrix<gslpp::complex>& VCKM_i, double mu_i,
            double md_i, double mc_i, double ms_i, double mt_i,
            double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i,
            double me_i, double mmu_i, double mtau_i,
            double mnu1_i, double mnu2_i, double mnu3_i);

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
     * @return the PMNS matrix
     */
    gslpp::matrix<gslpp::complex> GetUPMNS() const { return UPMNS; }

    /**
     * @brief set the PMNS matrix
     * @param UPMNS the PMNS matrix
     */
    void SetUPMNS(gslpp::matrix<gslpp::complex> UPMNS) { this->UPMNS = UPMNS; }

    /**
     * @return the CKM matrix
     */
    gslpp::matrix<gslpp::complex> GetVCKM() const { return VCKM; }

    /**
     * @brief set the CKM matrix
     * @param VCKM the CKM matrix
     */
    void SetVCKM(gslpp::matrix<gslpp::complex> VCKM) { this->VCKM = VCKM; }

    /**
     * @return the bottom mass
     */
    double GetMb() const { return mb; }

    /**
     * @brief set the bottom mass
     * @param mb the bottom mass mb(mb)
     */
    void SetMb(double mb) { this->mb = mb; }

    /**
     * @return the charm mass
     */
    double GetMc() const { return mc; }

    /**
     * @brief set the charm mass
     * @param mc the charm mass mc(mc)
     */
    void SetMc(double mc) { this->mc = mc; }

    /**
     * @return the down mass
     */
    double GetMd() const { return md; }

    /**
     * @brief set the down mass
     * @param md the down mass at 2 GeV
     */
     void SetMd(double md) { this->md = md; }

    /**
     * @return the electron mass
     */
    double GetMe() const { return me; }

    /**
     * @brief set the electron mass
     * @param me the electron mass
     */
     void SetMe(double me) { this->me = me; }

    /**
     * @return the muon mass
     */
    double GetMmu() const { return mmu; }

    /**
     * @brief set the muon mass
     * @param mmu the muon mass
     */
     void SetMmu(double mmu) { this->mmu = mmu; }

    /**
     * @return the lightest neutrino mass
     */
    double GetMnu1() const { return mnu1; }

    /**
     * @brief set the lightest neutrino mass
     * @param mnu1 the lightest neutrino mass
     */
     void SetMnu1(double mnu1) { this->mnu1 = mnu1; }

    /**
     * @return the middle neutrino mass
     */
    double GetMnu2() const { return mnu2; }

    /**
     * @brief set the middle neutrino mass
     * @param mnu2 the middle neutrino mass
     */
     void SetMnu2(double mnu2) { this->mnu2 = mnu2; }

    /**
     * @return the heaviest neutrino mass
     */
    double GetMnu3() const { return mnu3; }

    /**
     * @brief set the heaviest neutrino mass
     * @param mnu3 the heaviest neutrino mass
     */
     void SetMnu3(double mnu3) { this->mnu3 = mnu3; }

    /**
     * @return the strange mass
     */
    double GetMs() const { return ms; }

    /**
     * @brief set the strange mass
     * @param ms the strange mass at 2 GeV
     */
     void SetMs(double ms) { this->ms = ms; }

    /**
     * @return the top mass
     */
    double GetMt() const { return mt; }

    /**
     * @brief set the top mass
     * @param mt the top mass mt(mt)
     */
     void SetMt(double mt) { this->mt = mt; }

    /**
     * @return the tau mass
     */
    double GetMtau() const { return mtau; }

    /**
     * @brief set the tau mass
     * @param mtau the tau mass
     */
     void SetMtau(double mtau) { this->mtau = mtau; }

    /**
     * @return the up mass
     */
    double GetMu() const { return mu; }

     /**
     * @brief set the up mass
     * @param mu the up mass at 2 GeV
     */
    void SetMu(double mu) { this->mu = mu; }

private:
    gslpp::matrix<gslpp::complex> VCKM,UPMNS;
    double mu, md, mc, ms, mt, mb, me, mmu, mtau, mnu1, mnu2, mnu3;

};

#endif	/* STANDARDMODEL_H */
