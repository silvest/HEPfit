/* 
 * File:   MFV.h
 * Author: silvest
 *
 * Created on September 24, 2010, 10:53 AM
 */

#ifndef MFV_H
#define	MFV_H

#include <SUSY.h>

/**
 * @class MFV
 * @brief Minimal Flavour Violating SUSY model
 */
class MFV : public SUSY {
public:
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
    MFV(const gslpp::matrix<gslpp::complex>& VCKM_i,
            double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
            double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i,
            double me_i, double mmu_i, double mtau_i, double mnu1_i,
            double mnu2_i, double mnu3_i, double tanb_i, gslpp::complex muH_i);

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
     * @param SM_i reference to a StandardModel object
     */
   MFV(const SUSY& SUSY_i);

    /**
     * @brief MFV copy contructor
     * @param orig reference to a MFV object
     */
    MFV(const MFV& orig);

    /**
     * @brief MFV destructor
     */
    virtual ~MFV();

    void setParameters(double a1, double a2, double a3, gslpp::complex a4,
        gslpp::complex a5, double b1, double b2, gslpp::complex b3, double b5,
        double b6, gslpp::complex b7, gslpp::complex b8, double c4,
        gslpp::complex c7, gslpp::complex c10, gslpp::complex c11);
    /**
     * 
     * @return universal trilinear down-type squark coupling
     */
    double getAd() const { return Ad; }

    /**
     * @brief set universal trilinear down-type squark coupling
     * @param Ad universal trilinear down-type squark coupling
     */
    void setAd(double Ad) { this->Ad = Ad; }

    /**
     *
     * @return universal trilinear charged slepton coupling
     */
    double getAe() const { return Ae; }

    /**
     * @brief set universal trilinear charged slepton coupling
     * @param Ae universal trilinear charged slepton coupling
     */
    void setAe(double Ae) { this->Ae = Ae; }

    /**
     * 
     * @return universal trilinear sneutrino coupling
     */
    double getAn() const { return An; }

    /**
     * @brief set universal trilinear sneutrino coupling
     * @param An universal trilinear sneutrino coupling
     */
    void setAn(double An) { this->An = An; }

    /**
     *
     * @return universal trilinear up-type squark coupling
     */
    double getAu() const { return Au; }

    /**
     * @brief set universal trilinear up-type squark coupling
     * @param Au universal trilinear up-type squark coupling
     */
    void setAu(double Au) { this->Au = Au; }

    /**
     * 
     * @return right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     */
    double getMDtilde() const { return mDtilde; }

    /**
     * @brief set right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     * @param mDtilde right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     */
    void setMDtilde(double mDtilde) { this->mDtilde = mDtilde; }

    /**
     * 
     * @return right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     */
    double getMEtilde() const { return mEtilde; }

    /**
     * @brief set right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     * @param mEtilde right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     */
    void setMEtilde(double mEtilde) { this->mEtilde = mEtilde; }

    /**
     *
     * @return slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     */
    double getMLtilde() const { return mLtilde; }

    /**
     * @brief set slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     * @param mLtilde slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     */
    void setMLtilde(double mLtilde) { this->mLtilde = mLtilde; }

    /**
     *
     * @return right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     */
    double getMNtilde() const { return mNtilde; }

    /**
     * @brief set right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     * @param mNtilde right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     */
    void setMNtilde(double mNtilde) { this->mNtilde = mNtilde; }

    /**
     * 
     * @return squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     */
    double getMQtilde() const { return mQtilde; }

    /**
     * @brief set squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     * @param mQtilde squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     */
    void setMQtilde(double mQtilde) { this->mQtilde = mQtilde; }

    /**
     *
     * @return mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     */
    double getMUtilde() const { return mUtilde; }

    /**
     * @brief set mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     * @param mUtilde mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     */
    void setMUtilde(double mUtilde) { this->mUtilde = mUtilde; }

private:
    double mQtilde, mUtilde, mDtilde, Au, Ad, mLtilde, mEtilde, mNtilde, Ae, An;

};

#endif	/* MFV_H */
