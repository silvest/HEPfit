/* 
 * File:   MFV.h
 * Author: silvest
 *
 * Created on September 24, 2010, 10:53 AM
 */

#ifndef MFV_H
#define	MFV_H

#include <StandardModel.h>

/**
 * @class MFV
 * @brief Minimal Flavour Violating SUSY model
 */
class MFV: public StandardModel {
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
     * @param SM_i reference to a StandardModel object
     */
    MFV(double mQtilde_i, double mUtilde_i, double mDtilde_i,
            double Au_i, double Ad_i, double mLtilde_i, double mEtilde_i,
            double mNtilde_i, double Ae_i, double An_i, double m1_i,
            double m2_i, double m3_i, double muH_i, double tanb_i, double mHp_i,
            const gslpp::matrix<gslpp::complex>& VCKM_i,
            double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
            double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i,
            double me_i, double mmu_i, double mtau_i, double mnu1_i,
            double mnu2_i, double mnu3_i);

    /**
     * @brief MFV copy contructor
     * @param orig reference to a MFV object
     */
    MFV(const MFV& orig);

    /**
     * @brief MFV destructor
     */
    virtual ~MFV();

    /**
     * 
     * @return universal trilinear down-type squark coupling
     */
    double GetAd() const { return Ad; }

    /**
     * @brief set universal trilinear down-type squark coupling
     * @param Ad universal trilinear down-type squark coupling
     */
    void SetAd(double Ad) { this->Ad = Ad; }

    /**
     *
     * @return universal trilinear charged slepton coupling
     */
    double GetAe() const { return Ae; }

    /**
     * @brief set universal trilinear charged slepton coupling
     * @param Ae universal trilinear charged slepton coupling
     */
    void SetAe(double Ae) { this->Ae = Ae; }

    /**
     * 
     * @return universal trilinear sneutrino coupling
     */
    double GetAn() const { return An; }

    /**
     * @brief set universal trilinear sneutrino coupling
     * @param An universal trilinear sneutrino coupling
     */
    void SetAn(double An) { this->An = An; }

    /**
     *
     * @return universal trilinear up-type squark coupling
     */
    double GetAu() const { return Au; }

    /**
     * @brief set universal trilinear up-type squark coupling
     * @param Au universal trilinear up-type squark coupling
     */
    void SetAu(double Au) { this->Au = Au; }

    /**
     * 
     * @return bino soft DRbar mass @f$m_{1}(m_{1})@f$
     */
    double GetM1() const { return m1; }

    /**
     * @brief set bino soft DRbar mass @f$m_{1}(m_{1})@f$
     * @param m1 bino soft DRbar mass @f$m_{1}(m_{1})@f$
     */
    void SetM1(double m1) { this->m1 = m1; }

    /**
     *
     * @return wino soft DRbar mass @f$m_{2}(m_{2})@f$
     */
    double GetM2() const { return m2; }

    /**
     * @brief set wino soft DRbar mass @f$m_{2}(m_{2})@f$
     * @param m2 wino soft DRbar mass @f$m_{2}(m_{2})@f$
     */
    void SetM2(double m2) { this->m2 = m2; }

    /**
     * 
     * @return gluino soft DRbar mass @f$m_{3}(m_{3})@f$
     */
    double GetM3() const { return m3; }

    /**
     * @brief set gluino soft DRbar mass @f$m_{3}(m_{3})@f$
     * @param m3 gluino soft DRbar mass @f$m_{3}(m_{3})@f$
     */
    void SetM3(double m3) { this->m3 = m3; }

    /**
     * 
     * @return right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     */
    double GetMDtilde() const { return mDtilde; }

    /**
     * @brief set right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     * @param mDtilde right-handed down-type universal squark soft DRbar mass @f$\tilde{m}_{D}(\tilde{m}_{D})@f$
     */
    void SetMDtilde(double mDtilde) { this->mDtilde = mDtilde; }

    /**
     * 
     * @return right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     */
    double GetMEtilde() const { return mEtilde; }

    /**
     * @brief set right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     * @param mEtilde right-handed charged slepton universal soft DRbar mass @f$\tilde{m}_{E}(\tilde{m}_{E})@f$
     */
    void SetMEtilde(double mEtilde) { this->mEtilde = mEtilde; }

    /**
     * 
     * @return charged Higgs mass @f$m_{H^+}@f$
     */
    double GetMHp() const { return mHp; }

    /**
     * @brief set charged Higgs mass @f$m_{H^+}@f$
     * @param mHp charged Higgs mass @f$m_{H^+}@f$
     */
    void SetMHp(double mHp) { this->mHp = mHp; }

    /**
     *
     * @return slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     */
    double GetMLtilde() const { return mLtilde; }

    /**
     * @brief set slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     * @param mLtilde slepton doublet universal soft DRbar mass @f$\tilde{m}_{L}(\tilde{m}_{L})@f$
     */
    void SetMLtilde(double mLtilde) { this->mLtilde = mLtilde; }

    /**
     *
     * @return right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     */
    double GetMNtilde() const { return mNtilde; }

    /**
     * @brief set right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     * @param mNtilde right-handed sneutrino universal soft DRbar mass @f$\tilde{m}_{N}(\tilde{m}_{N})@f$
     */
    void SetMNtilde(double mNtilde) { this->mNtilde = mNtilde; }

    /**
     * 
     * @return squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     */
    double GetMQtilde() const { return mQtilde; }

    /**
     * @brief set squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     * @param mQtilde squark doublet universal soft DRbar mass @f$\tilde{m}_{Q}(\tilde{m}_{Q})@f$
     */
    void SetMQtilde(double mQtilde) { this->mQtilde = mQtilde; }

    /**
     *
     * @return mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     */
    double GetMUtilde() const { return mUtilde; }

    /**
     * @brief set mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     * @param mUtilde mUtilde_i right-handed up-type universal squark soft DRbar mass @f$\tilde{m}_{U}(\tilde{m}_{U})@f$
     */
    void SetMUtilde(double mUtilde) { this->mUtilde = mUtilde; }

    /**
     *
     * @return superpotential @f$\mu@f$ term
     */
    double GetMuH() const { return muH; }

    /**
     * @brief set superpotential @f$\mu@f$ term
     * @param muH superpotential @f$\mu@f$ term
     */
    void SetMuH(double muH) { this->muH = muH; }

    /**
     *
     * @return @f$\tan \beta @f$
     */
    double GetTanb() const { return tanb; }

    /**
     * @brief set @f$\tan \beta @f$
     * @param tanb @f$\tan \beta @f$
     */
    void SetTanb(double tanb) { this->tanb = tanb; }

private:
    double mQtilde, mUtilde, mDtilde, Au,  Ad,  mLtilde, mEtilde, mNtilde,
           Ae, An, m1, m2, m3,  muH, tanb, mHp;

};

#endif	/* MFV_H */

