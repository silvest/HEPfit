/* 
 * File:   SUSY.h
 * Author: marco
 *
 * Created on December 2, 2010, 3:32 PM
 */

#ifndef SUSY_H
#define	SUSY_H

#include <StandardModel.h>

/**
 * @class SUSY
 * @brief generic SUSY model
 */
class SUSY: public StandardModel {
public:
    /**
     * @brief SUSY constructor
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
//    SUSY(const gslpp::matrix<gslpp::complex>& VCKM_i,
//        double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
//        double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i,
//        double me_i, double mmu_i, double mtau_i, double mnu1_i,
//        double mnu2_i, double mnu3_i, double tanb_i, gslpp::complex muH_i);
    /**
     * @brief SUSY constructor
     * @param StandardModel_i reference to a StandardModel object
     */
    SUSY(const StandardModel& StandardModel_i, double tanb_i,
        gslpp::complex muH_i);
    /**
     * @brief SUSY constructor
     * @param orig reference to a SUSY object
     */
//    SUSY(const SUSY& orig);
    /**
     * @brief SUSY destructor
     */
    virtual ~SUSY();


    ///////////////////////////////////////////////////////////////////////////

    /**
     * 
     * @return the down-type VEV
     */
    double v1();

    /**
     *
     * @return the up-type VEV
     */
    double v2();


    ///////////////////////////////////////////////////////////////////////////

    /**
     * 
     * @return gluino mass
     */
    double getM3() const {
        return m3;
    }

    /**
     * @brief set the gluino mass
     * @param m3 gluino mass
     */
    void setM3(double m3) {
        this->m3 = m3;
    }

    /**
     *
     * @return pseudoscalar Higgs mass
     */
    double getMHa() const {
        return mHa;
    }

    /**
     * @brief set the pseudoscalar Higgs mass
     * @param mA pseudoscalar Higgs mass
     */
    void setMHa(double mHa) {
        this->mHa = mHa;
    }

    /**
     *
     * @return heavy Higgs mass
     */
    double getMHh() const {
        return mHh;
    }

    /**
     * @brief set the heavy Higgs mass
     * @param mH heavy Higgs mass
     */
    void setMHh(double mHh) {
        this->mHh = mHh;
    }

    /**
     *
     * @return charged Higgs mass
     */
    double getMHp() const {
        return mHp;
    }

    /**
     * @brief set the charged Higgs mass
     * @param mHp charged Higgs mass
     */
    void setMHp(double mHp) {
        this->mHp = mHp;
    }

    /**
     *
     * @return tan beta
     */
    double getTanb() const {
        return tanb;
    }

    /**
     * @brief set tan beta, sin beta and cos beta
     * @param tanb tan beta
     */
    void setTanb(double tanb);

    /**
     *
     * @return sin beta
     */
    double getSinb() const {
        return sinb;
    }

    /**
     * @brief set tan beta, sin beta and cos beta
     * @param sinb sin beta
     */
    void setSinb(double sinb);

    /**
     *
     * @return cos beta
     */
    double getCosb() const {
        return cosb;
    }

    /**
     * @brief set tan beta, sin beta and cos beta
     * @param cosb cos beta
     */
    void setCosb(double cosb);

    /**
     * 
     * @return superpotential Higgs mixing term 
     */
    gslpp::complex getMuH() const {
        return muH;
    }

    /**
     *
     * @param muH superpotential Higgs mixing term
     */
    void setMuH(gslpp::complex muH) {
        this->muH = muH;
    }

    /**
     *
     * @return chargino mass
     */
    gslpp::vector<double> getMch() const {
        return Mch;
    }

    /**
     * @brief set chargino mass
     * @param Mch the chargino mass
     */
    void setMch(gslpp::vector<double> Mch) {
        this->Mch = Mch;
    }

    /**
     *
     * @return neutralino mass
     */
    gslpp::vector<double> getMneu() const {
        return Mneu;
    }

    /**
     * @brief set the neutralino mass
     * @param Mneu the neutralino mass
     */
    void setMneu(gslpp::vector<double> Mneu) {
        this->Mneu = Mneu;
    }

    /**
     *
     * @return sdown mass squared
     */
    gslpp::vector<double> getMsd2() const {
        return Msd2;
    }

    /**
     * @brief set the sdown mass squared
     * @param Msd the sdown mass squared
     */
    void setMsd2(gslpp::vector<double> Msd2) {
        this->Msd2 = Msd2;
    }

    /**
     *
     * @return charged slepton mass squared
     */
    gslpp::vector<double> getMsl2() const {
        return Msl2;
    }

    /**
     * @brief set the slepton mass squared
     * @param Msl2 the slepton mass squared
     */
    void setMsl2(gslpp::vector<double> Msl2) {
        this->Msl2 = Msl2;
    }

    /**
     *
     * @return sneutrino mass squared
     */
    gslpp::vector<double> getMsn2() const {
        return Msn2;
    }

    /**
     * @brief set the sneutrino mass squared
     * @param Msn2 the sneutrino mass squared
     */
    void setMsn2(gslpp::vector<double> Msn2) {
        this->Msn2 = Msn2;
    }

    /**
     *
     * @return up squark mass squared
     */
    gslpp::vector<double> getMsu2() const {
        return Msu2;
    }

    /**
     * @brief set the up squark mass squared
     * @param Msu2 the up squark mass squared
     */
    void setMsu2(gslpp::vector<double> Msu2) {
        this->Msu2 = Msu2;
    }

    /**
     *
     * @return neutralino rotation matrix
     */
    gslpp::matrix<gslpp::complex> getN() const {
        return N;
    }

    /**
     * @brief set the neutralino rotation matrix
     * @param N the neutralino rotation matrix
     */
    void setN(gslpp::matrix<gslpp::complex> N) {
        this->N = N;
    }

    /**
     *
     * @return down squark rotation matrix
     */
    gslpp::matrix<gslpp::complex> getRd() const {
        return Rd;
    }

    /**
     * @brief set the down squark rotation matrix
     * @param Rd the down squark rotation matrix
     */
    void setRd(gslpp::matrix<gslpp::complex> Rd) {
        this->Rd = Rd;
    }

    /**
     *
     * @return charged slepton rotation matrix
     */
    gslpp::matrix<gslpp::complex> getRl() const {
        return Rl;
    }

    /**
     * @brief set the charged slepton rotation matrix
     * @param Rl the charged slepton rotation matrix
     */
    void setRe(gslpp::matrix<gslpp::complex> Rl) {
        this->Rl = Rl;
    }

    /**
     *
     * @return sneutrino rotation matrix
     */
    gslpp::matrix<gslpp::complex> getRn() const {
        return Rn;
    }

    /**
     * @brief set the sneutrino rotation matrix
     * @param Rn the sneutrino rotation matrix
     */
    void setRn(gslpp::matrix<gslpp::complex> Rn) {
        this->Rn = Rn;
    }

    /**
     *
     * @return up squark rotation matrix
     */
    gslpp::matrix<gslpp::complex> getRu() const {
        return Ru;
    }

    /**
     * @brief set the up squark rotation matrix
     * @param Ru the up squark rotation matrix
     */
    void setRu(gslpp::matrix<gslpp::complex> Ru) {
        this->Ru = Ru;
    }

    /**
     *
     * @return negative chargino rotation matrix
     */
    gslpp::matrix<gslpp::complex> getU() const {
        return U;
    }

    /**
     * @brief set the negative chargino rotation matrix
     * @param U the negative chargino rotation matrix
     */
    void setU(gslpp::matrix<gslpp::complex> U) {
        this->U = U;
    }

    /**
     *
     * @return positive chargino rotation matrix
     */
    gslpp::matrix<gslpp::complex> getV() const {
        return V;
    }

    /**
     * @brief set the positive chargino rotation matrix
     * @param V the positive chargino rotation matrix
     */
    void setV(gslpp::matrix<gslpp::complex> V) {
        this->V = V;
    }


    ///////////////////////////////////////////////////////////////////////////
    /* Functions for EW precision observables */

    /**
     * @return the W boson mass, including radiative corrections
     */
    double Mw() const;

    /**
     * @param[in] INDF fermion index [0-9] (see EWphysics::flavour_st_to_int())
     * @return the ratio of the effective vector coupling constants @f$g_Z^f=g_V^f/g_A^f@f$ for INDF
     */
    gslpp::complex gZf(const int INDF) const; // gZf = gVf/gAf

    /**
     * @param[in] INDF fermion index [0-9] (see EWphysics::flavour_st_to_int())
     * @return the weak form factor for INDF
     */
    gslpp::complex rhoZf(const int INDF) const;

    /**
     * @return the radiative-correction factor @f$\Delta r@f$
     */
    double Delta_r() const;


    ///////////////////////////////////////////////////////////////////////////

private:
    void setY(double tanb_i);

protected:
    gslpp::matrix<gslpp::complex> Ru, Rd, Rl, Rn, U, V, N;
    gslpp::vector<double> Msu2, Msd2, Msl2, Msn2, Mch, Mneu;
    double mHh, mHa, mHp, tanb, sinb, cosb, m3;
    gslpp::complex muH;
};

#endif	/* SUSY_H */

