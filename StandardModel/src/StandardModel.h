/*
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODEL_H
#define	STANDARDMODEL_H

#include <gslpp.h>
#include "QCD.h"
#include "CKM.h"
#include "WilsonCoefficient.h"
#include "StandardModelMatching.h"

using namespace gslpp;
class EWSM; // forward reference to EWSM class

/**
 * @class StandardModel
 * @ingroup StandardModel
 * @brief A class for the Standard Model.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is ....
 *
 * Model parameters:
 * \li \b ale:&nbsp; the fine-structure constant @f$\alpha(0)@f$,
 * \li \b dAle5Mz:&nbsp; the five-flavour hadronic contribution to the electromagnetic coupling,
 * \li \b GF:&nbsp; the Fermi constant in @f${\rm GeV}^{-2}@f$,
 * \li \b mHl:&nbsp; the Higgs mass in GeV,
 * \li \b delMw:&nbsp; the theoretical uncertainty in @f$M_W@f$ in GeV,
 * \li \b delSin2th_l:&nbsp; the theoretical uncertainty in @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$,
 * \li \b delRhoZ_nu:&nbsp;
 * \li \b delRhoZ_e:&nbsp;
 * \li \b delRhoZ_b:&nbsp;
 * \li \b muw:&nbsp;
 * \li \b mneutrino_1:&nbsp;
 * \li \b mneutrino_2:&nbsp;
 * \li \b mneutrino_3:&nbsp;
 * \li \b melectron:&nbsp; the electron mass in GeV,
 * \li \b mmu:&nbsp; the muon mass in GeV,
 * \li \b mtau:&nbsp; the tau mass in GeV,
 * \li \b lambda:&nbsp;
 * \li \b A:&nbsp;
 * \li \b rhob:&nbsp;
 * \li \b etab:&nbsp;
 * \li \b EpsK:&nbsp; the experimental value of @f$\varepsilon_{K}@f$,
 * \li \b phiEpsK:&nbsp; the experimental value of @f$\Delta M_{K}/(\Delta\Gamma_{K}/2)@f$,
 * \li \b DeltaMK:&nbsp; the experimental value of @f$\Delta M_{K}@f$ in GeV,
 * \li \b KbarEpsK:&nbsp;
 * \li \b Dmk:&nbsp; the SM contribution to @f$\Delta m_{K}@f$ in GeV,
 * \li \b SM_M12D:&nbsp; the SM amplitude of the @f$D^{0}-\bar{D}^{0}@f$ mixing,
 *
 * Flags:
 * \li \b WithoutNonUniversalVC:&nbsp;
 * \li \b ApproximateGqOverGb:&nbsp;
 * \li \b ApproximateGammaZ:&nbsp;
 * \li \b ApproximateSigmaH:&nbsp;
 * \li \b RhoZbFromGuOverGb:&nbsp;
 * \li \b RhoZbFromGdOverGb:&nbsp;
 * \li \b TestSubleadingTwoLoopEW:&nbsp;
 * \li \b EWCHMN:&nbsp;
 *
 */
class StandardModel: public QCD {
public:

    /**
     * An enum type for leptons.
     */
    enum lepton {
        NEUTRINO_1,
        ELECTRON, /**< Electron */
        NEUTRINO_2,
        MU, /**< Muon */
        NEUTRINO_3,
        TAU /**< Tau */
    };

    static const int NSMvars = 27;
    static const std::string SMvars[NSMvars];
    static const int NSMflags = 8;
    static const std::string SMflags[NSMflags];

    /**
     * A constructor.
     */
    StandardModel();

    virtual std::string ModelName() const
    {
        return "StandardModel";
    }

    virtual bool Init(const std::map<std::string, double>& DPars);

    virtual bool PreUpdate();

    /**
     * Updates the SM parameters found in the argument.
     * @param[in] DPars A map containing the parameters (all as double) to be updated.
     */
    virtual bool Update(const std::map<std::string, double>& DPars);

    virtual bool PostUpdate();

    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ///////////////////////////////////////////////////////////////////////////
    // Flags

    virtual bool setFlag(const std::string name, const bool& value);

    virtual bool CheckFlags() const;

    bool IsFlagWithoutNonUniversalVC() const
    {
        return FlagWithoutNonUniversalVC;
    }

    bool IsFlagApproximateGqOverGb() const
    {
        return FlagApproximateGqOverGb;
    }

    bool IsFlagApproximateGammaZ() const
    {
        return FlagApproximateGammaZ;
    }

    bool IsFlagApproximateSigmaH() const
    {
        return FlagApproximateSigmaH;
    }

    bool IsFlagRhoZbFromGuOverGb() const
    {
        return FlagRhoZbFromGuOverGb;
    }

    bool IsFlagRhoZbFromGdOverGb() const
    {
        return FlagRhoZbFromGdOverGb;
    }

    bool IsFlagTestSubleadingTwoLoopEW() const
    {
        return FlagTestSubleadingTwoLoopEW;
    }

    bool IsFlagEWCHMN() const
    {
        return FlagEWCHMN;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // Initialization and Matching

    virtual bool InitializeModel();

    virtual void setEWSMflags(EWSM& myEWSM);

    virtual StandardModelMatching* getMyMatching() const
    {
        return myStandardModelMatching;
    }

    virtual void setMyMatching(StandardModelMatching* myMatching)
    {
        this->myStandardModelMatching = myMatching;
    }

    virtual const double matchingScale() const
    {
        return muw;
    }


    ///////////////////////////////////////////////////////////////////////////
    // get and set methods for class members

    Particle getLeptons(const StandardModel::lepton p) const
    {
        return leptons[p];
    }

    /**
     * @return the electromagnetic coupling
     */
    double getAle() const
    {
        return ale;
    }

    /**
     * Sets the fine-structure constant.
     * @param[in] ale the fine-structure constant.
     */
    void setAle(double ale)
    {
        this->ale = ale;
    }

    /**
     * @return @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$.
     */
    double getDAle5Mz() const
    {
        return dAle5Mz;
    }

    /**
     * Sets @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$.
     * @param[in] dAle5Mz @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$.
     */
    void setDAle5Mz(double dAle5Mz)
    {
        this->dAle5Mz = dAle5Mz;
    }

    /**
     * @return The Fermi constant.
     */
    double getGF() const
    {
        return GF;
    }

    /**
     * Sets the Fermi constant.
     * @param[in] GF the Fermi constant.
     */
    void setGF(double GF)
    {
        this->GF = GF;
    }

    /**
     * @return The Higgs mass.
     */
    double getMHl() const
    {
        return mHl;
    }

    /**
     * Sets the Higgs mass.
     * @param[in] mHl The Higgs mass.
     */
    void setMHl(double mHl)
    {
        this->mHl = mHl;
    }

    /**
     * @return Theoretical uncertainty in the approximate formula for @f$M_W@f$.
     */
    double getDelMw() const
    {
        return delMw;
    }

    /**
     * @return Theoretical uncertainty in the approximate formula for the leptonic weak mixing angle.
     */
    double getDelSin2th_l() const
    {
        return delSin2th_l;
    }

    /**
     * @return Theoretical uncertainty in the total width of Z boson.
     */
    double getDelGammaZ() const
    {
        return delGammaZ;
    }

    double getDelRhoZ_nu() const
    {
        return delRhoZ_nu;
    }

    double getDelRhoZ_e() const
    {
        return delRhoZ_e;
    }

    double getDelRhoZ_b() const
    {
        return delRhoZ_b;
    }

    /**
     * @return The CKM matrix.
     */
    matrix<complex> getVCKM() const
    {
        return VCKM;
    }

    /**
     * Sets the CKM matrix.
     * @param[in] VCKM The CKM matrix.
     */
    void setVCKM(matrix<complex> VCKM)
    {
        this->VCKM = VCKM;
    }

    CKM getCKM() const
    {
        return myCKM;
    }

    double getLambda() const 
    {
        return lambda;
    }

    double getA() const
    {
        return A;
    }

    double getEtab() const
    {
        return etab;
    }
    double getRhob() const
    {
        return rhob;
    }

    /**
     * @return The PMNS matrix.
     */
    matrix<complex> getUPMNS() const
    {
        return UPMNS;
    }

    /**
     * Sets the PMNS matrix.
     * @param[in] UPMNS The PMNS matrix.
     */
    void setUPMNS(matrix<complex> UPMNS)
    {
        this->UPMNS = UPMNS;
    }

    /**
     * @return The up Yukawa matrix.
     */
    matrix<complex> getYu() const
    {
        return Yu;
    }

    /**
     * @return The down Yukawa matrix.
     */
    matrix<complex> getYd() const
    {
        return Yd;
    }

    /**
     * @return The neutrino Yukawa matrix.
     */
    matrix<complex> getYn() const
    {
        return Yn;
    }

    /**
     * @return The charged lepton Yukawa matrix.
     */
    matrix<complex> getYe() const
    {
        return Ye;
    }

    /**
     * Sets the Standard Model contribution to @f$ \Delta m_{K} @f$.
     * @param[in] Dmk The Standard Model contribution to @f$ \Delta m_{K} @f$.
     */
    void setDmk(double Dmk)
    {
        this->Dmk = Dmk;
    }

    /**
     * @return The Standard Model contribution to @f$ \Delta m_{K} @f$.
     */
    double getDmk() const
    {
        return Dmk;
    }

    /**
     * @return The Standard Model amplitude of the @f$ D^{0} - \bar{D}^{0} @f$ mixing.
     */
    double getSM_M12D() const
    {
        return SM_M12D;
    }

    double getMuw() const
    {
        return muw;
    }

    double getKbarEpsK() const
    {
        return KbarEpsK;
    }

    /**
     * @return The experimental value of @f$ \Delta M_{K}/(\Delta\Gamma_{K}/2) @f$.
     */
    double getphiEpsK() const
    {
        return phiEpsK;
    }

    /**
     * @return The experimental value of @f$ \Delta M_{K} @f$.
     */
    double getDeltaMK() const
    {
        return DeltaMK;
    }

    /**
     * @return The experimental value of @f$ \varepsilon_{K} @f$.
     */
    double getEpsK() const {
        return EpsK;
    }


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @return The VEV.
     */
    virtual double v() const;

    /**
     * @return The W boson mass at tree level.
     */
    virtual double Mw_tree() const;

    /**
     * Computes the running electromagnetic coupling alpha(mu) in the on-shell
     * scheme, where the top-quark contribution is not included.
     * @param[in] mu A scale @f$\mu@f$ in GeV.
     * @param[in] order (=LO, FULLNLO)
     * @return @f$\alpha(\mu)@f$ in the on-shell scheme.
     */
    double ale_OS(const double mu, orders order=FULLNLO) const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return The W boson mass without weak corrections.
     */
    double Mw0() const;

    /**
     * @return @f$sin^2\theta_W@f$ without weak corrections.
     */
    double s02() const;

    /**
     * @return @f$\cos^2\theta_W@f$ without weak corrections.
     */
    double c02() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return A pointer to EWSM object.
     */
    EWSM* getEWSM() const
    {
        return myEWSM;
    }

    /**
     * @param[in] s An invariant mass squared.
     * @return The leptonic corrections to alpha at @f$M_Z@f$.
     */
    double DeltaAlphaLepton(const double s) const;

    /**
     * @return The sum of the leptonic and hadronic corrections to alpha at @f$M_Z@f$.
     */
    double DeltaAlphaL5q() const;

    /**
     * @return The total (leptonic+hadronic+top) corrections to alpha at @f$M_Z@f$.
     */
    double DeltaAlpha() const;

    /**
     * @return The electromagnetic coupling at @f$M_Z@f$, @f$\alpha(M_Z)@f$.
     */
    double alphaMz() const;

    /**
     * @return The W boson mass.
     */
    virtual double Mw() const;

    /**
     * @return @f$M_W^2/M_Z^2@f$.
     */
    virtual double cW2() const;

    /**
     * @return @f$1-M_W^2/M_Z^2@f$.
     */
    virtual double sW2() const;

    /**
     * @param[in] l lepton.
     * @return The effective coupling @f$\rho_Z^l@f$.
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark.
     * @return The effective coupling @f$\rho_Z^q@f$.
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l name of lepton.
     * @return The effective coupling @f$\kappa_Z^l@f$ in the SM.
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of quark.
     * @return The effective coupling @f$\kappa_Z^q@f$ in the SM.
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l lepton.
     * @return The effective vector coupling for neutral-current interactions @f$g_V^l@f$.
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark.
     * @return The effective vector coupling for neutral-current interactions @f$g_V^q@f$.
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @param[in] l lepton.
     * @return The effective axial-vector coupling for neutral-current interactions @f$g_A^l@f$.
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark.
     * @return The effective axial-vector coupling for neutral-current interactions @f$g_A^q@f$.
     */
    virtual complex gAq(const StandardModel::quark q) const;

    /**
     * @return The total width of the W boson.
     */
    virtual double GammaW() const;

    /**
     * @return SM contribution to @f$\epsilon_1@f$.
     */
    double epsilon1_SM() const;

    /**
     * @return SM contribution to @f$\epsilon_2@f$.
     */
    double epsilon2_SM() const;

    /**
     * @return SM contribution to @f$\epsilon_3@f$.
     */
    double epsilon3_SM() const;

    /**
     * @return SM contribution to @f$\epsilon_b@f$.
     */
    double epsilonb_SM() const;

    /**
     * @return @f$\epsilon_1@f$.
     */
    virtual double epsilon1() const;

    /**
     * @return @f$\epsilon_2@f$.
     */
    virtual double epsilon2() const;

    /**
     * @return @f$\epsilon_3@f$.
     */
    virtual double epsilon3() const;

    /**
     * @return @f$\epsilon_b@f$.
     */
    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////
    // CKM parameters

    // Angles
    double computeBeta() const;
    double computeGamma() const;
    double computeAlpha() const;
    double computeBetas() const;

    // Lambda_q
    complex computelamt() const;
    complex computelamc() const;
    complex computelamu() const;

    complex computelamt_d() const;
    complex computelamc_d() const;
    complex computelamu_d() const;

    complex computelamt_s() const;
    complex computelamc_s() const;
    complex computelamu_s() const;

    // Sides
    double getRt() const;
    double getRts() const;
    double getRb() const;


    ////////////////////////////////////////////////////////////////////////
protected:
    virtual void setParameter(const std::string, const double&);
    virtual void computeCKM();
    virtual void computeYukawas();

    matrix<complex> VCKM;
    matrix<complex> UPMNS;
    matrix<complex> Yu;
    matrix<complex> Yd;
    matrix<complex> Yn;
    matrix<complex> Ye;
    Particle leptons[6];
    EWSM* myEWSM;
    CKM myCKM;

    // model parameters
    double ale;
    double dAle5Mz;
    double GF;
    double mHl;
    double delMw;
    double delSin2th_l;
    double delGammaZ;
    double delRhoZ_nu;
    double delRhoZ_e;
    double delRhoZ_b;
    double muw;
    double lambda;
    double A;
    double rhob;
    double etab;
    double EpsK;
    double phiEpsK;
    double DeltaMK;
    double KbarEpsK;
    double Dmk;
    double SM_M12D;
    
    
    ////////////////////////////////////////////////////////////////////////    
private:
    bool FlagWithoutNonUniversalVC;
    bool FlagApproximateGqOverGb;
    bool FlagApproximateGammaZ;
    bool FlagApproximateSigmaH;
    bool FlagRhoZbFromGuOverGb;
    bool FlagRhoZbFromGdOverGb;
    bool FlagTestSubleadingTwoLoopEW;
    bool FlagEWCHMN;
    bool requireCKM;
    bool requireYe;
    bool requireYn;
    StandardModelMatching* myStandardModelMatching;

};

#endif	/* STANDARDMODEL_H */
