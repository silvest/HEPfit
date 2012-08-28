/*
 * File:   StandardModel.h
 * Author: silvest
 *
 * Created on November 30, 2010, 1:27 PM
 */

#ifndef STANDARDMODEL_H
#define	STANDARDMODEL_H

#include <gslpp.h>
#include "QCD.h"
#include "CKM.h"
#include "WilsonCoefficient.h"
#include "StandardModelMatching.h"

using namespace gslpp;
class EWSM; // forward reference to ESSM class


/**
 * @class StandardModel
 * @brief Standard Model Class
 */
class StandardModel: public QCD {
public:
    enum lepton {NEUTRINO_1,ELECTRON,NEUTRINO_2,MU,NEUTRINO_3,TAU};
    static const int NSMvars = 20;
    static const std::string SMvars[NSMvars];
    
    StandardModel(const bool bDebug_i=false);

    virtual std::string ModelName() const {
        return "StandardModel";
    }
    
    virtual StandardModelMatching* GetMyMatching() const {
        return myStandardModelMatching;
    }

    virtual void SetMyMatching(StandardModelMatching* myMatching) {
        this->myStandardModelMatching = myMatching;
    }


    /**
     * updates the SM parameters found in the argument
     * @param a map containing the parameters (all as double) to be updated
     * "GF"
     * "mneutrino_1"
     * "mneutrino_2"
     * "mneutrino_3"
     * "melectron"
     * "mmu"
     * "mtau"
     * "lambda"
     * "A"
     * "rhob"
     * "etab"
     * "ale"
     * "dAle5Mz"
     * "mHl"
     * "muw"
     * "phiEpsK"
     * "DeltaMK"
     * "KbarEpsK"
     * "Dmk"
     * "SM_M12D" 
     */
    virtual bool Update(const std::map<std::string, double>&);
    
    virtual bool PreUpdate();
     
    virtual bool PostUpdate();

    virtual bool Init(const std::map<std::string, double>&);
        
    virtual bool CheckParameters(const std::map<std::string, double>&);

    
    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * @return the VEV
     */
    double v() const;

    /**
     * @return the W boson mass at tree level
     */
    double Mw_tree() const;    
    
    /**
     * @return sin^2\theta_W without weak corrections, but with \alpha(Mz^2)
     */
    double s02() const;

    /**
     * @return cos^2\theta_W without weak corrections, but with \alpha(Mz^2)
     */
    double c02() const;
    
    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @return pointer to EWSM object
     */
    EWSM* getEWSM() const {
        return myEWSM;
    }
    
    /**
     * @return the leptonic corrections to alpha at Mz
     */
    double DeltaAlphaLepton() const;

    /**
     * @return the sum of the leptonic and hadronic corrections to alpha at Mz
     */
    double DeltaAlphaL5q() const;
    
    /**
     * @return the total (leptonic+hadronic+top) corrections to alpha at Mz
     */
    double DeltaAlpha() const;
    
    /**
     * @brief electromagnetic coupling alpha at Mz
     * @return alpha(Mz)
     */
    double alphaMz() const;
    
    /**
     * @return the W boson mass
     */
    virtual double Mw() const;

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const;
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const;
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l lepton
     * @return rho_Z^l
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;
    
    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q quark
     * @return rho_Z^q
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @brief effective coupling kappa_Z^l
     * @param[in] l name of lepton
     * @return kappa_Z^l in the SM
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief effective coupling kappa_Z^q
     * @param[in] q name of quark
     * @return kappa_Z^q in the SM
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;       
    
    /**
     * @brief vector effective coupling for neutral-current interactions
     * @param[in] l lepton
     * @return g_V^l
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief vector effective coupling for neutral-current interactions
     * @param[in] q quark
     * @return g_V^q
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @brief axial-vector effective coupling for neutral-current interactions
     * @param[in] l lepton
     * @return g_A^l
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief axial-vector effective coupling for neutral-current interactions
     * @param[in] q quark
     * @return g_A^q
     */
    virtual complex gAq(const StandardModel::quark q) const; 
    
    /**
     * @return the total width of the W boson
     */
    virtual double GammaW() const;
    
    /**
     * @return NP contribution to oblique parameter S
     */
    virtual double obliqueS() const {
        return 0.0;
    };
        
    /**
     * @return NP contribution to oblique parameter T
     */
    virtual double obliqueT() const {
        return 0.0;
    };
    
    /**
     * @return NP contribution to oblique parameter U
     */
    virtual double obliqueU() const {
        return 0.0;
    };
    
    
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
     * @brief set the value Standard Model contribution to \f$ \Delta m_{K} \f$
     * @param Dmk Standard Model contribution to \f$ \Delta m_{K} \f$
     */
    void setDmk(double Dmk) {
        this->Dmk = Dmk;
    }
    
    /**
     * 
     * @return the value of Standard Model contribution to \f$ \Delta m_{K} \f$
     */
    double getDmk() const {
        return Dmk;
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
     * @return the Standard Model amplitude of the \f$ D^{0} - \bar{D}^{0} \f$ mixing
     */
    double getSM_M12D() const {
        return SM_M12D;
    }

    virtual const double matchingScale() const ;

    ////////////////////////////////////////////////////////////////////////
    
    // Angles
    double getBeta() const;
    double getGamma() const;
    double getAlpha() const;
    double getBetas() const;

    // Lambda_q
    complex getlamt() const;
    complex getlamc() const;
    complex getlamu() const;

    complex getlamt_d() const;
    complex getlamc_d() const;
    complex getlamu_d() const;

    complex getlamt_s() const;
    complex getlamc_s() const;
    complex getlamu_s() const;

    // Sides
    double getRt() const;
    double getRts() const;
    double getRb() const;

    CKM getCKM() const {
        return myCKM;
    }
    
    ///////////////////////////////////////////////////////////////////////////

    Particle getLeptons(const int p) const {
        return leptons[p];
    }

    double getMuw() const {
        return muw;
    }
    
    double getKbarEpsK() const {
        return KbarEpsK;
    }
    
    double getphiEpsK() const {
        return phiEpsK;
    }
    
    double getDeltaMK() const {
        return DeltaMK;
    }

    double GetA() const {
        return A;
    }

    double GetEtab() const {
        return etab;
    }

    double GetLambda() const {
        return lambda;
    }

    double GetRhob() const {
        return rhob;
    }
    
    virtual bool SetFlag(const std::string, const bool&);
    
    double GetMHl() const {
        return mHl;
    }    
       
    
    virtual bool InitializeMatching();
    
protected:
    virtual void SetParameter(const std::string, const double&);
    
    
    
    matrix<complex> VCKM,UPMNS, Yu, Yd, Yn, Ye;
    double GF, ale, dAle5Mz, mHl, lambda, A, rhob, etab, Dmk;
    double muw, KbarEpsK, phiEpsK, DeltaMK, SM_M12D;
    CKM myCKM;
    Particle leptons[6];
    EWSM* myEWSM;
    

private:
    bool computeCKM, computeYe, computeYn;
    StandardModelMatching* myStandardModelMatching;
    
   
};

#endif	/* STANDARDMODEL_H */

