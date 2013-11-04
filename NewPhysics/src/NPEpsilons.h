/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_H
#define	NPEPSILONS_H

#include <StandardModel.h>
#include <EWepsilons.h>

/**
 * @class NPEpsilons
 * @brief A class for new physics with the @f$\varepsilon@f$ parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 *
 * Flags:
 * \li \b EWABC:&nbsp; use EW_ABC class based on the formulae in Eqs.(7)-(14)
 * of IJMP, A7, 1031-1058 (1998) by Altarelli et al.
 * \li \b EWABC2:&nbsp; use use the approximate formulae in Eqs.(16)-(20) of
 * IJMP, A7, 1031-1058 (1998) by Altarelli et al.
 */
class NPEpsilons : public StandardModel  {
public:
    static const int NEPSILONvars = 4;
    static const std::string EPSILONvars[NEPSILONvars];
    static const int NEPSILONflags = 6;
    static const std::string EPSILONflags[NEPSILONflags];
    
    /**
     * @brief NPEpsilons constructor
     */
    NPEpsilons();

    virtual std::string ModelName() const 
    {
        return "NPEpsilons";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();  
    virtual void SetEWSMflags(EWSM& myEWSM);    

    virtual bool SetFlag(const std::string, const bool&); 
    
    
    ////////////////////////////////////////////////////////////////////////

    bool IsFlagEWABC() const
    {
        return FlagEWABC;
    }

    bool IsFlagEWABC2() const
    {
        return FlagEWABC2;
    }

    
    ////////////////////////////////////////////////////////////////////////

    virtual double epsilon1() const 
    {
        if (FlagEpsilon1SM) 
            return epsilon1_SM();
        else
            return myEpsilon_1;
    }

    virtual double epsilon2() const 
    {
        if (FlagEpsilon2SM) 
            return epsilon2_SM();
        else
            return myEpsilon_2;
    }

    virtual double epsilon3() const 
    {
        if (FlagEpsilon3SM) 
            return epsilon3_SM();
        else
            return myEpsilon_3;
    }

    virtual double epsilonb() const 
    {
        if (FlagEpsilonbSM) 
            return epsilonb_SM();
        else
            return myEpsilon_b;
    }    

    
    ////////////////////////////////////////////////////////////////////////     
    
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

    
    ////////////////////////////////////////////////////////////////////////   

protected:    
    double myEpsilon_1, myEpsilon_2, myEpsilon_3, myEpsilon_b;
    virtual void setParameters(const std::string name, const double& value);
    
    
    ////////////////////////////////////////////////////////////////////////     
    
private:
    bool FlagEpsilon1SM, FlagEpsilon2SM, FlagEpsilon3SM, FlagEpsilonbSM;
    bool FlagEWABC, FlagEWABC2;
    EWepsilons* myEWepsilons;
    
};

#endif	/* NPEPSILONS_H */

