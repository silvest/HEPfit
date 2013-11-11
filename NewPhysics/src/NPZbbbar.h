/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPZBBBAR_H
#define	NPZBBBAR_H

#include "NPbase.h"

/**
 * @class NPZbbbar
 * @brief A class for new physics with non-standard @f$Zb\bar{b}@f$ couplings. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPZbbbar : public NPbase  {
public:
    static const int NZbbbarVars = 2;
    static const std::string ZbbbarVars[NZbbbarVars];
    static const int NZbbbarflags = 2;
    static const std::string Zbbbarflags[NZbbbarflags];

    NPZbbbar();

    virtual std::string ModelName() const 
    {
        return "NPZbbbar";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();  
    virtual void setEWSMflags(EWSM& myEWSM);    

    virtual bool setFlag(const std::string, const bool&); 
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    bool IsFlagNotLinearizedNP() const
    {
        return FlagNotLinearizedNP;
    }
    
    ////////////////////////////////////////////////////////////////////////    

    virtual double deltaGVl(StandardModel::lepton l) const;
    virtual double deltaGVq(StandardModel::quark q) const;
    virtual double deltaGAl(StandardModel::lepton l) const;
    virtual double deltaGAq(StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////    
    
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
     * @return epsilon_b
     */
    virtual double epsilonb() const;
    
    
    ////////////////////////////////////////////////////////////////////////   

protected:
    virtual void setParameter(const std::string name, const double& value);
    

    ////////////////////////////////////////////////////////////////////////   

private:
    /* These variables may be used as the deviations in the left-handed 
     * and right-handed couplings if the flag "NPZbbbarLR" is set to true.
     * Therefore, they should not be used directly. Instead, the functions
     * deltaGVq() and deltaGAq() have to be called. */
    double myDeltaGVb, myDeltaGAb;

    /*
     * If true,
     *    myDeltaGVb --> delta g_L^b
     *    myDeltaGAb --> delta g_R^b
     */
    bool FlagNPZbbbarLR;

    bool FlagNotLinearizedNP;

};

#endif	/* NPZBBBAR_H */

