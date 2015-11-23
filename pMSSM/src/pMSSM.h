/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PMSSM_H
#define	PMSSM_H

#include <SUSY.h>

/**
 * @addtogroup pMSSM
 * @brief A module for the phenomenological MSSM.
 * @{
 */

/**
 * @class pMSSM
 * @brief A class for the phenomenological MSSM. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Sets all soft masses and trilinear couplings for squarks and sleptons using simplifying assumptions.
 */
class pMSSM : public SUSY {
public:
    static const int NpMSSMvars = 13;
    static const std::string pMSSMvars[NpMSSMvars];

    pMSSM();

    virtual bool InitializeModel();
    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool PreUpdate();
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool PostUpdate();
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

protected:
    virtual void setParameter(const std::string, const double&);
    virtual void SetSoftTerms();
    double msQ12, msQ3, msU12, msU3, msD12, msD3, msL12, msL3, msE12, msE3, 
           AU, AD, AE;
};

/**
 * @}
 */

#endif	/* pMSSM_H */
