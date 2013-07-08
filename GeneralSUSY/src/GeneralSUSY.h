/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALSUSY_H
#define	GENERALSUSY_H

#include <SUSY.h>

/**
 * @addtogroup GeneralSUSY
 * @brief A project for the general MSSM.
 * @{
 */

/**
 * @class GeneralSUSY
 * @brief A class for the general MSSM. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details Sets all soft masses and trilinear couplings for squarks and sleptons.
 */
class GeneralSUSY : public SUSY {
public:
    static const int NGeneralSUSYvars = 126;
    static const std::string GeneralSUSYvars[NGeneralSUSYvars];

    GeneralSUSY();

    virtual std::string ModelName() const
    {
        return "GeneralSUSY";
    }

    virtual bool InitializeModel();
    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool PreUpdate();
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool PostUpdate();
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

protected:
    virtual void SetParameter(const std::string, const double&);
    virtual void SetSoftTerms();
    double msQ2_11r, msQ2_12r, msQ2_12i, msQ2_13r, msQ2_13i, msQ2_22r, msQ2_23r, msQ2_23i, msQ2_33r,
           msU2_11r, msU2_12r, msU2_12i, msU2_13r, msU2_13i, msU2_22r, msU2_23r, msU2_23i, msU2_33r,
           msD2_11r, msD2_12r, msD2_12i, msD2_13r, msD2_13i, msD2_22r, msD2_23r, msD2_23i, msD2_33r,
           msL2_11r, msL2_12r, msL2_12i, msL2_13r, msL2_13i, msL2_22r, msL2_23r, msL2_23i, msL2_33r,
           msE2_11r, msE2_12r, msE2_12i, msE2_13r, msE2_13i, msE2_22r, msE2_23r, msE2_23i, msE2_33r,
           msN2_11r, msN2_12r, msN2_12i, msN2_13r, msN2_13i, msN2_22r, msN2_23r, msN2_23i, msN2_33r,
           TU_11r, TU_12r, TU_13r, TU_21r, TU_22r, TU_23r, TU_31r, TU_32r, TU_33r,
           TU_11i, TU_12i, TU_13i, TU_21i, TU_22i, TU_23i, TU_31i, TU_32i, TU_33i,
           TD_11r, TD_12r, TD_13r, TD_21r, TD_22r, TD_23r, TD_31r, TD_32r, TD_33r,
           TD_11i, TD_12i, TD_13i, TD_21i, TD_22i, TD_23i, TD_31i, TD_32i, TD_33i,
           TE_11r, TE_12r, TE_13r, TE_21r, TE_22r, TE_23r, TE_31r, TE_32r, TE_33r,
           TE_11i, TE_12i, TE_13i, TE_21i, TE_22i, TE_23i, TE_31i, TE_32i, TE_33i,
           TN_11r, TN_12r, TN_13r, TN_21r, TN_22r, TN_23r, TN_31r, TN_32r, TN_33r,
           TN_11i, TN_12i, TN_13i, TN_21i, TN_22i, TN_23i, TN_31i, TN_32i, TN_33i;
};

/**
 * @}
 */

#endif	/* GeneralSUSY_H */
