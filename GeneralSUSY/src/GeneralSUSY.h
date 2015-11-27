/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALSUSY_H
#define	GENERALSUSY_H

#include "SUSY.h"

/**
 * @addtogroup GeneralSUSY
 * @brief A module for the general MSSM.
 */

/**
 * @class GeneralSUSY
 * @brief A class for the general MSSM. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Sets all soft masses and trilinear couplings for squarks and sleptons.
 */
class GeneralSUSY : public SUSY {
public:
    static const int NGeneralSUSYvars = 126;
    static const std::string GeneralSUSYvars[NGeneralSUSYvars];

    GeneralSUSY();

    virtual bool InitializeModel();
    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool PreUpdate();
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool PostUpdate();
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

protected:
    virtual void setParameter(const std::string, const double&);
    virtual void SetSoftTerms();
    double msQhat2_11r, msQhat2_12r, msQhat2_12i, msQhat2_13r, msQhat2_13i, msQhat2_22r, msQhat2_23r, msQhat2_23i, msQhat2_33r,
           msUhat2_11r, msUhat2_12r, msUhat2_12i, msUhat2_13r, msUhat2_13i, msUhat2_22r, msUhat2_23r, msUhat2_23i, msUhat2_33r,
           msDhat2_11r, msDhat2_12r, msDhat2_12i, msDhat2_13r, msDhat2_13i, msDhat2_22r, msDhat2_23r, msDhat2_23i, msDhat2_33r,
           msLhat2_11r, msLhat2_12r, msLhat2_12i, msLhat2_13r, msLhat2_13i, msLhat2_22r, msLhat2_23r, msLhat2_23i, msLhat2_33r,
           msEhat2_11r, msEhat2_12r, msEhat2_12i, msEhat2_13r, msEhat2_13i, msEhat2_22r, msEhat2_23r, msEhat2_23i, msEhat2_33r,
           msNhat2_11r, msNhat2_12r, msNhat2_12i, msNhat2_13r, msNhat2_13i, msNhat2_22r, msNhat2_23r, msNhat2_23i, msNhat2_33r,
           TUhat_11r, TUhat_12r, TUhat_13r, TUhat_21r, TUhat_22r, TUhat_23r, TUhat_31r, TUhat_32r, TUhat_33r,
           TUhat_11i, TUhat_12i, TUhat_13i, TUhat_21i, TUhat_22i, TUhat_23i, TUhat_31i, TUhat_32i, TUhat_33i,
           TDhat_11r, TDhat_12r, TDhat_13r, TDhat_21r, TDhat_22r, TDhat_23r, TDhat_31r, TDhat_32r, TDhat_33r,
           TDhat_11i, TDhat_12i, TDhat_13i, TDhat_21i, TDhat_22i, TDhat_23i, TDhat_31i, TDhat_32i, TDhat_33i,
           TEhat_11r, TEhat_12r, TEhat_13r, TEhat_21r, TEhat_22r, TEhat_23r, TEhat_31r, TEhat_32r, TEhat_33r,
           TEhat_11i, TEhat_12i, TEhat_13i, TEhat_21i, TEhat_22i, TEhat_23i, TEhat_31i, TEhat_32i, TEhat_33i,
           TNhat_11r, TNhat_12r, TNhat_13r, TNhat_21r, TNhat_22r, TNhat_23r, TNhat_31r, TNhat_32r, TNhat_33r,
           TNhat_11i, TNhat_12i, TNhat_13i, TNhat_21i, TNhat_22i, TNhat_23i, TNhat_31i, TNhat_32i, TNhat_33i;
};

/**
 * @}
 */

#endif	/* GeneralSUSY_H */
