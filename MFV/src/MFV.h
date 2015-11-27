/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MFV_H
#define	MFV_H

#include <SUSY.h>
#include "Xmatrices.h"

/**
 * @class MFV
 * @brief calculation of soft masses and trilinear couplings for squarks and sleptons in 
 * Minimal Flavour Violating SUSY models as in arXiv:0807.0801v2
 */
class MFV : public SUSY {
public:
    static const int NMFVvars = 33;
    static const std::string MFVvars[NMFVvars];

    MFV();

    virtual bool InitializeModel();
    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool PreUpdate();
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool PostUpdate();
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    
protected:
    virtual void setParameter(const std::string, const double&);
    virtual void SetSoftTerms();

    double a1, a2, a3, a6, a7, x1, x2, y1, y3, y6, y7;
    gslpp::complex a4, a5, a8, y2, y4, y5, w1, w2, w3, w4, w5;
    Xmatrices X;
};

#endif	/* MFV_H */
