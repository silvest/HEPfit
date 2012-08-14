/* 
 * File:   MFV.h
 * Author: silvest
 *
 * Created on September 24, 2010, 10:53 AM
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

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool PreUpdate();
    virtual bool PostUpdate();
    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

private:
    void SetSoftTerms(void);
    void SetParameter(const std::string, const double&);
    double a1, a2, a3, a6, a7, x1, x2, y1, y3, y6, y7;
    complex a4, a5, a8, y2, y4, y5, w1, w2, w3, w4, w5;
    Xmatrices X;
};

#endif	/* MFV_H */
