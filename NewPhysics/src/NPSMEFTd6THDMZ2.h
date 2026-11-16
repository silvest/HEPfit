/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/*
 * File:   NPSMEFTd6THDMZ2.h
 * Author: V. Miralles
 *
 * Created on November 10, 2026, 7:22 PM
 */

#ifndef NPSMEFTD6THDMZ2_H
#define NPSMEFTD6THDMZ2_H

#include "NPSMEFTd6ATHDM.h"

/**
 * @class NPSMEFTd6THDMZ2
 * @brief Tree-level matching of the THDM with a softly broken Z2 symmetry (types I, II, X, Y and inert)
 * onto the SMEFT, built on top of %NPSMEFTd6ATHDM.
 * @details The model uses the same parameters as %GeneralTHDMZ2: logtb = log10(tan(beta)), bma = beta-alpha
 * and m12_2, with the three heavy scalars degenerate, mH = mA = mH+ = sqY2, which is also the matching scale.
 * The quartic couplings lambda1,...,lambda5 of the Z2 basis are computed as in %GeneralTHDMZ2 and rotated
 * to the Higgs basis to get the Z's of %NPSMEFTd6ATHDM, while the Z2 symmetry fixes etau, etad, etae
 * as functions of tan(beta) and of the model type. The matching onto %NPSMEFTd6MFV is inherited unchanged.
 * The model type is chosen with the flag modelType = type1, type2, typeX, typeY or inert. With the flag
 * Decoupling = true, m12_2 is not an input and is fixed to mH^2 tan(beta)/(1+tan(beta)^2).
 * Since the Z2 symmetry is RG stable, the alignment holds at all scales.
 */
class NPSMEFTd6THDMZ2: public NPSMEFTd6ATHDM {
public:

    static const int NNPSMEFTd6THDMZ2Vars = 3;

    static const std::string NPSMEFTd6THDMZ2Vars[NNPSMEFTd6THDMZ2Vars];

    NPSMEFTd6THDMZ2();

    virtual bool PostUpdate();

    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string name, const bool value);

    virtual bool setFlagStr(const std::string name, const std::string value);

    double gettanb() const {
        return (tanb);
    }

    double getbma() const {
        return (bma);
    }

    double getcosbma() const {
        return (cosbma);
    }

    /**
     * @brief The value of m12_2 used in the computation, which is the derived one if Decoupling = true.
     */
    double getm12_2() const {
        return (m12_2);
    }

    std::string getZ2ModelType() const {
        return flag_model;
    }

protected:

    double logtb = 0.; ///< log10(tan(beta)), input parameter.
    double bma = M_PI/2.; ///< beta-alpha, input parameter.
    double m12_2 = 0.; ///< The soft Z2 breaking term m_12^2, input parameter unless Decoupling = true.
    double tanb = 1.; ///< tan(beta)
    double cosbma = 0.; ///< cos(beta-alpha), only to be used as an observable.
    double Z1 = 0.; ///< Z1 of the Higgs basis, only to be used as an observable.
    double lambda1_Z2 = 0.; ///< lambda1 of the Z2 basis, only to be used as an observable.
    double lambda2_Z2 = 0.; ///< lambda2 of the Z2 basis, only to be used as an observable.
    double lambda3_Z2 = 0.; ///< lambda3 of the Z2 basis, only to be used as an observable.
    double lambda4_Z2 = 0.; ///< lambda4 of the Z2 basis, only to be used as an observable.
    double lambda5_Z2 = 0.; ///< lambda5 of the Z2 basis, only to be used as an observable.

    virtual void setParameter(const std::string name, const double& value);

    /**
     * @brief Sets etau, etad and etae of %NPSMEFTd6ATHDM according to the model type.
     */
    void setZ2Yukawas();

    /**
     * @brief Computes lambda1,...,lambda5 of the Z2 basis from tan(beta), beta-alpha, m12_2 and the
     * heavy mass sqY2, as in %GeneralTHDMZ2, and rotates them to the Higgs basis (Haber-O'Neil
     * conventions, the same as in %GeneralTHDM) to set Z1,...,Z7.
     */
    void setZ2Potential();

    bool CheckModelType(const std::string modeltype) const;

private:

    std::string flag_model = ""; ///< The Z2 model type: type1, type2, typeX, typeY or inert.
    bool flag_decoupling = false; ///< If true, m12_2 = mH^2 tan(beta)/(1+tan(beta)^2).

};

#endif /* NPSMEFTD6THDMZ2_H */

