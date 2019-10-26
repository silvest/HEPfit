/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"

#ifndef NPDF2_H
#define NPDF2_H

/**
 * @class NPDF2
 * @brief Model for general constraints NP contributions to Delta F=2 processes.
 */
class NPDF2 : public StandardModel {
public:
    static const int NNPDF2vars = 6;

    static const std::string NPDF2vars[NNPDF2vars];
    
    /**
     * @brief NPDF2 constructor
     */
    NPDF2();
    
    /**
     * @brief A method to set the value of a parameter of %NPDF2.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    void setParameter(const std::string name, const double& value);
    
    /**
     * @brief A method to check if all the mandatory parameters for %NPDF2
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The ratio of the absolute value of the $B_d$ mixing amplitude over the Standard Model value
     * @return @f$\vert (M_{12}^{bd})_\mathrm{full}/(M_{12}^{bd})_\mathrm{SM}\vert@f$
     */
    double getCBd() const
    {
        return CBd;
    }

    /**
     * @brief The ratio of the absolute value of the $B_s$ mixing amplitude over the Standard Model value
     * @return @f$\vert (M_{12}^{bs})_\mathrm{full}/(M_{12}^{bs})_\mathrm{SM}\vert@f$
     */
    double getCBs() const
    {
        return CBs;
    }

    /**
     * @brief The ratio of the real part of the $K$ mixing amplitude over the Standard Model value
     * @return @f$(\mathrm{Re} M_{12}^{sd})_\mathrm{full}/(\mathrm{Re} M_{12}^{sd})_\mathrm{SM}\vert@f$
     */
    double getCDMK() const
    {
        return CDMK;
    }

    /**
     * @brief The ratio of the imaginary part of the $K$ mixing amplitude over the Standard Model value
     * @return @f$(\mathrm{Im} M_{12}^{sd})_\mathrm{full}/(\mathrm{Im} M_{12}^{sd})_\mathrm{SM}\vert@f$
     */
    double getCepsK() const
    {
        return CepsK;
    }

    /**
     * @brief Half the relative phase of the $B_s$ mixing amplitude w.r.t. the Standard Model one
     * @return @f$ 1/2 (\mathrm{arg}((M_{12}^{bs})_\mathrm{full})-\mathrm{arg}((M_{12}^{bs})_\mathrm{SM}))\vert@f$
     */
    double getPhiBs() const
    {
        return PhiBs;
    }

    /**
     * @brief Half the relative phase of the $B_d$ mixing amplitude w.r.t. the Standard Model one
     * @return @f$ 1/2 (\mathrm{arg}((M_{12}^{bd})_\mathrm{full})-\mathrm{arg}((M_{12}^{bd})_\mathrm{SM}))\vert@f$
     */
    double getPhiBd() const
    {
        return phiBd;
    }

    
private:

    double CepsK,CDMK,CBd,phiBd,CBs,PhiBs;///< Modifications to the Standard %Model couplings.
    
};

#endif /* NPDF2_H */

