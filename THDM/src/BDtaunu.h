/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BDTAUNU_H
#define	BDTAUNU_H

#include <ThObservable.h>

class THDM;

/**
 * @addtogroup THDM
 * @brief A module for the @f$Z_2@f$ symmetric Two-Higgs-Doublet models.
 * @details This module contains code necessary for analyzing theoretical
 * and experimental constraints on @f$Z_2@f$ symmetric Two-Higgs-Doublet models.
 * @{
 */

/**
 * @class BDtaunu
 * @ingroup THDM
 * @brief A class for @f$B \to D^{(*)} \tau \nu@f$ decays in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details We use the parametrized approach and the numerical values from @cite Lees:2013uzd.
 */
class BDtaunu : public ThObservable {
public:

    /**
     * @brief Constructor of the BDtaunu class.
     */
    BDtaunu(const StandardModel& SM_i);

    /**
     * @brief Destructor of the BDtaunu class.
     */
    virtual ~BDtaunu();

    /**
     * @brief Empty function.
     */
    double computeThValue();

    const THDM * myTHDM;

protected:

private:
};

/**
 * @class RBDtaunu
 * @ingroup THDM
 * @brief A class for @f$B \to D \tau \nu@f$ decays in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class RBDtaunu: public BDtaunu {
public:

    /**
     * @brief Constructor of the RBDtaunu class.
     */
    RBDtaunu(const StandardModel& SM_i);

    /**
     * @brief Calculates the ratio of the branching fractions of @f$B \to D \tau \nu@f$ and @f$B \to D \ell \nu@f$ decays.
     * @details The expression is Eq. (20) in @cite Lees:2013uzd.
     */
    double computeThValue();
    
private:
};

/**
 * @class RBDstartaunu
 * @ingroup THDM
 * @brief A class for @f$B \to D^* \tau \nu@f$ decays in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class RBDstartaunu: public BDtaunu {
public:
    
    /**
     * @brief Constructor of the RBDstartaunu class.
     */
    RBDstartaunu(const StandardModel& SM_i);
    
    /**
     * @brief Calculates the ratio of the branching fractions of @f$B \to D^* \tau \nu@f$ and @f$B \to D^* \ell \nu@f$ decays.
     */
    double computeThValue();
    
private:
};

/**
 * @class obsBDtaunu_SM
 * @ingroup THDM
 * @brief A class for the parametrized value of the ratio of the branching fractions of @f$B \to D \tau \nu@f$ and @f$B \to D \ell \nu@f$ decays in the Standard Model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class maps the @f${\cal R}(D)_{\rm SM}@f$ parameter to a virtual observable in order to include the correlations.
 */
class obsBDtaunu_SM: public BDtaunu {
public:
    
    /**
     * @brief Constructor of the obsBDtaunu_SM class.
     */
    obsBDtaunu_SM(const StandardModel& SM_i);

    /**
     * @brief Returns the parametrized value of the ratio of the branching fractions of @f$B \to D \tau \nu@f$ and @f$B \to D \ell \nu@f$ decays in the Standard Model.
     */
    double computeThValue();
    
private:
};

/**
 * @class obsBDtaunu_A
 * @ingroup THDM
 * @brief A class for the parameter @f$A_D@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class maps the @f$A_D@f$ parameter to a virtual observable in order to include the correlations.
 */
class obsBDtaunu_A: public BDtaunu {
public:

    /**
     * @brief Constructor of the obsBDtaunu_A class.
     */
    obsBDtaunu_A(const StandardModel& SM_i);

    /**
     * @brief Returns the @f$A_D@f$ parameter.
     */
    double computeThValue();

private:
};

/**
 * @class obsBDtaunu_B
 * @ingroup THDM
 * @brief A class for the parameter @f$B_D@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class maps the @f$B_D@f$ parameter to a virtual observable in order to include the correlations.
 */
class obsBDtaunu_B: public BDtaunu {
public:

    /**
     * @brief Constructor of the obsBDtaunu_B class.
     */
    obsBDtaunu_B(const StandardModel& SM_i);

    /**
     * @brief Returns the @f$B_D@f$ parameter.
     */
    double computeThValue();

private:
};

/**
 * @class obsBDstartaunu_SM
 * @ingroup THDM
 * @brief A class for the parametrized value of the ratio of the branching fractions of @f$B \to D^* \tau \nu@f$ and @f$B \to D^* \ell \nu@f$ decays in the Standard Model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class maps the @f${\cal R}(D^*)_{\rm SM}@f$ parameter to a virtual observable in order to include the correlations.
 */
class obsBDstartaunu_SM: public BDtaunu {
public:

    /**
     * @brief Constructor of the obsBDstartaunu_SM class.
     */
    obsBDstartaunu_SM(const StandardModel& SM_i);

    /**
     * @brief Returns the parametrized value of the ratio of the branching fractions of @f$B \to D^* \tau \nu@f$ and @f$B \to D^* \ell \nu@f$ decays in the Standard Model.
     */
    double computeThValue();

private:
};

/**
 * @class obsBDstartaunu_A
 * @ingroup THDM
 * @brief A class for the parameter @f$A_{D^*}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class maps the @f$A_{D^*}@f$ parameter to a virtual observable in order to include the correlations.
 */
class obsBDstartaunu_A: public BDtaunu {
public:

    /**
     * @brief Constructor of the obsBDstartaunu_A class.
     */
    obsBDstartaunu_A(const StandardModel& SM_i);

    /**
     * @brief Returns the @f$A_{D^*}@f$ parameter.
     */
    double computeThValue();

private:
};

/**
 * @class obsBDstartaunu_B
 * @ingroup THDM
 * @brief A class for the parameter @f$B_{D^*}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class maps the @f$B_{D^*}@f$ parameter to a virtual observable in order to include the correlations.
 */
class obsBDstartaunu_B: public BDtaunu {
public:

    /**
     * @brief Constructor of the obsBDstartaunu_B class.
     */
    obsBDstartaunu_B(const StandardModel& SM_i);

    /**
     * @brief Returns the @f$B_{D^*}@f$ parameter.
     */
    double computeThValue();
    
private:
};

/**
 * @}
 */

#endif	/* BDTAUNU_H */
