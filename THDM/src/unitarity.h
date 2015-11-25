/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef UNITARITY_H
#define	UNITARITY_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"

/**
 * @class unitarity
 * @ingroup THDM 
 * @brief An observable class for the requirement of tree level perturbative 
 * unitarity.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require the unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 * The 12 eigenvalues of the S-matrix can be found in @cite Akeroyd:2000wc,Ginzburg:2005dt.
 * They should be smaller than @f$8\pi@f$ in magnitude to preserve unitarity of the S-matrix,
 * however corresponding two-loop calculations in the SM show 
 * that a more reasonable upper bound on their absolute values would be @f$2\pi@f$
 * (see discussion in @cite Baglio:2014nea).
 */
class unitarity : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   unitarity(const StandardModel& SM_i);
     
    /**
     * @brief An empty constructor.
     * @return
     */
    double computeThValue();

    const THDM * myTHDM;
    
};

class unitarity1: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{21+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity1(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{21+}=\frac{1}{2}(\lambda_1+\lambda_2+\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_5^2})@f$
     * @details Corresponds to @f$c_+@f$ from equation (2.8) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity2: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{21-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity2(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{21-}=\frac{1}{2}(\lambda_1+\lambda_2-\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_5^2})@f$
     * @details Corresponds to @f$c_-@f$ from equation (2.8) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity3: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{01+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity3(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{01+}=\frac{1}{2}(\lambda_1+\lambda_2+\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_4^2})@f$
     * @details Corresponds to @f$b_+@f$ from equation (2.7) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity4: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{01-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity4(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{01-}=\frac{1}{2}(\lambda_1+\lambda_2-\sqrt{(\lambda_1-\lambda_2)^2+4\lambda_4^2})@f$
     * @details Corresponds to @f$b_-@f$ from equation (2.7) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity5: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{00+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity5(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{00+}=\frac{3}{2}(\lambda_1+\lambda_2)+\frac{1}{2}\sqrt{9(\lambda_1-\lambda_2)^2+4(2\lambda_3+\lambda_4)^2}@f$
     * @details Corresponds to @f$a_+@f$ from equation (2.6) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity6: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{00-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity6(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{00+}=\frac{3}{2}(\lambda_1+\lambda_2)-\frac{1}{2}\sqrt{9(\lambda_1-\lambda_2)^2+4(2\lambda_3+\lambda_4)^2}@f$
     * @details Corresponds to @f$a_-@f$ from equation (2.6) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity7: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{21}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity7(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{21}=\lambda_3+\lambda_4@f$
     * @details Corresponds to @f$f_1=f_2@f$ from equation (2.5) in @cite Akeroyd:2000wc
     */
////	LAod21 =lambada3+lambada4
    double computeThValue();
};

class unitarity8: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{20}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity8(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{20}=\lambda_3-\lambda_4@f$
     * @details Corresponds to @f$p_1@f$ from equation (2.9) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity9: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{01+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity9(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{01+}=\lambda_3+\lambda_5@f$
     * @details Corresponds to @f$f_-@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity10: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{01-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity10(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{01-}=\lambda_3-\lambda_5@f$
     * @details Corresponds to @f$e_2@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity11: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{00+}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity11(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{00+}=\lambda_3+2\lambda_4+3\lambda_5@f$
     * @details Corresponds to @f$f_+@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

class unitarity12: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{odd}_{00-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity12(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{odd}_{00-}=\lambda_3+2\lambda_4-3\lambda_5@f$
     * @details Corresponds to @f$e_1@f$ from equations (2.5) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

#endif	/* UNITARITY_H */

