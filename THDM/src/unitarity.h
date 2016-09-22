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
#include "THDMcache.h"

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

/**
 * @class unitarity1
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{21+}@f$.
 */
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

/**
 * @class unitarity2
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{21-}@f$.
 */
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

/**
 * @class unitarity3
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{01+}@f$.
 */
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

/**
 * @class unitarity4
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{01-}@f$.
 */
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

/**
 * @class unitarity5
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{00+}@f$.
 */
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

/**
 * @class unitarity6
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{even}_{00-}@f$.
 */
class unitarity6: public unitarity {
public:

    /**
     * @brief Constructor for @f$\Lambda^{even}_{00-}@f$ from equations (9) in @cite Ginzburg:2005dt.
     */
    unitarity6(const StandardModel& SM_i);

    /**
     * @return @f$\Lambda^{even}_{00-}=\frac{3}{2}(\lambda_1+\lambda_2)-\frac{1}{2}\sqrt{9(\lambda_1-\lambda_2)^2+4(2\lambda_3+\lambda_4)^2}@f$
     * @details Corresponds to @f$a_-@f$ from equation (2.6) in @cite Akeroyd:2000wc.
     */
    double computeThValue();
};

/**
 * @class unitarity7
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{21}@f$.
 */
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
    double computeThValue();
};

/**
 * @class unitarity8
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{20}@f$.
 */
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

/**
 * @class unitarity9
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{01+}@f$.
 */
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

/**
 * @class unitarity10
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{01-}@f$.
 */
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

/**
 * @class unitarity11
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{00+}@f$.
 */
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

/**
 * @class unitarity12
 * @ingroup THDM 
 * @brief Eigenvalue of the @f$\Phi \Phi \to \Phi \Phi@f$ scattering matrix.
 * @details @f$\Lambda^{odd}_{00-}@f$.
 */
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

/**
 * @class unitarityNLO1
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO1: public ThObservable {
public:

    /**
     * @brief unitarityNLO1 constructor.
     */
    unitarityNLO1(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO1@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO2
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO2: public ThObservable {
public:

    /**
     * @brief unitarityNLO2 constructor.
     */
    unitarityNLO2(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO2@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO3
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO3: public ThObservable {
public:

    /**
     * @brief unitarityNLO3 constructor.
     */
    unitarityNLO3(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO3@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO4
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO4: public ThObservable {
public:

    /**
     * @brief unitarityNLO4 constructor.
     */
    unitarityNLO4(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO4@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO5
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO5: public ThObservable {
public:

    /**
     * @brief unitarityNLO5 constructor.
     */
    unitarityNLO5(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO5@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO6
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO6: public ThObservable {
public:

    /**
     * @brief unitarityNLO6 constructor.
     */
    unitarityNLO6(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO6@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO7
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO7: public ThObservable {
public:

    /**
     * @brief unitarityNLO7 constructor.
     */
    unitarityNLO7(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO7@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO8
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO8: public ThObservable {
public:

    /**
     * @brief unitarityNLO8 constructor.
     */
    unitarityNLO8(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO8@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO9
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO9: public ThObservable {
public:

    /**
     * @brief unitarityNLO9 constructor.
     */
    unitarityNLO9(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO9@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO10
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO10: public ThObservable {
public:

    /**
     * @brief unitarityNLO10 constructor.
     */
    unitarityNLO10(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO10@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO11
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO11: public ThObservable {
public:

    /**
     * @brief unitarityNLO11 constructor.
     */
    unitarityNLO11(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO11@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO12
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO12: public ThObservable {
public:

    /**
     * @brief unitarityNLO12 constructor.
     */
    unitarityNLO12(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO12@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO13
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO13: public ThObservable {
public:

    /**
     * @brief unitarityNLO13 constructor.
     */
    unitarityNLO13(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO13@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO14
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO14: public ThObservable {
public:

    /**
     * @brief unitarityNLO14 constructor.
     */
    unitarityNLO14(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO14@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO15
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO15: public ThObservable {
public:

    /**
     * @brief unitarityNLO15 constructor.
     */
    unitarityNLO15(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO15@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO16
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO16: public ThObservable {
public:

    /**
     * @brief unitarityNLO16 constructor.
     */
    unitarityNLO16(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO16@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO17
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO17: public ThObservable {
public:

    /**
     * @brief unitarityNLO17 constructor.
     */
    unitarityNLO17(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO17@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO18
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO18: public ThObservable {
public:

    /**
     * @brief unitarityNLO18 constructor.
     */
    unitarityNLO18(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO18@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO19
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO19: public ThObservable {
public:

    /**
     * @brief unitarityNLO19 constructor.
     */
    unitarityNLO19(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO19@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO20
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO20: public ThObservable {
public:

    /**
     * @brief unitarityNLO20 constructor.
     */
    unitarityNLO20(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO20@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO21
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO21: public ThObservable {
public:

    /**
     * @brief unitarityNLO21 constructor.
     */
    unitarityNLO21(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO21@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO22
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO22: public ThObservable {
public:

    /**
     * @brief unitarityNLO22 constructor.
     */
    unitarityNLO22(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO22@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO23
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO23: public ThObservable {
public:

    /**
     * @brief unitarityNLO23 constructor.
     */
    unitarityNLO23(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO23@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO24
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO24: public ThObservable {
public:

    /**
     * @brief unitarityNLO24 constructor.
     */
    unitarityNLO24(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO24@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO25
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO25: public ThObservable {
public:

    /**
     * @brief unitarityNLO25 constructor.
     */
    unitarityNLO25(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO25@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLO26
 * @ingroup THDM
 * @brief .
 */
class unitarityNLO26: public ThObservable {
public:

    /**
     * @brief unitarityNLO26 constructor.
     */
    unitarityNLO26(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLO26@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev1
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev1: public ThObservable {
public:

    /**
     * @brief unitarityNLOev1 constructor.
     */
    unitarityNLOev1(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev1@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev2
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev2: public ThObservable {
public:

    /**
     * @brief unitarityNLOev2 constructor.
     */
    unitarityNLOev2(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev2@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev3
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev3: public ThObservable {
public:

    /**
     * @brief unitarityNLOev3 constructor.
     */
    unitarityNLOev3(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev3@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev4
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev4: public ThObservable {
public:

    /**
     * @brief unitarityNLOev4 constructor.
     */
    unitarityNLOev4(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev4@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev5
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev5: public ThObservable {
public:

    /**
     * @brief unitarityNLOev5 constructor.
     */
    unitarityNLOev5(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev5@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev6
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev6: public ThObservable {
public:

    /**
     * @brief unitarityNLOev6 constructor.
     */
    unitarityNLOev6(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev6@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev7
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev7: public ThObservable {
public:

    /**
     * @brief unitarityNLOev7 constructor.
     */
    unitarityNLOev7(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev7@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev8
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev8: public ThObservable {
public:

    /**
     * @brief unitarityNLOev8 constructor.
     */
    unitarityNLOev8(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev8@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev9
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev9: public ThObservable {
public:

    /**
     * @brief unitarityNLOev9 constructor.
     */
    unitarityNLOev9(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev9@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev10
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev10: public ThObservable {
public:

    /**
     * @brief unitarityNLOev10 constructor.
     */
    unitarityNLOev10(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev10@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev11
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev11: public ThObservable {
public:

    /**
     * @brief unitarityNLOev11 constructor.
     */
    unitarityNLOev11(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev11@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev12
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev12: public ThObservable {
public:

    /**
     * @brief unitarityNLOev12 constructor.
     */
    unitarityNLOev12(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev12@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev13
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev13: public ThObservable {
public:

    /**
     * @brief unitarityNLOev13 constructor.
     */
    unitarityNLOev13(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev13@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev14
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev14: public ThObservable {
public:

    /**
     * @brief unitarityNLOev14 constructor.
     */
    unitarityNLOev14(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev14@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev15
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev15: public ThObservable {
public:

    /**
     * @brief unitarityNLOev15 constructor.
     */
    unitarityNLOev15(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev15@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev16
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev16: public ThObservable {
public:

    /**
     * @brief unitarityNLOev16 constructor.
     */
    unitarityNLOev16(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev16@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev17
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev17: public ThObservable {
public:

    /**
     * @brief unitarityNLOev17 constructor.
     */
    unitarityNLOev17(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev17@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityNLOev18
 * @ingroup THDM
 * @brief .
 */
class unitarityNLOev18: public ThObservable {
public:

    /**
     * @brief unitarityNLOev18 constructor.
     */
    unitarityNLOev18(const StandardModel& SM_i);

    /**
     * @return @f$unitarityNLOev18@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp1
 * @ingroup THDM
 * @brief .
 */
class unitarityRp1: public ThObservable {
public:

    /**
     * @brief unitarityRp1 constructor.
     */
    unitarityRp1(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp1@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp2
 * @ingroup THDM
 * @brief .
 */
class unitarityRp2: public ThObservable {
public:

    /**
     * @brief unitarityRp2 constructor.
     */
    unitarityRp2(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp2@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp3
 * @ingroup THDM
 * @brief .
 */
class unitarityRp3: public ThObservable {
public:

    /**
     * @brief unitarityRp3 constructor.
     */
    unitarityRp3(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp3@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp4
 * @ingroup THDM
 * @brief .
 */
class unitarityRp4: public ThObservable {
public:

    /**
     * @brief unitarityRp4 constructor.
     */
    unitarityRp4(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp4@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp5
 * @ingroup THDM
 * @brief .
 */
class unitarityRp5: public ThObservable {
public:

    /**
     * @brief unitarityRp5 constructor.
     */
    unitarityRp5(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp5@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp6
 * @ingroup THDM
 * @brief .
 */
class unitarityRp6: public ThObservable {
public:

    /**
     * @brief unitarityRp6 constructor.
     */
    unitarityRp6(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp6@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp7
 * @ingroup THDM
 * @brief .
 */
class unitarityRp7: public ThObservable {
public:

    /**
     * @brief unitarityRp7 constructor.
     */
    unitarityRp7(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp7@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp8
 * @ingroup THDM
 * @brief .
 */
class unitarityRp8: public ThObservable {
public:

    /**
     * @brief unitarityRp8 constructor.
     */
    unitarityRp8(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp8@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp9
 * @ingroup THDM
 * @brief .
 */
class unitarityRp9: public ThObservable {
public:

    /**
     * @brief unitarityRp9 constructor.
     */
    unitarityRp9(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp9@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp10
 * @ingroup THDM
 * @brief .
 */
class unitarityRp10: public ThObservable {
public:

    /**
     * @brief unitarityRp10 constructor.
     */
    unitarityRp10(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp10@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp11
 * @ingroup THDM
 * @brief .
 */
class unitarityRp11: public ThObservable {
public:

    /**
     * @brief unitarityRp11 constructor.
     */
    unitarityRp11(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp11@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp12
 * @ingroup THDM
 * @brief .
 */
class unitarityRp12: public ThObservable {
public:

    /**
     * @brief unitarityRp12 constructor.
     */
    unitarityRp12(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp12@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp13
 * @ingroup THDM
 * @brief .
 */
class unitarityRp13: public ThObservable {
public:

    /**
     * @brief unitarityRp13 constructor.
     */
    unitarityRp13(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp13@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp14
 * @ingroup THDM
 * @brief .
 */
class unitarityRp14: public ThObservable {
public:

    /**
     * @brief unitarityRp14 constructor.
     */
    unitarityRp14(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp14@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp15
 * @ingroup THDM
 * @brief .
 */
class unitarityRp15: public ThObservable {
public:

    /**
     * @brief unitarityRp15 constructor.
     */
    unitarityRp15(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp15@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp16
 * @ingroup THDM
 * @brief .
 */
class unitarityRp16: public ThObservable {
public:

    /**
     * @brief unitarityRp16 constructor.
     */
    unitarityRp16(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp16@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp17
 * @ingroup THDM
 * @brief .
 */
class unitarityRp17: public ThObservable {
public:

    /**
     * @brief unitarityRp17 constructor.
     */
    unitarityRp17(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp17@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp18
 * @ingroup THDM
 * @brief .
 */
class unitarityRp18: public ThObservable {
public:

    /**
     * @brief unitarityRp18 constructor.
     */
    unitarityRp18(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp18@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp19
 * @ingroup THDM
 * @brief .
 */
class unitarityRp19: public ThObservable {
public:

    /**
     * @brief unitarityRp19 constructor.
     */
    unitarityRp19(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp19@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp20
 * @ingroup THDM
 * @brief .
 */
class unitarityRp20: public ThObservable {
public:

    /**
     * @brief unitarityRp20 constructor.
     */
    unitarityRp20(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp20@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp21
 * @ingroup THDM
 * @brief .
 */
class unitarityRp21: public ThObservable {
public:

    /**
     * @brief unitarityRp21 constructor.
     */
    unitarityRp21(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp21@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityRp22
 * @ingroup THDM
 * @brief .
 */
class unitarityRp22: public ThObservable {
public:

    /**
     * @brief unitarityRp22 constructor.
     */
    unitarityRp22(const StandardModel& SM_i);

    /**
     * @return @f$unitarityRp22@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR1
 * @ingroup THDM
 * @brief .
 */
class unitarityR1: public ThObservable {
public:

    /**
     * @brief unitarityR1 constructor.
     */
    unitarityR1(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR1@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR2
 * @ingroup THDM
 * @brief .
 */
class unitarityR2: public ThObservable {
public:

    /**
     * @brief unitarityR2 constructor.
     */
    unitarityR2(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR2@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR3
 * @ingroup THDM
 * @brief .
 */
class unitarityR3: public ThObservable {
public:

    /**
     * @brief unitarityR3 constructor.
     */
    unitarityR3(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR3@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR4
 * @ingroup THDM
 * @brief .
 */
class unitarityR4: public ThObservable {
public:

    /**
     * @brief unitarityR4 constructor.
     */
    unitarityR4(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR4@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR5
 * @ingroup THDM
 * @brief .
 */
class unitarityR5: public ThObservable {
public:

    /**
     * @brief unitarityR5 constructor.
     */
    unitarityR5(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR5@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR6
 * @ingroup THDM
 * @brief .
 */
class unitarityR6: public ThObservable {
public:

    /**
     * @brief unitarityR6 constructor.
     */
    unitarityR6(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR6@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR7
 * @ingroup THDM
 * @brief .
 */
class unitarityR7: public ThObservable {
public:

    /**
     * @brief unitarityR7 constructor.
     */
    unitarityR7(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR7@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR8
 * @ingroup THDM
 * @brief .
 */
class unitarityR8: public ThObservable {
public:

    /**
     * @brief unitarityR8 constructor.
     */
    unitarityR8(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR8@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR9
 * @ingroup THDM
 * @brief .
 */
class unitarityR9: public ThObservable {
public:

    /**
     * @brief unitarityR9 constructor.
     */
    unitarityR9(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR9@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR10
 * @ingroup THDM
 * @brief .
 */
class unitarityR10: public ThObservable {
public:

    /**
     * @brief unitarityR10 constructor.
     */
    unitarityR10(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR10@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR11
 * @ingroup THDM
 * @brief .
 */
class unitarityR11: public ThObservable {
public:

    /**
     * @brief unitarityR11 constructor.
     */
    unitarityR11(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR11@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR12
 * @ingroup THDM
 * @brief .
 */
class unitarityR12: public ThObservable {
public:

    /**
     * @brief unitarityR12 constructor.
     */
    unitarityR12(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR12@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR13
 * @ingroup THDM
 * @brief .
 */
class unitarityR13: public ThObservable {
public:

    /**
     * @brief unitarityR13 constructor.
     */
    unitarityR13(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR13@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR14
 * @ingroup THDM
 * @brief .
 */
class unitarityR14: public ThObservable {
public:

    /**
     * @brief unitarityR14 constructor.
     */
    unitarityR14(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR14@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR15
 * @ingroup THDM
 * @brief .
 */
class unitarityR15: public ThObservable {
public:

    /**
     * @brief unitarityR15 constructor.
     */
    unitarityR15(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR15@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR16
 * @ingroup THDM
 * @brief .
 */
class unitarityR16: public ThObservable {
public:

    /**
     * @brief unitarityR16 constructor.
     */
    unitarityR16(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR16@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR17
 * @ingroup THDM
 * @brief .
 */
class unitarityR17: public ThObservable {
public:

    /**
     * @brief unitarityR17 constructor.
     */
    unitarityR17(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR17@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR18
 * @ingroup THDM
 * @brief .
 */
class unitarityR18: public ThObservable {
public:

    /**
     * @brief unitarityR18 constructor.
     */
    unitarityR18(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR18@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR19
 * @ingroup THDM
 * @brief .
 */
class unitarityR19: public ThObservable {
public:

    /**
     * @brief unitarityR19 constructor.
     */
    unitarityR19(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR19@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR20
 * @ingroup THDM
 * @brief .
 */
class unitarityR20: public ThObservable {
public:

    /**
     * @brief unitarityR20 constructor.
     */
    unitarityR20(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR20@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR21
 * @ingroup THDM
 * @brief .
 */
class unitarityR21: public ThObservable {
public:

    /**
     * @brief unitarityR21 constructor.
     */
    unitarityR21(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR21@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitarityR22
 * @ingroup THDM
 * @brief .
 */
class unitarityR22: public ThObservable {
public:

    /**
     * @brief unitarityR22 constructor.
     */
    unitarityR22(const StandardModel& SM_i);

    /**
     * @return @f$unitarityR22@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00evenpRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya00evenpRe: public ThObservable {
public:

    /**
     * @brief unitaritya00evenpRe constructor.
     */
    unitaritya00evenpRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00evenpRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00evenpIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya00evenpIm: public ThObservable {
public:

    /**
     * @brief unitaritya00evenpIm constructor.
     */
    unitaritya00evenpIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00evenpIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00evenmRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya00evenmRe: public ThObservable {
public:

    /**
     * @brief unitaritya00evenmRe constructor.
     */
    unitaritya00evenmRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00evenmRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00evenmIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya00evenmIm: public ThObservable {
public:

    /**
     * @brief unitaritya00evenmIm constructor.
     */
    unitaritya00evenmIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00evenmIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00oddpRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya00oddpRe: public ThObservable {
public:

    /**
     * @brief unitaritya00oddpRe constructor.
     */
    unitaritya00oddpRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00oddpRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00oddpIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya00oddpIm: public ThObservable {
public:

    /**
     * @brief unitaritya00oddpIm constructor.
     */
    unitaritya00oddpIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00oddpIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00oddmRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya00oddmRe: public ThObservable {
public:

    /**
     * @brief unitaritya00oddmRe constructor.
     */
    unitaritya00oddmRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00oddmRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya00oddmIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya00oddmIm: public ThObservable {
public:

    /**
     * @brief unitaritya00oddmIm constructor.
     */
    unitaritya00oddmIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya00oddmIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01evenpRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya01evenpRe: public ThObservable {
public:

    /**
     * @brief unitaritya01evenpRe constructor.
     */
    unitaritya01evenpRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01evenpRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01evenpIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya01evenpIm: public ThObservable {
public:

    /**
     * @brief unitaritya01evenpIm constructor.
     */
    unitaritya01evenpIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01evenpIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01evenmRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya01evenmRe: public ThObservable {
public:

    /**
     * @brief unitaritya01evenmRe constructor.
     */
    unitaritya01evenmRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01evenmRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01evenmIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya01evenmIm: public ThObservable {
public:

    /**
     * @brief unitaritya01evenmIm constructor.
     */
    unitaritya01evenmIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01evenmIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01oddpRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya01oddpRe: public ThObservable {
public:

    /**
     * @brief unitaritya01oddpRe constructor.
     */
    unitaritya01oddpRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01oddpRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01oddpIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya01oddpIm: public ThObservable {
public:

    /**
     * @brief unitaritya01oddpIm constructor.
     */
    unitaritya01oddpIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01oddpIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01oddmRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya01oddmRe: public ThObservable {
public:

    /**
     * @brief unitaritya01oddmRe constructor.
     */
    unitaritya01oddmRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01oddmRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya01oddmIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya01oddmIm: public ThObservable {
public:

    /**
     * @brief unitaritya01oddmIm constructor.
     */
    unitaritya01oddmIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya01oddmIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya10oddRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya10oddRe: public ThObservable {
public:

    /**
     * @brief unitaritya10oddRe constructor.
     */
    unitaritya10oddRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya10oddRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya10oddIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya10oddIm: public ThObservable {
public:

    /**
     * @brief unitaritya10oddIm constructor.
     */
    unitaritya10oddIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya10oddIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya11evenpRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya11evenpRe: public ThObservable {
public:

    /**
     * @brief unitaritya11evenpRe constructor.
     */
    unitaritya11evenpRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya11evenpRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya11evenpIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya11evenpIm: public ThObservable {
public:

    /**
     * @brief unitaritya11evenpIm constructor.
     */
    unitaritya11evenpIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya11evenpIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya11evenmRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya11evenmRe: public ThObservable {
public:

    /**
     * @brief unitaritya11evenmRe constructor.
     */
    unitaritya11evenmRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya11evenmRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya11evenmIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya11evenmIm: public ThObservable {
public:

    /**
     * @brief unitaritya11evenmIm constructor.
     */
    unitaritya11evenmIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya11evenmIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya11oddRe
 * @ingroup THDM
 * @brief .
 */
class unitaritya11oddRe: public ThObservable {
public:

    /**
     * @brief unitaritya11oddRe constructor.
     */
    unitaritya11oddRe(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya11oddRe@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class unitaritya11oddIm
 * @ingroup THDM
 * @brief .
 */
class unitaritya11oddIm: public ThObservable {
public:

    /**
     * @brief unitaritya11oddIm constructor.
     */
    unitaritya11oddIm(const StandardModel& SM_i);

    /**
     * @return @f$unitaritya11oddIm@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

#endif	/* UNITARITY_H */

