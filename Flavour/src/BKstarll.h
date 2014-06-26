/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BKSTARLL_H
#define	BKSTARLL_H

#include <math.h>
#include "StandardModelMatching.h"
#include <StandardModel.h>
#include <ThObservable.h>


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the decay B -> K^*ll. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BKstarll : public ThObservable {
public:
    BKstarll(const StandardModel& SM_i, int lep_i);
    virtual ~BKstarll();
    virtual double computeThValue()=0;
        /*dummy variables*/
    double GF;            //Fermi constant
    double ale;           //alpha electromagnetic
    double Mm;            //muon mass
    double MB;            //B meson mass
    double Mb;            //b quark mass
    double Ms;            //s quark mass
    double MW;            //W boson mass
    complex lambda_t;      //Vckm factor
    double h_lambda;      //parameter that contains the contribution from the hadronic hamiltonian  
    double q2;            //q^2 of the decay
    double k2;            //square of the 3-momentum k
    double b;             //BF of the decay K^* -> K pi
    double C7,C9,C7p,C9p; //WC
    double C10,C10p;      //WC
    double CS,CSp;        //WC
    
    
    
    
    
    /**
    * @brief \f$ V_L \f$ ; NOT CODED, SHOULD BE ADDED AS INPUT FROM LATTICE
    * @param[in] i polarization
    * @return return the matrix element V_L(lambda)
    */
    double V_L(int i);

    
    /**
    * @brief \f$ V_L \f$ ; NOT CODED, SHOULD BE ADDED AS INPUT FROM LATTICE
    * @param[in] i polarization
    * @return return the matrix element V_R(lambda)
    */
    double V_R(int i);


    /**
    * @brief \f$ V_L \f$ ; NOT CODED, SHOULD BE ADDED AS INPUT FROM LATTICE
    * @param[in] i polarization
    * @return return the matrix element T_L(lambda)
    */
    double T_L(int i);


    /**
    * @brief \f$ V_L \f$ ; NOT CODED, SHOULD BE ADDED AS INPUT FROM LATTICE 
    * @param[in] i polarization
    * @return return the matrix element T_R(lambda)
    */
    double T_R(int i);


    /**
    * 
    * @brief \f$ V_L \f$ ; NOT CODED, SHOULD BE ADDED AS INPUT FROM LATTICE 
    * @return return the matrix element S_L
    */
    double S_L;


    /**
    * 
    * @brief \f$ V_L \f$ ; NOT CODED, SHOULD BE ADDED AS INPUT FROM LATTICE 
    * @return return the matrix element S_R
    */
    double S_R;
    
    
    /**
    * @brief \f$ N \f$ 
    * @return return the helicity amplitude normalization factor N
    */
    
    
    /**
    * @brief \f$ H_V \f$ 
    * @param[in] i polarization lambda
    * @return return the helicity amplitude H_V(lambda)
    */
    gslpp::complex H_V(int i);


    /**
    * @brief \f$ H_A \f$ 
    * @param[in] i polarization lambda
    * @return return the helicity amplitude H_A(lambda)
    */
    gslpp::complex H_A(int i);


    /**
    * @brief \f$ H_S \f$ 
    * @return return the helicity amplitude H_S
    */
    gslpp::complex H_S();


    /**
    * @brief \f$ H_P \f$ 
    * @return return the helicity amplitude H_P
    */
    gslpp::complex H_P();
    
    
    /**
    * @brief \f$ beta \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor beta used in the angular coefficients I_i
    */
    double beta (double q2);
    
    
    /**
    * @brief \f$ lambda \f$ 
    * @param[in] k2 square of the 3-momentum k of the recoiling meson in the B rest frame
    * @return return the factor lambda used in the angular coefficients I_i
    */
    double lambda(double k2);

    
    /**
    * @brief \f$ F \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] b BF of the decay K* -> K pi
    * @return return the factor F used in the angular coefficients I_i
    */
    double F(double q2, double b);
    
    
    /**
    * i values:
    * 0 = 1c
    * 1 = 1s
    * 2 = 2c
    * 3 = 2s
    * 4 = 3
    * 5 = 4
    * 6 = 5
    * 7 = 6s
    * 8 = 6c
    * 9 = 7
    * 10 = 8
    * 11 = 9
    */
    
    
    /**
    * @brief \f$ I_{i} \f$ 
    * @param[in] i index of the angular coefficient
    * @return return the angular coefficient I_i
    */
    double  I(int i);
    
    
    /**
    * @brief \f$ I_{i} \f$ 
    * @param[in] i index of the angular coefficient
    * @return return the angular coefficient \bar{I}_i
    */
    double I_bar(int i);
    
    
    /**
    * @brief \f$ Sigma_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @return return the CP avarage Sigma_i
    */
    double Sigma(int i);

    

private:
    const StandardModel& mySM;

};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable P_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_1 : public BKstarll{
public:
    
    /**
    * @brief \f$ P_{1} \f$ 
    */
    complex N();
    
    /**
    * @return return the clean observable P_1
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable P_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_2 : public BKstarll{
public:
    
    /**
    * @brief \f$ P_{2} \f$ 
    */
    P_1(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable P_2
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable P_3. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_3 : public BKstarll{
public:
    
    /**
    * @brief \f$ P_{3} \f$ 
    */
    P_2(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable P_3
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable P'_4. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_4Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ P'_{4} \f$ 
    */
    P_3(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable P'_4
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable P'_5. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_5Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ P'_{5} \f$ 
    */
    P_4Prime(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable P'_5
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable P'_6. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_6Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ P'_{6} \f$ 
    */
    P_5Prime(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable P'_6
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GammaPrime : public BKstarll{
public:
    
    /**
    * @brief \f$ Gamma' \f$ 
    */
    P_6Prime(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable Gamma'
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable F_L. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class F_L : public BKstarll{
public:
    
    /**
    * @brief \f$ F_L \f$ 
    */
    GammaPrime(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable F_L
    */
    double computeThValue ();
    
    const StandardModel& mySM;
    F_L(const StandardModel& SM_i, int lep_i);
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable M'_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class M_1Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ M'_1 \f$ 
    */
    int lep;
    
    /**
    * @return return the clean observable M'_1
    */
    double computeThValue ();
};


/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the clean observable M'_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class M_2Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ M'_2 \f$ 
    */
    M_1Prime(const StandardModel& SM_i, int lep_i);
    
    /**
    * @return return the clean observable M'_2
    */
    double computeThValue ();
};


#endif	/* BKSTARLL_H */

    M_2Prime(const StandardModel& SM_i, int lep_i);
