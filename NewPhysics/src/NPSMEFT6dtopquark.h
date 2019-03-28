/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFT6DTOPQUARK_H
#define NPSMEFT6DTOPQUARK_H


#include "gslpp.h"
#include "ThObservable.h"
#include "NPbase.h"
#include <stdexcept>

class NPSMEFT6dtopquark : public NPbase {
public:
    
    
/**
*　@brief The number of the model parameters in %NPSMEFT6dtopquarkVars. 
*/
static const int NNPSMEFT6dtopquarkVars = 9;
   
/*
@brief A string array containing the labels of the model parameters in NPSMEFT6dtopquark 
*/
    static const std::string NPSMEFT6dtopquarkVars[NNPSMEFT6dtopquarkVars];
    
    
   NPSMEFT6dtopquark();
   
   
       /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the SM value (FlagEpsilon1SM=true) or the SM plus new physics
     * value (FlagEpsilon1SM=false) of \f$\varepsilon_1\f$
     */
    //virtual double myC_phit() const;
    
    /**
     *
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    
    double getNPSMEFT6dtopquark_C_phit() const {
        return C_phit;
    }
    
    double getNPSMEFT6dtopquark_C_phiQ3() const {
        return C_phiQ3;
    }
    
    double getNPSMEFT6dtopquark_C_phiQ1() const {
        return C_phiQ1;
    }
    
    double getNPSMEFT6dtopquark_C_tW() const {
        return C_tW;
    }
    
    double getNPSMEFT6dtopquark_C_tB() const {
        return C_tB;
    }
    
    double getNPSMEFT6dtopquark_C_tphi() const {
        return C_tphi;
    }
    
    double getNPSMEFT6dtopquark_C_phib() const {
        return C_phib;
    }
    
    double getNPSMEFT6dtopquark_C_bW() const {
        return C_bW;
    }
    
    double getNPSMEFT6dtopquark_C_bB() const {
        return C_bB;
    }
    
    double getNPSMEFT6dtopquark_flag_Quadratic() const {
        return flag_Quadratic;
    }
    



    /////////////////////////////////////////////////////////////////
    protected:
    virtual void setParameter(const std::string name, const double& value);
        //mutable Matching<NPSMEFT6dtopquarkMatching,NPSMEFT6dtopquark> NPSMEFT6dtopquarkM;
    double C_phit; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQ3; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQ1; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tW; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tB; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tphi; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phib; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bW; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bB; ///< The dimension-6 operator coefficient \f$C_{G}\f$.

    virtual bool setFlag(const std::string, const bool);

        
////////////////////////////////////////////////////////////////////////
private:

    bool flag_Quadratic;



    

};


/**
 * @class C_phit
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_phit: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_phit(const StandardModel& SM_i);

    /**
     * @return The value of C_phit
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_phiQ3
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phiQ3
 */
class C_phiQ3: public ThObservable {
public:

    /**
     * @brief C_phiQ3 constructor.
     */
    C_phiQ3(const StandardModel& SM_i);

    /**
     * @return The value of C_phiQ3
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_phiQ1
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_phiQ1: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_phiQ1(const StandardModel& SM_i);

    /**
     * @return The value of C_phiQ1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



/**
 * @class C_tW
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tW
 */
class C_tW: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_tW(const StandardModel& SM_i);

    /**
     * @return The value of C_tW
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



/**
 * @class C_tB
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_tB: public ThObservable {
public:

    /**
     * @brief C_tB constructor.
     */
    C_tB(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



/**
 * @class C_tphi
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tphi
 */
class C_tphi: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_tphi(const StandardModel& SM_i);

    /**
     * @return The value of C_tphi
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



/**
 * @class C_phib
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_phib: public ThObservable {
public:

    /**
     * @brief C_phib constructor.
     */
    C_phib(const StandardModel& SM_i);

    /**
     * @return The value of C_phib
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_bW
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_bW: public ThObservable {
public:

    /**
     * @brief C_bW constructor.
     */
    C_bW(const StandardModel& SM_i);

    /**
     * @return The value of C_bW
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_bB
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_bB: public ThObservable {
public:

    /**
     * @brief C_bB constructor.
     */
    C_bB(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


//Observables from LEP1

class Rb_NPSMEFT6dtopquark : public ThObservable {
public:   

    Rb_NPSMEFT6dtopquark(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class AFBLR : public ThObservable {
public:   

    AFBLR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


//Observables from LHC Run 2 (also prospects at High lumniosity)

class sigmattZ : public ThObservable {
public:   

    sigmattZ(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigmattA_1 : public ThObservable {
public:   

    sigmattA_1(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmattA_2 : public ThObservable {
public:   

    sigmattA_2(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmattH : public ThObservable {
public:   

    sigmattH(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmattW : public ThObservable {
public:   

    sigmattW(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmatq : public ThObservable {
public:   

    sigmatq(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmatW : public ThObservable {
public:   

    sigmatW(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmatqZ : public ThObservable {
public:   

    sigmatqZ(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class F0 : public ThObservable {
public:   

    F0(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class FL : public ThObservable {
public:   

    FL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


//Prospects of Linear Colders at 250 GeV
//250 bb observables


class sigma_250_bb_eLpR : public ThObservable {
public:   

    sigma_250_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_250_bb_eLpR : public ThObservable {
public:   

    a_250_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_250_bb_eRpL : public ThObservable {
public:   

    sigma_250_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_250_bb_eRpL : public ThObservable {
public:   

    a_250_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

//Prospects of Linear Colders at 500 GeV
//500 bb observables


class sigma_500_bb_eLpR : public ThObservable {
public:   

    sigma_500_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_500_bb_eLpR : public ThObservable {
public:   

    a_500_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_500_bb_eRpL : public ThObservable {
public:   

    sigma_500_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_500_bb_eRpL : public ThObservable {
public:   

    a_500_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




//500 tt observables

class sigma_500_tt_eLpR : public ThObservable {
public:   

    sigma_500_tt_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_500_tt_eLpR : public ThObservable {
public:   

    a_500_tt_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_500_tt_eRpL : public ThObservable {
public:   

    sigma_500_tt_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_500_tt_eRpL : public ThObservable {
public:   

    a_500_tt_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class pt_500_tt_eLpR : public ThObservable {
public:   

    pt_500_tt_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class pt_500_tt_eRpL : public ThObservable {
public:   

    pt_500_tt_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

    //OPTIMIZED OBSERVABLES


class op1 : public ThObservable {
public:   

    op1(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class op2 : public ThObservable {
public:   

    op2(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class op3 : public ThObservable {
public:   

    op3(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op4 : public ThObservable {
public:   

    op4(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



#endif /* NPSMEFT6DTOPQUARK_H */
