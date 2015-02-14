/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVLLOBSERVABLES_H
#define	MVLLOBSERVABLES_H

#include "MVll.h"
#include <StandardModel.h>
#include <ThObservable.h>


/**
 * @class P_1
 * @ingroup flavour
 * @brief A class for the clean observable P_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_1 : public ThObservable{
public:
    
    /**
    * @brief \f$ P_{1} \f$ 
    */
    P_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    
    /**
    * @return return the clean observable P_1
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    
};


/**
 * @class P_2
 * @ingroup flavour
 * @brief A class for the clean observable P_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_2 : public ThObservable {
public:
    
    /**
    * @brief \f$ P_{2} \f$ 
    */
    P_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable P_2
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    
};


/**
 * @class P_3
 * @ingroup flavour
 * @brief A class for the clean observable P_3. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_3 : public ThObservable{
public:
    
    /**
    * @brief \f$ P_{3} \f$ 
    */
    P_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable P_3
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class P_4Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_4. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_4Prime : public ThObservable{
public:
    
    /**
    * @brief \f$ P'_{4} \f$ 
    */
    P_4Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable P'_4
    */
    double computeThValue ();
    
protected:
    
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
};


/**
 * @class P_5Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_5. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_5Prime : public ThObservable{
public:
    
    /**
    * @brief \f$ P'_{5} \f$ 
    */

    P_5Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable P'_5
    */
    double computeThValue ();
    
protected:
    
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class P_6Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_6. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_6Prime : public ThObservable{
public:
    
    /**
    * @brief \f$ P'_{6} \f$ 
    */
    P_6Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @return return the clean observable P'_6
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class P_8Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_8. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_8Prime : public ThObservable{
public:
    
    /**
    * @brief \f$ P'_{8} \f$ 
    */
    P_8Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @return return the clean observable P'_8
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class GammaPrime
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GammaPrime : public ThObservable{
public:
    
    /**
    * @brief \f$ Gamma' \f$ 
    */
    GammaPrime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    double computeGammaPrime(double qmin, double qmax, StandardModel::lepton lep);
    
    /**
    * @return return the clean observable Gamma'
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class A_FB
 * @ingroup flavour
 * @brief A class for the clean observable A_{FB}. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class A_FB : public GammaPrime{
public:
    
    /**
    * @brief \f$ A_{FB} \f$ 
    */
    A_FB(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable A_FB
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class BF
 * @ingroup flavour
 * @brief A class for the Branching Fraction. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BR_MVll : public GammaPrime{
public:
    
    /**
    * @brief \f$BR_{B \to K^* l^+l^-}\f$ 
    */
    BR_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return the branching fraction of \f$B\to K^* l^+l^-\f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    
};


/**
 * @class F_L
 * @ingroup flavour
 * @brief A class for the clean observable F_L. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class F_L : public GammaPrime{
public:
    
    /**
    * @brief \f$ F_{L} \f$ 
    */
    F_L(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    double computeFL(double qmin, double qmax, StandardModel::lepton lep);

    
    /**
    * @return return the clean observable F_L
    */
    double computeThValue ();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class R_K^*
 * @ingroup flavour
 * @brief A class for the Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class R_MVll : public GammaPrime{
public:
    
    /**
    * @brief \f$ BR_{mu}/BR_{e} \f$ 
    */
    R_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @return the ratio between branching fractions of \f$ B\to K^* \mu^+ \mu^- \f$ and \f$ B\to K^* e^+ e^- \f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep1;
    StandardModel::lepton lep2;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
};


/**
 * @class R_K^*_L
 * @ingroup flavour
 * @brief A class for the longitudinal Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class RL_MVll : public F_L{
public:
    
    /**
    * @brief \f$ BR_{mu,L}/BR_{e,L} \f$ 
    */
    RL_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @return the ratio between longitudinal branching fractions of \f$ B\to K^* \mu^+ \mu^- \f$ and \f$ B\to K^* e^+ e^- \f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep1;
    StandardModel::lepton lep2;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
};


/**
 * @class R_K^*_T
 * @ingroup flavour
 * @brief A class for the transverse Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class RT_MVll : public F_L{
public:
    
    /**
    * @brief \f$ BR_{mu/,TBR_{e,T} \f$ 
    */
    RT_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @return the ratio between transverse branching fractions of \f$ B\to K^* \mu^+ \mu^- \f$ and \f$ B\to K^* e^+ e^- \f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep1;
    StandardModel::lepton lep2;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
};


/**
 * @class R_6
 * @ingroup flavour
 * @brief A class for the Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class R_6 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Sigma_{6,mu} / \Sigma_{6,e} \f$ 
    */
    R_6(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @return the ratio between Sigma_6 of \f$ B\to K^* \mu^+ \mu^- \f$ and \f$ B\to K^* e^+ e^- \f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep1;
    StandardModel::lepton lep2;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
};


/**
 * @class ACP
 * @ingroup flavour
 * @brief A class for the clean observable ACP. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ACP_MVll : public GammaPrime{
public:
    
    /**
    * @brief \f$ A_{CP} \f$ 
    */
    ACP_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable ACP
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class P3CP
 * @ingroup flavour
 * @brief A class for the clean observable P3CP. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P3CP : public ThObservable{
public:
    
    /**
    * @brief \f$ P_3^{CP} \f$ 
    */
    P3CP(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable P3CP
    */
    double computeThValue ();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class M_1Prime
 * @ingroup flavour
 * @brief A class for the clean observable M'_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class M_1Prime : public ThObservable{
public:
    
    /**
    * @brief \f$ M'_{1} \f$ 
    */
    M_1Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable M'_1
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class M_2Prime
 * @ingroup flavour
 * @brief A class for the clean observable M'_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class M_2Prime : public ThObservable{
public:
    
    /**
    * @brief \f$ M'_{2} \f$ 
    */
    M_2Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the clean observable M'_2
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class S_3
 * @ingroup flavour
 * @brief A class for the observable S_3. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class S_3 : public GammaPrime{
public:
    
    /**
    * @brief \f$ S_{3} \f$ 
    */
    S_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable S_3
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class S_4
 * @ingroup flavour
 * @brief A class for the observable S_4. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class S_4: public GammaPrime{
public:
    
    /**
    * @brief \f$ S_{4} \f$ 
    */
    S_4(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable S_4
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class S_5
 * @ingroup flavour
 * @brief A class for the observable S_5. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class S_5 : public GammaPrime{
public:
    
    /**
    * @brief \f$ S_{5} \f$ 
    */
    S_5(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable S_5
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class S_7
 * @ingroup flavour
 * @brief A class for the observable S_7. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class S_7 : public GammaPrime{
public:
    
    /**
    * @brief \f$ S_{7} \f$ 
    */
    S_7(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable S_7
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class S_8
 * @ingroup flavour
 * @brief A class for the observable S_8. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class S_8 : public GammaPrime{
public:
    
    /**
    * @brief \f$ S_{8} \f$ 
    */
    S_8(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable S_8
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class S_9
 * @ingroup flavour
 * @brief A class for the observable S_9. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class S_9 : public GammaPrime{
public:
    
    /**
    * @brief \f$ S_{9} \f$ 
    */
    S_9(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable S_9
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class A_6
 * @ingroup flavour
 * @brief A class for the observable A_6. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class A_6 : public GammaPrime{
public:
    
    /**
    * @brief \f$ A_{6} \f$ 
    */
    A_6(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable A_6
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class A_9
 * @ingroup flavour
 * @brief A class for the observable A_9. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class A_9 : public GammaPrime{
public:
    
    /**
    * @brief \f$ A_{9} \f$ 
    */
    A_9(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the observable A_9
    */
    double computeThValue ();
   
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class Vp_Vm
 * @ingroup flavour
 * @brief A class for the ratio V_+/V_-. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Vp_Vm : public ThObservable{
public:
    
    /**
    * @brief \f$ V_{+}/V_{-} \f$ 
    */
    Vp_Vm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the ratio V_+/V_-
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class Tp_Tm
 * @ingroup flavour
 * @brief A class for the ratio T_+/T_-. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Tp_Tm : public ThObservable{
public:
    
    /**
    * @brief \f$ T_{+}/T_{-} \f$ 
    */
    Tp_Tm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the ratio T_+/T_-
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class V0_m_T0_V0_p_T0
 * @ingroup flavour
 * @brief A class for the ratio (V_0 - T_0) / (V_0 + T_0). 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class V0_m_T0_V0_p_T0 : public ThObservable{
public:
    
    /**
    * @brief \f$ (V_{0} - T_{0}) / (V_{0} + T_{0}) \f$ 
    */
    V0_m_T0_V0_p_T0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @return return the ratio (V_0 - T_0) / (V_0 + T_0)
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class Vm_m_Tm_Vm_p_Tm
 * @ingroup flavour
 * @brief A class for the ratio (V_- - T_-) / (V_- + T_-). 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Vm_m_Tm_Vm_p_Tm : public ThObservable{
public:
    
    /**
    * @brief \f$ (V_{-} - T_{-}) / (V_{-} + T_{-}) \f$ 
    */
    Vm_m_Tm_Vm_p_Tm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @return return the ratio (V_- - T_-) / (V_- + T_-)
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};


/**
 * @class V0_m_S_V0_p_S
 * @ingroup flavour
 * @brief A class for the ratio (V_0 - S) / (V_0 + S). 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class V0_m_S_V0_p_S : public ThObservable{
public:
    
    /**
    * @brief \f$ (V_{0} - S) / (V_{0} + S) \f$ 
    */
    V0_m_S_V0_p_S(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @return return the ratio (V_0 - S) / (V_0 + S)
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;

};



/**
 * @class Delta_C9_1
 * @ingroup flavour
 * @brief A class for the correction to C9 Delta_C9_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Delta_C9_1 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Delta_{C9}^1 \f$ 
    */
    Delta_C9_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @return return the factor Delta_C9_1
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    unsigned int typ;
};



/**
 * @class Delta_C9_2
 * @ingroup flavour
 * @brief A class for the correction to C9 Delta_C9_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Delta_C9_2 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Delta_{C9}^2 \f$ 
    */
    Delta_C9_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @return return the factor Delta_C9_2
    */
    double computeThValue();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    unsigned int typ;
};
    
    
    
    /**
 * @class Delta_C9_3
 * @ingroup flavour
 * @brief A class for the correction to C9 Delta_C9_3. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Delta_C9_3 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Delta_{C9}^3 \f$ 
    */
    Delta_C9_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @return return the factor Delta_C9_3
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    unsigned int typ;
};

/**
 * @class Delta_C7_1
 * @ingroup flavour
 * @brief A class for the correction to C7 Delta_C7_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Delta_C7_1 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Delta_{C7}^1 \f$ 
    */
    Delta_C7_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @return return the factor Delta_C7_1
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    unsigned int typ;
};



/**
 * @class Delta_C7_2
 * @ingroup flavour
 * @brief A class for the correction to C7 Delta_C7_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Delta_C7_2 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Delta_{C7}^2 \f$ 
    */
    Delta_C7_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @return return the factor Delta_C7_2
    */
    double computeThValue();
    
private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    unsigned int typ;
};
    
    
    
/**
 * @class Delta_C7_3
 * @ingroup flavour
 * @brief A class for the correction to C7 Delta_C7_3. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Delta_C7_3 : public ThObservable{
public:
    
    /**
    * @brief \f$ \Delta_{C7}^3 \f$ 
    */
    Delta_C7_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @return return the factor Delta_C7_3
    */
    double computeThValue();

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    unsigned int typ;
};


#endif	/* MVLLOBSERVABLES_H */

    