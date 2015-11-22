/* 
 * Copyright (C) 2014 HEPfit Collaboration
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
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_1>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_1>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_1>=\frac {<\Sigma_3>}{2<\Sigma_{2s}>}\,.
 * @f]
 */
class P_1 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    
    /**
    * @brief The binned observable @f$<P_1>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_1>@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class P_2
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_2>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_2>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_2>=\frac {<\Sigma_6>}{8<\Sigma_{2s}>}\,.
 * @f]
 */
class P_2 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<P_2>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_2>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class P_3
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_3>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_3>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_3>=- \frac {<\Sigma_9>}{4<\Sigma_{2s}>}\,.
 * @f]
 */
class P_3 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<P_3>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_3>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class P_4Prime
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_4'>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_4'>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_4'>=\frac {<\Sigma_4>}{\sqrt{-<\Sigma_{2s}><\Sigma_{2c}>}}\,.
 * @f]
 */
class P_4Prime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_4Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<P_4'>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_4'>@f$
    */
    double computeThValue ();
    
protected:
    
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};


/**
 * @class P_5Prime
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_5'>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_5'>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_5'>=\frac {<\Sigma_5>}{2\sqrt{-<\Sigma_{2s}><\Sigma_{2c}>}}\,.
 * @f]
 */
class P_5Prime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_5Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<P_5'>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_5'>@f$
    */
    double computeThValue ();
    
protected:
    
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class P_6Prime
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_6'>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_6'>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_6'>=- \frac {<\Sigma_7>}{2 \sqrt{-<\Sigma_{2s}><\Sigma_{2c}>}}\,.
 * @f]
 */
class P_6Prime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_6Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The binned observable @f$<P_6'>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_6'>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class P_8Prime
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_8'>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_8'>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <P_8'>=- \frac {<\Sigma_8>}{\sqrt{-<\Sigma_{2s}><\Sigma_{2c}>}}\,.
 * @f]
 */
class P_8Prime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_8Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The binned observable @f$<P_8'>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_8'>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class GammaPrime
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<\Gamma'>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<\Gamma'>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <\Gamma'>=- \frac {1}{4} (<(3\Sigma_{1c} - \Sigma_{2c}) + 2(3\Sigma_{1s} - \Sigma_{2s})>)\,.
 * @f]
 */
class GammaPrime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    GammaPrime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
     * @brief A method to compute the binned observable @f$<\Gamma'>@f$ in @f$M \to V l^+l^-@f$ in a generic bin.
     * @param[in] qmin minimal value of the bin
     * @param[in] qmax maximal value of the bin
     * @param[in] lep final leptons of the decay
     * @return @f$<\Gamma'>_{[qmin,qmax]}@f$
     */
    double computeGammaPrime(double qmin, double qmax, StandardModel::lepton lep);
    
    /**
    * @brief The binned observable @f$<\Gamma'>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<\Gamma'>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_FB
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{FB}>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<A_{FB}>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <A_{FB}>=- \frac {3<\Sigma_{6s}>}{4<\Gamma'>} \,.
 * @f]
 */
class A_FB : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_FB(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{FB}>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<A_{FB}>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class BR_MVll
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<BR>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<BR>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class, and the meson width @f$W_M@f$:
 * @f[
 * <BR>= \frac {<\Gamma'>}{W_M} \,.
 * @f]
 */
class BR_MVll : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    BR_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<BR>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<BR>@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class F_L
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<F_L>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<F_L>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class, and the meson width @f$M_W@f$:
 * @f[
 * <BR>= \frac {3<\Sigma_{1c}>-<\Sigma_{2c}>}{4<\Gamma'>} \,.
 * @f]
 */
class F_L : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    F_L(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
     * @brief A method to compute the binned observable @f$<F_L>@f$ in @f$M \to V l^+l^-@f$ in a generic bin.
     * @param[in] qmin minimal value of the bin
     * @param[in] qmax maximal value of the bin
     * @param[in] lep final leptons of the decay
     * @return @f$<F_L>_{[qmin,qmax]}@f$
     */
    double computeFL(double qmin, double qmax, StandardModel::lepton lep);

    /**
    * @brief The binned observable @f$<F_L>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<F_L>@f$
    */
    double computeThValue ();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class R_MVll
 * @ingroup Flavour
 * @brief A class for the binned observable ratio @f$\frac {<BR>_{M \to V l_1^+l_1^-}}{<BR>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable ratio 
 * @f$\frac {<BR>_{M \to V l_1^+l_1^-}}{<BR>_{M \to V l_2^+l_2^-}}@f$ 
 * in @f$M \to V l^+l^-@f$ in terms of the binned observables @f$<\Gamma'>_{M \to V l_i^+l_i^-}@f$:
 * @f[
 * <R>_{M \to V l^+ l^-}= \frac {<\Gamma'>_{M \to V l_1^+l_1^-}}{<\Gamma'>_{M \to V l_2^+l_2^-}} \,.
 * @f]
 */
class R_MVll : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_1 first final leptons of the decay
     * @param[in] lep_1 second final leptons of the decay
     */
    R_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @brief The binned observable ratio @f$\frac {<BR>_{M \to V l_1^+l_1^-}}{<BR>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<R>_{M \to V l^+ l^-}@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep1; /**< First final leptons type. */
    StandardModel::lepton lep2; /**< Second final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};


/**
 * @class RL_MVll
 * @ingroup Flavour
 * @brief A class for the binned observable ratio @f$\frac {<BR_L>_{M \to V l_1^+l_1^-}}{<BR_L>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable longitudinal ratio 
 * @f$\frac {<BR_L>_{M \to V l_1^+l_1^-}}{<BR_L>_{M \to V l_2^+l_2^-}}@f$ 
 * in @f$M \to V l^+l^-@f$ in terms of the binned observables @f$<\Gamma'>_{M \to V l_i^+l_i^-}@f$
 * and @f$<{F_L}>_{M \to V l_i^+l_i^-}@f$:
 * @f[
 * <{R_L}>_{M \to V l^+ l^-}= \frac {<\Gamma'>_{M \to V l_1^+l_1^-}<{F_L}>_{M \to V l_1^+l_1^-}}{<\Gamma'>_{M \to V l_2^+l_2^-}<{F_L}>_{M \to V l_2^+l_2^-}} \,.
 * @f]
 */
class RL_MVll : public F_L{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_1 first final leptons of the decay
     * @param[in] lep_1 second final leptons of the decay
     */
    RL_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @brief The binned observable ratio @f$\frac {<BR_L>_{M \to V l_1^+l_1^-}}{<BR_L>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<R_L>_{M \to V l^+ l^-}@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep1; /**< First final leptons type. */
    StandardModel::lepton lep2; /**< Second final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};


/**
 * @class RT_MVll
 * @ingroup Flavour
 * @brief A class for the binned observable ratio @f$\frac {<BR_T>_{M \to V l_1^+l_1^-}}{<BR_T>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable transverse ratio 
 * @f$\frac {<BR_T>_{M \to V l_1^+l_1^-}}{<BR_T>_{M \to V l_2^+l_2^-}}@f$ 
 * in @f$M \to V l^+l^-@f$ in terms of the binned observables @f$<\Gamma'>_{M \to V l_i^+l_i^-}@f$
 * and @f$<{F_L}>_{M \to V l_i^+l_i^-}@f$:
 * @f[
 * <{R_T}>_{M \to V l^+ l^-}= \frac {<\Gamma'>_{M \to V l_1^+l_1^-}(1-<{F_L}>_{M \to V l_1^+l_1^-})}{<\Gamma'>_{M \to V l_2^+l_2^-}(1-<{F_L}>_{M \to V l_2^+l_2^-})} \,.
 * @f]
 */
class RT_MVll : public F_L{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_1 first final leptons of the decay
     * @param[in] lep_1 second final leptons of the decay
     */
    RT_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @brief The binned observable ratio @f$\frac {<BR_T>_{M \to V l_1^+l_1^-}}{<BR_T>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<R_T>_{M \to V l^+ l^-}@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep1; /**< First final leptons type. */
    StandardModel::lepton lep2; /**< Second final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};


/**
 * @class R_6
 * @ingroup Flavour
 * @brief A class for the binned observable ratio @f$\frac {<\Sigma_6>_{M \to V l_1^+l_1^-}}{<\Sigma_6>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable ratio 
 * @f$\frac {<\Sigma_6>_{M \to V l_1^+l_1^-}}{<\Sigma_6>_{M \to V l_2^+l_2^-}}@f$ 
 * in @f$M \to V l^+l^-@f$ in terms of the binned observables @f$<{\Sigma_6}>_{M \to V l_i^+l_i^-}@f$,
 * computed in the MVll class:
 * @f[
 * <{R_6}>_{M \to V l^+ l^-}= \frac {<\Sigma_6>_{M \to V l_1^+l_1^-}}{<\Sigma_6>_{M \to V l_2^+l_2^-}} \,.
 * @f]
 */
class R_6 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_1 first final leptons of the decay
     * @param[in] lep_1 second final leptons of the decay
     */
    R_6(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @brief The binned observable ratio @f$\frac {<BR_6>_{M \to V l_1^+l_1^-}}{<BR_6>_{M \to V l_2^+l_2^-}}@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<R_6>_{M \to V l^+ l^-}@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep1; /**< First final leptons type. */
    StandardModel::lepton lep2; /**< Second final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};


/**
 * @class ACP_MVll
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{CP}>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<A_{CP}>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP asymmetry helicity coefficients 
 * @f$<\Delta_i>@f$, computed in the MVll class, and the @f$\Gamma'@f$:
 * @f[
 * <A_{CP}>= -\frac {3<\Delta_{1c}> - <\Delta_{2c}> + 2(3<\Delta_{1s}> - <\Delta_{2s}>)}{4<\Gamma'>} \,.
 * @f]
 */
class ACP_MVll : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    ACP_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{CP}>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<A_{CP}>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class P3CP
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<P_3^{CP}>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<P_3^{CP}>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP asymmetry helicity coefficients 
 * @f$<\Delta_i>@f$ and the CP average helicity coefficients @f$<\Sigma_i>@f$,
 *  computed in the MVll class:
 * @f[
 * <P_3^{CP}>= -\frac {<\Delta_{9}>}{4<\Sigma_{2s}>} \,.
 * @f]
 */
class P3CP : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P3CP(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<P_3^{CP}>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<P_3^{CP}>@f$
    */
    double computeThValue ();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class M_1Prime
 * @ingroup Flavour
 * @brief A class for the observable @f$M_1'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$M_1'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the helicity amplitudes 
 * @f$H_{A,V}^{+,-}@f$ computed in the MVll class:
 * @f[
 * M_1'= \frac {|H_V^+|^2 + |H_V^-|^2 - |H_A^+|^2 - |H_A^-|^2}{2(|H_V^+|^2 + |H_V^-|^2 + |H_A^+|^2 + |H_A^-|^2)} \,.
 * @f]
 */
class M_1Prime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    M_1Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$M_1'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$M_1'@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class M_2Prime
 * @ingroup Flavour
 * @brief A class for the observable @f$M_2'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$M_2'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the helicity amplitudes 
 * @f$H_{A,V}^0,\,H_{P,S}@f$ computed in the MVll class:
 * @f[
 * M_2'= \frac {q^2/(2m_l^2)^2(|H_P|^2 + \beta|H_S|^2) + (|H_V^0|^2 - |H_A^0|^2)}{|H_V^0|^2 + |H_A^0|^2} \,.
 * @f]
 */
class M_2Prime : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    M_2Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$M_2'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$M_2'@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_3
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<S_3>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<S_3>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <S_3>= \frac {<\Sigma_{3}>}{<\Gamma'>} \,.
 * @f]
 */
class S_3 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<S_3>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<S_3>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_4
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<S_4>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<S_4>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <S_4>= \frac {<\Sigma_{4}>}{<\Gamma'>} \,.
 * @f]
 */
class S_4: public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_4(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<S_4>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<S_4>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_5
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<S_5>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<S_5>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <S_5>= \frac {<\Sigma_{5}>}{<\Gamma'>} \,.
 * @f]
 */
class S_5 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_5(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<S_5>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<S_5>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_7
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<S_7>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<S_7>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <S_7>= \frac {<\Sigma_{7}>}{<\Gamma'>} \,.
 * @f]
 */
class S_7 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_7(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<S_7>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<S_7>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_8
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<S_8>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<S_8>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <S_8>= \frac {<\Sigma_{8}>}{<\Gamma'>} \,.
 * @f]
 */
class S_8 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_8(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<S_8>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<S_8>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_9
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<S_9>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<S_9>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MVll class:
 * @f[
 * <S_9>= \frac {<\Sigma_{9}>}{<\Gamma'>} \,.
 * @f]
 */
class S_9 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_9(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<S_9>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<S_9>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_6
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_6>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<A_6>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Delta_i>@f$, computed in the MVll class:
 * @f[
 * <A_6>= \frac {<\Delta_{6}>}{<\Gamma'>} \,.
 * @f]
 */
class A_6 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_6(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<A_6>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<A_6>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_9
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_9>@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<A_9>@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the binned CP average helicity coefficients 
 * @f$<\Delta_i>@f$, computed in the MVll class:
 * @f[
 * <A_9>= \frac {<\Delta_{9}>}{<\Gamma'>} \,.
 * @f]
 */
class A_9 : public GammaPrime{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_9(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$<A_9>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<A_9>@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class V0
 * @ingroup Flavour
 * @brief A class for the form factor @f$<V_0>@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<V_0>@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class V0 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    V0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The averaged form factor @f$<V_0>@f$ in @f$B \to K^*@f$.
    * @return @f$<V_0>@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class Vp
 * @ingroup Flavour
 * @brief A class for the form factor @f$<V_+>@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<V_+>@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class Vp : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    Vp(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The averaged form factor @f$<V_+>@f$ in @f$B \to K^*@f$.
    * @return @f$<V_+>@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class Vm
 * @ingroup Flavour
 * @brief A class for the form factor @f$<V_->@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<V_->@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class Vm : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    Vm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The averaged form factor @f$<V_->@f$ in @f$B \to K^*@f$.
    * @return @f$<V_->@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class T0
 * @ingroup Flavour
 * @brief A class for the form factor @f$<T_0>@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<T_0>@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class T0 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    T0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The averaged form factor @f$<T_0>@f$ in @f$B \to K^*@f$.
    * @return @f$<T_0>@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class Tp
 * @ingroup Flavour
 * @brief A class for the form factor @f$<T_+>@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<T_+>@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class Tp : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    Tp(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The averaged form factor @f$<T_+>@f$ in @f$B \to K^*@f$.
    * @return @f$<T_+>@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class Tm
 * @ingroup Flavour
 * @brief A class for the form factor @f$<T_->@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<T_->@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class Tm : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    Tm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The averaged form factor @f$<T_->@f$ in @f$B \to K^*@f$.
    * @return @f$<T_->@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S
 * @ingroup Flavour
 * @brief A class for the form factor @f$<S>@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the averaged form factor @f$<S>@f$ in 
 * @f$B \to K^*@f$ using LCSR results at low @f$q^2@f$ and lattice results at high @f$q^2@f$
 */
class S : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The averaged form factor @f$<S>@f$ in @f$B \to K^*@f$.
    * @return @f$<S>@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};



/**
 * @class gtilde_1
 * @ingroup Flavour
 * @brief A class for the correction @f$\tilde{g}_1@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the correction @f$\tilde{g}_1@f$ in 
 * @f$B \to K^*@f$ 
 */
class gtilde_1 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    gtilde_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @brief The averaged correction @f$\tilde{g}_1@f$ in @f$B \to K^*@f$.
    * @return @f$\tilde{g}_1@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    unsigned int typ;/**< Observable type: 1 for averaged real part, 2 for averaged 
                      * immaginary part, 3 for averaged absolute part, 4 for averaged
                      * argument part, 5 for real part, 6 for immaginary part,
                      * 7 for absolute part, 8 for argument part. */
};



/**
 * @class gtilde_2
 * @ingroup Flavour
 * @brief A class for the correction @f$\tilde{g}_2@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the correction @f$\tilde{g}_2@f$ in 
 * @f$B \to K^*@f$ 
 */
class gtilde_2 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    gtilde_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @brief The averaged correction @f$\tilde{g}_1@f$ in @f$B \to K^*@f$.
    * @return @f$\tilde{g}_1@f$
    */
    double computeThValue();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    unsigned int typ;/**< Observable type: 1 for averaged real part, 2 for averaged 
                      * immaginary part, 3 for averaged absolute part, 4 for averaged
                      * argument part, 5 for real part, 6 for immaginary part,
                      * 7 for absolute part, 8 for argument part. */
};
    
    
    
 /**
 * @class gtilde_3
 * @ingroup Flavour
 * @brief A class for the correction @f$\tilde{g}_3@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the correction @f$\tilde{g}_3@f$ in 
 * @f$B \to K^*@f$
 */
class gtilde_3 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    gtilde_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @brief The averaged correction @f$\tilde{g}_1@f$ in @f$B \to K^*@f$.
    * @return @f$\tilde{g}_1@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    unsigned int typ;/**< Observable type: 1 for averaged real part, 2 for averaged 
                      * immaginary part, 3 for averaged absolute part, 4 for averaged
                      * argument part, 5 for real part, 6 for immaginary part,
                      * 7 for absolute part, 8 for argument part. */
};

/**
 * @class h_0
 * @ingroup Flavour
 * @brief A class for the correction @f$h_0@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the correction @f$h_0@f$ in 
 * @f$B \to K^*@f$
 */
class h_0 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    h_0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @brief The correction @f$h_0@f$ in @f$B \to K^*@f$.
    * @return @f$h_0@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    unsigned int typ;/**< Observable type: 1 for averaged real part, 2 for averaged 
                      * immaginary part, 3 for averaged absolute part, 4 for averaged
                      * argument part. */
};



/**
 * @class h_p
 * @ingroup Flavour
 * @brief A class for the correction @f$h_+@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the correction @f$h_+@f$ in 
 * @f$B \to K^*@f$
 */
class h_p : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    h_p(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @brief The correction @f$h_+@f$ in @f$B \to K^*@f$.
    * @return @f$h_+@f$
    */
    double computeThValue();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    unsigned int typ;/**< Observable type: 1 for averaged real part, 2 for averaged 
                      * immaginary part, 3 for averaged absolute part, 4 for averaged
                      * argument part. */
};
    
    
    
/**
 * @class h_m
 * @ingroup Flavour
 * @brief A class for the correction @f$h_-@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the correction @f$h_-@f$ in 
 * @f$B \to K^*@f$
 */
class h_m : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    h_m(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i);

    /**
    * @brief The correction @f$h_-@f$ in @f$B \to K^*@f$.
    * @return @f$h_-@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    unsigned int typ;/**< Observable type: 1 for averaged real part, 2 for averaged 
                      * immaginary part, 3 for averaged absolute part, 4 for averaged
                      * argument part. */
};



/**
 * @class hp0_hm0
 * @ingroup Flavour
 * @brief A class for the absolute value of the ratio @f$h_+^{(0)}/h_-^{(0)}@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute value of the ratio @f$h_+^{(0)}/h_-^{(0)}@f$ in 
 * @f$B \to K^*@f$
 */
class hp0_hm0 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    hp0_hm0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The absolute value of the ratio @f$h_+^{(0)}/h_-^{(0)}@f$ in @f$B \to K^*@f$.
    * @return @f$h_+^{(0)}/h_-^{(0)}@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class hm0_h00
 * @ingroup Flavour
 * @brief A class for the absolute value of the ratio @f$h_-^{(0)}/h_0^{(0)}@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute value of the ratio @f$h_-^{(0)}/h_0^{(0)}@f$ in 
 * @f$B \to K^*@f$
 */
class hm0_h00 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] typ_i observable type
     */
    hm0_h00(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The absolute value of the ratio @f$h_-^{(0)}/h_0^{(0)}@f$ in @f$B \to K^*@f$.
    * @return @f$h_-^{(0)}/h_0^{(0)}@f$
    */
    double computeThValue();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};

/***********************************************************************************************************************************
FUNCTIONAL
***********************************************************************************************************************************/
        
/**
 * @class P_1f
 * @ingroup Flavour
 * @brief A class for the observable @f$P_1@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_1@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_1=\frac {\Sigma_3}{2\Sigma_{2s}}\,.
 * @f]
 */
class P_1f : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_1f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    
    /**
    * @brief The observable @f$P_1@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_1@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class P_2f
 * @ingroup Flavour
 * @brief A class for the binned observable @f$P_2@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_2@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_2=\frac {\Sigma_6}{8\Sigma_{2s}}\,.
 * @f]
 */
class P_2f : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_2f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$P_2@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_2@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class P_3f
 * @ingroup Flavour
 * @brief A class for the observable @f$P_3@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_3@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_3=- \frac {\Sigma_9}{4\Sigma_{2s}}\,.
 * @f]
 */
class P_3f : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_3f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$P_3@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_3@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class P_4Primef
 * @ingroup Flavour
 * @brief A class for the observable @f$P_4'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_4'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_4'=\frac {\Sigma_4}{\sqrt{-\Sigma_{2s}\Sigma_{2c}}}\,.
 * @f]
 */
class P_4Primef : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_4Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$P_4'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_4'@f$
    */
    double computeThValue ();
    
protected:
    
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};


/**
 * @class P_5Primef
 * @ingroup Flavour
 * @brief A class for the observable @f$P_5'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_5'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_5'=\frac {\Sigma_5}{2\sqrt{-\Sigma_{2s}\Sigma_{2c}}}\,.
 * @f]
 */
class P_5Primef : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_5Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$P_5'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_5'@f$
    */
    double computeThValue ();
    
protected:
    
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class P_6Primef
 * @ingroup Flavour
 * @brief A class for the observable @f$P_6'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_6'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_6'=- \frac {\Sigma_7}{2 \sqrt{-\Sigma_{2s}\Sigma_{2c}}}\,.
 * @f]
 */
class P_6Primef : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_6Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The observable @f$P_6'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_6'@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class P_8Primef
 * @ingroup Flavour
 * @brief A class for the observable @f$P_8'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$P_8'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * P_8'=- \frac {\Sigma_8}{\sqrt{-\Sigma_{2s}\Sigma_{2c}}}\,.
 * @f]
 */
class P_8Primef : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    P_8Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);

    /**
    * @brief The observable @f$P_8'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$P_8'@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class GammaPrimef
 * @ingroup Flavour
 * @brief A class for the observable @f$\Gamma'@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$\Gamma'@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * \Gamma'=- \frac {1}{4} ((3\Sigma_{1c} - \Sigma_{2c}) + 2(3\Sigma_{1s} - \Sigma_{2s}))\,.
 * @f]
 */
class GammaPrimef : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    GammaPrimef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
     * @brief A method to compute the observable @f$\Gamma'@f$ in @f$M \to V l^+l^-@f$ at a generic @f$q^2@f$.
     * @param[in] qmin @f$q^2@f$ value
     * @param[in] lep final leptons of the decay
     * @return @f$\Gamma'(q^2)@f$
     */
    double computeGammaPrimef(double qmin, StandardModel::lepton lep);
    
    /**
    * @brief The observable @f$\Gamma'@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$\Gamma'@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class A_FBf
 * @ingroup Flavour
 * @brief A class for the observable @f$A_{FB}@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$A_{FB}@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * A_{FB}=- \frac {3\Sigma_{6s}}{4\Gamma'} \,.
 * @f]
 */
class A_FBf : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    A_FBf(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$A_{FB}@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$A_{FB}@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class BR_MVllf
 * @ingroup Flavour
 * @brief A class for the observable @f$BR@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$BR@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class, and the meson width @f$W_M@f$:
 * @f[
 * BR= \frac {\Gamma'}{W_M} \,.
 * @f]
 */
class BRf_MVll : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    BRf_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$BR@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$BR@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
    
};


/**
 * @class F_Lf
 * @ingroup Flavour
 * @brief A class for the observable @f$F_L@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$F_L@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP averaged helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class, and the meson width @f$M_W@f$:
 * @f[
 * BR= \frac {3\Sigma_{1c}-\Sigma_{2c}}{4\Gamma'} \,.
 * @f]
 */
class F_Lf : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    F_Lf(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
     * @brief A method to compute the observable @f$F_L@f$ in @f$M \to V l^+l^-@f$ at a generic @f$ q^2@f$.
     * @param[in] qmin generic value for @f$ q^2@f$
     * @param[in] lep final leptons of the decay
     * @return @f$F_L(q^2)@f$
     */
    double computeFLf(double qmin, StandardModel::lepton lep);

    /**
    * @brief The observable @f$F_L@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$F_L@f$
    */
    double computeThValue ();

private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};      

/**
 * @class S_3f
 * @ingroup Flavour
 * @brief A class for the observable @f$S_3@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$S_3@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP average helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * S_3= \frac {\Sigma_{3}}{\Gamma'} \,.
 * @f]
 */
class S_3f : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_3f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$S_3@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$S_3@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_4f
 * @ingroup Flavour
 * @brief A class for the observable @f$S_4@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$S_4@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP average helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * S_4= \frac {\Sigma_{4}}{\Gamma'} \,.
 * @f]
 */
class S_4f: public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_4f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$S_4@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$S_4@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_5f
 * @ingroup Flavour
 * @brief A class for the observable @f$S_5@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$S_5@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP average helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * S_5= \frac {\Sigma_{5}}{\Gamma'} \,.
 * @f]
 */
class S_5f : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_5f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$S_5@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$S_5@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_7f
 * @ingroup Flavour
 * @brief A class for the observable @f$S_7@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$S_7@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP average helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * S_7= \frac {\Sigma_{7}}{\Gamma'} \,.
 * @f]
 */
class S_7f : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_7f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$S_7@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$S_7@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_8f
 * @ingroup Flavour
 * @brief A class for the observable @f$S_8@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$S_8@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP average helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * S_8= \frac {\Sigma_{8}}{\Gamma'} \,.
 * @f]
 */
class S_8f : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_8f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$S_8@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$S_8@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};


/**
 * @class S_9f
 * @ingroup Flavour
 * @brief A class for the observable @f$S_9@f$ in @f$M \to V l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the observable @f$S_9@f$ in 
 * @f$M \to V l^+l^-@f$ in terms of the CP average helicity coefficients 
 * @f$\Sigma_i@f$, computed in the MVll class:
 * @f[
 * S_9= \frac {\Sigma_{9}}{\Gamma'} \,.
 * @f]
 */
class S_9f : public GammaPrimef{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    S_9f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The observable @f$S_9@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$S_9@f$
    */
    double computeThValue ();
   
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */

};
#endif	/* MVLLOBSERVABLES_H */

    