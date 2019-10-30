/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GEORGIMACHACEK_H
#define	GEORGIMACHACEK_H

//#include "StandardModel.h"
#include "GMMatching.h"
#include "NPbase.h"

/**
 * @class GeorgiMachacek
 * @ingroup GeorgiMachacek
 * @brief A base class for the GeorgiMachacek model. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The Georgi-Machacek model extends the Standard Model by two scalar triplets.
 * Among the implemented theoretical constraints are positivity of the scalar
 * potential and the unitarity conditions.
 * As experimental constraints, Higgs signal strengths and direct searches for 
 * neutral, singly and doubly charged scalars
 * are available.
 *
 * 
 * @anchor GeorgiMachacekParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %GeorgiMachacek are summarized below.
 * 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * 
 * 
 * <tr>
 *   <td class="mod_name">%vDelta</td>
 *   <td class="mod_symb">@f$ v_{\Delta}  @f$</td>
 *   <td class="mod_desc"> The vacuum expectation value of the electroweak-triplet </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%alpha</td>
 *   <td class="mod_symb">@f$ \alpha  @f$</td>
 *   <td class="mod_desc"> Rotation angle which diagonalizes the singlet subspace  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mHh or %mHhsq</td>
 *   <td class="mod_symb">@f$m_{H_1}@f$</td>
 *   <td class="mod_desc">The mass of the "non-125 GeV" CP-even Higgs state.
 *   mHhsq (not mHh) used when  flag %use_sq_masses is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mA or %mAsq</td>
 *   <td class="mod_symb">@f$m_{H_3}@f$</td>
 *   <td class="mod_desc">The mass of the triplet Higgs states.
 *   mAsq (not mA) used when  flag %use_sq_masses is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mH5 or %mH5sq</td>
 *   <td class="mod_symb">@f$m_{H_5}@f$</td>
 *   <td class="mod_desc">The mass of the quintet Higgs states.</td>
 *   mH5sq (not mH5) used when  flag %use_sq_masses is set to true
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mu1</td>
 *   <td class="mod_symb">@f$\mu_{1}@f$</td>
 *   <td class="mod_desc">The @f$\mu_{1}@f$ potential parameter according to arxiv:1511.00865 @cite Chiang:2015amq </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mu2</td>
 *   <td class="mod_symb">@f$\mu_{2}@f$</td>
 *   <td class="mod_desc">The @f$\mu_{2}@f$ potential parameter according to arxiv:1511.00865 @cite Chiang:2015amq </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Q_GM</td>
 *   <td class="mod_symb">@f$Q_{\text{GM}}@f$</td>
 *   <td class="mod_desc">The GM scale.</td>
 * </tr>
 * </table>
 * 
 * 
 * @anchor GeorgiMachacekFlags
 * <h3>%Model flags</h3>
 *
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%use_sq_masses</td>
 *   <td class="mod_valu">true/false</td>
 *   <td class="mod_desc">This flag switched the parameters from mHh, mA, mH5
 *   to mHhsq, mAsq, mH5sq when set to true. </td>
 * </tr>
 * </table>
 *
 * 
 * 
 */

class GMcache; //forward reference to GMcache class



/**
 * @class GeorgiMachacek
 * @ingroup GeorgiMachacek
 * @brief A base class for the Georgi-Machacek model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class GeorgiMachacek: public NPbase {
public:

    static const int NGMvars = 8;
    //The parameters of the Higgs potential for Georgi Machacek model according to 1511.00865v1
    //We choose the physical basis: vDelta, alpha, mHh=mH1, mA=mH3, mH5, Mu1, Mu2
    static std::string GMvars[NGMvars];
    
    /**
     * @brief GeorgiMachacek constructor
     */
    GeorgiMachacek();
    
    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %GeorgiMachacek parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %GeorgiMachacek
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %GeorgiMachacek.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %GeorgiMachacek.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %GeorgiMachacek
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A get method to access the member reference of type GeorgiMachacekMatching.
     * @return a reference to a GeorgiMachacekMatching object
     */
    virtual GMMatching& getMatching() const
    {
        return GMM.getObj();
    }

    /**
     * @brief A method get the GeorgiMachacekCache
     * @return a object of the type %GeorgiMachacekCache
     */
    GMcache* getMyGMCache() const
    {
        return myGMcache;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Flags

//    virtual bool setFlagStr(const std::string name, const std::string value);
    virtual bool setFlag(const std::string, const bool);

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to get \f$\log(\tan \beta)$\f$
     * @return \f$\log(\tan \beta)$\f$
     */
    double getvDelta() const {
        return vDelta;
    }

    /**
     * @brief A method to get \f$\alpha$\f$
     * @return \f$\alpha$\f$
     */
    double getalpha() const {
        return alpha;
    }

    /**
     * @brief A method to get @f$\cos \alpha@f$
     * @return @f$\cos \alpha@f$
     */
    double getcosa() const{
        return cos(alpha);
    }

    /**
     * @brief A method to get @f$\sin \alpha@f$
     * @return @f$\sin \alpha@f$
     */
    double getsina() const{
        return sin(alpha);
    }

    /**
     * @brief A method to get the squared mass of the lighter singlet Higgs
     * @return squared mass of the lighter singlet Higgs
     */
    double getmHl2() const {
        if(flag_use_sq_masses) {
            if(mHhsq < mHl2) {
                return mHhsq;
            }
            else
            {
                return mHl2;
            }
        }
        else
        {
            if(mHh*mHh < mHl2) {
                return mHh*mHh;
            }
            else
            {
                return mHl2;
            }
        }
    }

    /**
     * @brief A method to get the squared mass of the heavier singlet Higgs
     * @return squared mass of the heavier singlet Higgs
     */
    double getmHh2() const {
        if(flag_use_sq_masses) {
            if(mHhsq < 0.) {
                return 0.;
            }
            else if(mHhsq < mHl2) {
                return mHl2;
            }
            else
            {
                return mHhsq;
            }
        }
        else
        {
            if(mHh*mHh < mHl2) {
                return mHl2;
            }
            else
            {
                return mHh*mHh;
            }
        }
    }

    /**
     * @brief A method to get the mass of the heavier singlet Higgs
     * @return mass of the heavier singlet Higgs
     */
    double getmHh() const {
        if(flag_use_sq_masses) {
            if(mHhsq < 0.) {
                return 0.;
            }
            else if(mHhsq < mHl2) {
                return sqrt(mHl2);
            }
            else
            {
                return sqrt(mHhsq);
            }
        }
        else
        {
            if(mHh*mHh < mHl2) {
                return sqrt(mHl2);
            }
            else
            {
                return mHh;
            }
        }
    }

    /**
     * @brief A method to get the squared mass of the singlet Higgs input
     * @return squared mass of the singlet Higgs input
     */
    double getinputmHh2() const {
        if(flag_use_sq_masses) {
            if(mHhsq < 0.) {
                return 0.;
            }
                return mHhsq;
        }
        else
        {
                return mHh*mHh;
        }
    }

    /**
     * @brief A method to get the squared triplet mass
     * @return squared triplet mass
     */
    double getmAsq() const {
        if(flag_use_sq_masses) {
            return mAsq;
        }
        else
        {
            return mA*mA;
        }
    }

    /**
     * @brief A method to get the triplet mass
     * @return triplet mass
     */
    double getmA() const {
        if(flag_use_sq_masses) {
            if(mA < 0.) {
                return 0.;
            }
            else
            {
                return sqrt(mAsq);
            }
        }
        else
        {
                return mA;
        }
    }

    /**
     * @brief A method to get the squared quintet mass
     * @return squared quintet mass
     */
    double getmH5sq() const {
        if(flag_use_sq_masses) {
            return mH5sq;
        }
        else
        {
            return mH5*mH5;
        }
    }

    /**
     * @brief A method to get the quintet mass
     * @return quintet mass
     */
    double getmH5() const {
        if(flag_use_sq_masses) {
            if(mH5 < 0.) {
                return 0.;
            }
            else
            {
                return sqrt(mH5sq);
            }
        }
        else
        {
                return mH5;
        }
    }

    /**
     * @brief A method to get the massive parameter of the scalar potential \f$\mu_1$\f
     * @return massive parameter of the scalar potential \f$\mu_1$\f
     */
    double getMu1() const {
        return Mu1;
    }

    /**
     * @brief A method to get the massive parameter of the scalar potential \f$\mu_2$\f
     * @return massive parameter of the scalar potential \f$\mu_2$\f
     */
    double getMu2() const {
        return Mu2;
    }

    /**
     * @brief A method to get the Georgi-Machacek scale
     * @return Georgi-Machacek scale
     */
    double getQ_GM() const {
        return Q_GM;
    }

    virtual double Mw() const;

    virtual double muggH(const double sqrt_s) const;
    virtual double muVBF(const double sqrt_s) const;
    virtual double mueeWBF(const double sqrt_s) const;
    virtual double muWH(const double sqrt_s) const;
    virtual double muZH(const double sqrt_s) const;
    virtual double mueeZH(const double sqrt_s) const;
    virtual double muVH(const double sqrt_s) const;
    virtual double muVBFpVH(const double sqrt_s) const;
    virtual double muttH(const double sqrt_s) const;
    virtual double GammaTotal() const;
    virtual double BrHggRatio() const;
    virtual double BrHWWRatio() const;
    virtual double BrHZZRatio() const;
    virtual double BrHZgaRatio() const;
    virtual double BrHgagaRatio() const;
    virtual double BrHmumuRatio() const;
    virtual double BrHtautauRatio() const;
    virtual double BrHccRatio() const;
    virtual double BrHbbRatio() const;
    virtual double muggHgaga(const double sqrt_s) const;
    virtual double muVBFHgaga(const double sqrt_s) const;
    virtual double muVHgaga(const double sqrt_s) const;
    virtual double muttHgaga(const double sqrt_s) const;
    virtual double muggHZZ(const double sqrt_s) const;
    virtual double muVBFHZZ(const double sqrt_s) const;
    virtual double muVHZZ(const double sqrt_s) const;
    virtual double muttHZZ(const double sqrt_s) const;
    virtual double muggHWW(const double sqrt_s) const;
    virtual double muVBFHWW(const double sqrt_s) const;
    virtual double muVHWW(const double sqrt_s) const;
    virtual double muttHWW(const double sqrt_s) const;
    virtual double muggHtautau(const double sqrt_s) const;
    virtual double muVBFHtautau(const double sqrt_s) const;
    virtual double muVHtautau(const double sqrt_s) const;
    virtual double muttHtautau(const double sqrt_s) const;
    virtual double muggHbb(const double sqrt_s) const;
    virtual double muVBFHbb(const double sqrt_s) const;
    virtual double muVHbb(const double sqrt_s) const;
    virtual double muttHbb(const double sqrt_s) const;
    virtual double muppHmumu(const double sqrt_s) const;
    virtual double muppHZga(const double sqrt_s) const;
    virtual double computeGammaTotalRatio() const;


protected: 
    /**
     * @brief A method to set the value of a parameter of %GeorgiMachacek.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);

    mutable Matching<GMMatching,GeorgiMachacek> GMM; ///< An object of type Matching.

private:

    GMcache* myGMcache;

    double vDelta, alpha, mHh, mA, mH5, mHhsq, mAsq, mH5sq, Mu1, Mu2, Q_GM;
    double mHl2;
    bool flag_use_sq_masses;
//    double sign(const double x) const;

};

/**
 * @}
 */

#endif	/* GEORGIMACHACEK_H */
