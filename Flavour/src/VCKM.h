/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VCKM_H
#define VCKM_H

#include "ThObservable.h"

/**
* @class VCKM
* @ingroup Flavour
* @brief A class for the CKM elements @f$V_{ij} @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM element @f$V_{ij} @f$.
*/
class VCKM : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] obsFlag_1 the first index of the CKM element
     * @param[in] obsFlag_2 the second index of the CKM element
     */
    VCKM(const StandardModel& SM_i, unsigned int obsFlag_1,  unsigned int obsFlag_2);
    
    /**
     * @brief Destructor.
     */
    virtual ~VCKM();

    /**
     * @return The CKM element @f$V_{ij} @f$
     */
    double computeThValue();
    
private:
    unsigned int obs_1;
    unsigned int obs_2;
        
};

/**
* @class CKM_Alpha
* @ingroup Flavour
* @brief A class for the CKM angle @f$\alpha @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\alpha @f$ in degrees.
*/
class CKM_Alpha : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Alpha(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\alpha @f$ in degrees
     */
    double computeThValue ();
      
};

/**
* @class CKM_Beta
* @ingroup Flavour
* @brief A class for the CKM angle @f$\beta @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\beta @f$ in degrees.
*/
class CKM_Beta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Beta(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\beta @f$ in degrees
     */
    double computeThValue ();
      
};

/**
* @class CKM_Betas
* @ingroup Flavour
* @brief A class for the CKM angle @f$\beta_s @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\beta_s @f$ in degrees.
*/
class CKM_Betas : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Betas(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\beta_s @f$ in degrees
     */
    double computeThValue ();
      
};

/**
* @class CKM_Gamma
* @ingroup Flavour
* @brief A class for the CKM angle @f$\gamma @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\gamma @f$ in degrees.
*/
class CKM_Gamma : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Gamma(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\gamma @f$ in degrees
     */
    double computeThValue ();
      
};

/**
* @class CKM_2BpG
* @ingroup Flavour
* @brief A class for the sum of CKM angles @f$2 \beta + \gamma @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the sum of CKM angles 
* @f$2 \beta + \gamma @f$ in degrees.
*/
class CKM_2BpG : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_2BpG(const StandardModel& SM_i);
    
    /**
     * @return The sum of CKM angles @f$2 \beta + \gamma @f$ in degrees
     */
    double computeThValue ();
      
};

/**
* @class CKM_S2Beta
* @ingroup Flavour
* @brief A class for @f$\sin 2 \beta @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute @f$\sin 2 \beta @f$, with @f$ \beta @f$
* being a CKM angle.
*/
class CKM_S2Beta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_S2Beta(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin 2 \beta @f$
     */
    double computeThValue ();
      
};

/**
* @class CKM_C2Beta
* @ingroup Flavour
* @brief A class for @f$\cos 2 \beta @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute @f$\cos 2 \beta @f$, with @f$ \beta @f$
* being a CKM angle.
*/
class CKM_C2Beta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_C2Beta(const StandardModel& SM_i);
    
    /**
     * @return @f$\cos 2 \beta @f$
     */
    double computeThValue ();
      
};

/**
* @class CKM_SinTheta12
* @ingroup Flavour
* @brief A class for @f$\sin \theta_{12} @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute @f$\sin \theta_{12} @f$, with @f$\theta_{12} @f$
* being one of the angles of the PDG parameterization of the CKM matrix.
*/
class CKM_SinTheta12 : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_SinTheta12(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin \theta_{12} @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

/**
* @class CKM_SinTheta13
* @ingroup Flavour
* @brief A class for @f$\sin \theta_{13} @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute @f$\sin \theta_{13} @f$, with @f$\theta_{13} @f$
* being one of the angles of the PDG parameterization of the CKM matrix.
*/
class CKM_SinTheta13 : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_SinTheta13(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin \theta_{13} @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

/**
* @class CKM_SinTheta23
* @ingroup Flavour
* @brief A class for @f$\sin \theta_{23} @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute @f$\sin \theta_{23} @f$, with @f$\theta_{23} @f$
* being one of the angles of the PDG parameterization of the CKM matrix.
*/
class CKM_SinTheta23 : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_SinTheta23(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin \theta_{13} @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

/**
* @class CKM_Delta
* @ingroup Flavour
* @brief A class for @f$\delta @f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute @f$\delta @f$, 
* the phase of the PDG parameterization of the CKM matrix.
*/
class CKM_Delta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Delta(const StandardModel& SM_i);
    
    /**
     * @return the phase @f$\delta @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

/**
* @class J_CP
* @ingroup Flavour
* @brief A class for the Jarlskog determinant. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the Jarlskog determinant.
*/
class J_CP : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    J_CP(const StandardModel& SM_i);
    
    /**
     * @return the Jarlskog determinant
     */
    double computeThValue ();
      
};

/**
* @class CKM_Rt
* @ingroup Flavour
* @brief A class for the CKM parameters ratio @f$R_t@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM parameters ratio
* @f$R_t=|(V_{td} V_{tb}^*)/(V_{cd}V_{cb}^*)|@f$.
*/
class CKM_Rt : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Rt(const StandardModel& SM_i);
    
    /**
     * @brief @f$R_t=|(V_{td} V_{tb}^*)/(V_{cd}V_{cb}^*)|@f$.
     * @return @f$R_t@f$
     */
    double computeThValue ();
      
};

/**
* @class CKM_Rt_dms
* @ingroup Flavour
* @brief An auxiliary class for the CKM UT plot from the constraint from @f$\Delta M_s@f$.
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to produce the UT plot with the constraints from
* @f$\Delta M_s@f$ and @f$\Delta M_d@f$.
*/
class CKM_Rt_dms : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Rt_dms(const StandardModel& SM_i);

    /**
     * @brief @f$R_t=|(V_{td} V_{tb}^*)/(V_{cd}V_{cb}^*)|@f$.
     * @return @f$R_t * \sqrt{1+\Lambda^2)*(1-2 \rho)}@f$
     */
    double computeThValue ();

};

/**
* @class CKM_Rts
* @ingroup Flavour
* @brief A class for the CKM parameters ratio @f$R_{ts}@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM parameters ratio
* @f$R_{ts}=|(V_{ts} V_{tb}^*)/(V_{cs}V_{cb}^*)|@f$.
*/
class CKM_Rts : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Rts(const StandardModel& SM_i);
    
    /**
     * @brief @f$R_{ts}=|(V_{ts} V_{tb}^*)/(V_{cs}V_{cb}^*)|@f$.
     * @return @f$R_{ts}@f$
     */
    double computeThValue ();
      
};

/**
* @class CKM_Rb
* @ingroup Flavour
* @brief A class for the CKM parameters ratio @f$R_b@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM parameters ratio
* @f$R_{b}=|(V_{ud} V_{ub}^*)/(V_{ud}V_{ub}^*)|@f$.
*/
class CKM_Rb : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_Rb(const StandardModel& SM_i);
    
    /**
     * @brief @f$R_{b}=|(V_{ud} V_{ub}^*)/(V_{ud}V_{ub}^*)|@f$.
     * @return @f$R_{b}@f$
     */
    double computeThValue ();
      
};

/**
* @class CKM_VtdoVts
* @ingroup Flavour
* @brief A class for the CKM parameters ratio @f$|V_{td}/V_{ts}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM parameters ratio
* @f$|V_{td}/V_{ts}|@f$.
*/
class CKM_VtdoVts : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_VtdoVts(const StandardModel& SM_i);
    
    /**
     * @brief @f$|V_{td}/V_{ts}|@f$.
     * @return @f$|V_{td}/V_{ts}|@f$
     */
    double computeThValue ();
      
};


/**
* @class Abslam_t
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_t|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_t| = \vert V_{td} V_{ts}^* \vert@f$.
*/
class Abslam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_t(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{td} V_{ts}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_c
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_c|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_c| = \vert V_{cd} V_{cs}^* \vert@f$.
*/
class Abslam_c : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_c(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{cd} V_{cs}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_u
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_u|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_u| = \vert V_{ud} V_{us}^* \vert@f$.
*/
class Abslam_u : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_u(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{ud} V_{us}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_td
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_{td}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_{td}| = \vert V_{td} V_{tb}^* \vert@f$.
*/
class Abslam_td : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_td(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{td} V_{tb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_cd
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_{cd}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_{cd}| = \vert V_{cd} V_{cb}^* \vert@f$.
*/
class Abslam_cd : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_cd(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{cd} V_{cb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_ud
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_{ud}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_{ud}| = \vert V_{ud} V_{ub}^* \vert@f$.
*/
class Abslam_ud : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_ud(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{ud} V_{ub}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_ts
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_{ts}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_{ts}| = \vert V_{ts} V_{tb}^* \vert@f$.
*/
class Abslam_ts : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_ts(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{ts} V_{tb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_cs
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_{cs}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_{cs}| = \vert V_{cs} V_{cb}^* \vert@f$.
*/
class Abslam_cs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_cs(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{cs} V_{cb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Abslam_us
* @ingroup Flavour
* @brief A class for the absolute value of the CKM parameters combination @f$|\lambda_{us}|@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$|\lambda_{us}| = \vert V_{us} V_{ub}^* \vert@f$.
*/
class Abslam_us : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Abslam_us(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{us} V_{ub}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_t
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_t)@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_t) = \mathrm{Im}(V_{td} V_{ts}^*)@f$.
*/
class Imlam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_t(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{td} V_{ts}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_c
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_c)@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_c) = \mathrm{Im}(V_{cd} V_{cs}^*)@f$.
*/
class Imlam_c : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_c(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{cd} V_{cs}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_u
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_u)@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_u) = \mathrm{Im}(V_{ud} V_{us}^*)@f$.
*/
class Imlam_u : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_u(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{ud} V_{us}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_td
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_{td})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_{td}) = \mathrm{Im}(V_{td} V_{tb}^*)@f$.
*/
class Imlam_td : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_td(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{td} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_cd
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_{cd})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_{cd}) = \mathrm{Im}(V_{cd} V_{cb}^*)@f$.
*/
class Imlam_cd : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_cd(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{cd} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_ud
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_{ud})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_{ud}) = \mathrm{Im}(V_{ud} V_{ub}^*)@f$.
*/
class Imlam_ud : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_ud(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{ud} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_ts
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_{ts})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_{ts}) = \mathrm{Im}(V_{ts} V_{tb}^*)@f$.
*/
class Imlam_ts : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_ts(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{ts} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_cs
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_{cs})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_{cs}) = \mathrm{Im}(V_{cs} V_{cb}^*)@f$.
*/
class Imlam_cs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_cs(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{cs} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Imlam_us
* @ingroup Flavour
* @brief A class for the imaginary part of the CKM parameters combination @f$\mathrm{Im}(\lambda_{us})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Im}(\lambda_{us}) = \mathrm{Im}(V_{us} V_{ub}^*)@f$.
*/
class Imlam_us : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Imlam_us(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{us} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_t
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{t})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{t}) = \mathrm{Im}(V_{td} V_{ts}^*)@f$.
*/
class Relam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_t(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{td} V_{ts}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_c
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{c})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{c}) = \mathrm{Im}(V_{cd} V_{cs}^*)@f$.
*/
class Relam_c : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_c(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{cd} V_{cs}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_u
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{u})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{u}) = \mathrm{Im}(V_{ud} V_{us}^*)@f$.
*/
class Relam_u : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_u(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{ud} V_{us}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_td
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{td})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{td}) = \mathrm{Im}(V_{td} V_{tb}^*)@f$.
*/
class Relam_td : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_td(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{td} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_cd
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{cd})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{cd}) = \mathrm{Im}(V_{cd} V_{cb}^*)@f$.
*/
class Relam_cd : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_cd(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{cd} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_ud
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{ud})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{ud}) = \mathrm{Im}(V_{ud} V_{ub}^*)@f$.
*/
class Relam_ud : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_ud(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{ud} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_ts
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{ts})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{ts}) = \mathrm{Im}(V_{ts} V_{tb}^*)@f$.
*/
class Relam_ts : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_ts(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{ts} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_cs
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{cs})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{cs}) = \mathrm{Im}(V_{cs} V_{cb}^*)@f$.
*/
class Relam_cs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_cs(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{cs} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class Relam_us
* @ingroup Flavour
* @brief A class for the real part of the CKM parameters combination @f$\mathrm{Re}(\lambda_{us})@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the absolute value of the CKM parameters 
* combination @f$\mathrm{Re}(\lambda_{us}) = \mathrm{Im}(V_{us} V_{ub}^*)@f$.
*/
class Relam_us : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relam_us(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{us} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class CKM_rho
* @ingroup Flavour
* @brief A class for the CKM parameter @f$\rho@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM parameter @f$\rho@f$. 
*/
class CKM_rho : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_rho(const StandardModel& SM_i);
    
    /**
     * @return @f$\rho @f$
     */
    double computeThValue ();
    
private:
    
};

/**
* @class CKM_eta
* @ingroup Flavour
* @brief A class for the CKM parameter @f$\eta@f$. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM parameter @f$\eta@f$. 
*/
class CKM_eta : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CKM_eta(const StandardModel& SM_i);
    
    /**
     * @return @f$\eta @f$
     */
    double computeThValue ();
    
private:
    
};


#endif	/* VCKM_H */

