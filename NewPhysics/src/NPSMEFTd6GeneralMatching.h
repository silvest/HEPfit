/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFTD6GENERALMATCHING_H
#define NPSMEFTD6GENERALMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"
#include <TRandom3.h>

class NPSMEFTd6General;

/**
 * @class NPSMEFTd6GeneralMatching
 * @ingroup NewPhysics
 * @brief A class for the matching in the NPSMEFTd6_General model at the scale @f$ \mu_W @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  This class, after update, contains the SMEFT coefficients at the scale @f$ \mu_W @f$ defined in the SMEFT model
 */
class NPSMEFTd6GeneralMatching : public StandardModelMatching {
public:
    NPSMEFTd6GeneralMatching(const NPSMEFTd6General & NPSMEFTd6General_i);

    virtual ~NPSMEFTd6GeneralMatching();

    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */

    void updateLEFTGeneralParameters();

    /**
     * 
     * @brief \f$ \Delta S = 2 \f$
     * @return the vector of 8 wilson coefficients: SM + SMEFT
     */
    virtual std::vector<WilsonCoefficient>& CMdk2();

    /**
     * 
     * @brief \f$ \Delta C = 2 \f$
     * @return the vector of 8 wilson coefficients: SM + SMEFT
     */
    virtual std::vector<WilsonCoefficient>& CMdd2();

    /**
     * 
     * @brief \f$ \Delta S = 2 \f$
     * @return the vector of 8 wilson coefficients: SM + SMEFT
     */
    virtual std::vector<WilsonCoefficient>& CMdbd2();

    /**
     * 
     * @brief \f$ \Delta S = 2 \f$
     * @return the vector of 8 wilson coefficients: SM + SMEFT
     */
    virtual std::vector<WilsonCoefficient>& CMdbs2();

    /**
     * 
     * @return Wilson coefficients for \f$ \bar{d}_i u_j \bar{\nu} \ell_k \f$ operators in the JMS basis ordered as CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij
     */
    virtual std::vector<WilsonCoefficient>& CMdiujleptonknu(int i, int j, int k) ;

    /**
     * 
     * @return Wilson coefficients for \f$ K_{L} \rightarrow \pi \nu \nu \f$
     */
    virtual std::vector<WilsonCoefficient>& CMkpnn();
    

    //dimension 6 four-fermion operators involving all left-handed fields

    /**
     * @brief Return CnunuVLL
     * @return \f$ C_{\nu \nu}^{V,LL} \f$ 
     */    
    const gslpp::complex getCnunuVLL(int i, int j, int k, int l) const;

    /**
     * @brief Return CeeVLL
     * @return \f$ C_{ee}^{V,LL} \f$ 
     */    
    const gslpp::complex getCeeVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnueVLL
     * @return \f$ C_{\nu e}^{V,LL} \f$ 
     */    
    const gslpp::complex getCnueVLL(int i, int j, int k, int l) const;

    /**
     * @brief Return CnuuVLL
     * @return \f$ C_{\nu u}^{V,LL} \f$ 
     */    
    const gslpp::complex getCnuuVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnudVLL
     * @return \f$ C_{\nu d}^{V,LL} \f$ 
     */    
    const gslpp::complex getCnudVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CeuVLL
     * @return \f$ C_{eu}^{V,LL} \f$ 
     */    
    const gslpp::complex getCeuVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CedVLL
     * @return \f$ C_{ed}^{V,LL} \f$ 
     */    
    const gslpp::complex getCedVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnueduVLL
     * @return \f$ C_{\nu e d u}^{V,LL} \f$ 
     */    
    const gslpp::complex getCnueduVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CuuVLL
     * @return \f$ C_{uu}^{V,LL} \f$ 
     */    
    const gslpp::complex getCuuVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CddVLL
     * @return \f$ C_{dd}^{V,LL} \f$ 
     */    
    const gslpp::complex getCddVLL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudV1LL
     * @return \f$ C_{ud}^{V1,LL} \f$ 
     */    
    const gslpp::complex getCudV1LL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudV8LL
     * @return \f$ C_{ud}^{V8,LL} \f$ 
     */    
    const gslpp::complex getCudV8LL(int i, int j, int k, int l) const;
    
    //dimension 6 four-fermion operators involving all right-handed fields
    
    /**
     * @brief Return CeeVRR
     * @return \f$ C_{ee}^{V,RR} \f$ 
     */    
    const gslpp::complex getCeeVRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CeuVRR
     * @return \f$ C_{eu}^{V,RR} \f$ 
     */    
    const gslpp::complex getCeuVRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CedVRR
     * @return \f$ C_{ed}^{V,RR} \f$ 
     */    
    const gslpp::complex getCedVRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CuuVRR
     * @return \f$ C_{uu}^{V,RR} \f$ 
     */    
    const gslpp::complex getCuuVRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CddVRR
     * @return \f$ C_{dd}^{V,RR} \f$ 
     */    
    const gslpp::complex getCddVRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudV1RR
     * @return \f$ C_{ud}^{V1,RR} \f$ 
     */    
    const gslpp::complex getCudV1RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudV8RR
     * @return \f$ C_{ud}^{V8,RR} \f$ 
     */    
    const gslpp::complex getCudV8RR(int i, int j, int k, int l) const;
    
    //dimension 6 four-fermion operators involving a left-handed vector current and a right-handed vector current

    /**
     * @brief Return CnueVLR
     * @return \f$ C_{\nu e}^{V,LR} \f$ 
     */    
    const gslpp::complex getCnueVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CeeVLR
     * @return \f$ C_{e e}^{V,LR} \f$ 
     */    
    const gslpp::complex getCeeVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnuuVLR
     * @return \f$ C_{\nu u}^{V,LR} \f$ 
     */    
    const gslpp::complex getCnuuVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnudVLR
     * @return \f$ C_{\nu d}^{V,LR} \f$ 
     */    
    const gslpp::complex getCnudVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CeuVLR
     * @return \f$ C_{eu}^{V,LR} \f$ 
     */    
    const gslpp::complex getCeuVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CedVLR
     * @return \f$ C_{ed}^{V,LR} \f$ 
     */    
    const gslpp::complex getCedVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CueVLR
     * @return \f$ C_{ue}^{V,LR} \f$ 
     */    
    const gslpp::complex getCueVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CdeVLR
     * @return \f$ C_{de}^{V,LR} \f$ 
     */    
    const gslpp::complex getCdeVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnueduVLR
     * @return \f$ C_{\nu e d u}^{V,LR} \f$ 
     */    
    const gslpp::complex getCnueduVLR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CuuV1LR
     * @return \f$ C_{uu}^{V1,LR} \f$ 
     */    
    const gslpp::complex getCuuV1LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CuuV8LR
     * @return \f$ C_{uu}^{V8,LR} \f$ 
     */    
    const gslpp::complex getCuuV8LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudV1LR
     * @return \f$ C_{ud}^{V1,LR} \f$ 
     */    
    const gslpp::complex getCudV1LR(int i, int j, int k, int l) const; 
    
    /**
     * @brief Return CudV8LR
     * @return \f$ C_{ud}^{V8,LR} \f$ 
     */    
    const gslpp::complex getCudV8LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CduV1LR
     * @return \f$ C_{du}^{V1,LR} \f$ 
     */    
    const gslpp::complex getCduV1LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CduV8LR
     * @return \f$ C_{du}^{V8,LR} \f$ 
     */    
    const gslpp::complex getCduV8LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CddV1LR
     * @return \f$ C_{dd}^{V1,LR} \f$ 
     */    
    const gslpp::complex getCddV1LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CddV8LR
     * @return \f$ C_{dd}^{V8,LR} \f$ 
     */    
    const gslpp::complex getCddV8LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudduV1LR
     * @return \f$ C_{\uddu}^{V1,LR} \f$ 
     */    
    const gslpp::complex getCudduV1LR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudduV8LR
     * @return \f$ C_{uddu}^{V8,LR} \f$ 
     */    
    const gslpp::complex getCudduV8LR(int i, int j, int k, int l) const;
    
    
    //dimension 6 four-fermion operators involving two right-handed scalar densities or tensor currents
    
    /**
     * @brief Return CeeSRR
     * @return \f$ C_{ee}^{S,RR} \f$ 
     */    
    const gslpp::complex getCeeSRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CeuSRR
     * @return \f$ C_{eu}^{S,RR} \f$ 
     */    
    const gslpp::complex getCeuSRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CeuTRR
     * @return \f$ C_{eu}^{T,RR} \f$ 
     */    
    const gslpp::complex getCeuTRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CedSRR
     * @return \f$ C_{ed}^{S,RR} \f$ 
     */    
    const gslpp::complex getCedSRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CedTRR
     * @return \f$ C_{ed}^{T,RR} \f$ 
     */    
    const gslpp::complex getCedTRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnueduSRR
     * @return \f$ C_{\nu e d u}^{S,RR} \f$ 
     */    
    const gslpp::complex getCnueduSRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnueduTRR
     * @return \f$ C_{\nu e d u}^{T,RR} \f$ 
     */    
    const gslpp::complex getCnueduTRR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CuuS1RR
     * @return \f$ C_{uu}^{S1,RR} \f$ 
     */    
    const gslpp::complex getCuuS1RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CuuS8RR
     * @return \f$ C_{uu}^{S8,RR} \f$ 
     */    
    const gslpp::complex getCuuS8RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudS1RR
     * @return \f$ C_{ud}^{S1,RR} \f$ 
     */    
    const gslpp::complex getCudS1RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudS8RR
     * @return \f$ C_{ud}^{S8,RR} \f$ 
     */    
    const gslpp::complex getCudS8RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CddS1RR
     * @return \f$ C_{dd}^{S1,RR} \f$ 
     */    
    const gslpp::complex getCddS1RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CddS8RR
     * @return \f$ C_{dd}^{S8,RR} \f$ 
     */    
    const gslpp::complex getCddS8RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudduS1RR
     * @return \f$ C_{uddu}^{S1,RR} \f$ 
     */    
    const gslpp::complex getCudduS1RR(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CudduS8RR
     * @return \f$ C_{uddu}^{S8,RR} \f$ 
     */    
    const gslpp::complex getCudduS8RR(int i, int j, int k, int l) const;

    
    //dimension 6 four-fermion operators involving a right-handed and a left-handed scalar density, plus hermitian conjugates

    /**
     * @brief Return CeuSRL
     * @return \f$ C_{e u}^{S,RL} \f$ 
     */    
    const gslpp::complex getCeuSRL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CedSRL
     * @return \f$ C_{e d}^{S,RL} \f$ 
     */    
    const gslpp::complex getCedSRL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CnueduSRL
     * @return \f$ C_{\nu e d u}^{S,RL} \f$ 
     */    
    const gslpp::complex getCnueduSRL(int i, int j, int k, int l) const;
    
    /**
     * @brief Return CdGLR (chromomagnetic dipole operator)
     * @return \f$ C_{d G}^{LR}(i,j) \f$ 
     */
    const gslpp::complex getCdG(int i, int j) const;
 
    /**
     * @brief Return CdgLR (electric dipole operator)
     * @return \f$ C_{d \gamma}^{LR}(i,j) \f$ 
     */
    const gslpp::complex getCdg(int i, int j) const;

    //Fermion rotation matrices to mass-eigenstate basis
    
    /**
     * @brief Return VuL
     * @return \f$ V^{u}_L \f$ 
     */    
    const gslpp::matrix<gslpp::complex> getVuL() const;
    
    /**
     * @brief Return VuR
     * @return \f$ V^{u}_R \f$ 
     */    
    const gslpp::matrix<gslpp::complex> getVuR() const;
    
    /**
     * @brief Return VdL
     * @return \f$ V^{d}_L \f$ 
     */    
    const gslpp::matrix<gslpp::complex> getVdL() const;
    
    /**
     * @brief Return VdR
     * @return \f$ V^{d}_R \f$ 
     */    
    const gslpp::matrix<gslpp::complex> getVdR() const;
    
    /**
     * @brief Return VeL
     * @return \f$ V^{e}_L \f$ 
     */    
    const gslpp::matrix<gslpp::complex> getVeL() const;
    
    /**
     * @brief Return VeR
     * @return \f$ V^{e}_R \f$ 
     */    
    const gslpp::matrix<gslpp::complex> getVeR() const;

protected:

    //utility zeroes
    const std::array<std::array<gslpp::complex, 3>, 3> zero33 {};
    const std::array<std::array<gslpp::complex, 2>, 2> zero22 {};
    const std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> zero3333 {};
    const std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> zero3322 {};
    const std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> zero2233 {};
    const std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 2> zero2332 {};
    const std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 3> zero3332 {};
    const std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> zero2222 {};

    //WET coefficients following Manohar
    //dimension 5 L-conserving dipole operators
    std::array<std::array<gslpp::complex, 3>, 3> Ceg = {}; ///< The real part of the dimension-5 operator coefficient \f$(C_{e\gamma})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<gslpp::complex, 2>, 2> Cug = {}; ///< The real part of the dimension-5 operator coefficient \f$(C_{u\gamma})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<gslpp::complex, 3>, 3> Cdg = {}; ///< The real part of the dimension-5 operator coefficient \f$(C_{d\gamma})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<gslpp::complex, 2>, 2> CuG = {}; ///< The real part of the dimension-5 operator coefficient \f$(C_{uG})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<gslpp::complex, 3>, 3> CdG = {}; ///< The real part of the dimension-5 operator coefficient \f$(C_{dG})_{ij}(\Lambda_{\rm{EW}})\f$.

    //dimension 6 involving the gluon field strength

    double CG = 0.; ///< The dimension-6 operator coefficient \f$C_{G}(\Lambda_{\rm{EW}})\f$.
    double CGtilde = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{G}}(\Lambda_{\rm{EW}})\f$.

    //dimension 6 four-fermion operators involving all left-handed fields

    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CnunuVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu\nu}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CeeVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{ee}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CnueVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu e}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CnuuVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu u}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CnudVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu d}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CeuVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{e u}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CedVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{e d}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 3> CnueduVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu e d u}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> CuuVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{uu}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CddVLL = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{V,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudV1LL = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{V1,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudV8LL = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{V8,LL})_{ijkl}(\Lambda_{\rm{EW}})\f$.

    //dimension 6 four-fermion operators involving all right-handed fields

    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CeeVRR = {}; ///< The dimension-6 operator coefficient \f$(C_{ee}^{V,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CeuVRR = {}; ///< The dimension-6 operator coefficient \f$(C_{e u}^{V,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CedVRR = {}; ///< The dimension-6 operator coefficient \f$(C_{e d}^{V,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> CuuVRR = {}; ///< The dimension-6 operator coefficient \f$(C_{uu}^{V,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CddVRR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{V,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudV1RR = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{V1,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudV8RR = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{V8,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.

    //dimension 6 four-fermion operators involving a left-handed vector current and a right-handed vector current

    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CnueVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu e}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CeeVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{ee}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CnuuVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu u}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CnudVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu d}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CeuVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{e u}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CedVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{e d}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CueVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{u e}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CdeVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{d e}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 3> CnueduVLR = {}; ///< The dimension-6 operator coefficient \f$(C_{\nu e d u}^{V,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> CuuV1LR = {}; ///< The dimension-6 operator coefficient \f$(C_{uu}^{V1,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> CuuV8LR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{V8,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudV1LR = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{V1,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudV8LR = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{V8,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CduV1LR = {}; ///< The dimension-6 operator coefficient \f$(C_{du}^{V1,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CduV8LR = {}; ///< The dimension-6 operator coefficient \f$(C_{du}^{V8,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CddV1LR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{V1,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CddV8LR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{V8,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 2> CudduV1LR = {}; ///< The dimension-6 operator coefficient \f$(C_{uddu}^{V1,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 2> CudduV8LR = {}; ///< The dimension-6 operator coefficient \f$(C_{uddu}^{V8,LR})_{ijkl}(\Lambda_{\rm{EW}})\f$.

    //dimension 6 four-fermion operators involving two right-handed scalar densities or tensor currents

    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CeeSRR = {}; ///< The dimension-6 operator coefficient \f$(C_{ee}^{S,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CeuSRR = {}; ///< The dimension-6 operator coefficient \f$(C_{eu}^{S,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CeuTRR = {}; ///< The dimension-6 operator coefficient \f$(C_{eu}^{T,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CedSRR = {}; ///< The dimension-6 operator coefficient \f$(C_{ed}^{S,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CedTRR = {}; ///< The dimension-6 operator coefficient \f$(C_{ed}^{T,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 3> CnueduSRR = {}; ///< The dimension-6 operator coefficient \f$(C_{nedu}^{S,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 3> CnueduTRR = {}; ///< The dimension-6 operator coefficient \f$(C_{nedu}^{T,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> CuuS1RR = {}; ///< The dimension-6 operator coefficient \f$(C_{uu}^{S1,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 2>, 2> CuuS8RR = {}; ///< The dimension-6 operator coefficient \f$(C_{uu}^{S8,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudS1RR = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{S1,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 2>, 2> CudS8RR = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{S8,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CddS1RR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{S1,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CddS8RR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd}^{S8,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 2> CudduS1RR = {}; ///< The dimension-6 operator coefficient \f$(C_{uddu}^{S1,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 2> CudduS8RR = {}; ///< The dimension-6 operator coefficient \f$(C_{uddu}^{S8,RR})_{ijkl}(\Lambda_{\rm{EW}})\f$.

    //dimension 6 four-fermion operators involving a right-handed and a left-handed scalar density, plus hermitian conjugates

    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 2>, 3>, 3> CeuSRL = {}; ///< The dimension-6 operator coefficient \f$(C_{eu}^{S,RL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 3>, 3>, 3>, 3> CedSRL = {}; ///< The dimension-6 operator coefficient \f$(C_{ed}^{S,RL})_{ijkl}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<gslpp::complex, 2>, 3>, 3>, 3> CnueduSRL = {}; ///< The dimension-6 operator coefficient \f$(C_{nedu}^{S,RL})_{ijkl}(\Lambda_{\rm{EW}})\f$.

private:
    const NPSMEFTd6General & mySMEFT;
    double LambdaNP2;
    double v2;
    double v;
    gslpp::matrix<gslpp::complex> VuL, VuLd, VuR, VuRd, VdL, VdLd, VdR, VdRd, VeL, VeLd, VeR, VeRd, MU, MD;
    WilsonCoefficient mcd2, mcd1, mcbd, mcbs, mck2, mculeptonnu, mckpnn;
    //TRandom3 myrnd;

};

#endif /* NPSMEFTD6GENERALMATCHING_H */


