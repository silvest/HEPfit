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
    virtual  std::vector<WilsonCoefficient>& CMdiujleptonknu(int i, int j, int k) ;
    


protected:

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

    //dimension 6 four-fermion operators involving a left-handed current and a right-handed current

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
    gslpp::matrix<gslpp::complex> VuL, VuR, VdL, VdR, VeL, VeR;
    WilsonCoefficient mcd2, mcd1, mcbd, mcbs, mck2, mculeptonnu;
    TRandom3 myrnd;

};

#endif /* NPSMEFTD6GENERALMATCHING_H */


