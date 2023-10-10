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
    
    double getC2B() const
    {
        return C2B;
    }

    double getC2BS() const
    {
        return C2BS;
    }

    double getC2W() const
    {
        return C2W;
    }

    double getC2WS() const
    {
        return C2WS;
    }

    double getCDB() const
    {
        return CDB;
    }

    double getCDHB() const
    {
        return CDHB;
    }

    double getCDHW() const
    {
        return CDHW;
    }

    double getCDW() const
    {
        return CDW;
    }

    double getCG() const
    {
        return CG;
    }

    double getCGtilde() const
    {
        return CGtilde;
    }

    double getCH() const
    {
        return CH;
    }

    double getCHB() const
    {
        return CHB;
    }

    double getCHBtilde() const
    {
        return CHBtilde;
    }

    double getCHD() const
    {
        return CHD;
    }

    double getCHG() const
    {
        return CHG;
    }

    double getCHGtilde() const
    {
        return CHGtilde;
    }

    double getCHW() const
    {
        return CHW;
    }

    double getCHWB() const
    {
        return CHWB;
    }

    double getCHWtilde() const
    {
        return CHWtilde;
    }

    double getCHWtildeB() const
    {
        return CHWtildeB;
    }

    double getCHbox() const
    {
        return CHbox;
    }

    double getCT() const
    {
        return CT;
    }

    double getCW() const
    {
        return CW;
    }

    double getCWtilde() const
    {
        return CWtilde;
    }
    
    std::array<std::array<double, 3>, 3> getCHdI() const
    {
        return CHdI;
    }

    std::array<std::array<double, 3>, 3> getCHdR() const
    {
        return CHdR;
    }

    std::array<std::array<double, 3>, 3> getCHeI() const
    {
        return CHeI;
    }

    std::array<std::array<double, 3>, 3> getCHeR() const
    {
        return CHeR;
    }

    std::array<std::array<double, 3>, 3> getCHl1I() const
    {
        return CHl1I;
    }

    std::array<std::array<double, 3>, 3> getCHl1R() const
    {
        return CHl1R;
    }

    std::array<std::array<double, 3>, 3> getCHl3I() const
    {
        return CHl3I;
    }

    std::array<std::array<double, 3>, 3> getCHl3R() const
    {
        return CHl3R;
    }

    std::array<std::array<double, 3>, 3> getCHq1I() const
    {
        return CHq1I;
    }

    std::array<std::array<double, 3>, 3> getCHq1R() const
    {
        return CHq1R;
    }

    std::array<std::array<double, 3>, 3> getCHq3I() const
    {
        return CHq3I;
    }

    std::array<std::array<double, 3>, 3> getCHq3R() const
    {
        return CHq3R;
    }

    std::array<std::array<double, 3>, 3> getCHuI() const
    {
        return CHuI;
    }

    std::array<std::array<double, 3>, 3> getCHuR() const
    {
        return CHuR;
    }

    std::array<std::array<double, 3>, 3> getCHudI() const
    {
        return CHudI;
    }

    std::array<std::array<double, 3>, 3> getCHudR() const
    {
        return CHudR;
    }

    std::array<std::array<double, 3>, 3> getCdBI() const
    {
        return CdBI;
    }

    std::array<std::array<double, 3>, 3> getCdBR() const
    {
        return CdBR;
    }

    std::array<std::array<double, 3>, 3> getCdGI() const
    {
        return CdGI;
    }

    std::array<std::array<double, 3>, 3> getCdGR() const
    {
        return CdGR;
    }

    std::array<std::array<double, 3>, 3> getCdHI() const
    {
        return CdHI;
    }

    std::array<std::array<double, 3>, 3> getCdHR() const
    {
        return CdHR;
    }

    std::array<std::array<double, 3>, 3> getCdWI() const
    {
        return CdWI;
    }

    std::array<std::array<double, 3>, 3> getCdWR() const
    {
        return CdWR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCddI() const
    {
        return CddI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCddR() const
    {
        return CddR;
    }

    std::array<std::array<double, 3>, 3> getCeBI() const
    {
        return CeBI;
    }

    std::array<std::array<double, 3>, 3> getCeBR() const
    {
        return CeBR;
    }

    std::array<std::array<double, 3>, 3> getCeHI() const
    {
        return CeHI;
    }

    std::array<std::array<double, 3>, 3> getCeHR() const
    {
        return CeHR;
    }

    std::array<std::array<double, 3>, 3> getCeWI() const
    {
        return CeWI;
    }

    std::array<std::array<double, 3>, 3> getCeWR() const
    {
        return CeWR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCedI() const
    {
        return CedI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCedR() const
    {
        return CedR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCeeI() const
    {
        return CeeI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCeeR() const
    {
        return CeeR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCeuI() const
    {
        return CeuI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCeuR() const
    {
        return CeuR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCldI() const
    {
        return CldI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCldR() const
    {
        return CldR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCleI() const
    {
        return CleI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCleR() const
    {
        return CleR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCledqI() const
    {
        return CledqI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCledqR() const
    {
        return CledqR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClequ1I() const
    {
        return Clequ1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClequ1R() const
    {
        return Clequ1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClequ3I() const
    {
        return Clequ3I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClequ3R() const
    {
        return Clequ3R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCllI() const
    {
        return CllI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCllR() const
    {
        return CllR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClq1I() const
    {
        return Clq1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClq1R() const
    {
        return Clq1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClq3I() const
    {
        return Clq3I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getClq3R() const
    {
        return Clq3R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCluI() const
    {
        return CluI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCluR() const
    {
        return CluR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqd1I() const
    {
        return Cqd1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqd1R() const
    {
        return Cqd1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqd8I() const
    {
        return Cqd8I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqd8R() const
    {
        return Cqd8R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqeI() const
    {
        return CqeI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqeR() const
    {
        return CqeR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqq1I() const
    {
        return Cqq1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqq1R() const
    {
        return Cqq1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqq3I() const
    {
        return Cqq3I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqq3R() const
    {
        return Cqq3R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqu1I() const
    {
        return Cqu1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqu1R() const
    {
        return Cqu1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqu8I() const
    {
        return Cqu8I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCqu8R() const
    {
        return Cqu8R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCquqd1I() const
    {
        return Cquqd1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCquqd1R() const
    {
        return Cquqd1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCquqd8I() const
    {
        return Cquqd8I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCquqd8R() const
    {
        return Cquqd8R;
    }

    std::array<std::array<double, 3>, 3> getCuBI() const
    {
        return CuBI;
    }

    std::array<std::array<double, 3>, 3> getCuBR() const
    {
        return CuBR;
    }

    std::array<std::array<double, 3>, 3> getCuGI() const
    {
        return CuGI;
    }

    std::array<std::array<double, 3>, 3> getCuGR() const
    {
        return CuGR;
    }

    std::array<std::array<double, 3>, 3> getCuHI() const
    {
        return CuHI;
    }

    std::array<std::array<double, 3>, 3> getCuHR() const
    {
        return CuHR;
    }

    std::array<std::array<double, 3>, 3> getCuWI() const
    {
        return CuWI;
    }

    std::array<std::array<double, 3>, 3> getCuWR() const
    {
        return CuWR;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCud1I() const
    {
        return Cud1I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCud1R() const
    {
        return Cud1R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCud8I() const
    {
        return Cud8I;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCud8R() const
    {
        return Cud8R;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCuuI() const
    {
        return CuuI;
    }

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> getCuuR() const
    {
        return CuuR;
    }

protected:

    double CG = 0.; ///< The dimension-6 operator coefficient \f$C_{G}(\Lambda_{\rm{EW}})\f$.
    double CW = 0.; ///< The dimension-6 operator coefficient \f$C_{W}(\Lambda_{\rm{EW}})\f$.
    double CHG = 0.; ///< The dimension-6 operator coefficient \f$C_{HG}(\Lambda_{\rm{EW}})\f$.
    double CHW = 0.; ///< The dimension-6 operator coefficient \f$C_{HW}(\Lambda_{\rm{EW}})\f$.
    double CHB = 0.; ///< The dimension-6 operator coefficient \f$C_{HB}(\Lambda_{\rm{EW}})\f$.
    double CHWB = 0.; ///< The dimension-6 operator coefficient \f$C_{HWB}(\Lambda_{\rm{EW}})\f$.
    double CHD = 0.; ///< The dimension-6 operator coefficient \f$C_{HD}(\Lambda_{\rm{EW}})\f$.
    double CHbox = 0.; ///< The dimension-6 operator coefficient \f$C_{H\Box}(\Lambda_{\rm{EW}})\f$.
    double CH = 0.; ///< The dimension-6 operator coefficient \f$C_{H}(\Lambda_{\rm{EW}})\f$.
    double CGtilde = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{G}}(\Lambda_{\rm{EW}})\f$.
    double CWtilde = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{W}}(\Lambda_{\rm{EW}})\f$.
    double CHGtilde = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{G}}(\Lambda_{\rm{EW}})\f$.
    double CHWtilde = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{W}}(\Lambda_{\rm{EW}})\f$.
    double CHBtilde = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{B}}(\Lambda_{\rm{EW}})\f$.
    double CHWtildeB = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{W}B}(\Lambda_{\rm{EW}})\f$.


    ////////////////////////////////////////////////////////////////////////////////////////////////
    //These operators should be written in terms of those of the Warsaw basis
    double C2B = 0.; ///< The dimension-6 operator coefficient \f$C_{2W}(\Lambda_{\rm{EW}})\f$.
    double C2W = 0.; ///< The dimension-6 operator coefficient \f$C_{2B}(\Lambda_{\rm{EW}})\f$.
    double C2BS = 0.; ///< The dimension-6 operator coefficient \f$C_{2W}^{SILH}(\Lambda_{\rm{EW}})\f$.
    double C2WS = 0.; ///< The dimension-6 operator coefficient \f$C_{2B}^{SILH}(\Lambda_{\rm{EW}})\f$.
    double CDHB = 0.; ///< The dimension-6 operator coefficient \f$C_{DHB}(\Lambda_{\rm{EW}})\f$.
    double CDHW = 0.; ///< The dimension-6 operator coefficient \f$C_{DHW}(\Lambda_{\rm{EW}})\f$.
    double CDB = 0.; ///< The dimension-6 operator coefficient \f$C_{DB}(\Lambda_{\rm{EW}})\f$.
    double CDW = 0.; ///< The dimension-6 operator coefficient \f$C_{DW}(\Lambda_{\rm{EW}})\f$.
    double CT = 0.; ///< The dimension-6 operator coefficient \f$C_{T}(\Lambda_{\rm{EW}})\f$.
    ////////////////////////////////////////////////////////////////////////////////////////////////


    std::array<std::array<double, 3>,3> CHl1R = {}; ///< The dimension-6 operator coefficient \f$(C_{Hl}^{(1)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHl1I = {}; ///< The dimension-6 operator coefficient \f$(C_{Hl}^{(1)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHl3R = {}; ///< The dimension-6 operator coefficient \f$(C_{Hl}^{(3)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHl3I = {}; ///< The dimension-6 operator coefficient \f$(C_{Hl}^{(3)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHeR = {}; ///< The dimension-6 operator coefficient \f$(C_{He})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHeI = {}; ///< The dimension-6 operator coefficient \f$(C_{He})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHq1R = {}; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(1)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHq1I = {}; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(1)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHq3R = {}; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(3)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHq3I = {}; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(3)})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHuR = {}; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHuI = {}; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHdR = {}; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHdI = {}; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHudR = {}; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CHudI = {}; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CeHR = {}; ///< The dimension-6 operator coefficient \f$(C_{eH})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CeHI = {}; ///< The dimension-6 operator coefficient \f$(C_{eH})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuHR = {}; ///< The dimension-6 operator coefficient \f$(C_{uH})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuHI = {}; ///< The dimension-6 operator coefficient \f$(C_{uH})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdHR = {}; ///< The dimension-6 operator coefficient \f$(C_{dH})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdHI = {}; ///< The dimension-6 operator coefficient \f$(C_{dH})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuGR = {}; ///< The dimension-6 operator coefficient \f$(C_{uG})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuGI = {}; ///< The dimension-6 operator coefficient \f$(C_{uG})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuWR = {}; ///< The dimension-6 operator coefficient \f$(C_{uW})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuWI = {}; ///< The dimension-6 operator coefficient \f$(C_{uW})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuBR = {}; ///< The dimension-6 operator coefficient \f$(C_{uB})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CuBI = {}; ///< The dimension-6 operator coefficient \f$(C_{uB})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdGR = {}; ///< The dimension-6 operator coefficient \f$(C_{dG})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdGI = {}; ///< The dimension-6 operator coefficient \f$(C_{dG})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdWR = {}; ///< The dimension-6 operator coefficient \f$(C_{dW})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdWI = {}; ///< The dimension-6 operator coefficient \f$(C_{dW})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdBR = {}; ///< The dimension-6 operator coefficient \f$(C_{dB})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CdBI = {}; ///< The dimension-6 operator coefficient \f$(C_{dB})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CeWR = {}; ///< The dimension-6 operator coefficient \f$(C_{eW})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CeWI = {}; ///< The dimension-6 operator coefficient \f$(C_{eW})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CeBR = {}; ///< The dimension-6 operator coefficient \f$(C_{eB})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<double, 3>,3> CeBI = {}; ///< The dimension-6 operator coefficient \f$(C_{eB})_{ij}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CllR = {}; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CllI = {}; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clq1R = {}; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clq1I = {}; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clq3R = {}; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clq3I = {}; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CeeR = {}; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CeeI = {}; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CeuR = {}; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CeuI = {}; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CedR = {}; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CedI = {}; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CleR = {}; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CleI = {}; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CluR = {}; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CluI = {}; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CldR = {}; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CldI = {}; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CqeR = {}; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CqeI = {}; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CledqR = {}; ///< The dimension-6 operator coefficient \f$(C_{ledq})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CledqI = {}; ///< The dimension-6 operator coefficient \f$(C_{ledq})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqq1R = {}; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqq1I = {}; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqq3R = {}; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqq3I = {}; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CuuR = {}; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CuuI = {}; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CddR = {}; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> CddI = {}; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cud1R = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cud1I = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cud8R = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cud8I = {}; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqu1R = {}; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqu1I = {}; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqu8R = {}; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqu8I = {}; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqd1R = {}; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqd1I = {}; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqd8R = {}; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cqd8I = {}; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cquqd1R = {}; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cquqd1I = {}; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cquqd8R = {}; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Cquqd8I = {}; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clequ1R = {}; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clequ1I = {}; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(1)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clequ3R = {}; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(3)})_{ijkm}(\Lambda_{\rm{EW}})\f$.
    std::array<std::array<std::array<std::array<double, 3>,3>,3>,3> Clequ3I = {}; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(3)})_{ijkm}(\Lambda_{\rm{EW}})\f$.


private:
    const NPSMEFTd6General & mySMEFT;
    double LambdaNP2;
    double v2;
    double v;
    gslpp::matrix<gslpp::complex> VuL,VuR,VdL,VdR,VeL,VeR;
    
};

#endif /* NPSMEFTD6GENERALMATCHING_H */


