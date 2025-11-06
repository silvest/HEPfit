#ifndef BCLNU_H
#define BCLNU_H

class StandardModel;
#include <gsl/gsl_integration.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResultPtr.h>
#include <gsl/gsl_spline.h>
#include <memory>

class BClnu {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay (generic one containing b-quark)
     * @param[in] meson_j final meson of the decay (generic one containing c-quark)
     * @param[in] lep_i final leptons of the decay
     */
    BClnu(const StandardModel& SM_i, QCD::meson meson_i);

    /**
     * @brief Destructor.
     */
    virtual ~BClnu();

    /**
     * @brief Initialization of the parameters (later on we'll have the update)
     */
    std::vector<std::string> initializeBClnuParameters();

    /**
     * @brief Getter, will be useful later in moments computation
     */

    double getAmplsqfactor();

    double getBR();

    double computembkin();

    /**
     * @brief Functions for the Q2 moments computation
     * @param[in] i is the number of the moment
     * @param[in] u is a kinematic quantity defined in the moments functions that we'll see later on
     * @param[in] r is the ratio between the c - b quark masses
     * @param[in] dQ2 is the number of derivatives wrt q2
     * @param[in] dr is the number of derivatives wrt r
     * @param[in] isPERP is a flag for the choice of HQE basis
     */
    double Q2X(int i, double u, double r, int dQ2, int dr);

    double Q2XG(int i, double u, double r, int dQ2, int dr);

    double Q2XD(int i, double u, double r, int dQ2, int dr, int isPERP);

    /**
     * @brief Functions for the El moments computation
     * @param[in] i is the number of the moment
     * @param[in] elcuthat is the cut energy in units of b - quark masses
     * @param[in] r is the ratio between the c - b quark masses
     * @param[in] dEl is the number of derivatives wrt El
     * @param[in] dr is the number of derivatives wrt r
     */

    double ElX(int i, double elcuthat, double r, int dEl, int dr);

    double ElXpi(int i, double elcuthat, double r, int dEl, int dr);

    double ElXG(int i, double elcuthat, double r, int dEl, int dr);

    double ElXD(int i, double elcuthat, double r, int dEl, int dr);

    double ElXLS(int i, double elcuthat, double r, int dEl, int dr);

    /**
     * @brief Functions for the Mx moments computation
     * @param[in] i is the first moment number
     * @param[in] j is the second moment number. We need two of them since they're mixed moments by definition
     * @param[in] elcuthat is the cut energy in units of b - quark masses
     * @param[in] r is the ratio between the c - b quark masses
     * @param[in] dEl is the number of derivatives wrt El
     * @param[in] dr is the number of derivatives wrt r
     */

    double MxX(int i, int j, double elcuthat, double r, int dEl, int dr);

	double MxXpi(int i, int j, double elcuthat, double r, int dEl, int dr);

	double MxXG(int i, int j, double elcuthat, double r, int dEl, int dr);

	double MxXD(int i, int j, double elcuthat, double r, int dEl, int dr);

	double MxXLS(int i, int j, double elcuthat, double r, int dEl, int dr);

    /**
     * @brief Functions for the total rate computation
     * @param[in] r is the ratio between the c - b quark masses
     * @param[in] dr is the number of derivatives wrt r
     * @param[in] isPERP is a flag for the choice of HQE basis
     */

    double TrX(double r, int dr);

    double TrXG(double r, int dr);

    double TrXD(double r, int dr, int isPERP);

    double TrX1(double r, int dr);

    double TrX2(double r, int dr);

    double TrX3(double r, int dr);

    double TrXG1(double r, int dr);

    double TrXD1(double r, int dr);

    /**
     * @brief These functions are for the contributions due to the scheme change
     */

    double deltambkin(int order);

    double deltamcMS(int order);

    double RhoDPert(int order);

    double MuPiPert(int order);

    /**
     * @brief These are the functions to compute NLO and NLOpw corrections
     */

    double nloX1Q2(int moment, const std::string& cNP, double Q2cut, double r);

    double nloX1Q2DerivativeQ2(int moment, double Q2cut, double r);

    double nloX1Q2Derivativer(int moment, double Q2cut, double r);

    double nloX1El(int moment, const std::string& cNP, double Elcut, double r);

    double nloX1ElDerivativeEl(int moment, double Elcut, double r);

    double nloX1ElDerivativer(int moment, double Elcut, double r);

    double nloX1mix(int moment, int moment2, const std::string& cNP, double Elcut, double r);

    double nloX1mixDerivativeEl(int moment, int moment2, double Elcut, double r);

    double nloX1mixDerivativer(int moment, int moment2, double Elcut, double r);

    double nlopwX1Q2MuG(int moment, double Q2cut, double r);

	double nlopwX1Q2RhoD(int moment, double Q2cut, double r);

	double nlopwX1ElMuG(int moment, double Elcut, double r);

	double nlopwX1ElMuPi(int moment, double Elcut, double r);

	double nlopwX1ElRhoD(int moment, double Elcut, double r);

	double nlopwX1ElRhoLS(int moment, double Elcut, double r);

	double nlopwX1mixMuG(int moment, int moment2, double Elcut, double r);

	double nlopwX1mixMuPi(int moment, int moment2, double Elcut, double r);

	double nlopwX1mixRhoD(int moment, int moment2, double Elcut, double r);

	double nlopwX1mixRhoLS(int moment, int moment2, double Elcut, double r);

    double nnloX2Q2(int moment, double Q2cut, double r);

    // double nnloX2El(int moment, double Elcut, double r);

    double nnloX2ElnonBLM(int moment, double Elcut, double r);

    double nnlofitElNNLOnonBLM(int moment, double Elcut, double r);

    double nnloX2ElBLM(int moment, double Elcut, double r);

    // double nnloX2mix(int moment, int moment2, double Elcut, double r);

    double nnloX2mixnonBLM(int moment, int moment2, double Elcut, double r);

    double nnloX2mixBLM(int moment, int moment2, double Elcut, double r);

    double nnlofitMXNNLOnonBLM(int moment, double Elcut, double r);

    /**
     * @brief Functions for the moments computation
     */

    double Q2moment_1(double q2cut);

    double Q2moment_2(double q2cut);

    double Q2moment_3(double q2cut);

    double Elmoment_1(double elcut);

    double Elmoment_2(double elcut);

    double Elmoment_3(double elcut);

    double Mxmoment_1(double elcut);

    double Mxmoment_2(double elcut);

    double Mxmoment_3(double elcut);

    double XGamma();

    double XDeltaBR(double elcut);

    /**
     * @brief Function for the grid reading
     */

    std::vector<std::vector<std::vector<double>>> parse3DArray(const std::string& filename);

    std::vector<std::vector<std::vector<std::vector<double>>>> parse4DArray(const std::string& filename);

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> parse5DArray(const std::string& filename);


    /**
     * @brief These are useful functions for the NLO and NLOpw
     */

    std::vector<double> ChebPoints(double a, double b, int n);

    std::vector<double> ChebCoefficients(double a, double b, int n, const std::vector<double>& fpoints);

    double ChebPolynomial(const std::vector<double>& c, double a, double b, double s);

private:
    const StandardModel& mySM; /**< Model type */
    QCD::meson meson1; /**< Initial meson type */
    std::vector<std::string> bclnuParameters; /**< String of mandatory BClnu parameters */
    //bool NPanalysis; /**< A flag to switch to BSM analysis */

    double GF; /**< Fermi constant */
    double MB; /**< Mass of the B-meson */
    double r; /**< Dimensionless b - c mass ratio */
    double mcmc; /**< mc(mc) */
    double mbmb; /**< mb(mb) */
    double mbOS; /**< b-quark mass in the on-shell scheme */
    double mcOS; /**< c-quark mass in the on-shell scheme */
    double mbkin; /**< b-quark mass in the kinetic scheme */
    double MBhat; /**< Normalized mass of the B-meson */
    double mcMS; /**< c-quark mass in the @f$\bar{MS}@f$ */
    double alphas; /**< @f$\alpha_{s}(\mu_b)@f$ */
    double api; /**< @f$\alpha_{s}(\mu_b) / \pi @f$ */
    // gslpp::complex Vcb; /**< CKM factor of the decay */
    double BRXclnu; /**< Branching ratio as a parameter for the model */
    double beta0;
    double Aew;
    double amplsq_factor;   /**< Overall |A|^2 factor for gamma computations without |Vcb|^2 */

    double MBhat2;
    double MBhat3;

    double r2;
    double r3;
    double r4;
    double r5;
    double r6;
    double r7;
    double r8;
    double r9;
    double r10;
    double r11;
    double r12;
    double r13;
    double r14;
    double r15;
    double r16;
    double r17;
    double r18;
    double r19;
    double r20;
    double r21;
    double r22;
    double r23;
    double r24;
    double r25;
    double r26;
    double r27;
    double r28;
    double r29;
    double r30;
    double r31;
    double r32;
    double r33;
    double r34;
    double r35;
    double r36;
    double r38;
    double r40;
    double r42;
    double r44;

    double opr2;
    double omr2;
    double m1pr2;

    double Li2r;
    double Li2mr;

    double lr;
    double lr2;
    double l1pr;
    double l1mr;

    double delta;
    double delta2;
    double delta4;
    double delta5;
    double delta6;
    double delta7;
    double delta8;
    double delta9;
    double delta10;
    double delta11;
    double delta12;
    double delta13;
    double delta14;
    double delta15;
    double delta16;
    double delta17;
    double delta18;
    double delta19;
    double delta20;
    double delta21;
    double delta22;
    double delta23;
    double delta24;
    double delta25;
    double delta26;
    double delta27;
    double delta28;
    double delta29;
    double delta30;
    double delta31;
    double delta32;
    double delta33;
    double delta34;
    double delta35;
    double delta36;
    double delta37;
    double delta38;
    double delta39;
    double delta40;
    double delta41;
    double delta42;
    double delta43;
    double delta44;
    double delta45;
    double delta46;

    double l2d;
    double l2d2;
    double l2d3;


    /**
     * @brief Flag for the basis
     */
    int flagPERP;

    /**
     * @brief Wilson coefficients (still to be understood)
     */
    // double CVL; /**<Wilson coeffients @f$C_{V_L}@f$*/
    // double CVR; /**<Wilson coeffients @f$C_{V_R}@f$*/
    // double CSL; /**<Wilson coeffients @f$C_{S_L}'@f$*/
    // double CSR; /**<Wilson coeffients @f$C_{S_R}'@f$*/
    // double CT; /**<Wilson coeffients @f$C_{T}@f$*/

    /**
     * @brief Quantities useful for New Physics computations
     * It's always about Wilson Coefficients definitions
     */
    /*double CVL2;
    double CVR2;
    double CSL2;
    double CSR2;
    double CT2;
    double CSLT;
    double CSRT;
    double CSLSR;
    double CVLVR;*/

    /**
     * @brief HQE parameters for the functions
     */
    double MuPi; /**< HQE parameter of order 2 in mb powers */
    double MuG; /**< HQE parameter of order 2 in mb powers*/
    double RhoD; /**< HQE parameter of order 3 in mb powers */
    double RhoLS; /**< HQE parameter of order 3 in mb powers */

    /**
     * @brief Physical parameters for the functions
     */
    double scale_mbkin; /**< Renormalization scale for mb */
    double scale_mcMS; /**< Renormalization scale for mc */
    double scale_alphas; /**< Renormalization scale for alpha strong */
    double scale_muG; /**< Renormalization scale for HQE parameter */
    double scale_rhoD; /**< Renormalization scale for HQE parameter */
    double scale_rhoLS; /**< Renormalization scale for HQE parameter */

    /**
     * @brief Cache parameters
     */
     double cached_MuPi; /**< Cached value for MuPi */
     double cached_MuG; /**< Cached value for MuG */
     double cached_RhoD; /**< Cached value for RhoD */
     double cached_RhoLS; /**< Cached value for RhoLS */
     double cached_scale_mbkin; /**< Cached value for scale_mbkin */
     double cached_scale_mcMS; /**< Cached value for scale_mcMS */
     double cached_scale_muG; /**< Cached value for scale_muG */
     double cached_scale_rhoD; /**< Cached value for scale_rhoD */
     double cached_scale_rhoLS; /**< Cached value for scale_rhoLS */
    //  gslpp::complex cached_Vcb; /**< Cached value for Vcb */
     double cached_BRXclnu; /**< Cached value for BRXclnu */
     double cached_GF; /**< Cached value for GF */
     double cached_MB; /**< Cached value for MB */
     double cached_mcmc; /**< Cached value for mcmc */
     double cached_mbmb; /**< Cached value for mbmb */
     double cached_mbkin; /**< Cached value for mbkin */
     double cached_mcMS; /**< Cached value for mcMS */

     std::map<double, unsigned int > Q2X0ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X1ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X2ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X3ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X0ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X1ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X2ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X3ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X0ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X1ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X2ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X3ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X0ur20Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X1ur20Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X2ur20Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X3ur20Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X0ur11Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X1ur11Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X2ur11Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X3ur11Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X0ur02Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X1ur02Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X2ur02Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2X3ur02Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD0ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD1ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD2ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD3ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD0ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD1ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD2ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD3ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD0ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD1ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD2ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XD3ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG0ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG1ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG2ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG3ur00Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG0ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG1ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG2ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG3ur10Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG0ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG1ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG2ur01Cached;/**< Cache variable */
     std::map<double, unsigned int > Q2XG3ur01Cached;/**< Cache variable */

     std::map<double, unsigned int > X1Q2SM0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2SM1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2SM2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2SM3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2RhoD0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2RhoD1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2RhoD2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2RhoD3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2MuG0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2MuG1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2MuG2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2MuG3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2DerQ0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2DerQ1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2DerQ2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2DerQ3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2Derr0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2Derr1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2Derr2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1Q2Derr3Cached;/**< Cache variable */

     std::map<double, unsigned int > X2Q20Cached;/**< Cache variable */
     std::map<double, unsigned int > X2Q21Cached;/**< Cache variable */
     std::map<double, unsigned int > X2Q22Cached;/**< Cache variable */
     std::map<double, unsigned int > X2Q23Cached;/**< Cache variable */

     std::map<double, double > Q2X0ur00;/**< Cache variable */
     std::map<double, double > Q2X1ur00;/**< Cache variable */
     std::map<double, double > Q2X2ur00;/**< Cache variable */
     std::map<double, double > Q2X3ur00;/**< Cache variable */
     std::map<double, double > Q2X0ur10;/**< Cache variable */
     std::map<double, double > Q2X1ur10;/**< Cache variable */
     std::map<double, double > Q2X2ur10;/**< Cache variable */
     std::map<double, double > Q2X3ur10;/**< Cache variable */
     std::map<double, double > Q2X0ur01;/**< Cache variable */
     std::map<double, double > Q2X1ur01;/**< Cache variable */
     std::map<double, double > Q2X2ur01;/**< Cache variable */
     std::map<double, double > Q2X3ur01;/**< Cache variable */
     std::map<double, double > Q2X0ur20;/**< Cache variable */
     std::map<double, double > Q2X1ur20;/**< Cache variable */
     std::map<double, double > Q2X2ur20;/**< Cache variable */
     std::map<double, double > Q2X3ur20;/**< Cache variable */
     std::map<double, double > Q2X0ur11;/**< Cache variable */
     std::map<double, double > Q2X1ur11;/**< Cache variable */
     std::map<double, double > Q2X2ur11;/**< Cache variable */
     std::map<double, double > Q2X3ur11;/**< Cache variable */
     std::map<double, double > Q2X0ur02;/**< Cache variable */
     std::map<double, double > Q2X1ur02;/**< Cache variable */
     std::map<double, double > Q2X2ur02;/**< Cache variable */
     std::map<double, double > Q2X3ur02;/**< Cache variable */
     std::map<double, double > Q2XD0ur00;/**< Cache variable */
     std::map<double, double > Q2XD1ur00;/**< Cache variable */
     std::map<double, double > Q2XD2ur00;/**< Cache variable */
     std::map<double, double > Q2XD3ur00;/**< Cache variable */
     std::map<double, double > Q2XD0ur10;/**< Cache variable */
     std::map<double, double > Q2XD1ur10;/**< Cache variable */
     std::map<double, double > Q2XD2ur10;/**< Cache variable */
     std::map<double, double > Q2XD3ur10;/**< Cache variable */
     std::map<double, double > Q2XD0ur01;/**< Cache variable */
     std::map<double, double > Q2XD1ur01;/**< Cache variable */
     std::map<double, double > Q2XD2ur01;/**< Cache variable */
     std::map<double, double > Q2XD3ur01;/**< Cache variable */
     std::map<double, double > Q2XG0ur00;/**< Cache variable */
     std::map<double, double > Q2XG1ur00;/**< Cache variable */
     std::map<double, double > Q2XG2ur00;/**< Cache variable */
     std::map<double, double > Q2XG3ur00;/**< Cache variable */
     std::map<double, double > Q2XG0ur10;/**< Cache variable */
     std::map<double, double > Q2XG1ur10;/**< Cache variable */
     std::map<double, double > Q2XG2ur10;/**< Cache variable */
     std::map<double, double > Q2XG3ur10;/**< Cache variable */
     std::map<double, double > Q2XG0ur01;/**< Cache variable */
     std::map<double, double > Q2XG1ur01;/**< Cache variable */
     std::map<double, double > Q2XG2ur01;/**< Cache variable */
     std::map<double, double > Q2XG3ur01;/**< Cache variable */

     std::map<double, double > X1Q2SM0;/**< Cache variable */
     std::map<double, double > X1Q2SM1;/**< Cache variable */
     std::map<double, double > X1Q2SM2;/**< Cache variable */
     std::map<double, double > X1Q2SM3;/**< Cache variable */
     std::map<double, double > X1Q2RhoD0;/**< Cache variable */
     std::map<double, double > X1Q2RhoD1;/**< Cache variable */
     std::map<double, double > X1Q2RhoD2;/**< Cache variable */
     std::map<double, double > X1Q2RhoD3;/**< Cache variable */
     std::map<double, double > X1Q2MuG0;/**< Cache variable */
     std::map<double, double > X1Q2MuG1;/**< Cache variable */
     std::map<double, double > X1Q2MuG2;/**< Cache variable */
     std::map<double, double > X1Q2MuG3;/**< Cache variable */
     std::map<double, double > X1Q2DerQ0;/**< Cache variable */
     std::map<double, double > X1Q2DerQ1;/**< Cache variable */
     std::map<double, double > X1Q2DerQ2;/**< Cache variable */
     std::map<double, double > X1Q2DerQ3;/**< Cache variable */
     std::map<double, double > X1Q2Derr0;/**< Cache variable */
     std::map<double, double > X1Q2Derr1;/**< Cache variable */
     std::map<double, double > X1Q2Derr2;/**< Cache variable */
     std::map<double, double > X1Q2Derr3;/**< Cache variable */

     std::map<double, double > X2Q20;/**< Cache variable */
     std::map<double, double > X2Q21;/**< Cache variable */
     std::map<double, double > X2Q22;/**< Cache variable */
     std::map<double, double > X2Q23;/**< Cache variable */

     std::map<double, double > Q2X0ur00_2;/**< Cache variable */
     std::map<double, double > Q2X0ur00_3;/**< Cache variable */
     std::map<double, double > Q2X0ur00_4;/**< Cache variable */
     std::map<double, double > Q2X0ur00_5;/**< Cache variable */
     std::map<double, double > Q2X1ur00_2;/**< Cache variable */
     std::map<double, double > Q2X1ur00_3;/**< Cache variable */

     std::map<double, unsigned int > ElX0er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX1er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX2er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX3er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX0er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX1er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX2er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX3er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX0er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX1er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX2er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX3er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX0er20Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX1er20Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX2er20Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX3er20Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX0er11Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX1er11Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX2er11Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX3er11Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX0er02Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX1er02Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX2er02Cached;/**< Cache variable */
     std::map<double, unsigned int > ElX3er02Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD0er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD1er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD2er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD3er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD0er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD1er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD2er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD3er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD0er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD1er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD2er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXD3er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi0er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi1er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi2er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi3er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi0er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi1er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi2er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi3er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi0er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi1er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi2er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXpi3er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS0er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS1er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS2er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS3er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS0er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS1er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS2er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS3er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS0er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS1er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS2er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXLS3er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG0er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG1er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG2er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG3er00Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG0er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG1er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG2er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG3er10Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG0er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG1er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG2er01Cached;/**< Cache variable */
     std::map<double, unsigned int > ElXG3er01Cached;/**< Cache variable */

     std::map<double, unsigned int > X1ElSM0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElSM1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElSM2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElSM3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoD0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoD1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoD2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoD3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoLS0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoLS1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoLS2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElRhoLS3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuPi0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuPi1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuPi2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuPi3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuG0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuG1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuG2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElMuG3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerr0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerr1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerr2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerr3Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerEl0Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerEl1Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerEl2Cached;/**< Cache variable */
     std::map<double, unsigned int > X1ElDerEl3Cached;/**< Cache variable */

     std::map<double, unsigned int > X2ElBLM0Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElBLM1Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElBLM2Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElBLM3Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElnonBLM0Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElnonBLM1Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElnonBLM2Cached;/**< Cache variable */
     std::map<double, unsigned int > X2ElnonBLM3Cached;/**< Cache variable */
     std::map<double, unsigned int > fitElNNLOnonBLM0Cached;/**< Cache variable */
     std::map<double, unsigned int > fitElNNLOnonBLM1Cached;/**< Cache variable */
     std::map<double, unsigned int > fitElNNLOnonBLM2Cached;/**< Cache variable */
     std::map<double, unsigned int > fitElNNLOnonBLM3Cached;/**< Cache variable */

     std::map<double, double > ElX0er00;/**< Cache variable */
     std::map<double, double > ElX1er00;/**< Cache variable */
     std::map<double, double > ElX2er00;/**< Cache variable */
     std::map<double, double > ElX3er00;/**< Cache variable */
     std::map<double, double > ElX0er10;/**< Cache variable */
     std::map<double, double > ElX1er10;/**< Cache variable */
     std::map<double, double > ElX2er10;/**< Cache variable */
     std::map<double, double > ElX3er10;/**< Cache variable */
     std::map<double, double > ElX0er01;/**< Cache variable */
     std::map<double, double > ElX1er01;/**< Cache variable */
     std::map<double, double > ElX2er01;/**< Cache variable */
     std::map<double, double > ElX3er01;/**< Cache variable */
     std::map<double, double > ElX0er20;/**< Cache variable */
     std::map<double, double > ElX1er20;/**< Cache variable */
     std::map<double, double > ElX2er20;/**< Cache variable */
     std::map<double, double > ElX3er20;/**< Cache variable */
     std::map<double, double > ElX0er11;/**< Cache variable */
     std::map<double, double > ElX1er11;/**< Cache variable */
     std::map<double, double > ElX2er11;/**< Cache variable */
     std::map<double, double > ElX3er11;/**< Cache variable */
     std::map<double, double > ElX0er02;/**< Cache variable */
     std::map<double, double > ElX1er02;/**< Cache variable */
     std::map<double, double > ElX2er02;/**< Cache variable */
     std::map<double, double > ElX3er02;/**< Cache variable */
     std::map<double, double > ElXD0er00;/**< Cache variable */
     std::map<double, double > ElXD1er00;/**< Cache variable */
     std::map<double, double > ElXD2er00;/**< Cache variable */
     std::map<double, double > ElXD3er00;/**< Cache variable */
     std::map<double, double > ElXD0er10;/**< Cache variable */
     std::map<double, double > ElXD1er10;/**< Cache variable */
     std::map<double, double > ElXD2er10;/**< Cache variable */
     std::map<double, double > ElXD3er10;/**< Cache variable */
     std::map<double, double > ElXD0er01;/**< Cache variable */
     std::map<double, double > ElXD1er01;/**< Cache variable */
     std::map<double, double > ElXD2er01;/**< Cache variable */
     std::map<double, double > ElXD3er01;/**< Cache variable */
     std::map<double, double > ElXpi0er00;/**< Cache variable */
     std::map<double, double > ElXpi1er00;/**< Cache variable */
     std::map<double, double > ElXpi2er00;/**< Cache variable */
     std::map<double, double > ElXpi3er00;/**< Cache variable */
     std::map<double, double > ElXpi0er10;/**< Cache variable */
     std::map<double, double > ElXpi1er10;/**< Cache variable */
     std::map<double, double > ElXpi2er10;/**< Cache variable */
     std::map<double, double > ElXpi3er10;/**< Cache variable */
     std::map<double, double > ElXpi0er01;/**< Cache variable */
     std::map<double, double > ElXpi1er01;/**< Cache variable */
     std::map<double, double > ElXpi2er01;/**< Cache variable */
     std::map<double, double > ElXpi3er01;/**< Cache variable */
     std::map<double, double > ElXLS0er00;/**< Cache variable */
     std::map<double, double > ElXLS1er00;/**< Cache variable */
     std::map<double, double > ElXLS2er00;/**< Cache variable */
     std::map<double, double > ElXLS3er00;/**< Cache variable */
     std::map<double, double > ElXLS0er10;/**< Cache variable */
     std::map<double, double > ElXLS1er10;/**< Cache variable */
     std::map<double, double > ElXLS2er10;/**< Cache variable */
     std::map<double, double > ElXLS3er10;/**< Cache variable */
     std::map<double, double > ElXLS0er01;/**< Cache variable */
     std::map<double, double > ElXLS1er01;/**< Cache variable */
     std::map<double, double > ElXLS2er01;/**< Cache variable */
     std::map<double, double > ElXLS3er01;/**< Cache variable */
     std::map<double, double > ElXG0er00;/**< Cache variable */
     std::map<double, double > ElXG1er00;/**< Cache variable */
     std::map<double, double > ElXG2er00;/**< Cache variable */
     std::map<double, double > ElXG3er00;/**< Cache variable */
     std::map<double, double > ElXG0er10;/**< Cache variable */
     std::map<double, double > ElXG1er10;/**< Cache variable */
     std::map<double, double > ElXG2er10;/**< Cache variable */
     std::map<double, double > ElXG3er10;/**< Cache variable */
     std::map<double, double > ElXG0er01;/**< Cache variable */
     std::map<double, double > ElXG1er01;/**< Cache variable */
     std::map<double, double > ElXG2er01;/**< Cache variable */
     std::map<double, double > ElXG3er01;/**< Cache variable */

     std::map<double, double > X1ElSM0;/**< Cache variable */
     std::map<double, double > X1ElSM1;/**< Cache variable */
     std::map<double, double > X1ElSM2;/**< Cache variable */
     std::map<double, double > X1ElSM3;/**< Cache variable */
     std::map<double, double > X1ElRhoD0;/**< Cache variable */
     std::map<double, double > X1ElRhoD1;/**< Cache variable */
     std::map<double, double > X1ElRhoD2;/**< Cache variable */
     std::map<double, double > X1ElRhoD3;/**< Cache variable */
     std::map<double, double > X1ElRhoLS0;/**< Cache variable */
     std::map<double, double > X1ElRhoLS1;/**< Cache variable */
     std::map<double, double > X1ElRhoLS2;/**< Cache variable */
     std::map<double, double > X1ElRhoLS3;/**< Cache variable */
     std::map<double, double > X1ElMuPi0;/**< Cache variable */
     std::map<double, double > X1ElMuPi1;/**< Cache variable */
     std::map<double, double > X1ElMuPi2;/**< Cache variable */
     std::map<double, double > X1ElMuPi3;/**< Cache variable */
     std::map<double, double > X1ElMuG0;/**< Cache variable */
     std::map<double, double > X1ElMuG1;/**< Cache variable */
     std::map<double, double > X1ElMuG2;/**< Cache variable */
     std::map<double, double > X1ElMuG3;/**< Cache variable */
     std::map<double, double > X1ElDerr0;/**< Cache variable */
     std::map<double, double > X1ElDerr1;/**< Cache variable */
     std::map<double, double > X1ElDerr2;/**< Cache variable */
     std::map<double, double > X1ElDerr3;/**< Cache variable */
     std::map<double, double > X1ElDerEl0;/**< Cache variable */
     std::map<double, double > X1ElDerEl1;/**< Cache variable */
     std::map<double, double > X1ElDerEl2;/**< Cache variable */
     std::map<double, double > X1ElDerEl3;/**< Cache variable */

     std::map<double, double > X2ElBLM0;/**< Cache variable */
     std::map<double, double > X2ElBLM1;/**< Cache variable */
     std::map<double, double > X2ElBLM2;/**< Cache variable */
     std::map<double, double > X2ElBLM3;/**< Cache variable */
     std::map<double, double > X2ElnonBLM0;/**< Cache variable */
     std::map<double, double > X2ElnonBLM1;/**< Cache variable */
     std::map<double, double > X2ElnonBLM2;/**< Cache variable */
     std::map<double, double > X2ElnonBLM3;/**< Cache variable */
     std::map<double, double > fitElNNLOnonBLM0;/**< Cache variable */
     std::map<double, double > fitElNNLOnonBLM1;/**< Cache variable */
     std::map<double, double > fitElNNLOnonBLM2;/**< Cache variable */
     std::map<double, double > fitElNNLOnonBLM3;/**< Cache variable */

     std::map<double, double > ElX0er00_2;/**< Cache variable */
     std::map<double, double > ElX0er00_3;/**< Cache variable */
     std::map<double, double > ElX0er00_4;/**< Cache variable */
     std::map<double, double > ElX0er00_5;/**< Cache variable */
     std::map<double, double > ElX1er00_2;/**< Cache variable */
     std::map<double, double > ElX1er00_3;/**< Cache variable */

     std::map<double, unsigned int > MxX00er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX01er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX02er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX03er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX10er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX11er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX12er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX20er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX21er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX30er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX00er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX01er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX02er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX03er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX10er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX11er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX12er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX20er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX21er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX30er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX00er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX01er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX02er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX03er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX10er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX11er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX12er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX20er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX21er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX30er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX00er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX01er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX02er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX03er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX10er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX11er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX12er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX20er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX21er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX30er20Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX00er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX01er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX02er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX03er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX10er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX11er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX12er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX20er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX21er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX30er11Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX00er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX01er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX02er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX03er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX10er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX11er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX12er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX20er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX21er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxX30er02Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD00er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD01er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD02er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD03er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD10er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD11er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD12er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD20er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD21er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD30er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD00er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD01er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD02er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD03er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD10er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD11er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD12er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD20er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD21er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD30er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD00er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD01er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD02er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD03er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD10er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD11er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD12er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD20er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD21er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXD30er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi00er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi01er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi02er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi03er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi10er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi11er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi12er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi20er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi21er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi30er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi00er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi01er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi02er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi03er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi10er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi11er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi12er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi20er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi21er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi30er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi00er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi01er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi02er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi03er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi10er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi11er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi12er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi20er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi21er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXpi30er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS00er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS01er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS02er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS03er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS10er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS11er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS12er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS20er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS21er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS30er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS00er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS01er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS02er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS03er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS10er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS11er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS12er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS20er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS21er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS30er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS00er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS01er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS02er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS03er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS10er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS11er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS12er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS20er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS21er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXLS30er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG00er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG01er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG02er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG03er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG10er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG11er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG12er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG20er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG21er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG30er00Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG00er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG01er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG02er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG03er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG10er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG11er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG12er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG20er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG21er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG30er10Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG00er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG01er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG02er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG03er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG10er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG11er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG12er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG20er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG21er01Cached;/**< Cache variable */
     std::map<double, unsigned int > MxXG30er01Cached;/**< Cache variable */

     std::map<double, unsigned int > X1mixSM00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixSM30Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoD30Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixRhoLS30Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuPi30Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixMuG30Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerEl30Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr00Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr01Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr02Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr03Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr10Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr11Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr12Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr20Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr21Cached;/**< Cache variable */
     std::map<double, unsigned int > X1mixDerr30Cached;/**< Cache variable */

     std::map<double, unsigned int > X2mixBLM00Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM01Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM02Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM03Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM10Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM11Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM12Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM20Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM21Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixBLM30Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM00Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM01Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM02Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM03Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM10Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM11Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM12Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM20Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM21Cached;/**< Cache variable */
     std::map<double, unsigned int > X2mixnonBLM30Cached;/**< Cache variable */
     std::map<double, unsigned int > fitMXNNLOnonBLM0Cached;/**< Cache variable */
     std::map<double, unsigned int > fitMXNNLOnonBLM1Cached;/**< Cache variable */
     std::map<double, unsigned int > fitMXNNLOnonBLM2Cached;/**< Cache variable */
     std::map<double, unsigned int > fitMXNNLOnonBLM3Cached;/**< Cache variable */

     std::map<double, double > MxX00er00;/**< Cache variable */
     std::map<double, double > MxX01er00;/**< Cache variable */
     std::map<double, double > MxX02er00;/**< Cache variable */
     std::map<double, double > MxX03er00;/**< Cache variable */
     std::map<double, double > MxX10er00;/**< Cache variable */
     std::map<double, double > MxX11er00;/**< Cache variable */
     std::map<double, double > MxX12er00;/**< Cache variable */
     std::map<double, double > MxX20er00;/**< Cache variable */
     std::map<double, double > MxX21er00;/**< Cache variable */
     std::map<double, double > MxX30er00;/**< Cache variable */
     std::map<double, double > MxX00er10;/**< Cache variable */
     std::map<double, double > MxX01er10;/**< Cache variable */
     std::map<double, double > MxX02er10;/**< Cache variable */
     std::map<double, double > MxX03er10;/**< Cache variable */
     std::map<double, double > MxX10er10;/**< Cache variable */
     std::map<double, double > MxX11er10;/**< Cache variable */
     std::map<double, double > MxX12er10;/**< Cache variable */
     std::map<double, double > MxX20er10;/**< Cache variable */
     std::map<double, double > MxX21er10;/**< Cache variable */
     std::map<double, double > MxX30er10;/**< Cache variable */
     std::map<double, double > MxX00er01;/**< Cache variable */
     std::map<double, double > MxX01er01;/**< Cache variable */
     std::map<double, double > MxX02er01;/**< Cache variable */
     std::map<double, double > MxX03er01;/**< Cache variable */
     std::map<double, double > MxX10er01;/**< Cache variable */
     std::map<double, double > MxX11er01;/**< Cache variable */
     std::map<double, double > MxX12er01;/**< Cache variable */
     std::map<double, double > MxX20er01;/**< Cache variable */
     std::map<double, double > MxX21er01;/**< Cache variable */
     std::map<double, double > MxX30er01;/**< Cache variable */
     std::map<double, double > MxX00er20;/**< Cache variable */
     std::map<double, double > MxX01er20;/**< Cache variable */
     std::map<double, double > MxX02er20;/**< Cache variable */
     std::map<double, double > MxX03er20;/**< Cache variable */
     std::map<double, double > MxX10er20;/**< Cache variable */
     std::map<double, double > MxX11er20;/**< Cache variable */
     std::map<double, double > MxX12er20;/**< Cache variable */
     std::map<double, double > MxX20er20;/**< Cache variable */
     std::map<double, double > MxX21er20;/**< Cache variable */
     std::map<double, double > MxX30er20;/**< Cache variable */
     std::map<double, double > MxX00er11;/**< Cache variable */
     std::map<double, double > MxX01er11;/**< Cache variable */
     std::map<double, double > MxX02er11;/**< Cache variable */
     std::map<double, double > MxX03er11;/**< Cache variable */
     std::map<double, double > MxX10er11;/**< Cache variable */
     std::map<double, double > MxX11er11;/**< Cache variable */
     std::map<double, double > MxX12er11;/**< Cache variable */
     std::map<double, double > MxX20er11;/**< Cache variable */
     std::map<double, double > MxX21er11;/**< Cache variable */
     std::map<double, double > MxX30er11;/**< Cache variable */
     std::map<double, double > MxX00er02;/**< Cache variable */
     std::map<double, double > MxX01er02;/**< Cache variable */
     std::map<double, double > MxX02er02;/**< Cache variable */
     std::map<double, double > MxX03er02;/**< Cache variable */
     std::map<double, double > MxX10er02;/**< Cache variable */
     std::map<double, double > MxX11er02;/**< Cache variable */
     std::map<double, double > MxX12er02;/**< Cache variable */
     std::map<double, double > MxX20er02;/**< Cache variable */
     std::map<double, double > MxX21er02;/**< Cache variable */
     std::map<double, double > MxX30er02;/**< Cache variable */
     std::map<double, double > MxXD00er00;/**< Cache variable */
     std::map<double, double > MxXD01er00;/**< Cache variable */
     std::map<double, double > MxXD02er00;/**< Cache variable */
     std::map<double, double > MxXD03er00;/**< Cache variable */
     std::map<double, double > MxXD10er00;/**< Cache variable */
     std::map<double, double > MxXD11er00;/**< Cache variable */
     std::map<double, double > MxXD12er00;/**< Cache variable */
     std::map<double, double > MxXD20er00;/**< Cache variable */
     std::map<double, double > MxXD21er00;/**< Cache variable */
     std::map<double, double > MxXD30er00;/**< Cache variable */
     std::map<double, double > MxXD00er10;/**< Cache variable */
     std::map<double, double > MxXD01er10;/**< Cache variable */
     std::map<double, double > MxXD02er10;/**< Cache variable */
     std::map<double, double > MxXD03er10;/**< Cache variable */
     std::map<double, double > MxXD10er10;/**< Cache variable */
     std::map<double, double > MxXD11er10;/**< Cache variable */
     std::map<double, double > MxXD12er10;/**< Cache variable */
     std::map<double, double > MxXD20er10;/**< Cache variable */
     std::map<double, double > MxXD21er10;/**< Cache variable */
     std::map<double, double > MxXD30er10;/**< Cache variable */
     std::map<double, double > MxXD00er01;/**< Cache variable */
     std::map<double, double > MxXD01er01;/**< Cache variable */
     std::map<double, double > MxXD02er01;/**< Cache variable */
     std::map<double, double > MxXD03er01;/**< Cache variable */
     std::map<double, double > MxXD10er01;/**< Cache variable */
     std::map<double, double > MxXD11er01;/**< Cache variable */
     std::map<double, double > MxXD12er01;/**< Cache variable */
     std::map<double, double > MxXD20er01;/**< Cache variable */
     std::map<double, double > MxXD21er01;/**< Cache variable */
     std::map<double, double > MxXD30er01;/**< Cache variable */
     std::map<double, double > MxXpi00er00;/**< Cache variable */
     std::map<double, double > MxXpi01er00;/**< Cache variable */
     std::map<double, double > MxXpi02er00;/**< Cache variable */
     std::map<double, double > MxXpi03er00;/**< Cache variable */
     std::map<double, double > MxXpi10er00;/**< Cache variable */
     std::map<double, double > MxXpi11er00;/**< Cache variable */
     std::map<double, double > MxXpi12er00;/**< Cache variable */
     std::map<double, double > MxXpi20er00;/**< Cache variable */
     std::map<double, double > MxXpi21er00;/**< Cache variable */
     std::map<double, double > MxXpi30er00;/**< Cache variable */
     std::map<double, double > MxXpi00er10;/**< Cache variable */
     std::map<double, double > MxXpi01er10;/**< Cache variable */
     std::map<double, double > MxXpi02er10;/**< Cache variable */
     std::map<double, double > MxXpi03er10;/**< Cache variable */
     std::map<double, double > MxXpi10er10;/**< Cache variable */
     std::map<double, double > MxXpi11er10;/**< Cache variable */
     std::map<double, double > MxXpi12er10;/**< Cache variable */
     std::map<double, double > MxXpi20er10;/**< Cache variable */
     std::map<double, double > MxXpi21er10;/**< Cache variable */
     std::map<double, double > MxXpi30er10;/**< Cache variable */
     std::map<double, double > MxXpi00er01;/**< Cache variable */
     std::map<double, double > MxXpi01er01;/**< Cache variable */
     std::map<double, double > MxXpi02er01;/**< Cache variable */
     std::map<double, double > MxXpi03er01;/**< Cache variable */
     std::map<double, double > MxXpi10er01;/**< Cache variable */
     std::map<double, double > MxXpi11er01;/**< Cache variable */
     std::map<double, double > MxXpi12er01;/**< Cache variable */
     std::map<double, double > MxXpi20er01;/**< Cache variable */
     std::map<double, double > MxXpi21er01;/**< Cache variable */
     std::map<double, double > MxXpi30er01;/**< Cache variable */
     std::map<double, double > MxXLS00er00;/**< Cache variable */
     std::map<double, double > MxXLS01er00;/**< Cache variable */
     std::map<double, double > MxXLS02er00;/**< Cache variable */
     std::map<double, double > MxXLS03er00;/**< Cache variable */
     std::map<double, double > MxXLS10er00;/**< Cache variable */
     std::map<double, double > MxXLS11er00;/**< Cache variable */
     std::map<double, double > MxXLS12er00;/**< Cache variable */
     std::map<double, double > MxXLS20er00;/**< Cache variable */
     std::map<double, double > MxXLS21er00;/**< Cache variable */
     std::map<double, double > MxXLS30er00;/**< Cache variable */
     std::map<double, double > MxXLS00er10;/**< Cache variable */
     std::map<double, double > MxXLS01er10;/**< Cache variable */
     std::map<double, double > MxXLS02er10;/**< Cache variable */
     std::map<double, double > MxXLS03er10;/**< Cache variable */
     std::map<double, double > MxXLS10er10;/**< Cache variable */
     std::map<double, double > MxXLS11er10;/**< Cache variable */
     std::map<double, double > MxXLS12er10;/**< Cache variable */
     std::map<double, double > MxXLS20er10;/**< Cache variable */
     std::map<double, double > MxXLS21er10;/**< Cache variable */
     std::map<double, double > MxXLS30er10;/**< Cache variable */
     std::map<double, double > MxXLS00er01;/**< Cache variable */
     std::map<double, double > MxXLS01er01;/**< Cache variable */
     std::map<double, double > MxXLS02er01;/**< Cache variable */
     std::map<double, double > MxXLS03er01;/**< Cache variable */
     std::map<double, double > MxXLS10er01;/**< Cache variable */
     std::map<double, double > MxXLS11er01;/**< Cache variable */
     std::map<double, double > MxXLS12er01;/**< Cache variable */
     std::map<double, double > MxXLS20er01;/**< Cache variable */
     std::map<double, double > MxXLS21er01;/**< Cache variable */
     std::map<double, double > MxXLS30er01;/**< Cache variable */
     std::map<double, double > MxXG00er00;/**< Cache variable */
     std::map<double, double > MxXG01er00;/**< Cache variable */
     std::map<double, double > MxXG02er00;/**< Cache variable */
     std::map<double, double > MxXG03er00;/**< Cache variable */
     std::map<double, double > MxXG10er00;/**< Cache variable */
     std::map<double, double > MxXG11er00;/**< Cache variable */
     std::map<double, double > MxXG12er00;/**< Cache variable */
     std::map<double, double > MxXG20er00;/**< Cache variable */
     std::map<double, double > MxXG21er00;/**< Cache variable */
     std::map<double, double > MxXG30er00;/**< Cache variable */
     std::map<double, double > MxXG00er10;/**< Cache variable */
     std::map<double, double > MxXG01er10;/**< Cache variable */
     std::map<double, double > MxXG02er10;/**< Cache variable */
     std::map<double, double > MxXG03er10;/**< Cache variable */
     std::map<double, double > MxXG10er10;/**< Cache variable */
     std::map<double, double > MxXG11er10;/**< Cache variable */
     std::map<double, double > MxXG12er10;/**< Cache variable */
     std::map<double, double > MxXG20er10;/**< Cache variable */
     std::map<double, double > MxXG21er10;/**< Cache variable */
     std::map<double, double > MxXG30er10;/**< Cache variable */
     std::map<double, double > MxXG00er01;/**< Cache variable */
     std::map<double, double > MxXG01er01;/**< Cache variable */
     std::map<double, double > MxXG02er01;/**< Cache variable */
     std::map<double, double > MxXG03er01;/**< Cache variable */
     std::map<double, double > MxXG10er01;/**< Cache variable */
     std::map<double, double > MxXG11er01;/**< Cache variable */
     std::map<double, double > MxXG12er01;/**< Cache variable */
     std::map<double, double > MxXG20er01;/**< Cache variable */
     std::map<double, double > MxXG21er01;/**< Cache variable */
     std::map<double, double > MxXG30er01;/**< Cache variable */

     std::map<double, double > X1mixSM00;/**< Cache variable */
     std::map<double, double > X1mixSM01;/**< Cache variable */
     std::map<double, double > X1mixSM02;/**< Cache variable */
     std::map<double, double > X1mixSM03;/**< Cache variable */
     std::map<double, double > X1mixSM10;/**< Cache variable */
     std::map<double, double > X1mixSM11;/**< Cache variable */
     std::map<double, double > X1mixSM12;/**< Cache variable */
     std::map<double, double > X1mixSM20;/**< Cache variable */
     std::map<double, double > X1mixSM21;/**< Cache variable */
     std::map<double, double > X1mixSM30;/**< Cache variable */
     std::map<double, double > X1mixRhoD00;/**< Cache variable */
     std::map<double, double > X1mixRhoD01;/**< Cache variable */
     std::map<double, double > X1mixRhoD02;/**< Cache variable */
     std::map<double, double > X1mixRhoD03;/**< Cache variable */
     std::map<double, double > X1mixRhoD10;/**< Cache variable */
     std::map<double, double > X1mixRhoD11;/**< Cache variable */
     std::map<double, double > X1mixRhoD12;/**< Cache variable */
     std::map<double, double > X1mixRhoD20;/**< Cache variable */
     std::map<double, double > X1mixRhoD21;/**< Cache variable */
     std::map<double, double > X1mixRhoD30;/**< Cache variable */
     std::map<double, double > X1mixRhoLS00;/**< Cache variable */
     std::map<double, double > X1mixRhoLS01;/**< Cache variable */
     std::map<double, double > X1mixRhoLS02;/**< Cache variable */
     std::map<double, double > X1mixRhoLS03;/**< Cache variable */
     std::map<double, double > X1mixRhoLS10;/**< Cache variable */
     std::map<double, double > X1mixRhoLS11;/**< Cache variable */
     std::map<double, double > X1mixRhoLS12;/**< Cache variable */
     std::map<double, double > X1mixRhoLS20;/**< Cache variable */
     std::map<double, double > X1mixRhoLS21;/**< Cache variable */
     std::map<double, double > X1mixRhoLS30;/**< Cache variable */
     std::map<double, double > X1mixMuPi00;/**< Cache variable */
     std::map<double, double > X1mixMuPi01;/**< Cache variable */
     std::map<double, double > X1mixMuPi02;/**< Cache variable */
     std::map<double, double > X1mixMuPi03;/**< Cache variable */
     std::map<double, double > X1mixMuPi10;/**< Cache variable */
     std::map<double, double > X1mixMuPi11;/**< Cache variable */
     std::map<double, double > X1mixMuPi12;/**< Cache variable */
     std::map<double, double > X1mixMuPi20;/**< Cache variable */
     std::map<double, double > X1mixMuPi21;/**< Cache variable */
     std::map<double, double > X1mixMuPi30;/**< Cache variable */
     std::map<double, double > X1mixMuG00;/**< Cache variable */
     std::map<double, double > X1mixMuG01;/**< Cache variable */
     std::map<double, double > X1mixMuG02;/**< Cache variable */
     std::map<double, double > X1mixMuG03;/**< Cache variable */
     std::map<double, double > X1mixMuG10;/**< Cache variable */
     std::map<double, double > X1mixMuG11;/**< Cache variable */
     std::map<double, double > X1mixMuG12;/**< Cache variable */
     std::map<double, double > X1mixMuG20;/**< Cache variable */
     std::map<double, double > X1mixMuG21;/**< Cache variable */
     std::map<double, double > X1mixMuG30;/**< Cache variable */
     std::map<double, double > X1mixDerEl00;/**< Cache variable */
     std::map<double, double > X1mixDerEl01;/**< Cache variable */
     std::map<double, double > X1mixDerEl02;/**< Cache variable */
     std::map<double, double > X1mixDerEl03;/**< Cache variable */
     std::map<double, double > X1mixDerEl10;/**< Cache variable */
     std::map<double, double > X1mixDerEl11;/**< Cache variable */
     std::map<double, double > X1mixDerEl12;/**< Cache variable */
     std::map<double, double > X1mixDerEl20;/**< Cache variable */
     std::map<double, double > X1mixDerEl21;/**< Cache variable */
     std::map<double, double > X1mixDerEl30;/**< Cache variable */
     std::map<double, double > X1mixDerr00;/**< Cache variable */
     std::map<double, double > X1mixDerr01;/**< Cache variable */
     std::map<double, double > X1mixDerr02;/**< Cache variable */
     std::map<double, double > X1mixDerr03;/**< Cache variable */
     std::map<double, double > X1mixDerr10;/**< Cache variable */
     std::map<double, double > X1mixDerr11;/**< Cache variable */
     std::map<double, double > X1mixDerr12;/**< Cache variable */
     std::map<double, double > X1mixDerr20;/**< Cache variable */
     std::map<double, double > X1mixDerr21;/**< Cache variable */
     std::map<double, double > X1mixDerr30;/**< Cache variable */

     std::map<double, double > X2mixBLM00;/**< Cache variable */
     std::map<double, double > X2mixBLM01;/**< Cache variable */
     std::map<double, double > X2mixBLM02;/**< Cache variable */
     std::map<double, double > X2mixBLM03;/**< Cache variable */
     std::map<double, double > X2mixBLM10;/**< Cache variable */
     std::map<double, double > X2mixBLM11;/**< Cache variable */
     std::map<double, double > X2mixBLM12;/**< Cache variable */
     std::map<double, double > X2mixBLM20;/**< Cache variable */
     std::map<double, double > X2mixBLM21;/**< Cache variable */
     std::map<double, double > X2mixBLM30;/**< Cache variable */
     std::map<double, double > X2mixnonBLM00;/**< Cache variable */
     std::map<double, double > X2mixnonBLM01;/**< Cache variable */
     std::map<double, double > X2mixnonBLM02;/**< Cache variable */
     std::map<double, double > X2mixnonBLM03;/**< Cache variable */
     std::map<double, double > X2mixnonBLM10;/**< Cache variable */
     std::map<double, double > X2mixnonBLM11;/**< Cache variable */
     std::map<double, double > X2mixnonBLM12;/**< Cache variable */
     std::map<double, double > X2mixnonBLM20;/**< Cache variable */
     std::map<double, double > X2mixnonBLM21;/**< Cache variable */
     std::map<double, double > X2mixnonBLM30;/**< Cache variable */
     std::map<double, double > fitMXNNLOnonBLM0;/**< Cache variable */
     std::map<double, double > fitMXNNLOnonBLM1;/**< Cache variable */
     std::map<double, double > fitMXNNLOnonBLM2;/**< Cache variable */
     std::map<double, double > fitMXNNLOnonBLM3;/**< Cache variable */

     std::map<double, double > MxX00er00_2;/**< Cache variable */
     std::map<double, double > MxX00er00_3;/**< Cache variable */
     std::map<double, double > MxX00er00_4;/**< Cache variable */
     std::map<double, double > MxX00er00_5;/**< Cache variable */
     std::map<double, double > MxX01er00_2;/**< Cache variable */
     std::map<double, double > MxX01er00_3;/**< Cache variable */
     std::map<double, double > MxX10er00_2;/**< Cache variable */
     std::map<double, double > MxX10er00_3;/**< Cache variable */

     unsigned int MuPi_updated; /**< Cache flag for MuPi */
     unsigned int MuG_updated; /**< Cache flag for MuG */
     unsigned int RhoD_updated; /**< Cache flag for RhoD */
     unsigned int RhoLS_updated; /**< Cache flag for RhoLS */
     unsigned int scale_mbkin_updated; /**< Cache flag for scale_mbkin */
     unsigned int scale_mcMS_updated; /**< Cache flag for scale_mcMS */
     unsigned int scale_muG_updated; /**< Cache flag for scale_muG */
     unsigned int scale_rhoD_updated; /**< Cache flag for scale_rhoD */
     unsigned int scale_rhoLS_updated; /**< Cache flag for scale_rhoLS */
    //  unsigned int Vcb_updated; /**< Cache flag for Vcb */
     unsigned int BRXclnu_updated; /**< Cache flag for BRXclnu */
     unsigned int GF_updated; /**< Cache flag for GF */
     unsigned int MB_updated; /**< Cache flag for MB */
     unsigned int mcmc_updated; /**< Cache flag for mcmc */
     unsigned int mbmb_updated; /**< Cache flag for mbmb */

     unsigned int mcOS_updated; /**< Cache flag for mcOS */
     unsigned int mbOS_updated; /**< Cache flag for mbOS */
     unsigned int mbkin_updated; /**< Cache flag for mbkin */
     unsigned int MBhat_updated; /**< Cache flag for MBhat */
     unsigned int alphas_updated; /**< Cache flag for alphas */
     unsigned int mcMS_updated; /**< Cache flag for mcMS */
     unsigned int r_updated; /**< Cache flag for r */
     unsigned int amplsq_factor_updated; /**< Cache flag for amplsq_factor */

    std::vector<std::vector<std::vector<std::vector<double>>>> grid_np_Q2moments;
    std::vector<std::vector<std::vector<double>>> grid_DerQ2_Q2moments;
    std::vector<std::vector<std::vector<double>>> grid_Derr_Q2moments;
    std::vector<std::vector<std::vector<double>>> grid_NLO_MuG_Q2moments;
    std::vector<std::vector<std::vector<double>>> grid_NLO_RhoD_Q2moments;
    std::vector<std::vector<std::vector<double>>> grid_NNLO_Q2moments;

    std::vector<std::vector<std::vector<std::vector<double>>>> grid_np_Elmoments;
    std::vector<std::vector<std::vector<double>>> grid_DerEl_Elmoments;
    std::vector<std::vector<std::vector<double>>> grid_Derr_Elmoments;
    std::vector<std::vector<std::vector<double>>> grid_BLM_Elmoments;
    std::vector<std::vector<std::vector<double>>> grid_NLO_MuG_Elmoments;
    std::vector<std::vector<std::vector<double>>> grid_NLO_MuPi_Elmoments;

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> grid_np_mixmoments;
    std::vector<std::vector<std::vector<std::vector<double>>>> grid_DerEl_mixmoments;
    std::vector<std::vector<std::vector<std::vector<double>>>> grid_Derr_mixmoments;
    std::vector<std::vector<std::vector<std::vector<double>>>> grid_BLM_mixmoments;
    std::vector<std::vector<std::vector<std::vector<double>>>> grid_NLO_MuG_mixmoments;
    std::vector<std::vector<std::vector<std::vector<double>>>> grid_NLO_MuPi_mixmoments;

    /**
     * @brief The update parameter method for BClnu.
     */
    void updateParameters();

    /**
     * @brief The caching method for BClnu.
     */
    void checkCache();

    /**
     * @brief The caching method for BClnu.
     */
    void initializeQ2AuxiliaryFunctions(double q2cut);

    /**
     * @brief The caching method for BClnu.
     */
    void initializeElAuxiliaryFunctions(double Elcuthat);

    /**
     * @brief The caching method for BClnu.
     */
    void initializeMxAuxiliaryFunctions(double Elcuthat);
};

#endif /* BCLNU_H */
