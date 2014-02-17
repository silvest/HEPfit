/*
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef QCD_H
#define	QCD_H

#include "Model.h"
#include "Meson.h"
#include "OrderScheme.h"
#define MEPS 1.e-10 // mass precision

/**
 * @class QCD
 * @ingroup StandardModel
 * @brief A class for parametrs related to QCD, hadrons and quarks.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is a Model class that assigns and updates parameters
 * related to and derived from QCD. A complete list of parameters in the QCD
 * class can be found below. Thi s class includes, but is not limited to,
 * the running of the strong coupling constant (Full NNLO), running of the quark
 * masses and conversions between pole mass and \f$\overline{\mathrm{MS}}\f$ mass. All hadronization
 * parameters like the bag parameters for the mesons and their decay constants are
 * assigned and updated by this class.
 *
 * Model parameters:
 *
 * @anchor QCDInitialization
 * <h3>Initialization</h3>
 *
 * The constructor QCD() sets the charge and isospin of the quarks. It also sets the
 * mass scale of the light quarks UP, DOWN and STRANGE to 2 \f$GeV\f$. The cache is initialized
 * too along with the computation of \f$\zeta(2)\f$ and \f$\zeta(3)\f$.
 *
 * The initializations and updates of the model parameters and flags are explained
 * below.
 *
 *
 * @anchor QCDParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of QCD are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%AlsMz</td>
 *   <td class="mod_symb">@f$\alpha_s(M_Z)@f$</td>
 *   <td class="mod_desc">The strong coupling constant at the Z-boson mass.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mz</td>
 *   <td class="mod_symb">@f$M_Z@f$</td>
 *   <td class="mod_desc">The mass of the \f$Z\f$ boson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mup</td>
 *   <td class="mod_symb">@f$m_{u}@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ mass of the up quark at 2 GeV, \f$m_u(2\,\mathrm{GeV})\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mdown</td>
 *   <td class="mod_symb">@f$m_{d}@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ mass of the down quark at 2 GeV, \f$m_d(2\,\mathrm{GeV})\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mcharm</td>
 *   <td class="mod_symb">@f$m_{c}@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ scale-invariant mass of the charm quark, \f$m_c(m_c)\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mstrange</td>
 *   <td class="mod_symb">@f$m_{s}@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ mass of the strange quark at 2 GeV , \f$m_s(2\,\mathrm{GeV})\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mtop</td>
 *   <td class="mod_symb">@f$m_{t}@f$</td>
 *   <td class="mod_desc">The pole mass of the top quark in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mbottom</td>
 *   <td class="mod_symb">@f$m_{b}@f$</td>
 *   <td class="mod_desc">the \f$\overline{\mathrm{MS}}\f$ scale-invariant mass of the bottom quark, \f$m_b(m_b)\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mut</td>
 *   <td class="mod_symb">@f$\mu_t@f$</td>
 *   <td class="mod_desc">the threshold between six- and five-flavour theory in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mub</td>
 *   <td class="mod_symb">@f$\mu_b@f$</td>
 *   <td class="mod_desc">the threshold between five- and four-flavour theory in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%muc</td>
 *   <td class="mod_symb">@f$\mu_c@f$</td>
 *   <td class="mod_desc">the threshold between four- and three-flavour theory in GeV.</td>
 * </tr>
 * </table>
 *
 * The parameters below are associated with flavour observables
 * <table class="model">
 * <tr>
 *   <td class="mod_name">%MBd</td>
 *   <td class="mod_symb">@f$M_{B_d}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_d \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBd</td>
 *   <td class="mod_symb">@f$\tau_{B_d}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_d \f$ meson in \f$ ps^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBs</td>
 *   <td class="mod_symb">@f$M_{B_s}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_s \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBs</td>
 *   <td class="mod_symb">@f$\tau_{B_s}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_d \f$ meson in \f$ ps^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBp</td>
 *   <td class="mod_symb">@f$M_{B^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MK0</td>
 *   <td class="mod_symb">@f$M_{K^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MKp</td>
 *   <td class="mod_symb">@f$M_{K^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MD</td>
 *   <td class="mod_symb">@f$M_{D^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ D^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKl</td>
 *   <td class="mod_symb">@f$\tau_{K_L}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K_L \f$ meson in \f$ ps^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKp</td>
 *   <td class="mod_symb">@f$\tau{K^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K^\pm \f$ meson in \f$ ps^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBs</td>
 *   <td class="mod_symb">@f$F_{B_s}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ B_s \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBsoFBd</td>
 *   <td class="mod_symb">@f$F_{B_d}/F_{B_d}@f$</td>
 *   <td class="mod_desc">The ratio \f$ F_{B_s}/F_{B_d} \f$ necessary to compute \f$ F_{B_s} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FD</td>
 *   <td class="mod_symb">@f$F_{D^0}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ D^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsoBBd</td>
 *   <td class="mod_symb">@f$B_{B_s}/B_{B_d}@f$</td>
 *   <td class="mod_desc">The ratio \f$ B_{B_s}/B_{B_d} \f$ necessary to compute \f$ B_{B_s} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBs1 - %BBs5</td>
 *   <td class="mod_symb">@f$B^1_{B_s} - B^5_{B_s}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsscale</td>
 *   <td class="mod_symb">@f$\mu_{B_{s}}@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ B_s \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ B_s \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BD1 - %BD5</td>
 *   <td class="mod_symb">@f$B^1_{D} - B^5_{D}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5\f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BDscale</td>
 *   <td class="mod_symb">@f$\mu_D@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ D_0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BDscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ D_0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK1 - %BK5</td>
 *   <td class="mod_symb">@f$B^1_{K} - B^5_{K}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5\f$ in \f$ \Delta s = 2 \f$ processes in \f$ K^0 \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKscale</td>
 *   <td class="mod_symb">@f$\mu_K@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ K^0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ K^0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK(1/2)1 - %BK(1/2)10</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKd_scale</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKd_scheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK(3/2)1 - %BK(3/2)10</td>
 *   <td class="mod_symb"></td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ReA0_Kd</td>
 *   <td class="mod_symb">@f${\cal Re}(A_0(K\to\pi\pi))@f$</td>
 *   <td class="mod_desc">The experimental value of the real part of the amplitude for \f$K^0\to\pi\pi\f$ with \f$\Delta I=0\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ReA2_Kd</td>
 *   <td class="mod_symb">@f${\cal Re}(A_2(K\to\pi\pi))@f$</td>
 *   <td class="mod_desc">the experimental value of the real part of the amplitude for \f$K^0\to\pi\pi\f$ with \f$\Delta I=2\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Omega_eta_etap</td>
 *   <td class="mod_symb">@f$\Omega_{\eta/\eta'}@f$</td>
 *   <td class="mod_desc">The isospin breaking contribution in \f$K^0\to\pi\pi\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_Kp_P0enu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}{K^+\to\pi^0e^+\nu}@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$K^+\to\pi^0e^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_Kp_munu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(K^+\to\mu^+\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$K^+\to\mu^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_B_Xcenu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(B\to X_ce\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$B\to X_c e\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaP_cu</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The long-distance correction to the charm contribution of \f$K^+\to\pi^+\nu\bar{\nu}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kl</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">the isospin breaking corrections between @f$K_L\to\pi^0\nu\bar{\nu}@f$ and \f$K^+\to\pi^0 e^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kp</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The isospin breaking corrections between @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and \f$K^+\to\pi^0 e^+\nu\f$.</td>
 * </tr>
 * </table>
 *
 * The set of the model parameters are initialized and updated with the methods
 * Init() and Update(), respectively, where the former calls the latter.
 * In Update(), the methods PreUpdate() and PostUpdate() are called to run all
 * the procedures that are need to be executed before and after the model parameters
 * are updated. The \f$\overline{\mathrm{MS}}\f$ mass for the top quark is computed and the scale set
 * in PostUpdate() with the updated parameters. Inside the Update() method, the
 * individual model parameter is assigned with the protected member function
 * setParameter().
 *
 */

class QCD: public Model {
public:
    
    /**
     * @brief An enum type for mesons.
     */
    enum meson
    {
        B_D, /**< @f$B_d@f$ meson */
        B_S, /**< @f$B_s@f$ meson */
        B_P, /**< @f$B^\pm@f$ meson */
        K_0, /**< @f$K^0@f$ meson */
        K_P, /**< @f$K^\pm@f$ meson */
        D_0, /**< @f$D^0@f$ meson */
        P_0, /**< @f$\pi^0@f$ meson */
        P_P, /**< @f$\pi^\pm@f$ meson */
        MESON_END /**< The size of this enum. */
    };
    
    /**
     * @brief An enum type for quarks.
     */
    enum quark
    {
        UP, /**< Up quark */
        DOWN, /**< Down quark */
        CHARM, /**< Charm quark */
        STRANGE, /**< Strange quark */
        TOP, /**< Top quark */
        BOTTOM /**< Bottom quark */
    };
    
    static const int NQCDvars = 78; /**< The number of model parameters in QCD. */
    
    /**
     * @brief An array containing the labels under which all QCD parameters are stored
     * in a vector of ModelParameter via InputParser::ReadParameters().
     */
    static const std::string QCDvars[NQCDvars];
    
    /**
     * @brief Constructor.
     */
    QCD();
    
    /**
     * @brief A method to fetch the name of %QCD.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "QCD";
    }
    
    /**
     * @brief Converts an object of the enum type "orders" to the corresponding string.
     * @param[in] order an object of the enum type "orders"
     * @return the string of the given "order"
     */
    std::string orderToString(const orders order) const;
    
    ////////////////////////////////////////////////////////////////////////
    // Parameters
    
    /**
     * @brief Initializes the QCD parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %QCD
     * @details This method resets the internal flags #requireYu, #requireYd,
     * #computeBd #computeFBd and #computemt before updating the model parameters with the method Update().
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %QCD.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %QCD.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated. This includes 
     * \li computing the decay constatnt \f$F_{B_D}\f$ from \f$F_{B_s}\f$
     * \li computing the bag parameters \f$B_{B_d}\f$ from \f$B_{B_s}\f$
     * \li computing the \f$\overline{\rm MS}\f$ mass of the top quark at the \f$\overline{\rm MS}\f$ mass,
     * \f$m_t^{\overline{\rm MS}}(m_t^{\overline{\rm MS}})\f$ and setting the scale at the same value.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %StandardModel
     * have been provided in the model configuration file.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    
    ////////////////////////////////////////////////////////////////////////
    // Flags
    
    /**
     * @brief A method to set a flag of %QCD.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to set a flag of %QCD.
     * @param[in] name name of a model flag
     * @param[in] value the string to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);
    
    /**
     * @brief A method to check the sanity of the set of model flags.
     * @return a boolean that is true if the set of model flags is sane
     */
    virtual bool CheckFlags() const;
    
    
    ////////////////////////////////////////////////////////////////////////
    // get and set methods for class members
    
    /**
     * @brief A get method to access a meson as an object of the type Meson.
     * @param[in] m the name of a meson
     * @return the object of the meson specified in the argument
     */
    Meson getMesons(const meson m) const
    {
        return mesons[m];
    }
    
    /**
     * @brief A get method to access a quark as an object of the type Particle.
     * @param[in] q The name of a quark.
     * @return the object of the quark found in the argument
     */
    Particle getQuarks(const quark q) const
    {
        return quarks[q];
    }
    
    /**
     * @brief A get method to access the value of \f$\alpha_s(M_Z)\f$
     * @return the strong coupling constant at @f$M_Z@f$, @f$\alpha_s(M_Z)@f$
     */
    double getAlsMz() const
    {
        return AlsMz;
    }
    
    /**
     * @brief Sets the strong coupling constant at @f$M_Z@f$, @f$\alpha_s(M_Z)@f$.
     * @param[in] AlsMz @f$\alpha_s(M_Z)@f$
     */
    void setAlsMz(double AlsMz)
    {
        this->AlsMz = AlsMz;
    }
    
    /**
     * @brief A get method to access the mass of the \f$Z\f$ boson \f$M_Z\f$
     * @return the @f$Z@f$-boson mass @f$M_Z@f$
     */
    double getMz() const
    {
        return Mz;
    }
    
    /**
     * @brief Sets the @f$Z@f$ boson mass @f$M_Z@f$.
     * @param[in] Mz @f$M_Z@f$ in GeV
     */
    void setMz(double Mz)
    {
        this->Mz = Mz;
    }
    
    /**
     * @brief A get method to access the number of colours \f$N_c\f$
     * @return the number of colours
     */
    double getNc() const
    {
        return Nc;
    }
    
    /**
     * @brief Sets the number of colours.
     * @param[in] Nc the number of colours
     */
    void setNc(double Nc)
    {
        this->Nc = Nc;
    }
    
    /**
     * @brief A get method to access he threshold between six- and five-flavour theory in GeV
     * @return the threshold \f$\mu_t\f$
     */
    double getMut() const
    {
        return mut;
    }
    
    /**
     * @brief Sets the threshold between six- and five-flavour theory.
     * @param[in] mut the threshold between six- and five-flavour theory in GeV \f$\mu_t\f$
     */
    void setMut(double mut)
    {
        this->mut = mut;
    }
    
    /**
     * @brief A get method to access he threshold between five- and four-flavour theory in GeV
     * @return the threshold \f$\mu_b\f$
     */
    double getMub() const
    {
        return mub;
    }
    
    /**
     * @brief Sets the threshold between five- and four-flavour theory.
     * @param[in] mub the threshold between five- and four-flavour theory in GeV \f$\mu_b\f$
     */
    void setMub(double mub)
    {
        this->mub = mub;
    }
    
    /**
     * @brief A get method to access he threshold between four- and three-flavour theory in GeV
     * @return the threshold \f$\mu_c\f$
     */
    double getMuc() const
    {
        return muc;
    }
    
    /**
     * @brief Set the threshold between four- and three-flavour theory.
     * @param[in] muc the threshold between four- and three-flavour theory in GeV \f$\mu_c\f$
     */
    void setMuc(double muc)
    {
        this->muc = muc;
    }
    
    /**
     * @brief A get method to access the pole mass of the top quark
     * @return the pole mass of the top quark \f$m_t^{pole}\f$
     */
    double getMtpole() const
    {
        return mtpole;
    }
    
    /**
     * @brief A get method to access the Casimir Fator of %QCD
     * @return the Casimir factor
     */
    double getCF() const
    {
        return CF;
    }
    
    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta b = 2\f$
     * process in the \f$B_d\f$ meson system.
     * @return The vector of bag parameters
     */
    BParameter getBBd() const
    {
        return BBd;
    }
    
    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta b = 2\f$
     * process in the \f$B_s\f$ meson system.
     * @return The vector of bag parameters
     */
    BParameter getBBs() const
    {
        return BBs;
    }
    
    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta c = 2\f$
     * process in the \f$D^0\f$ meson system
     * @return The vector of bag parameters
     */
    BParameter getBD() const
    {
        return BD;
    }
    
    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta s = 2\f$
     * process in the \f$K^0\f$ meson system
     * @return The vector of bag parameters
     */
    BParameter getBK() const
    {
        return BK;
    }
    
    /**
     * @return
     */
    BParameter getBKd1() const
    {
        return BKd1;
    }
    
    /**
     * @return
     */
    BParameter getBKd3() const
    {
        return BKd3;
    }
    
    /**
     * @return The experimental value of the real part of the amplitude for
     * @f$K^0\to\pi\pi@f$ with @f$\Delta I=0@f$.
     */
    double getReA0_Kd() const
    {
        return ReA0_Kd;
    }
    
    /**
     * @return The experimental value of the real part of the amplitude for
     * @f$K^0\to\pi\pi@f$ with @f$\Delta I=2@f$.
     */
    double getReA2_Kd() const
    {
        return ReA2_Kd;
    }
    
    /**
     * @return The isospin breaking contribution in @f$K^0\to\pi\pi@f$.
     */
    double getOmega_eta_etap() const
    {
        return Omega_eta_etap;
    }
    
    /**
     * @return The experimental value for the branching ratio of @f$K^+\to\pi^0e^+\nu@f$.
     */
    double getBr_Kp_P0enu() const
    {
        return Br_Kp_P0enu;
    }
    
    /**
     * @return The experimental value for the branching ratio of @f$K^+\to\mu^+\nu@f$.
     */
    double getBr_Kp_munu() const
    {
        return Br_Kp_munu;
    }
    
    /**
     * @return The experimental value for the branching ratio of @f$B\to X_c e\nu@f$.
     */
    double getBr_B_Xcenu() const
    {
        return Br_B_Xcenu;
    }
    
    /**
     * @return The long-distance correction to the charm contribution of @f$K^+\to\pi^+\nu\bar{\nu}@f$.
     *
     * References:
     * [<A HREF="http://inspirehep.net/record/678222?ln=en" target="blank">Isidori et al.(2005)</A>],
     * [<A HREF="http://inspirehep.net/record/712083?ln=en" target="blank">Buras et al.(2006)</A>]
     */
    double getDeltaP_cu() const
    {
        return DeltaP_cu;
    }
    
    /**
     * @return The isospin breaking corrections between
     * @f$K_L\to\pi^0\nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$.
     */
    double getIB_Kl() const
    {
        return IB_Kl;
    }
    
    /**
     * @return The isospin breaking corrections between
     * @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$.
     */
    double getIB_Kp() const
    {
        return IB_Kp;
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief For accessing the active flavour threshold scales.
     * @param[in] i the index referring to active flavour thresholds.
     * @return the threshold scale: 1.0E10 (i = 0), \f$\mu_t\f$ (i = 1),
     * \f$\mu_b\f$ (i = 2), \f$\mu_c\f$ (i = 3) and 0. (default).
     */
    double Thresholds(const int i) const;
    
    /**
     * @brief The active flavour threshold above the scale \f$\mu\f$
     * as defined in QCD::Thresholds().
     * @param[in] mu a scale \f$\mu\f$ in GeV
     * @return the higher active flavour threshold
     */
    double AboveTh(const double mu) const;
    
    /**
     * @brief The active flavour threshold below the scale \f$\mu\f$
     * as defined in QCD::Thresholds().
     * @param[in] mu a scale \f$\mu\f$ in GeV
     * @return the lower active flavour threshold
     */
    double BelowTh(const double mu) const;
    
    /**
     * @brief The number of active flavour at scale @f$\mu@f$.
     * @param[in] mu a scale @f$\mu@f$ in GeV
     * @return active N_f
     */
    double Nf(const double mu) const;
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The \f$\beta_0(n_f)\f$ coefficient for a certain number of flavours \f$n_f\f$
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return @f$\beta_0(n_f)@f$
     */
    double Beta0(const double nf) const;
    
    /**
     * @brief The \f$\beta_1(n_f)\f$ coefficient for a certain number of flavours \f$n_f\f$
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return @f$\beta_1(n_f)@f$
     */
    double Beta1(const double nf) const;
    
    /**
     * @brief The \f$\beta_2(n_f)\f$ coefficient for a certain number of flavours \f$n_f\f$
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return @f$\beta_2(n_f)@f$
     */
    double Beta2(const double nf) const;
    
    /**
     * @brief Computes the running strong coupling @f$\alpha_s(\mu)@f$ from @f$\alpha_s(\mu_i)@f$
     * in the @f$\overline{\mathrm{MS}}@f$ scheme, where it is forbidden to across
     * a flavour threshould in the RG running from @f$\mu_i@f$ to @f$\mu@f$.
     * @param[in] mu a scale @f$\mu@f$ in GeV.
     * @param[in] alsi the initial value for the coupling at the scale given below.
     * @param[in] mu_i the initial scale @f$\mu_i@f$ in GeV.
     * @param[in] order LO, NLO or FULLNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the strong coupling constant @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme.
     */
    double AlsWithInit(const double mu, const double alsi, const double mu_i,
                       const orders order) const;
    
    /**
     * @brief Computes the running strong coupling @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme with the use of @f$\Lambda_{\rm QCD}@f$.
     * @param[in] mu A scale @f$\mu@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the strong coupling constant @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double AlsWithLambda(const double mu, const orders order) const;
    
    /**
     * @brief Computes the running strong coupling @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme. In the cases of LO, NLO and FULLNNLO,
     * the coupling is computed with AlsWithInit(). On the other hand, in the
     * cases of NNLO and FULLNNLO, the coupling is computed with AlsWithLambda().
     * @param[in] mu a scale @f$\mu@f$ in GeV.
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the strong coupling constant @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double Als(const double mu, const orders order = FULLNLO) const;
    
    /**
     * @brief Computes @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV
     */
    double logLambda(const double nf, orders order) const;
    
    /**
     * temporary function waiting for the implementation of NNLO etact
     * @param mu
     * @return
     */
    double Als4(const double mu) const;
    
    /**
     * temporary function waiting for the implementation of NNLO etact
     * @param mu_f
     * @param mu_i
     * @param m
     * @return
     */
    double Mrun4(const double mu_f, const double mu_i, const double m) const;
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @param[in] nf The number of active flavours.
     * @return The @f$\gamma_0@f$ coefficient.
     */
    double Gamma0(const double nf) const;
    
    /**
     * @param[in] nf The number of active flavours.
     * @return The @f$\gamma_1@f$ coefficient.
     */
    double Gamma1(const double nf) const;
    
    /**
     * @param[in] nf The number of active flavours.
     * @return The @f$\gamma_2@f$ coefficient.
     */
    double Gamma2(const double nf) const;
    
    /**
     * @brief Computes a running quark mass @f$m(\mu)@f$ from @f$m(m)@f$.
     * @param[in] mu a scale @f$\mu@f$ in GeV
     * @param[in] m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the running quark mass @f$m(\mu)@f$ in GeV
     */
    double Mrun(const double mu, const double m, const orders order = FULLNLO) const;
    
    /**
     * @brief Runs a quark mass from @f$\mu_i@f$ to @f$\mu_f@f$.
     * @param[in] mu_f a scale @f$\mu_f@f$ in GeV
     * @param[in] mu_i a scale @f$\mu_i@f$ in GeV
     * @param[in] m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(\mu_i)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the running quark mass @f$m(\mu_f)@f$ in GeV
     */
    double Mrun(const double mu_f, const double mu_i, const double m,
                const orders order = FULLNLO) const;
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Converts the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ to the pole mass
     * @param[in] mbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the pole mass in GeV
     */
    double Mbar2Mp(const double mbar, const orders order = FULLNLO) const;
    
    /**
     * @brief Converts a quark pole mass to the corresponding @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$.
     * @param[in] mp the pole mass of the bottom or top quark in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     */
    double Mp2Mbar(const double mp, const orders order = FULLNLO) const;
    
    /**
     * @brief Converts a quark mass from the @f$\overline{\mathrm{MS}}@f$ scheme to
     * the @f$\overline{\mathrm{DR}}@f$ scheme.
     * @param[in] MSscale the scale at which the @f$\overline{\mathrm{MS}}@f$ mass is defined
     * @param[in] MSbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @return the @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$ in GeV
     */
    double MS2DRqmass(const double MSscale, const double MSbar) const;
    
    /**
     * @brief Converts a quark mass from the @f$\overline{\mathrm{MS}}@f$ scheme to
     * the @f$\overline{\mathrm{DR}}@f$ scheme.
     * @param[in] MSbar The @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV.
     * @return The @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$ in GeV.
     */
    double MS2DRqmass(const double MSbar) const;
    
    ////////////////////////////////////////////////////////////////////////
    
protected:
    
    /**
     * @brief A method to set the value of a parameter of %QCD.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
    
    double Nc; ///< The number of colours.
    double CF; ///< The Casimir factor in the \f$SU(N_c)\f$ gauge theory.
    double mtpole;  ///< The pole mass of the top quark.
    Particle quarks[6]; ///< The vector of all SM quarks.
    Meson mesons[MESON_END]; ///< The vector of defined mesons.
    bool requireYu; ///< Switch for generating the Yukawa couplings to the up-type quarks.
    bool requireYd; ///< Switch for generating the Yukawa couplings to the down-type quarks.
    BParameter BBs; ///< The bag parameters for \f$\Delta b=2\f$ processes for the \f$B_s\f$ meson system.
    BParameter BBd; ///< The bag parameters for \f$\Delta b=2\f$ processes for the \f$B_d\f$ meson system.
    BParameter BD; ///< The bag parameters for \f$\Delta c=2\f$ processes for the \f$D^0\f$ meson system.
    BParameter BK; ///< The bag parameters for \f$\Delta s=2\f$ processes for the \f$K^0\f$ meson system.
    BParameter BKd1;
    BParameter BKd3;
    
    
    // model parameters
    double AlsMz; /**< The strong coupling constant at the Z-boson mass, \f$\alpha_s(M_Z)\f$. */
    double Mz; /**< The mass of the \f$Z\f$ boson in \f$GeV\f$ */
    double mut; /**< The threshold between six- and five-flavour theory in \f$GeV\f$. */
    double mub; /**< The threshold between five- and four-flavour theory in \f$GeV\f$. */
    double muc; /**< The threshold between four- and three-flavour theory in GeV. */
    double ReA0_Kd; /**< */
    double ReA2_Kd; /**< */
    double Omega_eta_etap; /**< */
    double Br_Kp_P0enu; /**< */
    double IB_Kl; /**< */
    double IB_Kp; /**< */
    double DeltaP_cu; /**< */
    double Br_Kp_munu; /**< */
    double Br_B_Xcenu; /**< */
    double BBsoBBd; /**< The ratio \f$ B_{B_s}/B_{B_d} \f$ necessary to compute \f$ B_{B_s} \f$. */
    double FBsoFBd; /**< The ratio \f$ F_{B_s}/F_{B_d} \f$ necessary to compute \f$ F_{B_s} \f$. */
    
private:
    
    /**
     * @brief
     * @param mu
     * @param logLambda
     * @param order
     * @return
     */
    double AlsWithLambda(const double mu, const double logLambda, const orders order) const;
    
    /**
     * @brief
     * @param logLambda6
     * @param logLambda5_in
     * @return
     */
    double ZeroNf6NLO(double *logLambda6, double *logLambda5_in) const;
    
    /**
     * @brief
     * @param logLambda5
     * @param order
     * @return
     */
    double ZeroNf5(double *logLambda5, double *order) const;
    
    /**
     * @brief
     * @param logLambda4
     * @param logLambda5_in
     * @return
     */
    double ZeroNf4NLO(double *logLambda4, double *logLambda5_in) const;
    
    /**
     * @brief
     * @param logLambda3
     * @param logLambda4_in
     * @return
     */
    double ZeroNf3NLO(double *logLambda3, double *logLambda4_in) const;
    
    /**
     * @brief
     * @param order
     * @return
     */
    double logLambda5(orders order) const;
    
    /**
     * @brief
     * @param nfNEW
     * @param nfORG
     * @param logLambdaORG
     * @return
     */
    double logLambdaNLO(const double nfNEW, const double nfORG, const double logLambdaORG) const;
    
    /**
     * @brief
     * @param muMatching
     * @param mf
     * @param nfNEW
     * @param nfORG
     * @param logLambdaORG
     * @param order
     * @return
     */
    double logLambda(const double muMatching, const double mf,
                     const double nfNEW, const double nfORG,
                     const double logLambdaORG, orders order) const;
    
    /**
     * @brief
     * @param nf_f
     * @param nf_i
     * @return
     */
    double threCorrForMass(const double nf_f, const double nf_i) const;
    
    /**
     * @brief
     * @param mu_f
     * @param mu_i
     * @param m
     * @param order
     * @return
     */
    double MrunTMP(const double mu_f, const double mu_i, const double m, const orders order) const;
    
    /**
     * @brief
     * @param mu
     * @param params
     * @return
     */
    double Mp2MbarTMP(double *mu, double *params) const;
    
    
    double zeta2; ///< \f$\zeta(2)\f$ computed from the <a href="http://www.gnu.org/software/gsl/" target=blank>gsl libraries</a>.
    double zeta3; ///< \f$\zeta(3)\f$ computed from the <a href="http://www.gnu.org/software/gsl/" target=blank>gsl libraries</a>.
    bool computeFBd; ///< Swith for computing \f$F_{B_d}\f$ from \f$F_{B_s}\f$.
    bool computeBd; ///< Swith for computing \f$B_{B_d}\f$ from \f$B_{B_s}\f$.
    bool computemt; ///< Switch for computing the \f$\overline{\mathrm{MS}}\f$ mass of the top quark.
    static const int CacheSize = 5; ///< Defines the depth of the cache.
    mutable double als_cache[8][CacheSize]; ///< Cache for \f$\alpha_s\f$.
    mutable double logLambda5_cache[4][CacheSize]; ///<
    mutable double logLambdaNLO_cache[9][CacheSize]; ///<
    mutable double mrun_cache[10][CacheSize]; ///< Cache for running quark mass.
    mutable double mp2mbar_cache[5][CacheSize]; ///< Cache for pole mass to msbar mass conversion.
    
    /**
     * @brief
     * @param cache
     * @param n
     */
    void CacheShift(double cache[][5], int n) const;
    
    
};

#endif	/* QCD_H */
