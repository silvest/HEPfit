#ifndef RGESolver_h
#define RGESolver_h
//RGESolver.h

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
//#include "IndependentIndices.h"

#include <unordered_map>
//#include <functionalz>
#include <boost/function.hpp>
//#include <boost/bind/bind.hpp>

#include "gslpp.h"

/** 
 * @brief A class that performs renormalization group evolution in the context of the SMEFT
 * @details The class solves the Renormalization Group Equations (RGEs) numerically. A faster, approximate 
 * solution that neglects the scale dependence of the anomalous dimension matrix is also available.
 * Only operators up to dimension six that preserve  lepton and baryon numbers 
 * are considered. The operator basis is the Warsaw basis, 
 * defined in https://arxiv.org/abs/1008.4884.
 * <tt>RGESolver</tt> splits real and imaginary part of each complex parameter. @n @n 
 * The numerical integration is performed with an adaptive step-size routine 
 * (the explicit embedded Runge-Kutta-Fehlberg method), using the
 * tools in the GNU Scientific Library. 
 * See https://www.gnu.org/software/gsl/doc/html/ode-initval.html for all the details. @n
 * The accuracy level of the numerical integration can be tuned selecting the parameters 
 * @f$\epsilon_{rel}@f$ and @f$\epsilon_{abs}@f$ using the
 *  dedicated setter functions. @n @n 
 * All the SMEFT coefficients are set using the \ref SetCoefficient methods and 
 * accessed with the \ref GetCoefficient methods. 
 * There exist three different signatures for each method, 
 * depending on the number of flavour indices of the 
 * parameter (0,2,4). @n 
 * These two routines must be used also for the SM parameters 
 *  \f$g_1,g_2,g_3,\lambda,m_h^2,\f$ 
 *   \f$\mathrm{Re}(\mathcal{Y}_u),\mathrm{Im}(\mathcal{Y}_u),\f$ 
 *   \f$\mathrm{Re}(\mathcal{Y}_d),\mathrm{Im}(\mathcal{Y}_d),\f$ 
 *   \f$\mathrm{Re}(\mathcal{Y}_e),\mathrm{Im}(\mathcal{Y}_e)\f$ 
 *  (we follow https://arxiv.org/abs/1308.2627 for what concerns
 * the conventions in the Higgs' sector). @n 
 * The routines \ref GetCKMAngle, \ref GetCKMPhase, \ref GetCKMRealPart, \ref GetCKMImagPart
 * should be used when interested in the CKM parameters or elements. The usage of this method is 
 * recommended after methods such \ref GenerateSMInitialConditions or \ref EvolveToBasis that choose 
 * a specific flavour basis ("UP" or "DOWN"), in which 
 * cases the CKM matrix is updated. 
 * A complete list of the keys that must be used to 
 * correctly invoke setter/getter methods are given in 
 * tables \ref SM, \ref 0F, \ref 2F and \ref 4F. @n
 * A summary of the operators symmetry classes is given in table \ref Sym. @n  
 * We follow http://www.utfit.org/UTfit/Formalism for what concerns the conventions for the CKM matrix.
 



 * @author S. Di Noi, L. Silvestrini.
 * @copyright GNU General Public License
 */


//tables with all the names of the coefficients.

/**
 * 
 * 
 * 

<table>
<caption id="SM"> Standard Model parameters. The labels in the left column 
must be used with the GetCoefficient/SetCoefficient methods, the ones 
in the right column must be used with GetCKMAngle methods.  </caption>
<tr style="vertical-align:top">
<td>
<table>
<tr> <th> Parameter     <th> Name  
<tr><td> \f$g_1\f$     <td> `g1`
<tr><td> \f$g_2\f$     <td> `g2`
<tr><td> \f$g_3\f$     <td> `g3`      
<tr><td> \f$\lambda\f$  <td> `lambda`      
<tr><td> \f$m_h^2\f$ \f$[\mathrm{GeV}^2]\f$  <td> `mh2`      
<tr><td> \f$\mathrm{Re}(\mathcal{Y}_u)\f$         <td> `YuR`      
<tr><td> \f$\mathrm{Im}(\mathcal{Y}_u)\f$         <td> `YuI`      
<tr><td> \f$\mathrm{Re}(\mathcal{Y}_d)\f$         <td> `YdR`      
<tr><td> \f$\mathrm{Im}(\mathcal{Y}_d)\f$         <td> `YdI`      
<tr><td> \f$\mathrm{Re}(\mathcal{Y}_e)\f$         <td> `YeR`      
<tr><td> \f$\mathrm{Im}(\mathcal{Y}_e)\f$         <td> `YeI`
</table> <td> 
<table>
<tr> <th> Parameter     <th> Name        
<tr><td> \f$ \sin(\theta_{12})\f$         <td> `s12`      
<tr><td> \f$ \sin(\theta_{13})\f$         <td> `s13`      
<tr><td> \f$ \sin(\theta_{23})\f$         <td> `s23`  
</table>
</table>

 *
 * 
 
<table>
<caption id="0F"> Scalar (and real) SMEFT operators. </caption>
<tr style="vertical-align:top">
<td>
<table>
<tr> <th colspan="2"> Classes 1-3
<tr> <th> Coefficient     <th> Name        
<tr><td>\f$C_{G}\f$         <td> `CG`      
<tr><td>\f$C_{\tilde{G}}\f$ <td> `CGtilde` 
<tr><td>\f$C_{W}\f$         <td> `CW`      
<tr><td>\f$C_{\tilde{W}}\f$ <td> `CWtilde` 
<tr><td>\f$C_H\f$           <td> `CH`      
<tr><td>\f$C_{H \Box} \f$   <td> `CHbox`   
<tr><td>\f$C_{HD}\f$        <td> `CHD`    
</table> <td> 
<table>
<tr><th colspan="2"> Class 4   
<tr><th> Coefficient     <th> Name           
<tr><td>\f$C_{HG}\f$         <td> `CHG`      
<tr><td>\f$C_{H\tilde{G}}\f$ <td> `CHGtilde` 
<tr><td>\f$C_{HW}\f$         <td> `CHW`      
<tr><td>\f$C_{H\tilde{W}}\f$ <td> `CHWtilde` 
<tr><td>\f$C_{HB}\f$         <td> `CHB`      
<tr><td>\f$C_{H\tilde{B}}\f$ <td> `CHBtilde` 
<tr><td>\f$C_{HWB}\f$         <td> `CHWB`      
<tr><td>\f$C_{H\tilde{W}B}\f$ <td> `CHWtildeB` 
</table>
</table>
  
 * 
 * 
 * 
 * 
 * 
<table>
<caption id="2F"> 2F SMEFT operators.  </caption>
<tr style="vertical-align:top"> <td>
<table>
<tr><th colspan="3"> Class 5   
<tr> <th> Coefficient     <th> Name   <th> Symmetry        
<tr><td>\f$\mathrm{Re}(C_{eH})\f$ <td> `CeHR` <td> WC1    
<tr><td>\f$\mathrm{Im}(C_{eH})\f$ <td> `CeHI`  <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{uH})\f$ <td> `CuHR`  <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{uH})\f$ <td> `CuHI`  <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{dH})\f$ <td> `CdHR`  <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{dH})\f$ <td> `CdHI`  <td> WC1
</table> <td> 
<table>
<tr><th colspan="3"> Class 6   
<tr> <th> Coefficient     <th> Name  <th> Symmetry   
<tr><td>\f$\mathrm{Re}(C_{eW})\f$ <td> `CeWR` <td> WC1  
<tr><td>\f$\mathrm{Im}(C_{eW})\f$ <td> `CeWI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{eB})\f$ <td> `CeBR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{eB})\f$ <td> `CeBI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{uG})\f$ <td> `CuGR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{uG})\f$ <td> `CuGI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{uW})\f$ <td> `CuWR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{uW})\f$ <td> `CuWI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{uB})\f$ <td> `CuBR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{uB})\f$ <td> `CuBI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{dG})\f$ <td> `CdGR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{dG})\f$ <td> `CdGI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{dW})\f$ <td> `CdWR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{dW})\f$ <td> `CdWI` <td> WC1 
<tr><td>\f$\mathrm{Re}(C_{dB})\f$ <td> `CdBR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{dB})\f$ <td> `CdBI` <td> WC1
</table> <td>  
<table> 
<tr> <th colspan="3"> Class 7   
<tr> <th> Coefficient     <th> Name <th> Symmetry  
<tr><td>\f$\mathrm{Re}(C_{Hl1})\f$ <td> `CHl1R` <td> WC2R  
<tr><td>\f$\mathrm{Im}(C_{Hl1})\f$ <td> `CHl1I` <td> WC2I 
<tr><td>\f$\mathrm{Re}(C_{Hl3})\f$ <td> `CHl3R` <td> WC2R 
<tr><td>\f$\mathrm{Im}(C_{Hl3})\f$ <td> `CHl3I` <td> WC2I 
<tr><td>\f$\mathrm{Re}(C_{He})\f$ <td> `CHeR` <td> WC2R 
<tr><td>\f$\mathrm{Im}(C_{He})\f$ <td> `CHeI` <td> WC2I 
<tr><td>\f$\mathrm{Re}(C_{Hq1})\f$ <td> `CHq1R` <td> WC2R 
<tr><td>\f$\mathrm{Im}(C_{Hq1})\f$ <td> `CHq1I` <td> WC2I 
<tr><td>\f$\mathrm{Re}(C_{Hq3})\f$ <td> `CHq3R` <td> WC2R  
<tr><td>\f$\mathrm{Im}(C_{Hq3})\f$ <td> `CHq3I` <td> WC2I 
<tr><td>\f$\mathrm{Re}(C_{Hu})\f$ <td> `CHuR` <td> WC2R 
<tr><td>\f$\mathrm{Im}(C_{Hu})\f$ <td> `CHuI` <td> WC2I 
<tr><td>\f$\mathrm{Re}(C_{Hd})\f$ <td> `CHdR` <td> WC2R 
<tr><td>\f$\mathrm{Im}(C_{Hd})\f$ <td> `CHdI` <td> WC2I
<tr><td>\f$\mathrm{Re}(C_{Hud})\f$ <td> `CHudR` <td> WC1 
<tr><td>\f$\mathrm{Im}(C_{Hud})\f$ <td> `CHudI` <td> WC1  
</table>
</table>

 * 
 * 
 * 
 *  
 
<table>
<caption id="4F"> 4F SMEFT Operators.</caption>
<tr style="vertical-align:top"><td>
<table>
<tr> <th colspan="3"> Class 8 \f$(\bar{L}L)(\bar{L}L)\f$  
<tr> <th> Coefficient     <th> Name <th> Symmetry  
<tr><td>\f$\mathrm{Re}(C_{ll})\f$ <td> `CllR` <td> WC6R 
<tr><td>\f$\mathrm{Im}(C_{ll})\f$ <td> `CllI` <td> WC6I 	
<tr><td>\f$\mathrm{Re}(C_{qq1})\f$ <td> `Cqq1R` <td> WC6R 
<tr><td>\f$\mathrm{Im}(C_{qq1})\f$ <td> `Cqq1I` <td> WC6I 	
<tr><td>\f$\mathrm{Re}(C_{qq3})\f$ <td> `Cqq3R` <td> WC6R 
<tr><td>\f$\mathrm{Im}(C_{qq3})\f$ <td> `Cqq3I` <td> WC6I 
<tr><td>\f$\mathrm{Re}(C_{lq1})\f$ <td> `Clq1R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{lq1})\f$ <td> `Clq1I` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{lq3})\f$ <td> `Clq3R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{lq3})\f$ <td> `Clq3I` <td> WC7I             
<tr><th colspan=3> Class 8 \f$(\bar{L}R)(\bar{L}R)\f$  
<tr><th> Coefficient     <th> Name <th> Symmetry  
<tr><td>\f$\mathrm{Re}(C_{quqd1})\f$ <td> `Cquqd1R` <td> WC5 
<tr><td>\f$\mathrm{Im}(C_{quqd1})\f$ <td> `Cquqd1I` <td> WC5 
<tr><td>\f$\mathrm{Re}(C_{quqd8})\f$ <td> `Cquqd8R` <td> WC5 
<tr><td>\f$\mathrm{Im}(C_{quqd8})\f$ <td> `Cquqs8I` <td> WC5 
<tr><td>\f$\mathrm{Re}(C_{lequ1})\f$ <td> `Clequ1R` <td> WC5 
<tr><td>\f$\mathrm{Im}(C_{lequ1})\f$ <td> `Clequ1I` <td> WC5 	
<tr><td>\f$\mathrm{Re}(C_{lequ3})\f$ <td> `Clequ3R` <td> WC5 
<tr><td>\f$\mathrm{Im}(C_{lequ3})\f$ <td> `Clequ3I` <td> WC5		
</table> <td>
<table>
<tr> <th colspan="3"> Class 8 \f$(\bar{R}R)(\bar{R}R)\f$  
<tr> <th> Coefficient     <th> Name <th> Symmetry  
<tr><td>\f$\mathrm{Re}(C_{ee})\f$ <td> `CeeR` <td> WC8R 
<tr><td>\f$\mathrm{Im}(C_{ee})\f$ <td> `CeeI` <td> WC8I 
<tr><td>\f$\mathrm{Re}(C_{uu})\f$ <td> `CuuR` <td> WC6R 
<tr><td>\f$\mathrm{Im}(C_{uu})\f$ <td> `CuuI` <td> WC6I 
<tr><td>\f$\mathrm{Re}(C_{dd})\f$ <td> `CddR` <td> WC6R 
<tr><td>\f$\mathrm{Im}(C_{dd})\f$ <td> `CddI` <td> WC6I 	
<tr><td>\f$\mathrm{Re}(C_{eu})\f$ <td> `CeuR` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{eu})\f$ <td> `CeuI` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{ed})\f$ <td> `CedR` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{ed})\f$ <td> `CedI` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{ud1})\f$ <td> `Cud1R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{ud1})\f$ <td> `Cud1I` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{ud8})\f$ <td> `Cud8R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{ud8})\f$ <td> `Cud8I` <td> WC7I             
<tr><th colspan="3"> Class 8 \f$(\bar{L}R)(\bar{R}L)\f$  
<tr><th> Coefficient     <th> Name <th> Symmetry  
<tr><td>\f$\mathrm{Re}(C_{ledq})\f$ <td> `CledqR` <td> WC5 
<tr><td>\f$\mathrm{Im}(C_{ledq})\f$ <td> `CledqI` <td> WC5 
</table> <td> 
<table> 
<tr> <th colspan="3" >Class 8 \f$(\bar{L}L)(\bar{R}R)\f$  
<tr> <th> Coefficient     <th> Name <th> Symmetry  
<tr><td>\f$\mathrm{Re}(C_{le})\f$ <td> `CleR` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{le})\f$ <td> `CleI` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{lu})\f$ <td> `CluR` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{lu})\f$ <td> `CluI` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{ld})\f$ <td> `CldR` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{ld})\f$ <td> `CldI` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{qe})\f$ <td> `CqeR` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{qe})\f$ <td> `CqeI` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{qu1})\f$ <td> `Cqu1R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{qu1})\f$ <td> `Cqu1I` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{qu8})\f$ <td> `Cqu8R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{qu8})\f$ <td> `Cqu8I` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{qd1})\f$ <td> `Cqd1R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{qd1})\f$ <td> `Cqd1I` <td> WC7I 
<tr><td>\f$\mathrm{Re}(C_{qd8})\f$ <td> `Cqd8R` <td> WC7R 
<tr><td>\f$\mathrm{Im}(C_{qd8})\f$ <td> `Cqd8I` <td> WC7I  
</table>
</table>
 * 

<table>
<caption id="Sym">Symmetry categories for operators in the SMEFT. nF indicates the number of
flavour indices for each category.
</caption>

<tr> <th> Parameter     <th> Name  
<tr><td> 0    <td> 0F scalar object
<tr><td> WC1    <td> 2F generic real matrix
<tr><td> WC2R   <td> 2F Hermitian matrix (real part)      
<tr><td> WC2I         <td> 2F Hermitian matrix (imaginary part)     
<tr><td> WC5         <td> 4F generic real object 
<tr><td> WC6R        <td> 4F two identical \f$ \bar{\psi} \psi \f$ currents (real part)    
<tr><td> WC6I        <td> 4F two identical \f$ \bar{\psi} \psi \f$ currents (imaginary part)   
<tr><td> WC7R       <td> 4F two independent \f$ \bar{\psi} \psi \f$ currents (real part)    
<tr><td> WC7I         <td>4F two independent \f$ \bar{\psi} \psi \f$ currents (imaginary part)   
<tr><td> WC8R        <td> \f$ \mathcal{C}_{ee}\f$ (real part)   
<tr><td> WC8I         <td> \f$ \mathcal{C}_{ee}\f$ (imaginary part)    
</table> 
 * 
 * 
 * 
 <table>
<caption id="SMInput">  SM parameters used by default to generate SM initial conditions at an arbitrary
scale. The scale at which these parameters are given is \f$ \mu = 173.65\f$ GeV. We follow http://www.utfit.org/UTfit/Formalism
for what concerns the conventions for the CKM matrix.  
 </caption>
<tr style="vertical-align:top">
<td>
<table>
<tr> <th> Parameter     <th> Value  
<tr>    <td> \f$g_1\f$ <td> 0.3573
<tr><td>    \f$g_2\f$       <td> 0.6511
<tr><td>\f$g_3\f$         <td> 1.161     
<tr><td>\f$\lambda\f$         <td> 0.1297      
<tr><td>\f$m_h^2\f$ \f$[\mathrm{GeV}^2]\f$  <td> 15650      
<tr><td>\f$\sin(\theta_{12})\f$         <td> 0.225     
<tr><td>\f$\sin(\theta_{13})\f$         <td> 0.042      
<tr><td>\f$\sin(\theta_{23})\f$         <td> 0.003675    
<tr><td>\f$\delta\f$ [rad]       <td> 1.1676     
</table> 
<td>
<table>
<tr><th> Parameter     <th> Value [GeV]   
<tr><td>\f$m_u\f$   <td> 0.0012      
<tr><td>\f$m_c\f$   <td> 0.640      
<tr><td>\f$m_t\f$   <td> 162.0      
<tr><td>\f$m_d\f$  <td> 0.0027      
<tr><td>\f$m_s\f$   <td> 0.052      
<tr><td>\f$m_b\f$   <td> 2.75     
<tr><td>\f$m_{e}\f$   <td> 0.000511     
<tr><td>\f$m_{\mu}\f$  <td> 0.1057    
<tr><td>\f$m_{\tau}\f$   <td> 1.776      
</table>
</table> 
 * 
 * 
 * 
 */

class RGESolver {
public:
    /**
     * @brief The default constructor.
     * @details It initializes to 0 all the SMEFT coefficients. 
     */
    RGESolver();

    /**
     * @brief The default destructor.
     */
    ~RGESolver();


    /**@name Parameters related to the numeric integration. */

    /**
     * @brief Getter for the relative error used in the numerical integration
     */
    double epsrel() {
        return epsrel_;
    }

    /**
     * @brief Getter for the absolute error used in the numerical integration
     */
    double epsabs() {
        return epsabs_;
    }

    /**
     * @brief Setter for the relative error used in the numerical integration
     * (default value = 0.005)
     */
    void Setepsrel(double epsrel) {
        epsrel_ = epsrel;
    }

    /**
     * @brief Setter for the absolute error used in the numerical integration 
     * (default value = e-13)
     */
    void Setepsabs(double epsabs) {
        epsabs_ = epsabs;
    }


    /** @name Evolution */
    /**
     * @brief Performs the RGE evolution
     * @details RGEs are solved with the chosen method from @p muI to @p muF.
     * Currently, the available methods are "Numeric" and "Approximate". @n 
     * The method takes as initial values the current values of the parameters, 
     * set with the \ref SetCoefficient functions. After completing the evolution
     * the values of the parameters are updated and are accessible with the
     * \ref GetCoefficient functions.
     * @param method solution method
     * @param muI initial energy scale (in GeV)
     * @param muF final energy scale (in GeV)  
     */
    void Evolve(std::string method, double muI, double muF);

    /**
     * @brief Performs the RGE evolution and the back rotation 
     * on the coefficients with flavour indices.
     * @details After the evolution, the CKM matrix is computed. A flavour rotation is performed 
     * on the coefficients to go in the chosen basis.
     * @param method solution method
     * @param muI initial energy scale (in GeV)
     * @param muF final energy scale (in GeV)  
     * @param basis flavour basis after the evolution 
     * ("UP" or "DOWN"). 
     */
    void EvolveToBasis(std::string method, double muI,
            double muF, std::string basis);


    /**
     * @brief Generates the initial conditions 
     * for Standard Model's parameters (gauge couplings,
     * Yukawa coupling, quartic coupling and Higgs' boson mass) at the scale 
     * <tt>mu</tt> (in GeV), using one-loop pure SM beta functions. Default low-energy
     * input is used. 
     * 
     * @details The initial conditions are generated at the scale <tt>mu</tt> starting from the values 
     * at \f$\mu = 173.65 \f$ GeV in table \ref SMInput. At the scale <tt>mu</tt> the CKM matrix is computed.
     * @param mu Scale (in GeV) at which the initial conditions 
     * are generated
     * @param basis Flavour basis (
     * <tt>"UP"</tt> or <tt>"DOWN"</tt>)
     * @param method Method used by \ref RGESolver
     * to run the SM parameters to the scale <tt>mu</tt> 
     * (<tt>"Numeric"</tt> or <tt>"Approximate"</tt>) 
     */
    void GenerateSMInitialConditions(double mu, std::string basis, std::string method);


    /** 
     * @brief Generates the initial conditions 
     * for Standard Model's parameters (gauge couplings,
     * Yukawa coupling, quartic coupling and Higgs' boson mass) at the scale 
     * <tt>mu</tt> (in GeV), using one-loop pure SM beta functions. User-defined low energy 
     * input is used.
     * 
     * @details The initial conditions are generated at the scale <tt>muFin</tt> starting from the 
     * inserted parameters at the scale <tt>muIn</tt>. This method should be used with usual fermion hierarchy 
     * (smallest mass for the 1st generation and greatest mass for the 3rd withoud mass degeneracy 
     * for all up and down quarks and for charged leptons). 
     * The generation of the initial conditions is performed only if all the masses are non-negative and 
     * if \f$ \sin \theta_{ij} \in (0,1) \f$, \f$\delta \in (\pi,\pi]\f$.
     * We follow http://www.utfit.org/UTfit for what concerns the conventions for the CKM matrix.
     * At the scale <tt>mu</tt> the CKM matrix is computed.
     * @param muIn Low-energy input scale (in GeV)
     * @param muFin  Scale (in GeV) at which the initial conditions 
     * are generated
     * @param basis Flavour basis (
     * <tt>"UP"</tt> or <tt>"DOWN"</tt>)
     * @param method Method used by \ref RGESolver
     * to run the SM parameters to the scale <tt>mu</tt> 
     * (<tt>"Numeric"</tt> or <tt>"Approximate"</tt>) 
     * @param g1in \f$g_1 \f$ 
     * @param g2in \f$g_2 \f$ 
     * @param g3in \f$g_3 \f$ 
     * @param lambdain \f$ \lambda \f$
     * @param mh2in \f$ m_h^2 \f$ (in GeV^2) 
     * @param Muin Array containing the masses of the up-type quarks in GeV in the order \f$(m_u,m_c,m_t)\f$ 
     * @param Mdin Array containing the masses of the down-type quarks in GeV in the order \f$(m_d,m_s,m_b)\f$ 
     * @param Mein Array containing the masses of the charged leptons in GeV in the order \f$(m_{e},m_{\mu},m_{\tau})\f$
     * @param s12in The sine of the CKM matrix angle \f$\sin \theta_{12} \f$ 
     * @param s13in The sine of the CKM matrix angle \f$\sin \theta_{13} \f$ 
     * @param s23in The sine of the CKM matrix angle \f$\sin \theta_{23} \f$ 
     * @param deltain The CKM matrix phase \f$ \delta \f$ 

     */
    void GenerateSMInitialConditions(double muIn, double muFin, std::string basis, std::string method,
            double g1in, double g2in, double g3in, double lambdain, double mh2in,
            double Muin[3], double Mdin[3], double Mein[3],
            double s12in, double s13in, double s23in, double deltain);


    /**
     * @brief Same as \ref Evolve, but only for the SM parameters. 
     * The user should use this method instead of \ref Evolve when 
     * interested in pure SM running. Using this function is the same of 
     * using \ref Evolve with all the SMEFT coefficients set to 0, but it is faster since 
     * it does compute only the evolution for the SM parameters.
     * @param method solution method
     * @param muI initial energy scale (in GeV)
     * @param muF final energy scale (in GeV)
     */
    void EvolveSMOnly(std::string method, double muI, double muF);


    /** @name Input/output   */






    /**
     * @brief Getter function for the CKM matrix angles 
     *  \f$\theta_{12},\theta_{13},\theta_{23}\f$. 
     * @details This method should be called only after methods that choose a specific
     *  flavour basis (as \ref GenerateSMInitialConditions or
     * \ref EvolveToBasis ), otherwise the CKM matrix is not updated.

     * @param name of the angle (see table \ref SM)
     * @return The selected CKM angle. 
     */
    double GetCKMAngle(std::string name);

    /**
     * @brief Getter function for the CKM matrix (real part)
     * @details This method should be called only after methods that choose a specific
     *  flavour basis (as \ref GenerateSMInitialConditions or
     * \ref EvolveToBasis ), otherwise the CKM matrix is not updated.
     * @return The real part of the selected CKM matrix element. 
     */
    double GetCKMRealPart(int i, int j) {
        return (CKM(i, j).real());
    };

    /**
     * @brief Getter function for the CKM matrix (imaginary part)
     * @details This method should be called only after methods that choose a specific
     *  flavour basis (as \ref GenerateSMInitialConditions or
     * \ref EvolveToBasis ), otherwise the CKM matrix is not updated.
     * @return The imaginary part of the selected CKM matrix element. 
     */
    double GetCKMImagPart(int i, int j) {
        return (CKM(i, j).imag());
    };


    /**
     * @brief Getter function for the CKM matrix phase \f$\delta\f$.
     * @details This method should be called only after methods that choose a specific
     *  flavour basis (as \ref GenerateSMInitialConditions or
     * \ref EvolveToBasis ), otherwise the CKM matrix is not updated.
     * @return  The CKM matrix phase \f$\delta\f$. 
     */
    double GetCKMPhase();



    //Setters for 0F,2F,4F

    /**
     * @brief Setter function for scalar/0F parameters (no flavour indices).
     * @param name name of the parameter (see table \ref 0F)
     * @param val its value
     */
    void SetCoefficient(std::string name, double val);
    /**
     * @brief Setter function for 2F parameters (2 flavour indices).
     * @details If at least 
     * one of the inserted indices is outside the [0:2] range, 
     * an error message is printed and no assignation is performed.
     * @param name name of the parameter (see table \ref 2F)
     * @param val its value
     * @param i first flavour index
     * @param j second flavour index
     */
    void SetCoefficient(std::string name, double val, int i, int j);
    /**
     * @brief Setter function for 4F parameters (4 flavour indices).
     * @details If at least 
     * one of the inserted indices is outside the [0:2] range, 
     * an error message is printed and no assignation is performed.
     * 
     * @param name name of the parameter (see table \ref 4F)
     * @param val its value
     * @param i first flavour index
     * @param j second flavour index
     * @param k third flavour index
     * @param l fourth flavour index
     */
    void SetCoefficient(std::string name, double val, int i, int j,
            int k, int l);



    //Getters for 0F,2F,4F
    /**
     * @brief Getter function for scalar/0F parameters (no flavour indices).
     * @param name name of the parameter (see table \ref 0F)
     * @return the requested parameter 
     */
    double GetCoefficient(std::string name);
    /**
     * @brief Getter function for 2F parameters (2 flavour indices).
     * @details If at least 
     * one of the inserted indices is outside the [0:2] range, 
     * an error message is printed and the value 0 is returned.
     * @param name name of the parameter (see table \ref 2F)
     * @param i first flavour index
     * @param j second flavour index
     * @return the requested parameter
     */
    double GetCoefficient(std::string name, int i, int j);

    /**
     * @brief Getter function for 4F parameters (4 flavour indices).
     * one of the inserted indices is outside the [0:2] range, 
     * an error message is printed and the value 0 is returned.
     * @param name name of the parameter (see table \ref 4F)
     * @param i first flavour index
     * @param j second flavour index
     * @param k third flavour index
     * @param l fourth flavour index
     * @return the requested parameter 
     */
    double GetCoefficient(std::string name, int i, int j,
            int k, int l);


    /**
     * @brief Resets all the SMEFT coefficients to 0 and the 
     * SM parameters to their default value. 
     * \f$\epsilon_{\textrm{abs}}\f$ and \f$\epsilon_{\textrm{rel}}\f$ are reset to 
     * their default value (in the UP basis).
     * @details 
     */
    void Reset();

    /**
     * @brief Saves the current values of parameters in a file
     * @details Currently, only "SLHA" format is implemented
     * @param filename Name of the output file 
     * @param format Format of the output file
     */
    void SaveOutputFile(std::string filename,
            std::string format);











private:




    /** @name Flavour */
    /**
     * @brief Goes into the chosen basis. No back-rotation
     * for SMEFT coefficients is performed. 
     * @param basis : allowed options are "UP","DOWN"
     */
    void GoToBasis(std::string basis);

    /**
     * @brief Same as \ref GoToBasis, but the rotation is performed 
     * only on the Yukawa matrices.
     * @param basis : allowed options are "UP","DOWN"
     */
    void GoToBasisSMOnly(std::string basis);


    /**
     * @brief Extracts from the CKM matrix the 4 
     * physical parameters. 
     */
    void ExtractParametersFromCKM();
    /**
     * @brief Starting from the current values of the 
     * masses, Yukawa matrices are generated in the chosen basis. 
     * @param basis
     */
    void FromMassesToYukawas(std::string basis);
    /**
     * @brief Computes the CKM matrix with the 
     * current values of the angles and phase
     */
    void UpdateCKM();
    /**
     * @brief Inserts the initial values of the SM parameters in the array @p x 
     * @details Only used in @p EvolveSMOnly
     */
    /** @name Private functions to perform the evolution*/

    void InitSMOnly();
    /**
     * @brief Saves the evolved values of the SM parameters from @p x  
     * @details Only used in @p EvolveSMOnly
     */
    void UpdateSMOnly();


    /**
     * @brief Sets all the SM parameters (and the SMInputScale)
     * at the default value and all the SMEFT coefficients to 0.
     */
    void SetSMDefaultInput();

    /**
     * @brief Inserts the initial values of the parameters in the array @p x 
     * @details Only used in @p Evolve
     */
    void Init();
    /**
     * @brief Saves the evolved values of the coefficients from @p x 
     * @details Only used in @p Evolve
     */

    void Update();

    /**
     * @brief Computes the beta functions for the SMEFT.
     * @param logmu value of the logarithm of the energy scale at which the beta functions are computed 
     * @param y 1D array in which are stored the current values of the parameters
     * @param f 1D array in which the beta functions for each parameters are saved
     * @param params eventual additional parameters (not used)
     * @return @p GSL_SUCCESS
     */
    static int func(double logmu, const double y[],
            double f[], void* params);

    /**
     * @brief Computes the beta functions for the SM only.
     * @param logmu value of the logarithm of the energy scale at which the beta functions are computed 
     * @param y 1D array in which are stored the current values of the parameters
     * @param f 1D array in which the beta functions for each parameters are saved
     * @param params eventual additional parameters (not used)
     * @return @p GSL_SUCCESS
     */
    static int funcSMOnly(double logmu, const double y[],
            double f[], void* params);


    /**@name GSL Objects */
    ///@{ 

    /**
     * @brief Relative error used in the integrator
     * 
     * ().   */
    double epsrel_; // = ;

    /**
     * @brief Absolute error used in the integrator
     */
    double epsabs_;

    /**
     * @brief Resets epsabs to its default value
     */
    void Resetepsabs() {
        epsabs_ = 0.0000000000001;
    };

    /**
     * @brief Resets epsrel to its default value
     */
    void Resetepsrel() {
        epsrel_ = 0.005;
    };



    gsl_odeiv2_system sys = {func, NULL, 2558, NULL};
    /*gsl_odeiv2_driver * d =
            gsl_odeiv2_driver_alloc_y_new(&sys,
            gsl_odeiv2_step_rkf45,
            0.1, epsrel_, epsabs_);*/
    gsl_odeiv2_step * s
            = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 2558);
    gsl_odeiv2_evolve * e
            = gsl_odeiv2_evolve_alloc(2558);

    gsl_odeiv2_system sysSMOnly = {funcSMOnly, NULL, 59, NULL};

    gsl_odeiv2_step * sSMOnly
            = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 59);
    gsl_odeiv2_evolve * eSMOnly
            = gsl_odeiv2_evolve_alloc(59);




    /*
    gsl_odeiv2_system sys_ = {func, NULL, 3};
    gsl_odeiv2_step * s_ = gsl_odeiv2_step_alloc(
            gsl_odeiv2_step_rkf45, 2558);
    gsl_odeiv2_control * con_ = gsl_odeiv2_control_standard_new(
            epsabs_, epsrel_, 1, 1);
    gsl_odeiv2_evolve* evo_ = gsl_odeiv2_evolve_alloc(2558);


    gsl_odeiv2_system sysSMOnly_ = {funcSMOnly, NULL, 3};
    gsl_odeiv2_step * sSMOnly_ = gsl_odeiv2_step_alloc(
            gsl_odeiv2_step_rkf45, 59);
    gsl_odeiv2_control * conSMOnly_ = gsl_odeiv2_control_standard_new(
            epsabs_, epsrel_, 1, 1);
    gsl_odeiv2_evolve* evoSMOnly_ = gsl_odeiv2_evolve_alloc(59);*/


    /** @brief 1D array for the integration */
    double x[2558];

    /// @}



    //-----------------------------------------------------------------------------------------

    /** @name Independent entries
     * We follow https://arxiv.org/abs/2010.16341 tab. 15, 16 for the symmetry class
     *  of the operators. @n 
     * Notice that operators in classes WC1 and WC5 have no restrictions
     *  for neither real nor imaginary part, each having @f$N_G ^2@f$ (for WC1) 
     * or @f$N_G ^4@f$ (for WC5).  */

    ///@{

    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC2 */
    static const int DWC2R = 6;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC2 */
    static const int DWC2I = 3;
    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC6 */
    static const int DWC6R = 27;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC6 */
    static const int DWC6I = 18;
    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC7 */
    static const int DWC7R = 45;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC7 */
    static const int DWC7I = 36;
    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC8 */
    static const int DWC8R = 21;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC8 */
    static const int DWC8I = 15;

    ///@}


    /** @name Indices chosen as independent entries
     * the element <tt>WCn_indices[a][b]</tt> must be interpretated as :
     * the b-th index (there are 4 in 4F operators and 2 in 2F) of the a-th 
     * independent entry for the WCn category, where n = 1,2R,2I...
     */
    ///@{
    /** @brief Independent indices for WC2R */
    static const int WC2R_indices[DWC2R][2];
    /** @brief Independent indices for WC2I */
    static const int WC2I_indices[DWC2I][2];
    /** @brief Independent indices for WC6R */
    static const int WC6R_indices[DWC6R][4];
    /** @brief Independent indices for WC6I */
    static const int WC6I_indices[DWC6I][4];
    /** @brief Independent indices for WC7R */
    static const int WC7R_indices[DWC7R][4];
    /** @brief Independent indices for WC7I */
    static const int WC7I_indices[DWC7I][4];
    /** @brief Independent indices for WC8R */
    static const int WC8R_indices[DWC8R][4];
    /** @brief Independent indices for WC8I */
    static const int WC8I_indices[DWC8I][4];


    ///@}







    /**@name Number of operators for each class
     * See https://arxiv.org/abs/1308.2627 tab. 1 for the full list of operators. 
     * There are 8 different classes, depending on the field content. This classification 
     * is different from the symmetry classification WC1, WC2R, ... @n
     * For each class  @p N is the number of operators and 
     * @p E the number of independent entries. When E is not explicitely
     * defined is understood  <tt> E=N </tt> (no flavour structure)
     */

    ///@{


    /**
     * @brief Number of fermion flavours 
     */
    static const int NG = 3;
    /**
     * @brief Dimension of matrices in flavour space
     */
    static const int DF = 9;
    /**
     * @brief Independent entries of a @f$N_G \times N_G@f$ real symmetric matrix
     */
    static const int DFs = (NG*NG + NG) / 2;
    /**
     * @brief Independent entries of a @f$N_G \times N_G@f$ real anti-symmetric matrix 
     */
    static const int DFa = (NG*NG - NG) / 2;

    /** @brief Number of gauge couplings 
     * @ingroup Miscellaneous*/
    static const int Ngauge = 3;
    /** @brief Number of Higgs' sector parameters
     * @ingroup Miscellaneous */
    static const int Nh = 2;
    /** @brief Number of Yukawa matrices 
     * @ingroup Miscellaneous*/
    static const int Nyukawa = 3;
    /** @brief Number of real parameters for each Yukawa matrix
     *  */
    static const int Eyuk = (Nyukawa * 2 * DF);
    static const int N1 = 4;
    static const int N23 = 3;
    static const int N4 = 8;

    static const int N5 = 3;
    static const int E5 = (N5 * 2 * DF);

    static const int N6 = 8;
    static const int E6 = (N6 * 2 * DF);

    static const int N7 = 8;
    /** @brief Number of Hermitian operators in class 7*/
    static const int N7H = 7;
    /** @brief Number of non-Hermitian operators in class 7*/
    static const int N7NH = 1;
    static const int E7 = (N7H*(DWC2R + DWC2I) + N7NH * 2 * DF);

    static const int N8_LLLL = 5;
    static const int E8_LLLL = 2 * (DWC7R + DWC7I) + 3 * (DWC6R + DWC6I);

    static const int N8_RRRR = 7;
    static const int E8_RRRR = 4 * (DWC7R + DWC7I) + 2 * (DWC6R + DWC6I) + 1 * (DWC8R + DWC8I);

    static const int N8_LLRR = 8;
    static const int E8_LLRR = 8 * (DWC7R + DWC7I);

    static const int N8_LRRL = 1;
    static const int E8_LRRL = 2 * NG*NG*NG*NG*N8_LRRL;
    static const int N8_LRLR = 4;
    static const int E8_LRLR = 2 * NG*NG*NG*NG*N8_LRLR;

    ///@}



    /** @name  Fractions 
     * Recurring fractions defined in order to increase efficiency 
     */

    ///@{
    static const double TWO_THIRDS;
    static const double FOUR_THIRDS;
    static const double EIGHT_THIRDS;
    static const double ONE_THIRD;
    static const double ONE_SIXTH;
    static const double TEN_THIRDS;
    ///@}


    /** @name Miscellaneous parameters
     */

    ///@{
    /**
     * @brief Kroenecker delta in flavour space 
     */
    static const double delta[3][3];


    /** @brief Dimension of the system */
    static const int dim = (Ngauge + Nh + Eyuk + N1 + N23 + N4 + E5 + E6 + E7 + E8_LLLL + E8_RRRR + E8_LLRR + E8_LRRL + E8_LRLR);

    /**@brief Number of colors */
    static const double NC;
    /**@brief Number of colors squared */
    static const double NC2;

    /** @brief Leading-order @f$g_1@f$  beta function
     * ( with @f$g_1@f$ normalized as usual and not as in GUT theories)*/
    static const double b01;
    /** @brief Leading-order @f$g_2@f$  beta function */
    static const double b02;
    /** @brief Leading-order @f$g_3@f$  beta function */
    static const double b03;

    //Casimirs
    /** @brief @f$\mathbf{SU(2)}@f$ adjoint Casimir */
    static const double cA2;
    /** @brief @f$\mathbf{SU(3)}@f$ adjoint Casimir */
    static const double cA3;
    /** @brief @f$\mathbf{SU(2)}@f$ fundamental Casimir */
    static const double cF2;
    /** @brief @f$\mathbf{SU(3)}@f$ fundamental Casimir */
    static const double cF3;
    ///@}

    //Hypercharges (and product of hypercharges)
    //

    /** @name Hypercharges (and products of hypercharges)
     * We use @f$Q = T_3 + Y @f$, being @f$T_3@f$ the z-component of the weak isospin 
     * and @f$Q@f$ the electric charge. The other possibility 
     * is @f$Q = T_3 + \frac{Y}{2} @f$ as in https://arxiv.org/abs/1308.2627 
     */
    ///@{
    static const double Yh;
    static const double Yh2;

    static const double Yq;
    static const double Yq2;
    static const double Yl;
    static const double Yl2;

    static const double Yu;
    static const double Yu2;
    static const double Yd;
    static const double Yd2;
    static const double Ye;
    static const double Ye2;

    static const double YhYu;
    static const double YhYd;
    static const double YhYe;
    static const double YhYq;
    static const double YhYl;

    static const double YuYd;
    static const double YuYe;
    static const double YuYq;
    static const double YuYl;

    static const double YdYe;
    static const double YdYq;
    static const double YdYl;

    static const double YeYq;
    static const double YeYl;

    static const double YlYq;
    ///@}




    /** @name Standard Model parameters 
     */

    ///@{
    /** @brief @f$g_1@f$  */
    double g1, /** @brief @f$g_2@f$  */ g2, /**@brief @f$g_3@f$  */ g3,
    /** @brief @f$m_h ^2 @f$ (Higgs boson mass squared) 
      @details See https://arxiv.org/pdf/1308.2627.pdf for the normalization */ mh2,
    /** @brief @f$ \lambda @f$ (Higgs quartic coupling)
     *  @details See https://arxiv.org/pdf/1308.2627.pdf for the normalization  */ lambda;
    gslpp::matrix<double> yuR, yuI, ydR, ydI, yeR, yeI;
    ///@}
    /**
     * @brief The CKM matrix
     */
    gslpp::matrix<gslpp::complex> CKM = gslpp::matrix<gslpp::complex>(3, 3, 0.);

    /**
     * @brief the scale at which the method 
     * \ref GenerateSMInitialConditions
     * for SM parameters.
     */
    double InputScale_SM; //GeV

    //CKM parameters
    /**
     * @brief \f$\delta\f$  
     */
    double CKM_delta;

    /**
     * @brief \f$\cos \theta_{12}\f$  
     */
    double c12;
    /**
     * @brief \f$\sin \theta_{12}\f$  
     */
    double s12;
    /**
     * @brief \f$\cos \theta_{13}\f$  
     */
    double c13;
    /**
     * @brief \f$\sin \theta_{13}\f$  
     */
    double s13;
    /**
     * @brief \f$\cos \theta_{23}\f$  
     */
    double c23;
    /**
     * @brief \f$\sin \theta_{23}\f$  
     */
    double s23;

    //masses
    /**
     *  @brief \f$m_u\f$ (GeV)
     * */
    double mu;
    /**
     *  @brief \f$m_c\f$ (GeV)
     * */
    double mc;/**
     *  @brief \f$m_t\f$ (GeV)
     * */
    double mt;
    
    /**
     *  @brief \f$m_d\f$ (GeV)
     * */
    double md;
    /**
     *  @brief \f$m_s\f$ (GeV)
     * */
    double ms;/**
     *  @brief \f$m_b\f$ (GeV)
     * */
    double mb;
    
    /**
     *  @brief \f$m_e\f$ (GeV)
     * */
    double mel;
    /**
     *  @brief \f$m_{\mu}\f$ (GeV)
     * */
    double mmu;/**
     *  @brief \f$m_{\tau}\f$ (GeV)
     * */
    double mtau;
    
            /**@name SMEFT dimension-six operators 
             * By default, all SMEFT dimension-six operators' coefficients are set to 0. 
             * See https://arxiv.org/abs/1308.2627 tab. 1 for the full list of operators. 
             * Each member has a class from 1 to 8 depending on its field contents, as well 
             * as a flavour symmetry classification (WC1, WC2R, WC2I...).  
             */

            ///@{

            /**  @brief @f$ C_G@f$ (class 1, scalar) */ double cG = 0.;
    /**  @brief @f$ C_{\tilde{G}}@f$ (class 1, scalar)*/double cGT = 0.;
    /**  @brief @f$ C_W@f$ (class 1, scalar)*/ double cW = 0.;
    /**  @brief @f$ C_{\tilde{W}} @f$ (class 1, scalar)*/double cWT = 0.;

    /** @brief @f$C_H@f$ (class 2, scalar)*/double cH = 0.;

    /** @brief @f$C_{H \Box} @f$ (class 3, scalar)*/ double cHBOX = 0.;
    /** @brief @f$C_{H D} @f$ (class 3, scalar)*/ double cHD = 0.;

    /** @brief @f$C_{H G} @f$ (class 4, scalar)*/ double cHG = 0.;
    /** @brief @f$C_{H \tilde{G}} @f$ (class 4, scalar)*/ double cHGT = 0.;
    /** @brief @f$C_{H W} @f$ (class 4, scalar)*/ double cHW = 0.;
    /** @brief @f$C_{H \tilde{W}} @f$ (class 4, scalar)*/ double cHWT = 0.;
    /** @brief @f$C_{H B} @f$ (class 4, scalar)*/ double cHB = 0.;
    /** @brief @f$C_{H \tilde{B}} @f$ (class 4, scalar)*/ double cHBT = 0.;
    /** @brief @f$C_{H WB} @f$ (class 4, scalar)*/ double cHWB = 0.;
    /** @brief @f$C_{H \tilde{W} B} @f$ (class 4, scalar)*/double cHWBT = 0.;



    /** @brief @f$ \mathrm{Re} \left[ C_{eH } \right]@f$ (class 5, WC1) */
    double ceHR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{eH } \right]@f$ (class 5, WC1) */
    double ceHI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{uH } \right]@f$ (class 5, WC1) */
    double cuHR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{uH } \right]@f$ (class 5, WC1) */
    double cuHI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{dH } \right]@f$ (class 5, WC1) */
    double cdHR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{dH } \right]@f$ (class 5, WC1) */
    double cdHI[3 * 3] = {0.};


    /** @brief @f$ \mathrm{Re} \left[ C_{eW } \right]@f$ (class 6, WC1) */ double ceWR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{eW } \right]@f$ (class 6, WC1) */ double ceWI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{eB } \right]@f$ (class 6, WC1) */ double ceBR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{eB } \right]@f$ (class 6, WC1) */ double ceBI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{uG } \right]@f$ (class 6, WC1) */double cuGR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{uG } \right]@f$ (class 6, WC1) */double cuGI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{uW } \right]@f$ (class 6, WC1) */double cuWR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{uW } \right]@f$ (class 6, WC1) */double cuWI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{uB } \right]@f$ (class 6, WC1) */double cuBR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{uB } \right]@f$ (class 6, WC1) */double cuBI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{dG } \right]@f$ (class 6, WC1) */double cdGR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{dG } \right]@f$ (class 6, WC1) */double cdGI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{dW } \right]@f$ (class 6, WC1) */double cdWR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{dW } \right]@f$ (class 6, WC1) */double cdWI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{dB } \right]@f$ (class 6, WC1) */double cdBR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{dW } \right]@f$ (class 6, WC1) */double cdBI[3 * 3] = {0.};




    /** @brief @f$ \mathrm{Re} \left[ C_{Hl1} \right]@f$ (class 7, WC2R) */ double cHl1R[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hl1} \right]@f$ (class 7, WC2I) */ double cHl1I[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{Hl3} \right]@f$ (class 7, WC2R) */double cHl3R[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hl3} \right]@f$ (class 7, WC2I) */double cHl3I[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{He} \right]@f$ (class 7, WC2R) */double cHeR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{He} \right]@f$ (class 7, WC2I) */double cHeI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{Hq1} \right]@f$ (class 7, WC2R) */double cHq1R[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hq1} \right]@f$ (class 7, WC2I) */double cHq1I[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{Hq3} \right]@f$ (class 7, WC2R) */double cHq3R[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hq3} \right]@f$ (class 7, WC2I) */double cHq3I[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{Hu} \right]@f$ (class 7, WC2R) */double cHuR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hu} \right]@f$ (class 7, WC2I) */double cHuI[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{Hd} \right]@f$ (class 7, WC2R) */double cHdR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hd} \right]@f$ (class 7, WC2I) */double cHdI[3 * 3] = {0.};

    /** @brief @f$ \mathrm{Re} \left[ C_{Hud} \right]@f$ (class 7, WC1) */double cHudR[3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{Hud} \right]@f$ (class 7, WC1) */double cHudI[3 * 3] = {0.};




    /** @brief @f$ \mathrm{Re} \left[ C_{ll} \right]@f$ (class 8-[LL][LL], WC6R) */double cllR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ll} \right]@f$ (class 8-[LL][LL], WC6I) */double cllI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qq1} \right]@f$ (class 8-[LL][LL], WC6R) */double cqq1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qq1} \right]@f$ (class 8-[LL][LL], WC6I) */double cqq1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qq3} \right]@f$ (class 8-[LL][LL], WC6R) */double cqq3R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qq3} \right]@f$ (class 8-[LL][LL], WC6I) */double cqq3I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{lq1} \right]@f$ (class 8-[LL][LL], WC7R) */double clq1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{lq1} \right]@f$ (class 8-[LL][LL], WC7I) */double clq1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{lq3} \right]@f$ (class 8-[LL][LL], WC7R) */double clq3R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{lq3} \right]@f$ (class 8-[LL][LL], WC7I) */double clq3I[3 * 3 * 3 * 3] = {0.};


    /** @brief @f$ \mathrm{Re} \left[ C_{uu} \right]@f$ (class 8-[RR][RR], WC6R) */double cuuR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{uu} \right]@f$ (class 8-[RR][RR], WC6I) */double cuuI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{dd} \right]@f$ (class 8-[RR][RR], WC6R) */double cddR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{dd} \right]@f$ (class 8-[RR][RR], WC6I) */double cddI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{ee} \right]@f$ (class 8-[RR][RR], WC8R) */double ceeR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ee} \right]@f$ (class 8-[RR][RR], WC8I) */double ceeI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{eu} \right]@f$ (class 8-[RR][RR], WC7R) */double ceuR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{eu} \right]@f$ (class 8-[RR][RR], WC7I) */double ceuI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{ed} \right]@f$ (class 8-[RR][RR], WC7R) */double cedR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ed} \right]@f$ (class 8-[RR][RR], WC7I) */ double cedI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{ud1} \right]@f$ (class 8-[RR][RR], WC7R) */double cud1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ud1} \right]@f$ (class 8-[RR][RR], WC7I) */double cud1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{ud8} \right]@f$ (class 8-[RR][RR], WC7R) */double cud8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ud1} \right]@f$ (class 8-[RR][RR], WC7I) */double cud8I[3 * 3 * 3 * 3] = {0.};



    /** @brief @f$ \mathrm{Re} \left[ C_{le} \right]@f$ (class 8-[LL][RR], WC7R) */double cleR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{le} \right]@f$ (class 8-[LL][RR], WC7I) */double cleI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{lu} \right]@f$ (class 8-[LL][RR], WC7R) */double cluR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{lu} \right]@f$ (class 8-[LL][RR], WC7I) */double cluI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{ld} \right]@f$ (class 8-[LL][RR], WC7R) */double cldR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ld} \right]@f$ (class 8-[LL][RR], WC7I) */double cldI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qe} \right]@f$ (class 8-[LL][RR], WC7R) */double cqeR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qe} \right]@f$ (class 8-[LL][RR], WC7I) */double cqeI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qu1} \right]@f$ (class 8-[LL][RR], WC7R) */double cqu1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qu1} \right]@f$ (class 8-[LL][RR], WC7I) */double cqu1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qu8} \right]@f$ (class 8-[LL][RR], WC7R) */double cqu8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qu8} \right]@f$ (class 8-[LL][RR], WC7I) */double cqu8I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qd1} \right]@f$ (class 8-[LL][RR], WC7R) */double cqd1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qd1} \right]@f$ (class 8-[LL][RR], WC7I) */ double cqd1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{qd8} \right]@f$ (class 8-[LL][RR], WC7R) */double cqd8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{qd8} \right]@f$ (class 8-[LL][RR], WC7I) */ double cqd8I[3 * 3 * 3 * 3] = {0.};




    /** @brief @f$ \mathrm{Re} \left[ C_{ledq} \right]@f$ (class 8-[LR][RL], WC5) */ double cledqR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{ledq} \right]@f$ (class 8-[LR][RL], WC5) */double cledqI[3 * 3 * 3 * 3] = {0.};



    /** @brief @f$ \mathrm{Re} \left[ C_{lequ1} \right]@f$ (class 8-[LR][LR], WC5) */double clequ1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{lequ1} \right]@f$ (class 8-[LR][LR], WC5) */double clequ1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{lequ3} \right]@f$ (class 8-[LR][LR], WC5) */double clequ3R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{lequ3} \right]@f$ (class 8-[LR][LR], WC5) */double clequ3I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{quqd1} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{quqd1} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Re} \left[ C_{quqd8} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \mathrm{Im} \left[ C_{quqd8} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd8I[3 * 3 * 3 * 3] = {0.};
    ///@}

    //-----------------------------------------------------------------------------






    /** @name Maps for I/O
     * Maps that connects the coefficients with the appropriate getter/setter functions
     * (based on their symmetry properties)
     * 
     */

    ///@{
    std::unordered_map<std::string, double*> CKMAngles;
    std::unordered_map<std::string, double*> Operators0F;
    std::unordered_map<std::string, boost::function<void(int, int, double) >> Setter2F;
    std::unordered_map<std::string, boost::function<double(int, int) >> Getter2F;
    std::unordered_map<std::string, boost::function<void(int, int, int, int, double) >> Setter4F;
    std::unordered_map<std::string, boost::function<double(int, int, int, int) >> Getter4F;
    ///@}


    /** @name Output file
     * @brief Stuff for output on file
     */
    ///@{
    /**
     * @brief Prints on file the coefficient @p c
     * 
     * @details This function is used only in the function
     * \ref SaveOutputFile. 
     * Currently only "SLHA" printing for WC1,WC2R/I,
     * WC5,WC6R/I,WC7R/I,WC8R/I is implemented. 
     * 
     * @param c coefficient
     * @param name printed name 
     * @param sym symmetry category of the operator
     * @param format chosen format 
     * @param f pointer to file
     */

    static inline void Print(double* c, std::string name,
            std::string sym, std::string format,
            std::ofstream& f);
    ///@}



    /** @name Setters and Getters 
     * Setter and getter function for each symmetry class. 
     * 
     */


    ///@{

    static inline void Yukawa_set(gslpp::matrix<double> *y, int i, int j, double val);
    static inline double Yukawa(gslpp::matrix<double> *y, int i, int j);

    static inline void WC1_set(double * c, int i, int j, double val);
    static inline double WC1(double * c, int i, int j);

    static inline double WC2R(double * c, int i, int j);
    static inline void WC2R_set(double * c, int i, int j, double val);
    static inline double WC2I(double * c, int i, int j);
    static inline void WC2I_set(double * c, int i, int j, double val);

    static inline double WC3(double * c, int i, int j);
    static inline void WC3_set(double * c, int i, int j, double val);

    static inline void WC5_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC5(double * c, int i, int j, int k, int l);

    static inline double WC6R(double * c, int i, int j, int k, int l);
    static inline void WC6R_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC6I(double * c, int i, int j, int k, int l);
    static inline void WC6I_set(double * c, int i, int j, int k, int l, double val);

    static inline void WC7R_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC7R(double * c, int i, int j, int k, int l);
    static inline double WC7I(double * c, int i, int j, int k, int l);
    static inline void WC7I_set(double * c, int i, int j, int k, int l, double val);

    static inline double WC8R(double * c, int i, int j, int k, int l);
    static inline void WC8R_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC8I(double * c, int i, int j, int k, int l);
    static inline void WC8I_set(double * c, int i, int j, int k, int l, double val);
    ///@}

};


#endif
