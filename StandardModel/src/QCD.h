/* 
 * File:   QCD.h
 * Author: marco
 *
 * Created on February 17, 2011, 2:13 PM
 */

#ifndef QCD_H
#define	QCD_H

#include "Model.h"
#include "Particle.h"
#include <Parameters.h>

class QCD: public Model {
public:
    enum mesons {B_D, B_S, B_P, K_0, K_P}; // update quark if changed!!!!!!
    enum quark {UP=mesons(K_P)+1,DOWN,CHARM,STRANGE,TOP,BOTTOM};
                           // update StandardModel::lepton if changed!!!!!!

    static const int NQCDvars = 18;
    /**
     * array containing the labels under which all QCD parameters must be
     * stored in a Parameters object
     */
    static const std::string QCDvars[NQCDvars];
//    /**
//     * default constructor
//     */
//    QCD();
//    /**
//     * constructor with default thresholds
//     * @param AlsM_i @f$\alpha_s(M)@f$
//     * @param M_i the scale @f$M@f$ at which @f$\alpha_s(M)@f$ is given
//     * @param mu_i up quark mass at 2 GeV
//     * @param md_i down quark mass at 2 GeV
//     * @param mc_i charm quark mass mc(mc)
//     * @param ms_i strange quark mass at 2 GeV
//     * @param mt_i top quark mass mt(mt)
//     * @param mb_i bottom quark mass mb(mb)
//     */
//    QCD(double AlsM_i, double M_i, double mu_i, double md_i, double ms_i,
//    double mc_i, double mb_i, double mt_i);
//    /**
//     * constructor with explicit thresholds
//     * @param AlsM_i @f$\alpha_s(M)@f$
//     * @param M_i the scale @f$M@f$ at which @f$\alpha_s(M)@f$ is given
//     * @param mu_i up quark mass at 2 GeV
//     * @param md_i down quark mass at 2 GeV
//     * @param mc_i charm quark mass mc(mc)
//     * @param ms_i strange quark mass at 2 GeV
//     * @param mt_i top quark mass mt(mt)
//     * @param mb_i bottom quark mass mb(mb)
//     * @param mu1_i threshold between six- and five-flavour theory
//     * @param mu2_i threshold between five- and four-flavour theory
//     * @param mu3_i threshold between four- and three-flavour theory
//     */
//    QCD(double AlsM_i, double M_i, double mu_i, double md_i, double ms_i,
//    double mc_i, double mb_i, double mt_i, double mu1_i, double mu2_i, double mu3_i);
//    /**
//     * copy constructor
//     * @param orig a QCD object
//     */
//    QCD(const QCD& orig);
    /**
     * Constructor for QCD
     * @param a Parameters object that must contain all the labels appearing in QCDvars
     */
    QCD(const Parameters&);
    virtual ~QCD();
    /**
     * the @f$\beta_0@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_0@f$ coefficient
     */
    double beta0(double nf) const;
    /**
     * the @f$\beta_1@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_1@f$ coefficient
     */
    double beta1(double nf) const;
    /**
     * the number of active flavour at scale @f$\mu@f$
     * @param mu the scale @f$\mu@f$ in GeV
     * @return the number of active flavour at scale @f$\mu@f$
     */
    double Nf(double mu) const;
    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param lam @f$\Lambda_\mathrm{QCD}@f$ with @f$n_f@f$ active flavours in GeV
     * @param nf the number of active flavours @f$n_f@f$
     * @param le order (=0 for LO, =1 for LO+NLO, =-1 for NLO)
     * @return @f$\alpha_s@f$
     */
    double als(double mu, double lam, double nf, int le) const;
    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param nf the number of active flavours
     * @param alsi the initial condition @f$\alpha_s(m_i)@f$
     * @param mi the scale @f$m_i@f$ in GeV
     * @param le order (=0 for LO, =1 for LO+NLO, =-1 for NLO)
     * @return @f$\alpha_s@f$
     */
    double als(double mu, double nf, double alsi, double mi, int le) const;
    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param nfmu the number of active flavours at the scale @f$\mu@f$
     * @param le order (=0 for LO, =1 for LO+NLO, =-1 for NLO)
     * @return @f$\alpha_s@f$
     */
    double als(double mu, double nfmu, int le) const;
    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param le order (=0 for LO, =1 for LO+NLO, =-1 for NLO)
     * @return @f$\alpha_s@f$
     */
    double als(double mu, int le) const;
    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @return @f$\alpha_s@f$ at the LO+NLO
     */
    double als(double mu) const;
    /**
     * @f$\Lambda_\mathrm{QCD}@f$ with four active flavours in GeV
     * @param le order (=0 for LO, =1 for LO+NLO)
     * @return @f$\Lambda_\mathrm{QCD}@f$
     */
    double lambda4(int le) const;

    /**
     *
     * @return @f$\alpha_s(M)@f$
     */
    double getAlsM() const {
        return AlsM;
    }

    /**
     * set the initial condition @f$\alpha_s(M)@f$
     * @param AlsM the initial condition @f$\alpha_s(M)@f$
     */
    void setAlsM(double AlsM) {
        this->AlsM = AlsM;
    }

    /**
     *
     * @return the scale M at which the initial condition for @f$\alpha_s(M)@f$ is given
     */
    double getM() const {
        return M;
    }

    /**
     * set the scale M at which the initial condition for @f$\alpha_s(M)@f$ is given
     * @param M the scale M in GeV
     */
    void setM(double M) {
        this->M = M;
    }

    /**
     *
     * @return the number of colours
     */
    double getNc() const {
        return Nc;
    }

    /**
     * set the number of colours
     * @param Nc the number of colours
     */
    void setNc(double Nc) {
        this->Nc = Nc;
    }

    /**
     *
     * @return the threshold between six- and five-flavour theory in GeV
     */
    double getMu1() const {
        return mu1;
    }

    /**
     * set the threshold between six- and five-flavour theory
     * @param mu1 the threshold between six- and five-flavour theory in GeV
     */
    void setMu1(double mu1) {
        this->mu1 = mu1;
    }

    /**
     *
     * @return the threshold between five- and four-flavour theory in GeV
     */
    double getMu2() const {
        return mu2;
    }

    /**
     * set the threshold between five- and four-flavour theory
     * @param mu2 the threshold between five- and four-flavour theory in GeV
     */
    void setMu2(double mu2) {
        this->mu2 = mu2;
    }

    /**
     *
     * @return the threshold between four- and three-flavour theory in GeV
     */
    double getMu3() const {
        return mu3;
    }

    /**
     *
     * @param p particle index
     * @return mass of the particle in GeV
     */
    virtual double getMass(int p) const {
        return(particles[p].getMass());
    }

    /**
     * set the mass of a particle
     * @param p particle index
     * @param m mass of the particle in GeV
     */
    virtual void setMass(int p, double m) {
        std::cout << sizeof(particles)/sizeof(particles[0]) << "  " << p << std::endl;
        particles[p].setMass(m);
    }

    /**
     * set the threshold between four- and three-flavour theory
     * @param mu3 the threshold between four- and three-flavour theory in GeV
     */
    void setMu3(double mu3) {
        this->mu3 = mu3;
    }

    /**
     * the running quark mass @f$m(\mu)@f$
     * @param mu the scale @f$\mu@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param nf the number of active flavours
     * @param le the order (=0 for LO, =1 for LO+NLO, =-1 for NLO)
     * @return the running quark mass @f$m(\mu)@f$
     */

    double mrun(double mu, double m, double nf, int le) const;
    /**
     * the running quark mass @f$m(\mu)@f$ at the LO+NLO
     * @param mu the scale @f$\mu@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param nf the number of active flavours
     * @return the running quark mass @f$m(\mu)@f$
     */
    double mrun(double mu, double m, double nf) const;
    /**
     * convert the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ to the pole mass
     * @param mbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @return the pole mass in GeV
     */
    double mbar2mp(double mbar) const;
    /**
     * convert the pole mass to the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param mp the pole mass in GeV
     * @return the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     */
    double mp2mbar(double mp) const;

    /**
     * updates the QCD parameters found in the argument
     * @param a Parameters object containing the parameters to be updated
     */
    void update(const Parameters&);

protected:
    double Nc, AlsM, M, mu1, mu2, mu3;
    Particle particles[quark(BOTTOM)+1];

//    void init(double AlsM_i, double M_i, double mu_i, double md_i, double ms_i,
//    double mc_i, double mb_i, double mt_i, double mu1_i, double mu2_i,
//        double mu3_i);
private:
    double thresholds(int i) const;
    double aboveth(double mu) const;
    double belowth(double mu) const;
    mutable double als_cache[5][5], lambda4_cache[2][5], mp2mbar_cache[4][5];


    double zero(double *x, double *) const;
    double mp2mbara(double * mu, double * mp) const;
    void CacheShift(double cache[][5], int n) const;
};

#endif	/* QCD_H */
