/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSTUVWXY.h"
#include <stdexcept>


const std::string NPSTUVWXY::STUVWXYvars[NSTUVWXYvars]
        = {"obliqueShat", "obliqueThat", "obliqueUhat",
    "obliqueV", "obliqueW", "obliqueX", "obliqueY"};

NPSTUVWXY::NPSTUVWXY()
: NPbase(), myLEP2oblique(trueSM)
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueShat", boost::cref(myObliqueShat)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueThat", boost::cref(myObliqueThat)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueUhat", boost::cref(myObliqueUhat)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueV", boost::cref(myObliqueV)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueW", boost::cref(myObliqueW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueX", boost::cref(myObliqueX)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueY", boost::cref(myObliqueY)));

}

void NPSTUVWXY::setParameter(const std::string name, const double& value)
{
    if (name.compare("obliqueShat") == 0)
        myObliqueShat = value;
    else if (name.compare("obliqueThat") == 0)
        myObliqueThat = value;
    else if (name.compare("obliqueUhat") == 0)
        myObliqueUhat = value;
    else if (name.compare("obliqueV") == 0)
        myObliqueV = value;
    else if (name.compare("obliqueW") == 0)
        myObliqueW = value;
    else if (name.compare("obliqueX") == 0)
        myObliqueX = value;
    else if (name.compare("obliqueY") == 0)
        myObliqueY = value;
    else
        NPbase::setParameter(name, value);
}

bool NPSTUVWXY::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NSTUVWXYvars; i++) {
        if (DPars.find(STUVWXYvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTUVWXY parameter "
                    << STUVWXYvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::epsilon1() const
{
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();
    double eps1 = trueSM.epsilon1();
    eps1 += myObliqueThat - myObliqueW + 2.0 * sqrt(s2) / sqrt(c2) * myObliqueX
            - s2 / c2*myObliqueY;
    return eps1;
}

double NPSTUVWXY::epsilon2() const
{
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();
    double eps2 = trueSM.epsilon2();
    eps2 += myObliqueUhat - myObliqueV - myObliqueW + 2.0 * sqrt(s2) / sqrt(c2) * myObliqueX;
    return eps2;
}

double NPSTUVWXY::epsilon3() const
{
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();
    double eps3 = trueSM.epsilon3();
    eps3 += myObliqueShat - myObliqueW + myObliqueX / sqrt(s2) / sqrt(c2) - myObliqueY;
    return eps3;
}

double NPSTUVWXY::epsilonb() const
{
    double epsb = trueSM.epsilonb();
    return epsb;
}


////////////////////////////////////////////////////////////////////////     

double NPSTUVWXY::obliqueS() const
{
    double sW2_SM = trueSM.sW2();
    double sW_SM = sqrt(sW2_SM);
    double cW_SM = trueSM.cW2();
    return ( (myObliqueShat - myObliqueW + myObliqueX / (sW_SM * cW_SM) - myObliqueY)
            * 4.0 * sW2_SM / alphaMz());
}

double NPSTUVWXY::obliqueT() const
{
    double sW2_SM = trueSM.sW2();
    double sW_SM = sqrt(sW2_SM);
    double cW2_SM = trueSM.cW2();
    double cW_SM = sqrt(cW2_SM);
    return ( (myObliqueThat - myObliqueW + 2.0 * sW_SM / cW_SM * myObliqueX
            - sW2_SM / cW2_SM * myObliqueY) / alphaMz());
}

double NPSTUVWXY::obliqueU() const
{
    double sW2_SM = trueSM.sW2();
    double sW_SM = sqrt(sW2_SM);
    double cW_SM = sqrt(trueSM.cW2());
    return ( (-myObliqueUhat + myObliqueV + myObliqueW
            - 2.0 * sW_SM / cW_SM * myObliqueX)*4.0 * sW2_SM / alphaMz());
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::GammaW() const
{
    double Gamma_W = trueSM.GammaW();

    double Wbar = (obliqueV() - obliqueW()) / alphaMz();

    double alpha = StandardModel::alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    Gamma_W *= 1.0 - 3.0 * alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT()
            - (c2 - s2) * obliqueU() / 2.0 / s2 - 2.0 * (c2 - s2) * Wbar)
            - (1.0 + c2) / 2.0 / (c2 - s2) * DeltaGF();

    return Gamma_W;
}


double NPSTUVWXY::LEP2sigmaMu(const double s) const
{
   double sigma_mu;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   sigma_mu = trueSM.LEP2sigmaMu(s) + 
              myLEP2oblique.sigma_l_LEP2_NP(QCD::lepton(MU), s, trueSM.getLeptons(MU).getMass(), ObParam);
   
   return sigma_mu;
}

double NPSTUVWXY::LEP2sigmaTau(const double s) const
{
   double sigma_tau;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   sigma_tau = trueSM.LEP2sigmaTau(s) + 
              myLEP2oblique.sigma_l_LEP2_NP(QCD::lepton(TAU), s, trueSM.getLeptons(TAU).getMass(), ObParam);
   
   return sigma_tau;
}

double NPSTUVWXY::LEP2sigmaHadron(const double s) const
{
   double sigma_had;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   sigma_had = trueSM.LEP2sigmaHadron(s) + 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(UP), s, trueSM.getmq(QCD::quark(UP),sqrt(s)), ObParam)+ 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(DOWN), s, trueSM.getmq(QCD::quark(DOWN),sqrt(s)), ObParam)+ 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(STRANGE), s, trueSM.getmq(QCD::quark(STRANGE),sqrt(s)), ObParam)+ 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(CHARM), s, trueSM.getmq(QCD::quark(CHARM),sqrt(s)), ObParam)+ 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(BOTTOM), s, trueSM.getmq(QCD::quark(BOTTOM),sqrt(s)), ObParam);
   
   return sigma_had;
}

double NPSTUVWXY::LEP2sigmaCharm(const double s) const
{
   double sigma_charm;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   sigma_charm = trueSM.LEP2sigmaCharm(s) + 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(CHARM), s, trueSM.getmq(QCD::quark(CHARM),sqrt(s)), ObParam);
   
   return sigma_charm;
}

double NPSTUVWXY::LEP2sigmaBottom(const double s) const
{
   double sigma_bottom;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   sigma_bottom = trueSM.LEP2sigmaHadron(s) + 
              myLEP2oblique.sigma_q_LEP2_NP(QCD::quark(BOTTOM), s, trueSM.getmq(QCD::quark(BOTTOM),sqrt(s)), ObParam);
   
   return sigma_bottom;
}

double NPSTUVWXY::LEP2AFBmu(const double s) const
{
   double AFB_mu;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   AFB_mu = trueSM.LEP2sigmaMu(s) + 
              myLEP2oblique.AFB_l_LEP2_NP(QCD::lepton(MU), s, trueSM.getLeptons(MU).getMass(), ObParam);
   
   return AFB_mu;
}

double NPSTUVWXY::LEP2AFBtau(const double s) const
{
   double AFB_tau;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   AFB_tau = trueSM.LEP2sigmaTau(s) + 
              myLEP2oblique.AFB_l_LEP2_NP(QCD::lepton(TAU), s, trueSM.getLeptons(TAU).getMass(), ObParam);
   
   return AFB_tau;
}

double NPSTUVWXY::LEP2AFBcharm(const double s) const
{
   double AFB_charm;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   AFB_charm = trueSM.LEP2sigmaCharm(s) + 
              myLEP2oblique.AFB_q_LEP2_NP(QCD::quark(CHARM), s, trueSM.getmq(QCD::quark(CHARM),sqrt(s)), ObParam);
   
   return AFB_charm;
}

double NPSTUVWXY::LEP2AFBbottom(const double s) const
{
   double AFB_bottom;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   AFB_bottom = trueSM.LEP2sigmaHadron(s) + 
              myLEP2oblique.AFB_q_LEP2_NP(QCD::quark(BOTTOM), s, trueSM.getmq(QCD::quark(BOTTOM),sqrt(s)), ObParam);
   
   return AFB_bottom;
}


double NPSTUVWXY::LEP2Rcharm(const double s) const
{
   double R_charm;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   R_charm = trueSM.LEP2sigmaCharm(s) + 
              myLEP2oblique.R_q_LEP2_NP(QCD::quark(CHARM), s, trueSM.getmq(QCD::quark(CHARM),sqrt(s)), ObParam);
   
   return R_charm;
}

double NPSTUVWXY::LEP2Rbottom(const double s) const
{
   double R_bottom;
   double ObParam[7] = {obliqueShat(), obliqueThat(), obliqueUhat(),
                             obliqueV(), obliqueW(), obliqueX(), obliqueY()};
   
   R_bottom = trueSM.LEP2sigmaHadron(s) + 
              myLEP2oblique.R_q_LEP2_NP(QCD::quark(BOTTOM), s, trueSM.getmq(QCD::quark(BOTTOM),sqrt(s)), ObParam);
   
   return R_bottom;
}

