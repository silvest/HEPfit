/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2GIMR.h"


LEP2GIMR::LEP2GIMR(const StandardModel& SM_i)
: SM(SM_i)
{   
}



double LEP2GIMR::sigma_l_LEP2_GIMR(const QCD::lepton l, const double s,
                                   const double GIMRParam_i[]) const 
{
    double Mz = SM.getMz();
    double Nf = 1.0;
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double Ql = SM.getLeptons(l).getCharge();
    double GammaZ = SM.Gamma_Z();
    double deltaGammaZ = GIMRParam_i[delta_GammaZ];
    double alpha = SM.ale_OS(sqrt(s), FULLNLO);
    double alpha2 = alpha*alpha;
    double deltaGLl = GIMRParam_i[delta_gLf];
    double deltaGRl = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRl = gR_l(l);
    double GLl = gL_l(l);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
//    double sW = sqrt(GIMR.getTrueSM().sW2());
    double sW2 = SM.sW2();
    double cW = sqrt(1. - SM.cW2());
    double cW2 = cW*cW;
    double CLL = GIMRParam_i[C_LL];
    double CLR = GIMRParam_i[C_LR];
    double CRL = GIMRParam_i[C_RL];
    double CRR = GIMRParam_i[C_RR];
    complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
    complex chiZ = s/denom;
    complex deltachiZ = chiZ/denom*(complex(deltaMzsq,-Mz*deltaGammaZ-GammaZ*deltaMzsq/(2.0*Mz)));
    complex tmp = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    
    
    double ds = alpha/(6.0)*Nf*(CLL*(GLl*GLe/(cW2*sW2)+Ql*Qe)+CLR*(GLe*GRl/(cW2*sW2)+Ql*Qe)
              +CRL*(GLl*GRe/(cW2*sW2)+Qe*Ql)+CRR*(GRl*GRe/(cW2*sW2)+Qe*Ql))*chiZ.real()
              + alpha2*Nf*M_PI/(3*cW2*cW2*sW2*sW2*s)*(2.0*chiZ.abs2()*(GLe*(deltaGLl*GLl*GLe
              + deltaGLe*GLl*GLl+deltaGLe*GRl*GRl+deltaGRl*GLe*GRl) 
              + GRe*GRe*(deltaGLl*GLl+deltaGRl*GRl)+deltaGRe*GRe*(GLl*GLl+GRl*GRl))
              + tmp.real()*(GLl*GLl+GRl*GRl)*(GLe*GLe+GRe*GRe)
              + 2.0*sW2*cW2*Ql*Qe*(chiZ.real()*((deltaGLl+deltaGRl)*(GLe+GRe)
              + (deltaGLe+deltaGRe)*(GLl+GRl))+deltachiZ.real()*(GLl+GRl)*(GLe+GRe)));
    
    return ds;
   }


double LEP2GIMR::sigma_q_LEP2_GIMR(const QCD::quark q, const double s,
                                   const double GIMRParam_i[]) const
{
    double Mz = SM.getMz();
    double Nf = 3.0;
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double Qq = SM.getQuarks(q).getCharge();
    double GammaZ = SM.Gamma_Z();
    double deltaGammaZ = GIMRParam_i[delta_GammaZ];
    double alpha = SM.ale_OS(sqrt(s), FULLNLO);
    double alpha2 = alpha*alpha;
    double deltaGLq = GIMRParam_i[delta_gLf];
    double deltaGRq = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRq = gR_q(q);
    double GLq = gL_q(q);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
//    double sW = sqrt(GIMR.getTrueSM().sW2());
    double sW2 = SM.sW2();
    double cW = sqrt(1. - SM.cW2());
    double cW2 = cW*cW;
    double CLL = GIMRParam_i[C_LL];
    double CLR = GIMRParam_i[C_LR];
    double CRL = GIMRParam_i[C_RL];
    double CRR = GIMRParam_i[C_RR];
    complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
    complex chiZ = s/denom;
    complex deltachiZ = chiZ/denom*(complex(deltaMzsq,-Mz*deltaGammaZ-GammaZ*deltaMzsq/(2.0*Mz)));
    complex tmp = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    
    double ds = alpha/(6.0)*Nf*(CLL*(GLq*GLe/(cW2*sW2)+Qq*Qe)+CLR*(GLe*GRq/(cW2*sW2)+Qq*Qe)
              +CRL*(GLq*GRe/(cW2*sW2)+Qe*Qq)+CRR*(GRq*GRe/(cW2*sW2)+Qe*Qq))*chiZ.real()
              + alpha2*Nf*M_PI/(3*cW2*cW2*sW2*sW2*s)*(2.0*chiZ.abs2()*(GLe*(deltaGLq*GLq*GLe
              + deltaGLe*GLq*GLq+deltaGLe*GRq*GRq+deltaGRq*GLe*GRq) 
              + GRe*GRe*(deltaGLq*GLq+deltaGRq*GRq)+deltaGRe*GRe*(GLq*GLq+GRq*GRq))
              + tmp.real()*(GLq*GLq+GRq*GRq)*(GLe*GLe+GRe*GRe)
              + 2.0*sW2*cW2*Qq*Qe*(chiZ.real()*((deltaGLq+deltaGRq)*(GLe+GRe)
              + (deltaGLe+deltaGRe)*(GLq+GRq))+deltachiZ.real()*(GLq+GRq)*(GLe+GRe)));
    
    return ds;
    }

double LEP2GIMR::AFB_l_LEP2_GIMR(const QCD::lepton l, 
                                  const double s, const double Dim6Coef_i[]) const
{
    
    return 0;
    }


double LEP2GIMR::AFB_q_LEP2_GIMR(const QCD::quark q, 
                                  const double s, const double Dim6Coef_i[]) const
{
    
    return 0;
   }



double LEP2GIMR::gL_l(const QCD::lepton l) const
{
    double gA = SM.getLeptons(l).getIsospin();
    double gV = SM.getLeptons(l).getIsospin()-2.0*SM.getLeptons(l).getCharge()*SM.sW2();
    double gL = (gA + gV)/2.;

    return gL;
    }


double LEP2GIMR::gR_l(const QCD::lepton l) const
{
    double gA = SM.getLeptons(l).getIsospin();
    double gV = SM.getLeptons(l).getIsospin()-2.0*SM.getLeptons(l).getCharge()*SM.sW2();
    double gR = (gV - gA)/2.;

    return gR;
    }


double LEP2GIMR::gL_q(const QCD::quark q) const
{
    double gA = SM.getQuarks(q).getIsospin();
    double gV = SM.getQuarks(q).getIsospin()-2.0*SM.getQuarks(q).getCharge()*SM.sW2();
    double gL = (gA + gV)/2.;

    return gL;
    }


double LEP2GIMR::gR_q(const QCD::quark q) const
{
    double gA = SM.getQuarks(q).getIsospin();
    double gV = SM.getQuarks(q).getIsospin()-2.0*SM.getQuarks(q).getCharge()*SM.sW2();
    double gR = (gV - gA)/2.;

    return gR;
    }


double LEP2GIMR::sigmaF_l_LEP2_GIMR(const QCD::lepton l, const double s,
                                   const double GIMRParam_i[]) const 
{
    double Mz = SM.getMz();
    double Nf = 1.0;
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double Ql = SM.getLeptons(l).getCharge();
    double GammaZ = SM.Gamma_Z();
    double deltaGammaZ = GIMRParam_i[delta_GammaZ];
    double alpha = SM.ale_OS(sqrt(s), FULLNLO);
    double alpha2 = alpha*alpha;
    double deltaGLl = GIMRParam_i[delta_gLf];
    double deltaGRl = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRl = gR_l(l);
    double GLl = gL_l(l);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
//    double sW = sqrt(GIMR.getTrueSM().sW2());
    double sW2 = SM.sW2();
    double cW = sqrt(1. - SM.cW2());
    double cW2 = cW*cW;
    double CLL = GIMRParam_i[C_LL];
    double CLR = GIMRParam_i[C_LR];
    double CRL = GIMRParam_i[C_RL];
    double CRR = GIMRParam_i[C_RR];
    complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
    complex chiZ = s/denom;
    complex deltachiZ = chiZ/denom*(complex(deltaMzsq,-Mz*deltaGammaZ-GammaZ*deltaMzsq/(2.0*Mz)));
    complex tmp = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    
    
    double ds = alpha*Nf/48.0*chiZ.real()*(7.0*CLL*(GLl*GLe/cW2/sW2+Qe*Ql)+CLR*(GLe*GRl/cW2/sW2+Ql*Qe)
                +CRL*(GLl*GRe/cW2/sW2+Ql*Qe)+7.0*CRR*(GRl*GRe/cW2/sW2+Ql*Qe))
                +M_PI*alpha2*Nf/24./cW2/cW2/sW2/sW2/s*(2.*chiZ.abs2()*(GLe*(7.0*GLl*(deltaGLl*GLe
                +deltaGLe*GLl)+deltaGLe*GRl*GRl+deltaGRl*GLe*GRl))+GRe*GRe*(deltaGLl*GLl+7.0*deltaGRl*GRl)
                +deltaGRe*GRe*(GLl*GLl+7.0*GRl*GRl))+tmp.real()*(GLl*GLl*(7.0*GLe*GLe+GRe*GRe)+GRl*GRl*(GLe*GLe+7.0*GRe*GRe))
                +2.*Ql*Qe*cW2*sW2*(chiZ.real()*(deltaGLl*(7.0*GLe+GRe)+deltaGLe*(7.0*GLl+GRl)
                +deltaGRl*(GLe+7.0*GRe)+deltaGRe*(GLl+7.0*GRl))+deltachiZ.real()*(GLl*(7.0*GLe+GRe)+GRl*(GLe+7.0*GRe)));
    
    return ds;
   }

double LEP2GIMR::sigmaB_l_LEP2_GIMR(const QCD::lepton l, const double s,
                                   const double GIMRParam_i[]) const 
{
    double Mz = SM.getMz();
    double Nf = 1.0;
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double Ql = SM.getLeptons(l).getCharge();
    double GammaZ = SM.Gamma_Z();
    double deltaGammaZ = GIMRParam_i[delta_GammaZ];
    double alpha = SM.ale_OS(sqrt(s), FULLNLO);
    double alpha2 = alpha*alpha;
    double deltaGLl = GIMRParam_i[delta_gLf];
    double deltaGRl = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRl = gR_l(l);
    double GLl = gL_l(l);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
//    double sW = sqrt(GIMR.getTrueSM().sW2());
    double sW2 = SM.sW2();
    double cW = sqrt(1. - SM.cW2());
    double cW2 = cW*cW;
    double CLL = GIMRParam_i[C_LL];
    double CLR = GIMRParam_i[C_LR];
    double CRL = GIMRParam_i[C_RL];
    double CRR = GIMRParam_i[C_RR];
    complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
    complex chiZ = s/denom;
    complex deltachiZ = chiZ/denom*(complex(deltaMzsq,-Mz*deltaGammaZ-GammaZ*deltaMzsq/(2.0*Mz)));
    complex tmp = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    
    
    double ds = alpha*Nf/48.0*chiZ.real()*(CLL*(GLl*GLe/cW2/sW2+Qe*Ql)+7.0*CLR*(GLe*GRl/cW2/sW2+Ql*Qe)
                +7.0*CRL*(GLl*GRe/cW2/sW2+Ql*Qe)+CRR*(GRl*GRe/cW2/sW2+Ql*Qe))
                +M_PI*alpha2*Nf/24./cW2/cW2/sW2/sW2/s*(2.*chiZ.abs2()*(GLe*(GLl*(deltaGLl*GLe
                +deltaGLe*GLl)+7.0*deltaGLe*GRl*GRl+7.0*deltaGRl*GLe*GRl))+GRe*GRe*(7.0*deltaGLl*GLl+deltaGRl*GRl)
                +deltaGRe*GRe*(7.0*GLl*GLl+GRl*GRl))+tmp.real()*(GLl*GLl*(GLe*GLe+7.0*GRe*GRe)+GRl*GRl*(7.0*GLe*GLe+GRe*GRe))
                +2.*Ql*Qe*cW2*sW2*(chiZ.real()*(deltaGLl*(GLe+7.0*GRe)+deltaGLe*(GLl+7.0*GRl)
                +deltaGRl*(7.0*GLe+GRe)+deltaGRe*(7.0*GLl+GRl))+deltachiZ.real()*(GLl*(GLe+7.0*GRe)+GRl*(7.0*GLe+GRe)));
    
    return ds;
   }


double LEP2GIMR::sigmaF_q_LEP2_GIMR(const QCD::quark q, const double s,
                                   const double GIMRParam_i[]) const
{
    double Mz = SM.getMz();
    double Nf = 3.0;
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double Qq = SM.getQuarks(q).getCharge();
    double GammaZ = SM.Gamma_Z();
    double deltaGammaZ = GIMRParam_i[delta_GammaZ];
    double alpha = SM.ale_OS(sqrt(s), FULLNLO);
    double alpha2 = alpha*alpha;
    double deltaGLq = GIMRParam_i[delta_gLf];
    double deltaGRq = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRq = gR_q(q);
    double GLq = gL_q(q);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
//    double sW = sqrt(GIMR.getTrueSM().sW2());
    double sW2 = SM.sW2();
    double cW = sqrt(1. - SM.cW2());
    double cW2 = cW*cW;
    double CLL = GIMRParam_i[C_LL];
    double CLR = GIMRParam_i[C_LR];
    double CRL = GIMRParam_i[C_RL];
    double CRR = GIMRParam_i[C_RR];
    complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
    complex chiZ = s/denom;
    complex deltachiZ = chiZ/denom*(complex(deltaMzsq,-Mz*deltaGammaZ-GammaZ*deltaMzsq/(2.0*Mz)));
    complex tmp = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    
    double ds = alpha*Nf/48.0*chiZ.real()*(7.0*CLL*(GLq*GLe/cW2/sW2+Qe*Qq)+CLR*(GLe*GRq/cW2/sW2+Qq*Qe)
                +CRL*(GLq*GRe/cW2/sW2+Qq*Qe)+7.0*CRR*(GRq*GRe/cW2/sW2+Qq*Qe))
                +M_PI*alpha2*Nf/24./cW2/cW2/sW2/sW2/s*(2.*chiZ.abs2()*(GLe*(7.0*GLq*(deltaGLq*GLe
                +deltaGLe*GLq)+deltaGLe*GRq*GRq+deltaGRq*GLe*GRq))+GRe*GRe*(deltaGLq*GLq+7.0*deltaGRq*GRq)
                +deltaGRe*GRe*(GLq*GLq+7.0*GRq*GRq))+tmp.real()*(GLq*GLq*(7.0*GLe*GLe+GRe*GRe)+GRq*GRq*(GLe*GLe+7.0*GRe*GRe))
                +2.*Qq*Qe*cW2*sW2*(chiZ.real()*(deltaGLq*(7.0*GLe+GRe)+deltaGLe*(7.0*GLq+GRq)
                +deltaGRq*(GLe+7.0*GRe)+deltaGRe*(GLq+7.0*GRq))+deltachiZ.real()*(GLq*(7.0*GLe+GRe)+GRq*(GLe+7.0*GRe)));
    
    return ds;
    }


double LEP2GIMR::sigmaB_q_LEP2_GIMR(const QCD::quark q, const double s,
                                   const double GIMRParam_i[]) const
{
    double Mz = SM.getMz();
    double Nf = 3.0;
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double Qq = SM.getQuarks(q).getCharge();
    double GammaZ = SM.Gamma_Z();
    double deltaGammaZ = GIMRParam_i[delta_GammaZ];
    double alpha = SM.ale_OS(sqrt(s), FULLNLO);
    double alpha2 = alpha*alpha;
    double deltaGLq = GIMRParam_i[delta_gLf];
    double deltaGRq = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRq = gR_q(q);
    double GLq = gL_q(q);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
//    double sW = sqrt(GIMR.getTrueSM().sW2());
    double sW2 = SM.sW2();
    double cW = sqrt(1. - SM.cW2());
    double cW2 = cW*cW;
    double CLL = GIMRParam_i[C_LL];
    double CLR = GIMRParam_i[C_LR];
    double CRL = GIMRParam_i[C_RL];
    double CRR = GIMRParam_i[C_RR];
    complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
    complex chiZ = s/denom;
    complex deltachiZ = chiZ/denom*(complex(deltaMzsq,-Mz*deltaGammaZ-GammaZ*deltaMzsq/(2.0*Mz)));
    complex tmp = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    
    double ds = alpha*Nf/48.0*chiZ.real()*(CLL*(GLq*GLe/cW2/sW2+Qe*Qq)+7.0*CLR*(GLe*GRq/cW2/sW2+Qq*Qe)
                +7.0*CRL*(GLq*GRe/cW2/sW2+Qq*Qe)+CRR*(GRq*GRe/cW2/sW2+Qq*Qe))
                +M_PI*alpha2*Nf/24./cW2/cW2/sW2/sW2/s*(2.*chiZ.abs2()*(GLe*(GLq*(deltaGLq*GLe
                +deltaGLe*GLq)+7.0*deltaGLe*GRq*GRq+7.0*deltaGRq*GLe*GRq))+GRe*GRe*(7.0*deltaGLq*GLq+deltaGRq*GRq)
                +deltaGRe*GRe*(7.0*GLq*GLq+GRq*GRq))+tmp.real()*(GLq*GLq*(GLe*GLe+7.0*GRe*GRe)+GRq*GRq*(7.0*GLe*GLe+GRe*GRe))
                +2.*Qq*Qe*cW2*sW2*(chiZ.real()*(deltaGLq*(GLe+7.0*GRe)+deltaGLe*(GLq+7.0*GRq)
                +deltaGRq*(7.0*GLe+GRe)+deltaGRe*(7.0*GLq+GRq))+deltachiZ.real()*(GLq*(GLe+7.0*GRe)+GRq*(7.0*GLe+GRe)));
    
    return ds;
    }