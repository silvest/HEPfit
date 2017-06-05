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
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRl = gR_l(l);
    double GLl = gL_l(l);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
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
    complex chideltachi = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    double dA1l=deltaA1l(l,GIMRParam_i);
    double dA2l=deltaA2l(l,GIMRParam_i);
    double dB1l=deltaB1l(l,GIMRParam_i);
    double dB2l=deltaB2l(l,GIMRParam_i);
    
    
    double ds = alpha/(6.0)*Nf*(CLL*(GLl*GLe/(cW2*sW2)*chiZ.real()+Ql*Qe)
              + CLR*(GLe*GRl/(cW2*sW2)*chiZ.real()+Ql*Qe)
              + CRL*(GLl*GRe/(cW2*sW2)*chiZ.real()+Qe*Ql)
              + CRR*(GRl*GRe/(cW2*sW2)*chiZ.real()+Qe*Ql))
              + 2.*alpha2*Nf*M_PI*Ql*Qe/(3.*cW2*sW2*s)*(chiZ.real()*(dA1l+dA2l)
              + deltachiZ.real()*(GLl+GRl)*(GLe+GRe))
              + alpha2*Nf*M_PI/(3.*cW2*cW2*sW2*sW2*s)*(2.*chiZ.abs2()*(dB1l+dB2l)
              + chideltachi.real()*(GLl*GLl+GRl*GRl)*(GLe*GLe+GRe*GRe));
    
    return ds;
   }


double LEP2GIMR::sigmaFminusB_l_LEP2_GIMR(const QCD::lepton l, const double s,
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
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRl = gR_l(l);
    double GLl = gL_l(l);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
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
    complex chideltachi = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    double dA1l=deltaA1l(l,GIMRParam_i);
    double dA2l=deltaA2l(l,GIMRParam_i);
    double dB1l=deltaB1l(l,GIMRParam_i);
    double dB2l=deltaB2l(l,GIMRParam_i);
    
    
    double ds = alpha/(8.0)*Nf*(CLL*(GLl*GLe/(cW2*sW2)*chiZ.real()+Ql*Qe)
              - CLR*(GLe*GRl/(cW2*sW2)*chiZ.real()+Ql*Qe)
              - CRL*(GLl*GRe/(cW2*sW2)*chiZ.real()+Qe*Ql)
              + CRR*(GRl*GRe/(cW2*sW2)*chiZ.real()+Qe*Ql))
              + alpha2*Nf*M_PI*Ql*Qe/(2.*cW2*sW2*s)*(chiZ.real()*(dA1l-dA2l)
              + deltachiZ.real()*(GLl-GRl)*(GLe-GRe))
              + alpha2*Nf*M_PI/(4.*cW2*cW2*sW2*sW2*s)*(2.*chiZ.abs2()*(dB1l-dB2l)
              + chideltachi.real()*(GLl*GLl-GRl*GRl)*(GLe*GLe-GRe*GRe));
    
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
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRq = gR_q(q);
    double GLq = gL_q(q);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
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
    complex chideltachi = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    double dA1q=deltaA1q(q,GIMRParam_i);
    double dA2q=deltaA2q(q,GIMRParam_i);
    double dB1q=deltaB1q(q,GIMRParam_i);
    double dB2q=deltaB2q(q,GIMRParam_i);
    
    double ds = alpha/(6.0)*Nf*(CLL*(GLq*GLe/(cW2*sW2)*chiZ.real()+Qq*Qe)
              + CLR*(GLe*GRq/(cW2*sW2)*chiZ.real()+Qq*Qe)
              + CRL*(GLq*GRe/(cW2*sW2)*chiZ.real()+Qe*Qq)
              + CRR*(GRq*GRe/(cW2*sW2)*chiZ.real()+Qe*Qq))
              + 2.*alpha2*Nf*M_PI*Qq*Qe/(3.*cW2*sW2*s)*(chiZ.real()*(dA1q+dA2q)
              + deltachiZ.real()*(GLq+GRq)*(GLe+GRe))
              + alpha2*Nf*M_PI/(3.*cW2*cW2*sW2*sW2*s)*(2.*chiZ.abs2()*(dB1q+dB2q)
              + chideltachi.real()*(GLq*GLq+GRq*GRq)*(GLe*GLe+GRe*GRe));
    
    return ds;
    }

double LEP2GIMR::sigmaFminusB_q_LEP2_GIMR(const QCD::quark q, const double s,
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
    double deltaMzsq = GIMRParam_i[delta_Mz2];
    double GRq = gR_q(q);
    double GLq = gL_q(q);
    double GRe = gR_l(StandardModel::ELECTRON);
    double GLe = gL_l(StandardModel::ELECTRON);
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
    complex chideltachi = (chiZ*deltachiZ.conjugate()+chiZ.conjugate()*deltachiZ);
    double dA1q=deltaA1q(q,GIMRParam_i);
    double dA2q=deltaA2q(q,GIMRParam_i);
    double dB1q=deltaB1q(q,GIMRParam_i);
    double dB2q=deltaB2q(q,GIMRParam_i);
    
    double ds = alpha/(8.0)*Nf*(CLL*(GLq*GLe/(cW2*sW2)*chiZ.real()+Qq*Qe)
              - CLR*(GLe*GRq/(cW2*sW2)*chiZ.real()+Qq*Qe)
              - CRL*(GLq*GRe/(cW2*sW2)*chiZ.real()+Qe*Qq)
              + CRR*(GRq*GRe/(cW2*sW2)*chiZ.real()+Qe*Qq))
              + alpha2*Nf*M_PI*Qq*Qe/(4.*cW2*sW2*s)*(chiZ.real()*(dA1q-dA2q)
              + deltachiZ.real()*(GLq-GRq)*(GLe-GRe))
              + alpha2*Nf*M_PI/(4.*cW2*cW2*sW2*sW2*s)*(2.*chiZ.abs2()*(dB1q-dB2q)
              + chideltachi.real()*(GLq*GLq-GRq*GRq)*(GLe*GLe-GRe*GRe));
    
    return ds;
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


double LEP2GIMR::deltaA1q(const QCD::quark q, const double GIMRParam_i[]) const
{
    double A1q;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    A1q = deltaGLe * gL_q(q) + deltaGRe * gR_q(q) 
        + gL_l(StandardModel::ELECTRON) * deltaGLf 
        + gR_l(StandardModel::ELECTRON) * deltaGRf;
    
    return A1q;
}

double LEP2GIMR::deltaA2q(const QCD::quark q, const double GIMRParam_i[]) const
{
    double A2q;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    A2q = deltaGLe * gR_q(q) + deltaGRe * gL_q(q) 
        + gL_l(StandardModel::ELECTRON) * deltaGRf 
        + gR_l(StandardModel::ELECTRON) * deltaGLf;
    
    return A2q;
}

double LEP2GIMR::deltaB1q(const QCD::quark q, const double GIMRParam_i[]) const
{
    double B1q;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double gL_q2 = gL_q(q) * gL_q(q);
    double gR_q2 = gR_q(q) * gR_q(q);
    double gL_l2 = gL_l(StandardModel::ELECTRON) * gL_l(StandardModel::ELECTRON);
    double gR_l2 = gR_l(StandardModel::ELECTRON) * gR_l(StandardModel::ELECTRON);
    B1q = gL_l(StandardModel::ELECTRON) * gL_q2 * deltaGLe 
        + gR_l(StandardModel::ELECTRON) * gR_q2 * deltaGRe 
        + gL_l2 * gL_q(q) * deltaGLf + gR_l2 * gR_q(q) * deltaGRf;
    
    return B1q;
}

double LEP2GIMR::deltaB2q(const QCD::quark q, const double GIMRParam_i[]) const
{
    double B2q;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double gL_q2 = gL_q(q) * gL_q(q);
    double gR_q2 = gR_q(q) * gR_q(q);
    double gL_l2 = gL_l(StandardModel::ELECTRON) * gL_l(StandardModel::ELECTRON);
    double gR_l2 = gR_l(StandardModel::ELECTRON) * gR_l(StandardModel::ELECTRON);
    B2q = gL_l(StandardModel::ELECTRON) * gR_q2 * deltaGLe 
        + gR_l(StandardModel::ELECTRON) * gL_q2 * deltaGRe 
        + gR_l2 * gL_q(q) * deltaGLf + gL_l2 * gR_q(q) * deltaGRf;
    
    return B2q;    
}

double LEP2GIMR::deltaA1l(const QCD::lepton l, const double GIMRParam_i[]) const
{
    double A1l;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    A1l = deltaGLe * gL_l(l) + deltaGRe * gR_l(l) 
        + gL_l(l) * deltaGLf + gR_l(l) * deltaGRf;
    
    return A1l;
}

double LEP2GIMR::deltaA2l(const QCD::lepton l, const double GIMRParam_i[]) const
{
    double A2l;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    A2l = deltaGLe * gR_l(l) + deltaGRe * gL_l(l) 
        + gL_l(l) * deltaGRf + gR_l(l) * deltaGLf;
    
    return A2l;
}

double LEP2GIMR::deltaB1l(const QCD::lepton l, const double GIMRParam_i[]) const
{
    double B1l;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double gL_l2 = gL_l(l) * gL_l(l);
    double gR_l2 = gR_l(l) * gR_l(l);
    B1l = gL_l(l) * gL_l2 * deltaGLe + gR_l(l) * gR_l2 * deltaGRe 
        + gL_l2 * gL_l(l) * deltaGLf + gR_l2 * gR_l(l) * deltaGRf;
    
    return B1l;    
}

double LEP2GIMR::deltaB2l(const QCD::lepton l, const double GIMRParam_i[]) const
{
    double B2l;
    double deltaGLf = GIMRParam_i[delta_gLf];
    double deltaGRf = GIMRParam_i[delta_gRf];
    double deltaGLe = GIMRParam_i[delta_gLe];
    double deltaGRe = GIMRParam_i[delta_gRe];
    double gL_l2 = gL_l(l) * gL_l(l);
    double gR_l2 = gR_l(l) * gR_l(l);
    B2l = gL_l(l) * gR_l2 * deltaGLe + gR_l(l) * gL_l2 * deltaGRe 
        + gR_l2 * gL_l(l) * deltaGLf + gL_l2 * gR_l(l) * deltaGRf;
    
    return B2l;     
}
