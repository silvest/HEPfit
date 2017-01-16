/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Metastability.h"
#include "SUSY.h"
#include <limits>
#include <math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <boost/bind.hpp>

FindAction::FindAction(const StandardModel& SM_i) :
    ThObservable(SM_i),
    potentialcoefficientspar(35,0.),
    mySUSY(static_cast<const SUSY&> (SM_i))
{
    mySUSYScalarPotential=new SUSYScalarPotential(SM_i);
}

FindAction::~FindAction()
{
  delete mySUSYScalarPotential;
}

double FindAction::computeThValue()
{

    /* 1. Definition of all potential minima */

    // define scalar potential
    gslpp::vector<double> potentialcoefficients=mySUSYScalarPotential->coefficients();

    double Vmin0=mySUSYScalarPotential->potential(potentialcoefficients, 0.0, 0.0, 0.0);

//    std::cout << "Vmin0 = " << Vmin0 << std::endl;

    gslpp::vector<double> dV0=mySUSYScalarPotential->potentialderivative(potentialcoefficients, 0.0, 0.0, 0.0);

//    std::cout << "dV0 = " << dV0 << std::endl;

    // calculate all minima (Ayan)
    // calculate V and dV at origin (V0 and dV0)
    // check whether there is a charge-breaking minimum very close to the origin
//    double S_0=0.0;
    gslpp::vector<double> S(10,0.);
    //define toy minima (remove this as soon as the true minimum finder is available)
    gslpp::vector<double> minimavector(3,0.);
    minimavector(0)=-807192.888;
    minimavector(1)=807260.056;
    minimavector(2)=-803413.309;
//    minimavector(3)=4.0;
//    minimavector(4)=5.0;
//    minimavector(5)=6.0;
//    minimavector(6)=7.0;
//    minimavector(7)=8.0;
//    minimavector(8)=9.0;
    //if mod3 length of minima is not zero, throw error
    int lengthofminima=3;
    int NofMinima=lengthofminima/3;

    gslpp::vector<double> deeperminima(lengthofminima+NofMinima,0.), dV(3,0.);
//    gslpp::vector<double> deeperminima(3,0.), dV(3,0.);
    int i,n=0;
    double x1, x2, x3, Vmin;

    for(i=0;i<NofMinima;i++)
    {
        x1=minimavector(3*i);
        x2=minimavector(3*i+1);
        x3=minimavector(3*i+2);
        Vmin=mySUSYScalarPotential->potential(potentialcoefficients, x1, x2, x3);
    std::cout << "Vmin = " << Vmin << std::endl;
        if(Vmin>=Vmin0)
        {
            continue;
        }
        else
        {
            deeperminima(4*n)=minimavector(3*i);
            deeperminima(4*n+1)=minimavector(3*i+1);
            deeperminima(4*n+2)=minimavector(3*i+2);
            deeperminima(4*n+3)=Vmin;
            n++;
        }
    }

    std::cout << "2." << std::endl;

    /* 2. Calculation of the tunneling rate for each of the minima */

//    gslpp::vector<double> path;
    //n is now the number of deeper minima
    for(i=0;i<n;i++)
    {
        x1=deeperminima(4*i);
        x2=deeperminima(4*i+1);
        x3=deeperminima(4*i+2);
        Vmin=deeperminima(4*i+3);
        dV=mySUSYScalarPotential->potentialderivative(potentialcoefficients, x1, x2, x3);
//        path=splinepath(x1,x2,x3,Vmin0,dV0,Vmin,dV);
        //the following function is supposed to find the cubic spline of the potential along a straight line from one to the other minimum
                int steps=100;
                gslpp::vector<double> linearpath(steps,0.);
                gslpp::vector<double> V(steps,0.);
//                gslpp::vector<double> spath(steps,0.);
                //maybe the following will not work in extreme cases
                double distance=x1*x1+x2*x2+x3*x3;
                double stepsize=distance/double(steps-1);
                int k;
                for(k=0;k<steps;k++)
                {
                    linearpath(k)=stepsize*k;
                    V(k)=mySUSYScalarPotential->potential(potentialcoefficients, x1*linearpath(k), x2*linearpath(k), x3*linearpath(k));
//    gsl_interp_accel *acc 
//      = gsl_interp_accel_alloc ();
//    gsl_spline *spline 
//      = gsl_spline_alloc (gsl_interp_cspline, 10);
//
//    gsl_spline_init (spline, x, y, 10);
//
//    for (xi = x[0]; xi < x[9]; xi += 0.01)
//      {
//        yi = gsl_spline_eval (spline, xi, acc);
//        printf ("%g %g\n", xi, yi);
//      }
//    gsl_spline_free (spline);
//        gsl_interp_accel_free (acc);

                }


    /* 2.1 Calculation of the tunneling rate along the undeformed (straight) path */

                    std::cout << "2.1.1" << std::endl;

//        double xmin = 0.001;
//        double xmax = std::numeric_limits<double>::infinity();
        
    /* 2.1.1 Finding the edge of the tunneling barrier */

        double min = 1.0e-12;
        double pos1=0.0;  //relative position of the metastable minimum
        double pos2=1.0;  //relative position of the stable minimum
        double pos = 0.5;  //relative position between metastable and stable minimum (first guess)
        while(fabs(pos1-pos2) > min)
        {
            double Vpos = mySUSYScalarPotential->potential(potentialcoefficients, pos*x1, pos*x2, pos*x3);
            if(Vpos > 0.0)  //maybe this should be ">"
            {
                pos1 = pos;
            }
            else
            {
                pos2 = pos;
            }
            pos = 0.5*(pos1+pos2);
        }
        //pos is now the edge of the barrier

    /* 2.1.2 Finding the typical scale of the barrier (using gsl 1D minimization) */

                    std::cout << "2.1.2" << std::endl;

        int status;
        int iter = 0, max_iter = 100;
        const gsl_min_fminimizer_type *T;
        gsl_min_fminimizer *s;
        double m=0.5*pos;
        double a=0.0, b=pos;
                    std::cout << "pos =" << pos << std::endl;
        potentialcoefficientspar=potentialcoefficients;
        x1par=x1;
        x2par=x2;
        x3par=x3;
        gsl_function F = convertToGslFunctionS(boost::bind(&FindAction::invertedpotential, &(*this), _1));
//        F.function = -&invertedpotential; //pointer?
        T = gsl_min_fminimizer_brent;
        s = gsl_min_fminimizer_alloc(T);
        gsl_min_fminimizer_set(s, &F, m, a, b);
        do
        {
            iter++;
            status = gsl_min_fminimizer_iterate(s);
            m = gsl_min_fminimizer_x_minimum(s);
                    std::cout << "m =" << m << std::endl;
            a = gsl_min_fminimizer_x_lower(s);
            b = gsl_min_fminimizer_x_upper(s);
            status = gsl_min_test_interval(a, b, 0.0, 1.0e-3);
        } while(status == GSL_CONTINUE && iter<max_iter);
        gsl_min_fminimizer_free(s);
        double barriertop=m;
        if(barriertop<=0.0 || barriertop>=pos)
        {
            throw std::runtime_error("Error in Metastability.cpp: Potential barrier top outside the barrier range!");
        }
        double barrierheight=mySUSYScalarPotential->potential(potentialcoefficients, barriertop*x1, barriertop*x2, barriertop*x3)-Vmin0;

        if(barrierheight<=0.0)
        {
            throw std::runtime_error("Error in Metastability.cpp: No potential barrier!");
        }
        double rscale = barriertop/sqrt(6.0*barrierheight);

    /* 2.1.3 Finding a phi solution for the tunneling using an overshooting/undershooting algorithm  */

        double x = -log(pos);
//        double xincrease = 5.0;

        double rmin = 1.e-4*rscale;
        double rmax = 1.e4*rscale;
        double dr0 = rmin;
        double drmin = 0.01*rmin;

        double delta_phi = distance;
        double delta_phi_cutoff = 1.e-2 * delta_phi;
        double epsabs[2] = {fabs(delta_phi)*1.e-4 , fabs(delta_phi/rscale)*1.e-4};
        double epsfrac[2] = {1.e-4 , 1.e-4};

        double eps = distance*1.e-3;
        gslpp::vector<double> inconds(3,0.);
        do
        {
            double delta_phi0 = distance - exp(-x)*delta_phi;
            double sdp = delta_phi0/distance; //scaled_delta_phi0

            double dV_at_delta_phi0 = (mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0-2.0*eps)*x1/distance, (delta_phi0-2.0*eps)*x2/distance, (delta_phi0-2.0*eps)*x3/distance) 
                                       -8.0*mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0-eps)*x1/distance, (delta_phi0-eps)*x2/distance, (delta_phi0-eps)*x3/distance)
                                       +8.0*mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0+eps)*x1/distance, (delta_phi0+eps)*x2/distance, (delta_phi0+eps)*x3/distance) 
                                       -mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0+2.0*eps)*x1/distance, (delta_phi0+2.0*eps)*x2/distance, (delta_phi0+2.0*eps)*x3/distance) ) / (12.0*eps);

            double d2V_at_phi0 = (-mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0-2.0*eps)*x1/distance, (delta_phi0-2.0*eps)*x2/distance, (delta_phi0-2.0*eps)*x3/distance)
                                  +16.0*mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0-eps)*x1/distance, (delta_phi0-eps)*x2/distance, (delta_phi0-eps)*x3/distance)
                                  -30.0*mySUSYScalarPotential->potential(potentialcoefficients, sdp*x1, sdp*x2, sdp*x3)
                                  +16.0*mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0+eps)*x1/distance, (delta_phi0+eps)*x2/distance, (delta_phi0+eps)*x3/distance)
                                  -mySUSYScalarPotential->potential(potentialcoefficients, (delta_phi0+2.0*eps)*x1/distance, (delta_phi0+2.0*eps)*x2/distance, (delta_phi0+2.0*eps)*x3/distance) ) / (12.0*eps*eps);

            inconds = InitialConditions(delta_phi0, rmin, delta_phi_cutoff, distance, dV_at_delta_phi0, d2V_at_phi0);
            if(!std::isfinite(inconds(0)) || !std::isfinite(x))
            {
                break;                
            }

            //r_y = rf, yf1, yf2, convergence=1(converged), 2(undershoot), 3(overshoot)
            gslpp::vector<double> r_y(4,0.);
            r_y = integrateProfile(inconds(0), inconds(1), inconds(2), dr0, epsfrac, epsabs, drmin, rmax, distance);
//
//                        # Check for overshoot, undershoot
//                        if ctype == "converged":
//                            break
//                        elif ctype == "undershoot": # x is too low
//                            xmin = x
//                            x = x*xincrease if xmax == np.inf else .5*(xmin+xmax)
//                        elif ctype == "overshoot": # x is too high
//                            xmax = x
//                            x = .5*(xmin+xmax)
//                        # Check if we've reached xtol
//                        if (xmax-xmin) < xtol:
//                            break

        }
        while(true);

            /* 2.1.4 Integrate the barrier to get the tunneling action */

    /* 2.2 Calculation of the tunneling rate along the deformed path */

        /* 2.2.1 Deforming the path */
        /* 2.2.2 Calculating the tunneling rate */

    }//loop over the minima



    //          deform path between origin and i (old path is a line, new path is p_i)
        //define splines between origin and i
            //take origin and i
            //for i_spline < ?
            //add one knot, calculate its V and dV
        //split into parallel and orthogonal components
        //solve parallel tunneling
        //calculate N and minimize it wrt orthogonal directions
    //          calculate action integrating along p_i
        int rlength = 10;
    gslpp::vector<double> r(rlength,0.);
    gslpp::vector<double> phi(rlength,0.);
    gslpp::vector<double> dphi(rlength,0.);
//        double rlength = length of r;
        double VphiMin_i = deformedV(phi(rlength));
        //integrate.simps(integrand,r)
        double integral = 0.0;
        for(int j=1;j<rlength;j++)
        {
            integral += (r(j)-r(j-1))*(Simpsonintegrand(r(j-1),phi(j-1),dphi(j-1),VphiMin_i)
                                       +4.0*Simpsonintegrand((r(j)+r(j-1))/2.0,(phi(j)+phi(j-1))/2.0,(dphi(j)+dphi(j-1))/2.0,VphiMin_i)
                                       +Simpsonintegrand(r(j),phi(j),dphi(j),VphiMin_i))/(3.0*(double)rlength);
        }

//        double volume = r(0)*r(0)*r(0)*r(0)*M_PI*M_PI/2.0;
        
//        S.assign(j, integral+volume*(deformedV(phi(0))-VphiMin_i));
    //          check whether action is larger than S_(i-1)

    return 0.0;
}

double FindAction::Simpsonintegrand(double r, double phi, double dphi, double VphiMin_i)
{
    return (0.5*dphi*dphi + deformedV(phi)-VphiMin_i) * 2.0*M_PI*M_PI*r*r*r;
}

double FindAction::deformedV(double phi)
{
    return 0.0;
}
////
////gslpp::vector<double> FindAction::splinepath(double x1, double x2, double x3, double Vmin0, gslpp::vector<double> dV0, double Vmin, gslpp::vector<double> dV)
////{
////    int steps=100;
////    gslpp::vector<double> linearpath(steps,0.);
////    gslpp::vector<double> spath(steps,0.);
////    //maybe the following will not work in extreme cases
////    double distance=x1*x1+x2*x2+x3*x3;
////    double stepsize=distance/double(steps-1);
////    int i;
////    for(i=0;i<steps;i++)
////    {
////        linearpath(i)=stepsize*i;
////        Vmin=mySUSYScalarPotential->potential(potentialcoefficients, x1, x2, x3);
////    }
////    
////    return spath;
////}

gslpp::vector<double> FindAction::InitialConditions(double delta_phi0, double rmin, double delta_phi_cutoff, double distance, double dV_at_delta_phi0, double d2V_at_phi0)
{
    double phi0 = distance + delta_phi0;
    double dV = dV_at_delta_phi0;
    double d2V = d2V_at_phi0;
    gslpp::vector<double> exsol(2,0.);
//switch this on!    exsol=ExactSolution(rmin, phi0, dV, d2V);
    double phi_r0=exsol(0);
    double dphi_r0=exsol(1);

    gslpp::vector<double> inc(3,0.);
    inc(0)=rmin;
    inc(1)=phi_r0;
    inc(2)=dphi_r0;

    if(fabs(phi_r0)<fabs(delta_phi_cutoff) || dphi_r0*delta_phi0>0)
    {
        double r = rmin;
        double rlast=r;
        while(std::isfinite(r))
        {
            rlast = r;
            r *= 10;
//switch this on!            exsol=ExactSolution(r, phi0, dV, d2V);
            if(fabs(exsol(0))>fabs(delta_phi_cutoff))
            {
                break;
            }
        }

        rpar=r;
        phi0par=phi0;
        dVpar=dV;
        d2Vpar=d2V;
        delta_phi_cutoffpar=delta_phi_cutoff;

        gsl_function F = convertToGslFunctionS(boost::bind(&FindAction::func, &(*this), _1));
        //gsl root finding algorithm
                      int status;
                      int iter = 0, max_iter = 100;
                      const gsl_root_fsolver_type *T;
                      gsl_root_fsolver *s;
                      double r2 = 0;
                      double x_lo = rlast, x_hi = r;
                      T = gsl_root_fsolver_brent;
                      s = gsl_root_fsolver_alloc (T);
                      gsl_root_fsolver_set (s, &F, x_lo, x_hi);
                      do
                        {
                          iter++;
                          status = gsl_root_fsolver_iterate (s);
                          r2 = gsl_root_fsolver_root (s);
                          x_lo = gsl_root_fsolver_x_lower (s);
                          x_hi = gsl_root_fsolver_x_upper (s);
                          status = gsl_root_test_interval (x_lo, x_hi,
                                                           0, 0.001);
                        }
                      while (status == GSL_CONTINUE && iter < max_iter);
                      gsl_root_fsolver_free (s);
        //now r2 is the root
        exsol = ExactSolution(r2, phi0, dV, d2V);
        inc(1)=exsol(0);
        inc(2)=exsol(1);        
    }

    return inc;
}

//This function calculates the analytical solution for barriers with low height and curvature (?)
gslpp::vector<double> FindAction::ExactSolution(double r, double phi0, double dV, double d2V)
{
        gslpp::vector<double> phi_dphi(2,0.);
        double phi = 0.0;
        double dphi = 0.0;

        double curv = sqrt(fabs(d2V));
        double curv_r = curv*r;
        double nu = 1.0;
////        double gsl_sf_gamma(double x) = special.gamma # Gamma function
////        double gsl_sf_bessel_In (int n, double x)
        ////Function: double gsl_sf_bessel_Jn (int n, double x)
////        iv, jv = special.iv, special.jv # (modified) Bessel function
        if(curv_r<1.e-2)
        {
            double s = (d2V>0) ? 1.0 : -1.0;
            double add_to_phi;
            for(int k=1;k<=4;k++)
            {
                add_to_phi = pow(0.5*curv_r,2.0*k-2.0) * pow(s,k) / (gsl_sf_gamma(k+1)*gsl_sf_gamma(k+1+nu));
                phi += add_to_phi;
                dphi += add_to_phi * (2.0*k);
            }
            phi *= 0.25 * gsl_sf_gamma(nu+1) * r*r * dV * s;
            dphi *= 0.25 * gsl_sf_gamma(nu+1) * r * dV * s;
            phi += phi0;
        }
        else if(d2V>0)
        {
            phi = (gsl_sf_gamma(nu+1)*pow(0.5*curv_r,-nu) *gsl_sf_bessel_In(int(nu),curv_r)-1.0) * dV/d2V;
            dphi = -nu*(pow(0.5*curv_r,-nu) / r) * gsl_sf_bessel_In(int(nu), curv_r);
            dphi += pow(0.5*curv_r,-nu) * 0.5*curv * (gsl_sf_bessel_In(int(nu)-1, curv_r)+gsl_sf_bessel_In(int(nu)+1, curv_r));
            dphi *= gsl_sf_gamma(nu+1) * dV/d2V;
            phi += phi0;
        }
        else
        {
            phi = (gsl_sf_gamma(nu+1)*pow(0.5*curv_r,-nu) * gsl_sf_bessel_Jn(int(nu), curv_r)-1.0) * dV/d2V;
            dphi = -nu*(pow(0.5*curv_r,-nu) / r) * gsl_sf_bessel_Jn(int(nu), curv_r);
            dphi += pow(0.5*curv_r,-nu) * 0.5*curv * (gsl_sf_bessel_Jn(int(nu)-1, curv_r)-gsl_sf_bessel_Jn(int(nu)+1, curv_r));
            dphi *= gsl_sf_gamma(nu+1) * dV/d2V;
            phi += phi0;
        }
        phi_dphi(0)=phi;
        phi_dphi(1)=dphi;
        return phi_dphi;
}

gslpp::vector<double> FindAction::integrateProfile(double r0, double y01, double y02, double dr0, double epsfrac[2], double epsabs[2], double drmin, double rmax, double distance)
{
    double dr = dr0;
////  dY= y02, dV(y[0])-3.0*y02/r
    gslpp::vector<double> dydr0(2,0.);
    dydr0 = dY(y01, y02, r0);
////    ysign = np.sign(y0[0]-self.phi_metaMin) 
////            (positive means we're heading down, negative means heading up.)
    rmax += r0;

//    int i=1;
//    int convergence_type = 0;

    gslpp::vector<double> ret_values(4,0.);

    do
    {
        gslpp::vector<double> r_y_drnext(4,0.);
//        r_y_drnext = rkqs(y01, y02, dydr0, r0, dr, epsfrac, epsabs);

////        step_rkqs(0)=t1;
////        step_rkqs(1)=y[0];
////        step_rkqs(2)=y[1];
////        step_rkqs(3)=dr;

        double r1 = r_y_drnext(0);
        double y11 = y01+r_y_drnext(1);
        double y12 = y02+r_y_drnext(2);

        gslpp::vector<double> dydr1 = dY(y11,y12,r1);

        if (r1 > rmax)
        {
            throw std::runtime_error("r > rmax");
        }
        else if (dr < drmin)
        {
            throw std::runtime_error("dr < drmin");
        }
//        else if ( fabs(y11 - distance) < 3.0*epsabs[0] && fabs(y12) < 3.0*epsabs[1] )
//        {
//            ret_values(0)=r1;
//            ret_values(1)=y11;
//            ret_values(2)=y12;
//            ret_values(3)=1;    //converged
//            break;
//        }
        else if ( y12*(y01-distance) > 0 || (y11-distance)*(y01-distance) < 0 )
        {
//            double f=y0*(1-t)**3 + 3*y1*(1-t)*(1-t)*t + 3*y2*(1-t)*t*t + y1*t**3;
//(self, y0, dy0, y1, dy1)=            y0*mt**3 + 3*y1*mt*mt*t + 3*y2*mt*t*t + y3*t**3
//            ret_values(0)=r1;
//            ret_values(1)=y11;
//            ret_values(2)=y12;
            ret_values(3)=2;    //undershoot
            break;
            
        }
//        else if ( fabs(y11 - distance) < 3.0*epsabs && fabs(y12) < 3.0*epsabs )
//        {
            
//        }
        
//          elif( y1[1]*(y0[0]-self.phi_metaMin) > 0 or (y1[0]-self.phi_metaMin)*(y0[0]-self.phi_metaMin) < 0 ):
//                f = cubicInterpFunction(y0, dr*dydr0, y1, dr*dydr1)
//                if(y1[1]*(y0[0]-self.phi_metaMin) > 0):
//                    # Extrapolate to where dphi(r) = 0
//                    x = optimize.brentq(lambda x: f(x)[1], 0, 1 )
//                    convergence_type = "undershoot"
//                else:
//                    # Extrapolate to where phi(r) = phi_metaMin
//                    x = optimize.brentq(lambda x: f(x)[0]-self.phi_metaMin, 0,1)
//                    convergence_type = "overshoot"
//                r = r0 + dr*x
//                y = f(x)
//                break
//            # Advance the integration variables
//            r0,y0,dydr0 = r1,y1,dydr1
//            dr = drnext
    }
    while(true);
//        # Check convergence for a second time. 
//        # The extrapolation in overshoot/undershoot might have gotten us within
//        # the acceptable error.
//        if (abs(y - np.array([self.phi_metaMin,0])) < 3*epsabs).all():
//            convergence_type = "converged"
//        return self._integrateProfile_rval(r, y, convergence_type)

    return ret_values;
}

gslpp::vector<double> FindAction::dY(double y1, double y2, double r)
{
    gslpp::vector<double> values(2,0.);
//    double alpha=3.0; /*spatial dimensions*/
    values(0)=y2;
//    values(1)=(mySUSYScalarPotential->potential(potentialcoefficients, (y1-2.0*eps)*x1/distance, (y1-2.0*eps)*x2/distance, (y1-2.0*eps)*x3/distance) 
//                                       -8.0*mySUSYScalarPotential->potential(potentialcoefficients, (y1-eps)*x1/distance, (y1-eps)*x2/distance, (y1-eps)*x3/distance)
//                                       +8.0*mySUSYScalarPotential->potential(potentialcoefficients, (y1+eps)*x1/distance, (y1+eps)*x2/distance, (y1+eps)*x3/distance) 
//                                       -mySUSYScalarPotential->potential(potentialcoefficients, (y1+2.0*eps)*x1/distance, (y1+2.0*eps)*x2/distance, (y1+2.0*eps)*x3/distance) ) / (12.0*eps)
//              -alpha*y2/r;
    return values;
}

int dYfunc(double r, const double y[], double ODE[], void *flags)
{
    gslpp::vector<double> dYvalues(2,0.);
//    dYvalues=dY(y[0],y[1],r);
    ODE[0]=dYvalues(0);
    ODE[1]=dYvalues(1);
    return 0;
}

int dYJac(double r, const double y[], double *dfdy, double dfdt[], void *order)
{
    return 0;
}

//gslpp::vector<double> FindAction::rkqs(double y01, double y02, gslpp::vector<double> dydr0, double r0, double dr, double epsfrac[2], double epsabs[2])
//{
//    gslpp::vector<double> step_rkqs(4,0.);
//    double mu = 0;
//    gsl_odeiv2_system ODEsystem = { dYfunc, dYJac, 2, &mu };
////    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&ODEsystem, gsl_odeiv2_step_rkck, dr, epsabs, epsfrac);
//    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkck;
//    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (2);
//    //CHANGE THIS! (for the moment take only the 0 element)
//    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (epsabs[0], epsfrac[0]);
//    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 2);
//    double t = r0;
//    double t1 = r0+dr;
//    double y[2] = {y01, y02};
//    int status = gsl_odeiv2_evolve_apply(e, c, s, &ODEsystem, &r0, t1, &dr, y);
////    gsl_odeiv2_driver_free (d);
//    step_rkqs(0)=t1;
//    step_rkqs(1)=y[0];
//    step_rkqs(2)=y[1];
//    step_rkqs(3)=dr;
//    gsl_odeiv2_evolve_free (e);
//    gsl_odeiv2_control_free (c);
//    gsl_odeiv2_step_free (s);
//    return step_rkqs;

////    dt = dt_try
////    while True:
////        dy,yerr = _rkck(y,dydt,t,f,dt,args)
////        errmax = np.nan_to_num( np.max( np.min([abs(yerr/epsabs), 
////            abs(yerr)/((abs(y)+1e-300)*epsfrac)],axis=0) ) )
////        if(errmax < 1.0):
////            break # Step succeeded
////        dttemp = 0.9*dt*errmax**-.25
////        dt = max(dttemp,dt*.1) if dt > 0 else min(dttemp,dt*.1)
////        if(t+dt==t):
////            raise IntegrationError("Stepsize rounds down to zero.")
////    if errmax > 1.89e-4:
////        dtnext = 0.9 * dt * errmax**-.2
////    else:
////        dtnext = 5*dt
////    return _rkqs_rval(dy, dt, dtnext)

//}

//////
//////double rk(double InitialValues[], unsigned long int NumberOfODEs, double r1, double r2, int order)
//////{
//////    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
//////    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, NumberOfODEs);
//////
//////    //Define the absolute (A) and relative (R) error on y at each step.
//////    //The real error will be compared to the following error estimate:
//////    //  A + R * |y_i|
//////    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-6, 0.0);
//////
//////    //Allocate space for the evolutor
//////    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(NumberOfRGEs);
//////
//////    //Possibility to define a set of parameters which the RGE's depend on
////////    double order = 1;
//////
//////    //Definition of the RGE system (the Jacobian is not necessary for the RK4 method; it's an empty function here)
////////    gsl_odeiv2_system RGEsystem = {mySM.RGEs, Jacobian, mySM.NumberOfRGEs, &RGEparameters};
//////    gsl_odeiv2_system RGEsystem = {RGEs, Jacobian, NumberOfRGEs, &order};
//////
//////    //Set starting and end point as natural logarithmic scales (conversion from decadic log scale)
//////    double t1 = Q1*log(10.0);
//////    double t2 = Q2*log(10.0);
//////    double tNLOuni = NLOuniscale*log(10.0);
//////
//////    //Set initial step size
//////    double InitialStepSize = 1e-6;
//////
//////    //Run!
//////    while (t1 < t2)
//////    {
//////        int status = gsl_odeiv2_evolve_apply (e, c, s, &RGEsystem, &t1, t2, &InitialStepSize, InitialValues);
//////        if(status != GSL_SUCCESS) break;
//////
//////        //intermediate checks if appropriate
//////        if(RGEcheck(InitialValues,t1,Rpeps,tNLOuni) != 0) break;
//////    }
//////
//////    gsl_odeiv2_evolve_free (e);
//////    gsl_odeiv2_control_free (c);
//////    gsl_odeiv2_step_free (s);
//////
////////    for(int i=0;i<14;i++) std::cout<<InitialValues[i]<<std::endl;
//////
//////    //Return the decadic log scale at which the evolution stopped
//////    return t1/log(10.0);
//////}
