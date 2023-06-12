//RGESolver.cc
#include "RGESolver.h"
#include <boost/bind/bind.hpp>
//#include <chrono> //To measure execution times



#include "StaticMembers.cpp"

#include "SettersAndGetters.cpp"
#include "SMInput.cpp"

RGESolver::RGESolver() : yuR(3, 0.), yuI(3, 0.), ydR(3, 0.), ydI(3, 0.), yeR(3, 0.), yeI(3, 0.) {
  
  using namespace boost::placeholders;
  
  //Assigns default values for epsabs and epsrel
  
  Resetepsrel();
  Resetepsabs();
  SetSMDefaultInput();
  
   
  //-----SM SCALAR PARAMETERS I/O----------
  {
    Operators0F["g1"] = & g1;
    Operators0F["g2"] = & g2;
    Operators0F["g3"] = & g3;
    Operators0F["mh2"] = & mh2;
    Operators0F["lambda"] = & lambda;
    
    CKMAngles["s12"] = & s12;
    CKMAngles["s13"] = & s13;
    CKMAngles["s23"] = & s23;


  }


  //-----------0F OPERATORS I/O------------
  {
    //Class 1-3 operators 
    Operators0F["CG"] = & cG;
    Operators0F["CGtilde"] = & cGT;
    Operators0F["CW"] = & cW;
    Operators0F["CWtilde"] = & cWT;

    Operators0F["CH"] = & cH;

    Operators0F["CHbox"] = & cHBOX;
    Operators0F["CHD"] = & cHD;

    //Class 4 operators 
    Operators0F["CHG"] = & cHG;
    Operators0F["CHGtilde"] = & cHGT;
    Operators0F["CHW"] = & cHW;
    Operators0F["CHWtilde"] = & cHWT;
    Operators0F["CHB"] = & cHB;
    Operators0F["CHBtilde"] = & cHBT;
    Operators0F["CHWB"] = & cHWB;
    Operators0F["CHWtildeB"] = & cHWBT;


  }
  //-----------2F OPERATORS I/O------------
  Setter2F["YuR"] = boost::bind(&Yukawa_set, &yuR, _1, _2, _3);
  Setter2F["YuI"] = boost::bind(&Yukawa_set, &yuI, _1, _2, _3);
  Getter2F["YuR"] = boost::bind(&Yukawa, &yuR, _1, _2);
  Getter2F["YuI"] = boost::bind(&Yukawa, &yuI, _1, _2);

  Setter2F["YeR"] = boost::bind(&Yukawa_set, &yeR, _1, _2, _3);
  Setter2F["YeI"] = boost::bind(&Yukawa_set, &yeI, _1, _2, _3);
  Getter2F["YeR"] = boost::bind(&Yukawa, &yeR, _1, _2);
  Getter2F["YeI"] = boost::bind(&Yukawa, &yeI, _1, _2);

  Setter2F["YdR"] = boost::bind(&Yukawa_set, &ydR, _1, _2, _3);
  Setter2F["YdI"] = boost::bind(&Yukawa_set, &ydI, _1, _2, _3);
  Getter2F["YdR"] = boost::bind(&Yukawa, &ydR, _1, _2);
  Getter2F["YdI"] = boost::bind(&Yukawa, &ydI, _1, _2);
  {

    //Class 5 operators 

    Setter2F["CuHR"] = boost::bind(&WC1_set, cuHR, _1, _2, _3);
    Setter2F["CuHI"] = boost::bind(&WC1_set, cuHI, _1, _2, _3);
    Getter2F["CuHR"] = boost::bind(&WC1, cuHR, _1, _2);
    Getter2F["CuHI"] = boost::bind(&WC1, cuHI, _1, _2);

    Setter2F["CdHR"] = boost::bind(&WC1_set, cdHR, _1, _2, _3);
    Setter2F["CdHI"] = boost::bind(&WC1_set, cdHI, _1, _2, _3);
    Getter2F["CdHR"] = boost::bind(&WC1, cdHR, _1, _2);
    Getter2F["CdHI"] = boost::bind(&WC1, cdHI, _1, _2);

    Setter2F["CeHR"] = boost::bind(&WC1_set, ceHR, _1, _2, _3);
    Setter2F["CeHI"] = boost::bind(&WC1_set, ceHI, _1, _2, _3);
    Getter2F["CeHR"] = boost::bind(&WC1, ceHR, _1, _2);
    Getter2F["CeHI"] = boost::bind(&WC1, ceHI, _1, _2);

    //Class 6 operators 

    Setter2F["CeWR"] = boost::bind(&WC1_set, ceWR, _1, _2, _3);
    Setter2F["CeWI"] = boost::bind(&WC1_set, ceWI, _1, _2, _3);
    Getter2F["CeWR"] = boost::bind(&WC1, ceWR, _1, _2);
    Getter2F["CeWI"] = boost::bind(&WC1, ceWI, _1, _2);

    Setter2F["CeBR"] = boost::bind(&WC1_set, ceBR, _1, _2, _3);
    Setter2F["CeBI"] = boost::bind(&WC1_set, ceBI, _1, _2, _3);
    Getter2F["CeBR"] = boost::bind(&WC1, ceBR, _1, _2);
    Getter2F["CeBI"] = boost::bind(&WC1, ceBI, _1, _2);



    Setter2F["CuGR"] = boost::bind(&WC1_set, cuGR, _1, _2, _3);
    Setter2F["CuGI"] = boost::bind(&WC1_set, cuGI, _1, _2, _3);
    Getter2F["CuGR"] = boost::bind(&WC1, cuGR, _1, _2);
    Getter2F["CuGI"] = boost::bind(&WC1, cuGI, _1, _2);

    Setter2F["CuWR"] = boost::bind(&WC1_set, cuWR, _1, _2, _3);
    Setter2F["CuWI"] = boost::bind(&WC1_set, cuWI, _1, _2, _3);
    Getter2F["CuWR"] = boost::bind(&WC1, cuWR, _1, _2);
    Getter2F["CuWI"] = boost::bind(&WC1, cuWI, _1, _2);

    Setter2F["CuBR"] = boost::bind(&WC1_set, cuBR, _1, _2, _3);
    Setter2F["CuBI"] = boost::bind(&WC1_set, cuBI, _1, _2, _3);
    Getter2F["CuBR"] = boost::bind(&WC1, cuBR, _1, _2);
    Getter2F["CuBI"] = boost::bind(&WC1, cuBI, _1, _2);



    Setter2F["CdGR"] = boost::bind(&WC1_set, cdGR, _1, _2, _3);
    Setter2F["CdGI"] = boost::bind(&WC1_set, cdGI, _1, _2, _3);
    Getter2F["CdGR"] = boost::bind(&WC1, cdGR, _1, _2);
    Getter2F["CdGI"] = boost::bind(&WC1, cdGI, _1, _2);

    Setter2F["CdWR"] = boost::bind(&WC1_set, cdWR, _1, _2, _3);
    Setter2F["CdWI"] = boost::bind(&WC1_set, cdWI, _1, _2, _3);
    Getter2F["CdWR"] = boost::bind(&WC1, cdWR, _1, _2);
    Getter2F["CdWI"] = boost::bind(&WC1, cdWI, _1, _2);

    Setter2F["CdBR"] = boost::bind(&WC1_set, cdBR, _1, _2, _3);
    Setter2F["CdBI"] = boost::bind(&WC1_set, cdBI, _1, _2, _3);
    Getter2F["CdBR"] = boost::bind(&WC1, cdBR, _1, _2);
    Getter2F["CdBI"] = boost::bind(&WC1, cdBI, _1, _2);



    //Class 7 operators 

    Setter2F["CHl1R"] = boost::bind(&WC2R_set, cHl1R, _1, _2, _3);
    Getter2F["CHl1R"] = boost::bind(&WC2R, cHl1R, _1, _2);
    Setter2F["CHl1I"] = boost::bind(&WC2I_set, cHl1I, _1, _2, _3);
    Getter2F["CHl1I"] = boost::bind(&WC2I, cHl1I, _1, _2);

    Setter2F["CHl3R"] = boost::bind(&WC2R_set, cHl3R, _1, _2, _3);
    Getter2F["CHl3R"] = boost::bind(&WC2R, cHl3R, _1, _2);
    Setter2F["CHl3I"] = boost::bind(&WC2I_set, cHl3I, _1, _2, _3);
    Getter2F["CHl3I"] = boost::bind(&WC2I, cHl3I, _1, _2);

    Setter2F["CHeR"] = boost::bind(&WC2R_set, cHeR, _1, _2, _3);
    Getter2F["CHeR"] = boost::bind(&WC2R, cHeR, _1, _2);
    Setter2F["CHeI"] = boost::bind(&WC2I_set, cHeI, _1, _2, _3);
    Getter2F["CHeI"] = boost::bind(&WC2I, cHeI, _1, _2);

    Setter2F["CHq1R"] = boost::bind(&WC2R_set, cHq1R, _1, _2, _3);
    Getter2F["CHq1R"] = boost::bind(&WC2R, cHq1R, _1, _2);
    Setter2F["CHq1I"] = boost::bind(&WC2I_set, cHq1I, _1, _2, _3);
    Getter2F["CHq1I"] = boost::bind(&WC2I, cHq1I, _1, _2);

    Setter2F["CHq3R"] = boost::bind(&WC2R_set, cHq3R, _1, _2, _3);
    Getter2F["CHq3R"] = boost::bind(&WC2R, cHq3R, _1, _2);
    Setter2F["CHq3I"] = boost::bind(&WC2I_set, cHq3I, _1, _2, _3);
    Getter2F["CHq3I"] = boost::bind(&WC2I, cHq3I, _1, _2);

    Setter2F["CHuR"] = boost::bind(&WC2R_set, cHuR, _1, _2, _3);
    Getter2F["CHuR"] = boost::bind(&WC2R, cHuR, _1, _2);
    Setter2F["CHuI"] = boost::bind(&WC2I_set, cHuI, _1, _2, _3);
    Getter2F["CHuI"] = boost::bind(&WC2I, cHuI, _1, _2);

    Setter2F["CHdR"] = boost::bind(&WC2R_set, cHdR, _1, _2, _3);
    Getter2F["CHdR"] = boost::bind(&WC2R, cHdR, _1, _2);
    Setter2F["CHdI"] = boost::bind(&WC2I_set, cHdI, _1, _2, _3);
    Getter2F["CHdI"] = boost::bind(&WC2I, cHdI, _1, _2);

    Setter2F["CHudR"] = boost::bind(&WC1_set, cHudR, _1, _2, _3);
    Getter2F["CHudR"] = boost::bind(&WC1, cHudR, _1, _2);
    Setter2F["CHudI"] = boost::bind(&WC1_set, cHudI, _1, _2, _3);
    Getter2F["CHudI"] = boost::bind(&WC1, cHudI, _1, _2);
  }
  //-----------4F OPERATORS I/O------------
  {
    //Class 8LLLL operators        
    Setter4F["CllR"] = boost::bind(&WC6R_set, cllR, _1, _2, _3, _4, _5);
    Getter4F["CllR"] = boost::bind(&WC6R, cllR, _1, _2, _3, _4);
    Setter4F["CllI"] = boost::bind(&WC6I_set, cllI, _1, _2, _3, _4, _5);
    Getter4F["CllI"] = boost::bind(&WC6I, cllI, _1, _2, _3, _4);

    Setter4F["Cqq1R"] = boost::bind(&WC6R_set, cqq1R, _1, _2, _3, _4, _5);
    Getter4F["Cqq1R"] = boost::bind(&WC6R, cqq1R, _1, _2, _3, _4);
    Setter4F["Cqq1I"] = boost::bind(&WC6I_set, cqq1I, _1, _2, _3, _4, _5);
    Getter4F["Cqq1I"] = boost::bind(&WC6I, cqq1I, _1, _2, _3, _4);

    Setter4F["Cqq3R"] = boost::bind(&WC6R_set, cqq3R, _1, _2, _3, _4, _5);
    Getter4F["Cqq3R"] = boost::bind(&WC6R, cqq3R, _1, _2, _3, _4);
    Setter4F["Cqq3I"] = boost::bind(&WC6I_set, cqq3I, _1, _2, _3, _4, _5);
    Getter4F["Cqq3I"] = boost::bind(&WC6I, cqq3I, _1, _2, _3, _4);

    Setter4F["Clq1R"] = boost::bind(&WC7R_set, clq1R, _1, _2, _3, _4, _5);
    Getter4F["Clq1R"] = boost::bind(&WC7R, clq1R, _1, _2, _3, _4);
    Setter4F["Clq1I"] = boost::bind(&WC7I_set, clq1I, _1, _2, _3, _4, _5);
    Getter4F["Clq1I"] = boost::bind(&WC7I, clq1I, _1, _2, _3, _4);

    Setter4F["Clq3R"] = boost::bind(&WC7R_set, clq3R, _1, _2, _3, _4, _5);
    Getter4F["Clq3R"] = boost::bind(&WC7R, clq3R, _1, _2, _3, _4);
    Setter4F["Clq3I"] = boost::bind(&WC7I_set, clq3I, _1, _2, _3, _4, _5);
    Getter4F["Clq3I"] = boost::bind(&WC7I, clq3I, _1, _2, _3, _4);

    //Class 8RRRR operators   
    Setter4F["CeeR"] = boost::bind(&WC8R_set, ceeR, _1, _2, _3, _4, _5);
    Getter4F["CeeR"] = boost::bind(&WC8R, ceeR, _1, _2, _3, _4);
    Setter4F["CeeI"] = boost::bind(&WC8I_set, ceeI, _1, _2, _3, _4, _5);
    Getter4F["CeeI"] = boost::bind(&WC8I, ceeI, _1, _2, _3, _4);

    Setter4F["CuuR"] = boost::bind(&WC6R_set, cuuR, _1, _2, _3, _4, _5);
    Getter4F["CuuR"] = boost::bind(&WC6R, cuuR, _1, _2, _3, _4);
    Setter4F["CuuI"] = boost::bind(&WC6I_set, cuuI, _1, _2, _3, _4, _5);
    Getter4F["CuuI"] = boost::bind(&WC6I, cuuI, _1, _2, _3, _4);

    Setter4F["CddR"] = boost::bind(&WC6R_set, cddR, _1, _2, _3, _4, _5);
    Getter4F["CddR"] = boost::bind(&WC6R, cddR, _1, _2, _3, _4);
    Setter4F["CddI"] = boost::bind(&WC6I_set, cddI, _1, _2, _3, _4, _5);
    Getter4F["CddI"] = boost::bind(&WC6I, cddI, _1, _2, _3, _4);

    Setter4F["CeuR"] = boost::bind(&WC7R_set, ceuR, _1, _2, _3, _4, _5);
    Getter4F["CeuR"] = boost::bind(&WC7R, ceuR, _1, _2, _3, _4);
    Setter4F["CeuI"] = boost::bind(&WC7I_set, ceuI, _1, _2, _3, _4, _5);
    Getter4F["CeuI"] = boost::bind(&WC7I, ceuI, _1, _2, _3, _4);

    Setter4F["CedR"] = boost::bind(&WC7R_set, cedR, _1, _2, _3, _4, _5);
    Getter4F["CedR"] = boost::bind(&WC7R, cedR, _1, _2, _3, _4);
    Setter4F["CedI"] = boost::bind(&WC7I_set, cedI, _1, _2, _3, _4, _5);
    Getter4F["CedI"] = boost::bind(&WC7I, cedI, _1, _2, _3, _4);

    Setter4F["Cud1R"] = boost::bind(&WC7R_set, cud1R, _1, _2, _3, _4, _5);
    Getter4F["Cud1R"] = boost::bind(&WC7R, cud1R, _1, _2, _3, _4);
    Setter4F["Cud1I"] = boost::bind(&WC7I_set, cud1I, _1, _2, _3, _4, _5);
    Getter4F["Cud1I"] = boost::bind(&WC7I, cud1I, _1, _2, _3, _4);

    Setter4F["Cud8R"] = boost::bind(&WC7R_set, cud8R, _1, _2, _3, _4, _5);
    Getter4F["Cud8R"] = boost::bind(&WC7R, cud8R, _1, _2, _3, _4);
    Setter4F["Cud8I"] = boost::bind(&WC7I_set, cud8I, _1, _2, _3, _4, _5);
    Getter4F["Cud8I"] = boost::bind(&WC7I, cud8I, _1, _2, _3, _4);


    //Class 8_LLRR
    Setter4F["CleR"] = boost::bind(&WC7R_set, cleR, _1, _2, _3, _4, _5);
    Getter4F["CleR"] = boost::bind(&WC7R, cleR, _1, _2, _3, _4);
    Setter4F["CleI"] = boost::bind(&WC7I_set, cleI, _1, _2, _3, _4, _5);
    Getter4F["CleI"] = boost::bind(&WC7I, cleI, _1, _2, _3, _4);

    Setter4F["CluR"] = boost::bind(&WC7R_set, cluR, _1, _2, _3, _4, _5);
    Getter4F["CluR"] = boost::bind(&WC7R, cluR, _1, _2, _3, _4);
    Setter4F["CluI"] = boost::bind(&WC7I_set, cluI, _1, _2, _3, _4, _5);
    Getter4F["CluI"] = boost::bind(&WC7I, cluI, _1, _2, _3, _4);

    Setter4F["CldR"] = boost::bind(&WC7R_set, cldR, _1, _2, _3, _4, _5);
    Getter4F["CldR"] = boost::bind(&WC7R, cldR, _1, _2, _3, _4);
    Setter4F["CldI"] = boost::bind(&WC7I_set, cldI, _1, _2, _3, _4, _5);
    Getter4F["CldI"] = boost::bind(&WC7I, cldI, _1, _2, _3, _4);

    Setter4F["CqeR"] = boost::bind(&WC7R_set, cqeR, _1, _2, _3, _4, _5);
    Getter4F["CqeR"] = boost::bind(&WC7R, cqeR, _1, _2, _3, _4);
    Setter4F["CqeI"] = boost::bind(&WC7I_set, cqeI, _1, _2, _3, _4, _5);
    Getter4F["CqeI"] = boost::bind(&WC7I, cqeI, _1, _2, _3, _4);

    Setter4F["Cqu1R"] = boost::bind(&WC7R_set, cqu1R, _1, _2, _3, _4, _5);
    Getter4F["Cqu1R"] = boost::bind(&WC7R, cqu1R, _1, _2, _3, _4);
    Setter4F["Cqu1I"] = boost::bind(&WC7I_set, cqu1I, _1, _2, _3, _4, _5);
    Getter4F["Cqu1I"] = boost::bind(&WC7I, cqu1I, _1, _2, _3, _4);

    Setter4F["Cqu8R"] = boost::bind(&WC7R_set, cqu8R, _1, _2, _3, _4, _5);
    Getter4F["Cqu8R"] = boost::bind(&WC7R, cqu8R, _1, _2, _3, _4);
    Setter4F["Cqu8I"] = boost::bind(&WC7I_set, cqu8I, _1, _2, _3, _4, _5);
    Getter4F["Cqu8I"] = boost::bind(&WC7I, cqu8I, _1, _2, _3, _4);

    Setter4F["Cqd1R"] = boost::bind(&WC7R_set, cqd1R, _1, _2, _3, _4, _5);
    Getter4F["Cqd1R"] = boost::bind(&WC7R, cqd1R, _1, _2, _3, _4);
    Setter4F["Cqd1I"] = boost::bind(&WC7I_set, cqd1I, _1, _2, _3, _4, _5);
    Getter4F["Cqd1I"] = boost::bind(&WC7I, cqd1I, _1, _2, _3, _4);

    Setter4F["Cqd8R"] = boost::bind(&WC7R_set, cqd8R, _1, _2, _3, _4, _5);
    Getter4F["Cqd8R"] = boost::bind(&WC7R, cqd8R, _1, _2, _3, _4);
    Setter4F["Cqd8I"] = boost::bind(&WC7I_set, cqd8I, _1, _2, _3, _4, _5);
    Getter4F["Cqd8I"] = boost::bind(&WC7I, cqd8I, _1, _2, _3, _4);

    //Operators classes 8_LRRL and 8_LRLR

    Setter4F["CledqR"] = boost::bind(&WC5_set, cledqR, _1, _2, _3, _4, _5);
    Getter4F["CledqR"] = boost::bind(&WC5, cledqR, _1, _2, _3, _4);
    Setter4F["CledqI"] = boost::bind(&WC5_set, cledqI, _1, _2, _3, _4, _5);
    Getter4F["CledqI"] = boost::bind(&WC5, cledqI, _1, _2, _3, _4);

    Setter4F["Cquqd1R"] = boost::bind(&WC5_set, cquqd1R, _1, _2, _3, _4, _5);
    Getter4F["Cquqd1R"] = boost::bind(&WC5, cquqd1R, _1, _2, _3, _4);
    Setter4F["Cquqd1I"] = boost::bind(&WC5_set, cquqd1I, _1, _2, _3, _4, _5);
    Getter4F["Cquqd1I"] = boost::bind(&WC5, cquqd1I, _1, _2, _3, _4);

    Setter4F["Cquqd8R"] = boost::bind(&WC5_set, cquqd8R, _1, _2, _3, _4, _5);
    Getter4F["Cquqd8R"] = boost::bind(&WC5, cquqd8R, _1, _2, _3, _4);
    Setter4F["Cquqd8I"] = boost::bind(&WC5_set, cquqd8I, _1, _2, _3, _4, _5);
    Getter4F["Cquqd8I"] = boost::bind(&WC5, cquqd8I, _1, _2, _3, _4);

    Setter4F["Clequ1R"] = boost::bind(&WC5_set, clequ1R, _1, _2, _3, _4, _5);
    Getter4F["Clequ1R"] = boost::bind(&WC5, clequ1R, _1, _2, _3, _4);
    Setter4F["Clequ1I"] = boost::bind(&WC5_set, clequ1I, _1, _2, _3, _4, _5);
    Getter4F["Clequ1I"] = boost::bind(&WC5, clequ1I, _1, _2, _3, _4);

    Setter4F["Clequ3R"] = boost::bind(&WC5_set, clequ3R, _1, _2, _3, _4, _5);
    Getter4F["Clequ3R"] = boost::bind(&WC5, clequ3R, _1, _2, _3, _4);
    Setter4F["Clequ3I"] = boost::bind(&WC5_set, clequ3I, _1, _2, _3, _4, _5);
    Getter4F["Clequ3I"] = boost::bind(&WC5, clequ3I, _1, _2, _3, _4);
  }
}

RGESolver::~ RGESolver() {
  gsl_odeiv2_step_free(s);
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_step_free(sSMOnly);
  gsl_odeiv2_evolve_free(eSMOnly);
}

void RGESolver::Init() {

  int n = 3;
  int a, b, i, j, k, l;
  x[0] = g2;
  x[1] = g1;
  x[2] = g3;
  x[n++] = lambda;
  x[n++] = mh2;

  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      int count = 0;
      x[n + count++ * DF] = yuR(i, j);
      x[n + count++ * DF] = yuI(i, j);
      x[n + count++ * DF] = ydR(i, j);
      x[n + count++ * DF] = ydI(i, j);
      x[n + count++ * DF] = yeR(i, j);
      x[n++ + count++ * DF] = yeI(i, j);
    }
  }
  n += (2. * Nyukawa - 1.) * DF;


  x[n++] = cG;
  x[n++] = cGT;
  x[n++] = cW;
  x[n++] = cWT;
  x[n++] = cH;
  x[n++] = cHBOX;
  x[n++] = cHD;

  x[n++] = cHG;
  x[n++] = cHB;
  x[n++] = cHW;
  x[n++] = cHWB;
  x[n++] = cHGT;
  x[n++] = cHBT;
  x[n++] = cHWT;
  x[n++] = cHWBT;

  //Class 5
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      a = 0;
      x[n + 2 * a * DF] = WC1(cuHR, i, j);
      x[n + (2 * a++ + 1) * DF] = WC1(cuHI, i, j);
      x[n + 2 * a * DF] = WC1(cdHR, i, j);
      x[n + (2 * a++ + 1) * DF] = WC1(cdHI, i, j);
      x[n + 2 * a * DF] = WC1(ceHR, i, j);
      x[n++ + (2 * a++ + 1) * DF] = WC1(ceHI, i, j);
    }
  }
  n += (N5 * 2 - 1) * DF;

  //Class 6
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      a = 0;
      x[n + a++ * DF] = WC1(ceWR, i, j);
      x[n + a++ * DF] = WC1(ceWI, i, j);
      x[n + a++ * DF] = WC1(ceBR, i, j);
      x[n + a++ * DF] = WC1(ceBI, i, j);
      x[n + a++ * DF] = WC1(cuGR, i, j);
      x[n + a++ * DF] = WC1(cuGI, i, j);
      x[n + a++ * DF] = WC1(cuWR, i, j);
      x[n + a++ * DF] = WC1(cuWI, i, j);
      x[n + a++ * DF] = WC1(cuBR, i, j);
      x[n + a++ * DF] = WC1(cuBI, i, j);
      x[n + a++ * DF] = WC1(cdGR, i, j);
      x[n + a++ * DF] = WC1(cdGI, i, j);
      x[n + a++ * DF] = WC1(cdWR, i, j);
      x[n + a++ * DF] = WC1(cdWI, i, j);
      x[n + a++ * DF] = WC1(cdBR, i, j);
      x[n++ + a++ * DF] = WC1(cdBI, i, j);
    }
  }
  n += (N6 * 2 - 1) * DF;

  //Class 7 


  for (i = 0; i < DWC2R; i ++) 
    x[n++] = WC2R(cHl1R, WC2R_indices[i][0], WC2R_indices[i][1]);


  for (i = 0; i < DWC2I; i ++) 
    x[n++] = WC2I(cHl1I, WC2I_indices[i][0], WC2I_indices[i][1]);


  for (i = 0; i < DWC2R; i ++)
    x[n++] = WC2R(cHl3R, WC2R_indices[i][0], WC2R_indices[i][1]);
        

  for (i = 0; i < DWC2I; i ++)
    x[n++] = WC2I(cHl3I, WC2I_indices[i][0], WC2I_indices[i][1]);
        

  for (i = 0; i < DWC2R; i ++) 
    x[n++] = WC2R(cHeR, WC2R_indices[i][0], WC2R_indices[i][1]);
        

  for (i = 0; i < DWC2I; i ++) 
    x[n++] = WC2I(cHeI, WC2I_indices[i][0], WC2I_indices[i][1]);
        

  for (i = 0; i < DWC2R; i ++) 
    x[n++] = WC2R(cHq1R, WC2R_indices[i][0], WC2R_indices[i][1]);
        

  for (i = 0; i < DWC2I; i ++) 
    x[n++] = WC2I(cHq1I, WC2I_indices[i][0], WC2I_indices[i][1]);
        

  for (i = 0; i < DWC2R; i ++) 
    x[n++] = WC2R(cHq3R, WC2R_indices[i][0], WC2R_indices[i][1]);
        

  for (i = 0; i < DWC2I; i ++)
    x[n++] = WC2I(cHq3I, WC2I_indices[i][0], WC2I_indices[i][1]);
        

  for (i = 0; i < DWC2R; i ++) 
    x[n++] = WC2R(cHuR, WC2R_indices[i][0], WC2R_indices[i][1]);
        

  for (i = 0; i < DWC2I; i ++) 
    x[n++] = WC2I(cHuI, WC2I_indices[i][0], WC2I_indices[i][1]);
        

  for (i = 0; i < DWC2R; i ++)
    x[n++] = WC2R(cHdR, WC2R_indices[i][0], WC2R_indices[i][1]);
        

  for (i = 0; i < DWC2I; i ++) 
    x[n++] = WC2I(cHdI, WC2I_indices[i][0], WC2I_indices[i][1]);
        

  for (i = 0; i < NG; i ++) 
    for (j = 0; j < NG; j ++) 
      x[n++] = WC1(cHudR, i, j);

  for (i = 0; i < NG; i ++) 
    for (j = 0; j < NG; j ++) 
      x[n++] = WC1(cHudI, i, j);


  //Class 8_LLLL
  for (b = 0; b < DWC6R; b ++) 
    x[n++] = WC6R(cllR, WC6R_indices[b][0], WC6R_indices[b][1],
		  WC6R_indices[b][2], WC6R_indices[b][3]);
                
  for (b = 0; b < DWC6I; b ++)
    x[n++] = WC6I(cllI, WC6I_indices[b][0], WC6I_indices[b][1],
		  WC6I_indices[b][2], WC6I_indices[b][3]);
                

  for (b = 0; b < DWC6R; b ++)
    x[n++] = WC6R(cqq1R, WC6R_indices[b][0], WC6R_indices[b][1],
		  WC6R_indices[b][2], WC6R_indices[b][3]);
                
  for (b = 0; b < DWC6I; b ++) 
    x[n++] = WC6I(cqq1I, WC6I_indices[b][0], WC6I_indices[b][1],
		  WC6I_indices[b][2], WC6I_indices[b][3]);
                
  for (b = 0; b < DWC6R; b ++) 
    x[n++] = WC6R(cqq3R, WC6R_indices[b][0], WC6R_indices[b][1],
		  WC6R_indices[b][2], WC6R_indices[b][3]);
                
  for (b = 0; b < DWC6I; b ++) 
    x[n++] = WC6I(cqq3I, WC6I_indices[b][0], WC6I_indices[b][1],
		  WC6I_indices[b][2], WC6I_indices[b][3]);
                

  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(clq1R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(clq1I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(clq3R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(clq3I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  //Class 8_RRRR
    
  for (b = 0; b < DWC8R; b ++)
    x[n++] = WC8R(ceeR, WC8R_indices[b][0], WC8R_indices[b][1],
		  WC8R_indices[b][2], WC8R_indices[b][3]);
                
  for (b = 0; b < DWC8I; b ++) 
    x[n++] = WC8I(ceeI, WC8I_indices[b][0], WC8I_indices[b][1],
		  WC8I_indices[b][2], WC8I_indices[b][3]);
                
  for (b = 0; b < DWC6R; b ++) 
    x[n++] = WC6R(cuuR, WC6R_indices[b][0], WC6R_indices[b][1],
		  WC6R_indices[b][2], WC6R_indices[b][3]);
                
  for (b = 0; b < DWC6I; b ++) 
    x[n++] = WC6I(cuuI, WC6I_indices[b][0], WC6I_indices[b][1],
		  WC6I_indices[b][2], WC6I_indices[b][3]);
                
  for (b = 0; b < DWC6R; b ++) 
    x[n++] = WC6R(cddR, WC6R_indices[b][0], WC6R_indices[b][1],
		  WC6R_indices[b][2], WC6R_indices[b][3]);
                
  for (b = 0; b < DWC6I; b ++)
    x[n++] = WC6I(cddI, WC6I_indices[b][0], WC6I_indices[b][1],
		  WC6I_indices[b][2], WC6I_indices[b][3]);
                


  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(ceuR, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(ceuI, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cedR, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cedI, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cud1R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cud1I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cud8R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cud8I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
    

  //Class 8_LLRR
    
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cleR, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cleI, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cluR, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cluI, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cldR, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cldI, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                

  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cqeR, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cqeI, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                

  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cqu1R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cqu1I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                

  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cqu8R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cqu8I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                

  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cqd1R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cqd1I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                

  for (b = 0; b < DWC7R; b ++) 
    x[n++] = WC7R(cqd8R, WC7R_indices[b][0], WC7R_indices[b][1],
		  WC7R_indices[b][2], WC7R_indices[b][3]);
                
  for (b = 0; b < DWC7I; b ++) 
    x[n++] = WC7I(cqd8I, WC7I_indices[b][0], WC7I_indices[b][1],
		  WC7I_indices[b][2], WC7I_indices[b][3]);
                
    



  //Class 8_LRRL
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      for (k = 0; k < NG; k ++) {
	for (l = 0; l < NG; l ++) {
	  x[n] = WC5(cledqR, i, j, k, l);
	  x[n++ + NG * NG * NG * NG] = WC5(cledqI, i, j, k, l);
	}
      }
    }
  }
  n += NG * NG * NG*NG;
  //Class 8_LRLR
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      for (k = 0; k < NG; k ++) {
	for (l = 0; l < NG; l ++) {
	  a = 0;
	  x[n + a++ * NG * NG * NG * NG] = WC5(cquqd1R, i, j, k, l);
	  x[n + a++ * NG * NG * NG * NG] = WC5(cquqd1I, i, j, k, l);
	  x[n + a++ * NG * NG * NG * NG] = WC5(cquqd8R, i, j, k, l);
	  x[n + a++ * NG * NG * NG * NG] = WC5(cquqd8I, i, j, k, l);
	  x[n + a++ * NG * NG * NG * NG] = WC5(clequ1R, i, j, k, l);
	  x[n + a++ * NG * NG * NG * NG] = WC5(clequ1I, i, j, k, l);
	  x[n + a++ * NG * NG * NG * NG] = WC5(clequ3R, i, j, k, l);
	  x[n++ + a++ * NG * NG * NG * NG] = WC5(clequ3I, i, j, k, l);
	}
      }
    }
  }
  n += NG * NG * NG * NG * (2 * N8_LRLR - 1);
}

void RGESolver::Update() {
  //Here we return to original structures

  int n = 0;
  int a, i, j, k, l;
 
  g2 = x[0];
  g1 = x[1];
  g3 = x[2];
  n += 3;
  lambda = x[n++];
  mh2 = x[n++];

  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      a = 0;
      yuR(i, j) = x[n + a++ * DF];
      yuI(i, j) = x[n + a++ * DF];
      ydR(i, j) = x[n + a++ * DF];
      ydI(i, j) = x[n + a++ * DF];
      yeR(i, j) = x[n + a++ * DF];
      yeI(i, j) = x[n++ + a++ * DF];
    }
  }
  n += (2. * Nyukawa - 1.) * DF;

  cG = x[n++];
  cGT = x[n++];
  cW = x[n++];
  cWT = x[n++];
  cH = x[n++];
  cHBOX = x[n++];
  cHD = x[n++];

  cHG = x[n++];
  cHB = x[n++];
  cHW = x[n++];
  cHWB = x[n++];
  cHGT = x[n++];
  cHBT = x[n++];
  cHWT = x[n++];
  cHWBT = x[n++];

  //class 5
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      a = 0;
      WC1_set(cuHR, i, j, x[n + 2 * a * DF]);
      WC1_set(cuHI, i, j, x[n + (2 * a++ + 1) * DF]);
      WC1_set(cdHR, i, j, x[n + 2 * a * DF]);
      WC1_set(cdHI, i, j, x[n + (2 * a++ + 1) * DF]);
      WC1_set(ceHR, i, j, x[n + 2 * a * DF]);
      WC1_set(ceHI, i, j, x[n++ + (2 * a++ + 1) * DF]);
    }
  }
  n += (N5 * 2 - 1) * DF;


  //Class 6
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      a = 0;
      WC1_set(ceWR, i, j, x[n + a++ * DF]);
      WC1_set(ceWI, i, j, x[n + a++ * DF]);
      WC1_set(ceBR, i, j, x[n + a++ * DF]);
      WC1_set(ceBI, i, j, x[n + a++ * DF]);
      WC1_set(cuGR, i, j, x[n + a++ * DF]);
      WC1_set(cuGI, i, j, x[n + a++ * DF]);
      WC1_set(cuWR, i, j, x[n + a++ * DF]);
      WC1_set(cuWI, i, j, x[n + a++ * DF]);
      WC1_set(cuBR, i, j, x[n + a++ * DF]);
      WC1_set(cuBI, i, j, x[n + a++ * DF]);
      WC1_set(cdGR, i, j, x[n + a++ * DF]);
      WC1_set(cdGI, i, j, x[n + a++ * DF]);
      WC1_set(cdWR, i, j, x[n + a++ * DF]);
      WC1_set(cdWI, i, j, x[n + a++ * DF]);
      WC1_set(cdBR, i, j, x[n + a++ * DF]);
      WC1_set(cdBI, i, j, x[n++ + a++ * DF]);
    }
  }
  n += (N6 * 2 - 1) * DF;

  //class 7
    
  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHl1R, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHl1I, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHl3R, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHl3I, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHeR, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHeI, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                

  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHq1R, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHq1I, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHq3R, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHq3I, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                

  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHuR, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHuI, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                

  for (i = 0; i < DWC2R; i ++) 
    WC2R_set(cHdR, WC2R_indices[i][0], WC2R_indices[i][1], x[n++]);
                
  for (i = 0; i < DWC2I; i ++) 
    WC2I_set(cHdI, WC2I_indices[i][0], WC2I_indices[i][1], x[n++]);
                

  for (i = 0; i < NG; i ++) 
    for (j = 0; j < NG; j ++) 
      WC1_set(cHudR, i, j, x[n++]);
                        
        
  for (i = 0; i < NG; i ++) 
    for (j = 0; j < NG; j ++) 
      WC1_set(cHudI, i, j, x[n++]);
                        
        


    
  //class 8_LLLL
    
  for (a = 0; a < DWC6R; a ++) 
    WC6R_set(cllR, WC6R_indices[a][0], WC6R_indices[a][1],
	     WC6R_indices[a][2], WC6R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6I; a ++) 
    WC6I_set(cllI, WC6I_indices[a][0], WC6I_indices[a][1],
	     WC6I_indices[a][2], WC6I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6R; a ++) 
    WC6R_set(cqq1R, WC6R_indices[a][0], WC6R_indices[a][1],
	     WC6R_indices[a][2], WC6R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6I; a ++) 
    WC6I_set(cqq1I, WC6I_indices[a][0], WC6I_indices[a][1],
	     WC6I_indices[a][2], WC6I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6R; a ++) 
    WC6R_set(cqq3R, WC6R_indices[a][0], WC6R_indices[a][1],
	     WC6R_indices[a][2], WC6R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6I; a ++) 
    WC6I_set(cqq3I, WC6I_indices[a][0], WC6I_indices[a][1],
	     WC6I_indices[a][2], WC6I_indices[a][3], x[n++]);
                

  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(clq1R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(clq1I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(clq3R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(clq3I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                

    
  //Class 8_RRRR
    
  for (a = 0; a < DWC8R; a ++) 
    WC8R_set(ceeR, WC8R_indices[a][0], WC8R_indices[a][1],
	     WC8R_indices[a][2], WC8R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC8I; a ++) 
    WC8I_set(ceeI, WC8I_indices[a][0], WC8I_indices[a][1],
	     WC8I_indices[a][2], WC8I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6R; a ++) 
    WC6R_set(cuuR, WC6R_indices[a][0], WC6R_indices[a][1],
	     WC6R_indices[a][2], WC6R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6I; a ++) 
    WC6I_set(cuuI, WC6I_indices[a][0], WC6I_indices[a][1],
	     WC6I_indices[a][2], WC6I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6R; a ++) 
    WC6R_set(cddR, WC6R_indices[a][0], WC6R_indices[a][1],
	     WC6R_indices[a][2], WC6R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC6I; a ++) 
    WC6I_set(cddI, WC6I_indices[a][0], WC6I_indices[a][1],
	     WC6I_indices[a][2], WC6I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(ceuR, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(ceuI, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                


  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cedR, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cedI, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                

  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cud1R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cud1I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cud8R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cud8I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                


    

  //Class 8_LLRR
    
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cleR, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cleI, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cluR, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cluI, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                

  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cldR, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cldI, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cqeR, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cqeI, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                

  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cqu1R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cqu1I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cqu8R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cqu8I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                



  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cqd1R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cqd1I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7R; a ++) 
    WC7R_set(cqd8R, WC7R_indices[a][0], WC7R_indices[a][1],
	     WC7R_indices[a][2], WC7R_indices[a][3], x[n++]);
                
  for (a = 0; a < DWC7I; a ++) 
    WC7I_set(cqd8I, WC7I_indices[a][0], WC7I_indices[a][1],
	     WC7I_indices[a][2], WC7I_indices[a][3], x[n++]);
                

    

  //Class 8_LRRL

  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      for (k = 0; k < NG; k ++) {
	for (l = 0; l < NG; l ++) {
	  WC5_set(cledqR, i, j, k, l, x[n]);
	  WC5_set(cledqI, i, j, k, l, x[n++ + NG * NG * NG * NG]);
	}
      }
    }
  }

  n += DF*DF;

  //Class 8_LRLR
  for (i = 0; i < NG; i ++) {
    for (j = 0; j < NG; j ++) {
      for (k = 0; k < NG; k ++) {
	for (l = 0; l < NG; l ++) {
	  a = 0;
	  WC5_set(cquqd1R, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(cquqd1I, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(cquqd8R, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(cquqd8I, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(clequ1R, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(clequ1I, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(clequ3R, i, j, k, l, x[n + a++ * NG * NG * NG * NG]);
	  WC5_set(clequ3I, i, j, k, l, x[n++ + a++ * NG * NG * NG * NG]);
	}
      }
    }
  }
  n += NG * NG * NG * NG * (2 * N8_LRLR - 1);

}

void RGESolver::Evolve(std::string method, double muI, double muF) {


  if (method != "Numeric" && method != "Approximate") {
    std::cout << "WARNING : invalid method\n"
      "Available methods: Numeric, Approximate"
	      << std::endl;
  }
  if (muF != muI) {
    Init();
    //Initial conditions are inserted 
    //in the array x


    //Numeric solution 


    //Driver
    /* if (method == "Numeric") {
       double tI = log(muI);
       double tF = log(muF);
       double ttmp = tI;

       //Initial step 
       double stepIn = (tF - tI)*0.01;
       //std::cout << "step : " << step_ << std::endl;
       double er = 0.001;
       double ea = 0.0000000000000000000000001;

       //std::cout << "ok before driver reset" << std::endl;

       gsl_odeiv2_driver_reset_hstart(d, stepIn);
       //std::cout << "ok after driver reset" << std::endl;

       int status = gsl_odeiv2_driver_apply(d, &ttmp, tF, x);
       //std::cout << "ok after evolution" << std::endl;

       if (status != GSL_SUCCESS) {
       printf("error, return value=%d\n", status);
       //break;
       }
       }*/
    //Low level 
    if (method == "Numeric") {
      double tI = log(muI);
      double tF = log(muF);
      double ttmp = tI;

      gsl_odeiv2_control * c
	= gsl_odeiv2_control_y_new(epsabs_, epsrel_);

      //Initial step 
      double stepIn = (tF - tI)*0.01;
      if (tF > tI) {
	while (ttmp < tF) {
	  int status = gsl_odeiv2_evolve_apply(e, c, s,
					       &sys, &ttmp, tF, &stepIn, x);
	  if (status != GSL_SUCCESS) {
	    printf("error, return value=%d\n", status);
	    break;
	  }
	}
      }
      if (tF < tI) {
	while (ttmp > tF) {
	  int status = gsl_odeiv2_evolve_apply(e, c, s,
					       &sys, &ttmp, tF, &stepIn, x);
	  if (status != GSL_SUCCESS) {
	    printf("error, return value=%d\n", status);
	    break;
	  }
	}
      }

      gsl_odeiv2_evolve_reset(e);
      gsl_odeiv2_step_reset(s);
      gsl_odeiv2_control_free(c);
    }



    //Approximate solution 
    //-------------------------------------------
    if (method == "Approximate") {
      double beta[dim] = {0.};
      double Log_muF_over_muI = log(muF / muI);
      /*int status = */func(10., x, beta, NULL);
      for (int i = 0; i < dim; i ++) {
	x[i] += beta[i] * Log_muF_over_muI;
      }
    }
    //-------------------------------------------

    Update();
    //Evolved values from x are put 
    //back in the coefficients 
  }
}

void RGESolver::SaveOutputFile(std::string filename,
			       std::string format) {
  using namespace std;

  if (format != "SLHA") {
    cout << "ERROR : FORMAT " << format << " NOT AVAILABLE" << endl;
    //cout << "Use \"SLHA\" or \"JSON\" " << endl;
    return;
  }

  string space = " ";
  ofstream outf;
  outf.open(filename);
  int i, j;
  if (format == "SLHA") {
    //SM print 

    int n = 1;
    outf << "Block Standard Model" << endl;
    outf << n++ << " " << g1 << "\t# g1" << endl;
    outf << n++ << " " << g2 << "\t# g2" << endl;
    outf << n++ << " " << g3 << "\t# g3" << endl;
    outf << n++ << " " << lambda << "\t# lambda" << endl;
    outf << n++ << " " << mh2 << "\t# mh2 [GeV^2]" << endl;

    outf << "Block Re[Yu]" << endl;
    for (i = 0; i < 3; i ++) {
      for (j = 0; j < 3; j ++) {
	outf << i + 1 << " " << j + 1
	     << " " << yuR(i, j) << "\t"
	     << "# Re[Yu]"
	     << "(" << i + 1 << ","
	     << j + 1 << ")" << endl;
      }
    }
    outf << "Block Im[Yu]" << endl;
    for (i = 0; i < 3; i ++) {
      for (j = 0; j < 3; j ++) {
	outf << i + 1 << " " << j + 1
	     << " " << yuI(i, j) << "\t"
	     << "# Im[Yu]"
	     << "(" << i + 1 << ","
	     << j + 1 << ")" << endl;
      }
    }
    outf << "Block Re[Yd]" << endl;
    for (i = 0; i < 3; i ++) {
      for (j = 0; j < 3; j ++) {
	outf << i + 1 << " " << j + 1
	     << " " << ydR(i, j) << "\t"
	     << "# Re[Yd]"
	     << "(" << i + 1 << ","
	     << j + 1 << ")" << endl;
      }
    }
    outf << "Block Im[Yd]" << endl;
    for (i = 0; i < 3; i ++) {
      for (j = 0; j < 3; j ++) {
	outf << i + 1 << " " << j + 1
	     << " " << ydI(i, j) << "\t"
	     << "# Im[Yd]"
	     << "(" << i + 1 << ","
	     << j + 1 << ")" << endl;
      }
    }
    outf << "Block Re[Ye]" << endl;
    for (i = 0; i < 3; i ++) {
      for (j = 0; j < 3; j ++) {
	outf << i + 1 << " " << j + 1
	     << " " << yeR(i, j) << "\t"
	     << "# Re[Ye]"
	     << "(" << i + 1 << ","
	     << j + 1 << ")" << endl;
      }
    }
    outf << "Block Im[Ye]" << endl;
    for (i = 0; i < 3; i ++) {
      for (j = 0; j < 3; j ++) {
	outf << i + 1 << " " << j + 1
	     << " " << yeI(i, j) << "\t"
	     << "# Im[Ye]"
	     << "(" << i + 1 << ","
	     << j + 1 << ")" << endl;
      }
    }


    //Scalar SMEFT coefficients 

    outf << "Block SMEFT Class 1" << endl;
    n = 1;
    outf << n++ << " " << cG << "\t#CG" << endl;
    outf << n++ << " " << cGT << "\t#CGtilde" << endl;
    outf << n++ << " " << cW << "\t#CW" << endl;
    outf << n++ << " " << cWT << "\t#CWtilde" << endl;
    
    outf << "Block SMEFT Class 2" << endl;
    n = 1;
    outf << n++ << " " << cH << "\t#CH" << endl;
    
    outf << "Block SMEFT Class 3" << endl;
    n = 1;
    outf << n++ << " " << cHBOX << "\t#CHbox" << endl;
    outf << n++ << " " << cHD << "\t#CHD" << endl;
    
    outf << "Block SMEFT Class 4" << endl;
    n = 1;
    outf << n++ << " " << cHG << "\t#CHG" << endl;
    outf << n++ << " " << cHGT << "\t#CHGtilde" << endl;
    outf << n++ << " " << cHW << "\t#CHW" << endl;
    outf << n++ << " " << cHWT << "\t#CHWtilde" << endl;
    outf << n++ << " " << cHB << "\t#CHB" << endl;
    outf << n++ << " " << cHBT << "\t#CHBtilde" << endl;
    outf << n++ << " " << cHWB << "\t#CHWB" << endl;
    outf << n++ << " " << cHWBT << "\t#CHWtildeB" << endl;
    

    //Class 5
    Print(cuHR, "Re[CuH]", "WC1", format, outf);
    Print(cuHI, "Im[CuH]", "WC1", format, outf);
    Print(cdHR, "Re[CdH]", "WC1", format, outf);
    Print(cdHI, "Im[CdH]", "WC1", format, outf);
    Print(ceHR, "Re[CeH]", "WC1", format, outf);
    Print(ceHI, "Im[CeH]", "WC1", format, outf);

    //Class 6
    Print(ceWR, "Re[CeW]", "WC1", format, outf);
    Print(ceWI, "Im[CeW]", "WC1", format, outf);
    Print(ceBR, "Re[CeB]", "WC1", format, outf);
    Print(ceBI, "Im[CeB]", "WC1", format, outf);

    Print(cuGR, "Re[CuG]", "WC1", format, outf);
    Print(cuGI, "Im[CuG]", "WC1", format, outf);
    Print(cuWR, "Re[CuW]", "WC1", format, outf);
    Print(cuWI, "Im[CuW]", "WC1", format, outf);
    Print(cuBR, "Re[CuB]", "WC1", format, outf);
    Print(cuBI, "Im[CuB]", "WC1", format, outf);

    Print(cdGR, "Re[CdG]", "WC1", format, outf);
    Print(cdGI, "Im[CdG]", "WC1", format, outf);
    Print(cdWR, "Re[CdW]", "WC1", format, outf);
    Print(cdWI, "Im[CdW]", "WC1", format, outf);
    Print(cdBR, "Re[CdB]", "WC1", format, outf);
    Print(cdBI, "Im[CdB]", "WC1", format, outf);

    //Class 7
    Print(cHl1R, "Re[CHl1]", "WC2R", format, outf);
    Print(cHl1I, "Im[CHl1]", "WC2I", format, outf);
    Print(cHl3R, "Re[CHl3]", "WC2R", format, outf);
    Print(cHl3I, "Im[CHl3]", "WC2I", format, outf);
    Print(cHeR, "Re[CHe]", "WC2R", format, outf);
    Print(cHeI, "Im[CHe]", "WC2I", format, outf);
    Print(cHq1R, "Re[CHq1]", "WC2R", format, outf);
    Print(cHq1I, "Im[CHq1]", "WC2I", format, outf);
    Print(cHq3R, "Re[CHq3]", "WC2R", format, outf);
    Print(cHq3I, "Im[CHq3]", "WC2I", format, outf);
    Print(cHuR, "Re[CHu]", "WC2R", format, outf);
    Print(cHuI, "Im[CHu]", "WC2I", format, outf);
    Print(cHdR, "Re[CHd]", "WC2R", format, outf);
    Print(cHdI, "Im[CHd]", "WC2I", format, outf);
    Print(cHudR, "Re[CHud]", "WC1", format, outf);
    Print(cHudI, "Im[CHud]", "WC1", format, outf);


    //Class 8LLLL
    Print(cllR, "Re[Cll]", "WC6R", format, outf);
    Print(cllI, "Im[Cll]", "WC6I", format, outf);
    Print(cqq1R, "Re[Cqq1]", "WC6R", format, outf);
    Print(cqq1I, "Im[Cqq1]", "WC6I", format, outf);
    Print(cqq3R, "Re[Cqq3]", "WC6R", format, outf);
    Print(cqq3I, "Im[Cqq3]", "WC6I", format, outf);
    Print(clq1R, "Re[Clq1]", "WC7R", format, outf);
    Print(clq1I, "Im[Clq1]", "WC7I", format, outf);
    Print(clq3R, "Re[Clq3]", "WC7R", format, outf);
    Print(clq3I, "Im[Clq3]", "WC7I", format, outf);

    //Class 8RRRR
    Print(ceeR, "Re[Cee]", "WC8R", format, outf);
    Print(ceeI, "Im[Cee]", "WC8I", format, outf);
    Print(cuuR, "Re[Cuu]", "WC6R", format, outf);
    Print(cuuI, "Im[Cuu]", "WC6I", format, outf);
    Print(cddR, "Re[Cdd]", "WC6R", format, outf);
    Print(cddI, "Im[Cdd]", "WC6I", format, outf);
    Print(ceuR, "Re[Ceu]", "WC7R", format, outf);
    Print(ceuI, "Im[Ceu]", "WC7I", format, outf);
    Print(cedR, "Re[Ced]", "WC7R", format, outf);
    Print(cedI, "Im[Ced]", "WC7I", format, outf);
    Print(cud1R, "Re[Cud1]", "WC7R", format, outf);
    Print(cud1I, "Im[Cud1]", "WC7I", format, outf);
    Print(cud8R, "Re[Cud8]", "WC7R", format, outf);
    Print(cud8I, "Im[Cud8]", "WC7I", format, outf);

    //Class 8LLRR 
    Print(cleR, "Re[Cle]", "WC7R", format, outf);
    Print(cleI, "Im[Cle]", "WC7I", format, outf);
    Print(cluR, "Re[Clu]", "WC7R", format, outf);
    Print(cluI, "Im[Clu]", "WC7I", format, outf);
    Print(cldR, "Re[Cld]", "WC7R", format, outf);
    Print(cldI, "Im[Cld]", "WC7I", format, outf);
    Print(cqeR, "Re[Cqe]", "WC7R", format, outf);
    Print(cqeI, "Im[Cqe]", "WC7I", format, outf);
    Print(cqu1R, "Re[Cqu1]", "WC7R", format, outf);
    Print(cqu1I, "Im[Cqu1]", "WC7I", format, outf);
    Print(cqu8R, "Re[Cqu8]", "WC7R", format, outf);
    Print(cqu8I, "Im[Cqu8]", "WC7I", format, outf);
    Print(cqd1R, "Re[Cqd1]", "WC7R", format, outf);
    Print(cqd1I, "Im[Cqd1]", "WC7I", format, outf);
    Print(cqd8R, "Re[Cqd8]", "WC7R", format, outf);
    Print(cqd8I, "Im[Cqd8]", "WC7I", format, outf);


    //Class 8LRRL
    Print(cledqR, "Re[Cledq]", "WC5", format, outf);
    Print(cledqI, "Im[Cledq]", "WC5", format, outf);

    //Class 8LRLR
    Print(cquqd1R, "Re[Cquqd1]", "WC5", format, outf);
    Print(cquqd1I, "Im[Cquqd1]", "WC5", format, outf);
    Print(cquqd8R, "Re[Cquqd8]", "WC5", format, outf);
    Print(cquqd8I, "Im[Cquqd8]", "WC5", format, outf);
    Print(clequ1R, "Re[Clequ1]", "WC5", format, outf);
    Print(clequ1I, "Im[Clequ1]", "WC5", format, outf);
    Print(clequ3R, "Re[Clequ3]", "WC5", format, outf);
    Print(clequ3I, "Im[Clequ3]", "WC5", format, outf);
  }

  outf.close();

}

void RGESolver::Print(double* c, std::string name,
		      std::string sym,
		      std::string format,
		      std::ofstream & f) {
  using namespace std;

  int i, j, k, l;

  if (format == "SLHA") {
    f << "Block " << name << endl;
    if (sym == "WC1") {
      for (i = 0; i < 3; i ++) {
	for (j = 0; j < 3; j ++) {
	  f << i + 1 << " " << j + 1
	    << " " << WC1(c, i, j) << "\t"
	    << "# " << name
	    << "(" << i + 1 << ","
	    << j + 1 << ")" << endl;
	}
      }
    }

    if (sym == "WC2R") {
      for (int n = 0; n < DWC2R; n ++) {
	i = WC2R_indices[n][0];
	j = WC2R_indices[n][1];
	f << i + 1 << " " << j + 1
	  << " " << WC2R(c, i, j) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << ")" << endl;
      }
    }

    if (sym == "WC2I") {
      for (int n = 0; n < DWC2I; n ++) {
	i = WC2I_indices[n][0];
	j = WC2I_indices[n][1];
	f << i + 1 << " " << j + 1
	  << " " << WC2I(c, i, j) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << ")" << endl;
      }
    }

    if (sym == "WC5") {
      for (i = 0; i < 3; i ++) {
	for (j = 0; j < 3; j ++) {
	  for (k = 0; k < 3; k ++) {
	    for (l = 0; l < 3; l ++) {
	      f << i + 1 << " " << j + 1
		<< " " << k + 1 << " " << l + 1
		<< " " << WC5(c, i, j, k, l) << "\t"
		<< "# " << name
		<< "(" << i + 1 << ","
		<< j + 1 << "," << k + 1 << ","
		<< l + 1 << ")" << endl;
	    }
	  }
	}
      }
    }

    if (sym == "WC6R") {
      for (int n = 0; n < DWC6R; n ++) {
	i = WC6R_indices[n][0];
	j = WC6R_indices[n][1];
	k = WC6R_indices[n][2];
	l = WC6R_indices[n][3];
	f << i + 1 << " " << j + 1
	  << " " << k + 1 << " " << l + 1
	  << " " << WC6R(c, i, j, k, l) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << "," << k + 1 << ","
	  << l + 1 << ")" << endl;
      }
    }

    if (sym == "WC6I") {
      for (int n = 0; n < DWC6I; n ++) {
	i = WC6I_indices[n][0];
	j = WC6I_indices[n][1];
	k = WC6I_indices[n][2];
	l = WC6I_indices[n][3];
	f << i + 1 << " " << j + 1
	  << " " << k + 1 << " " << l + 1
	  << " " << WC6I(c, i, j, k, l) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << "," << k + 1 << ","
	  << l + 1 << ")" << endl;
      }
    }

    if (sym == "WC7R") {
      for (int n = 0; n < DWC7R; n ++) {
	i = WC7R_indices[n][0];
	j = WC7R_indices[n][1];
	k = WC7R_indices[n][2];
	l = WC7R_indices[n][3];
	f << i + 1 << " " << j + 1
	  << " " << k + 1 << " " << l + 1
	  << " " << WC7R(c, i, j, k, l) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << "," << k + 1 << ","
	  << l + 1 << ")" << endl;
      }
    }

    if (sym == "WC7I") {
      for (int n = 0; n < DWC7I; n ++) {
	i = WC7I_indices[n][0];
	j = WC7I_indices[n][1];
	k = WC7I_indices[n][2];
	l = WC7I_indices[n][3];
	f << i + 1 << " " << j + 1
	  << " " << k + 1 << " " << l + 1
	  << " " << WC7I(c, i, j, k, l) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << "," << k + 1 << ","
	  << l + 1 << ")" << endl;
      }
    }

    if (sym == "WC8R") {
      for (int n = 0; n < DWC8R; n ++) {
	i = WC8R_indices[n][0];
	j = WC8R_indices[n][1];
	k = WC8R_indices[n][2];
	l = WC8R_indices[n][3];
	f << i + 1 << " " << j + 1
	  << " " << k + 1 << " " << l + 1
	  << " " << WC8R(c, i, j, k, l) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << "," << k + 1 << ","
	  << l + 1 << ")" << endl;
      }
    }

    if (sym == "WC8I") {
      for (int n = 0; n < DWC8I; n ++) {
	i = WC8I_indices[n][0];
	j = WC8I_indices[n][1];
	k = WC8I_indices[n][2];
	l = WC8I_indices[n][3];
	f << i + 1 << " " << j + 1
	  << " " << k + 1 << " " << l + 1
	  << " " << WC8I(c, i, j, k, l) << "\t"
	  << "# " << name
	  << "(" << i + 1 << ","
	  << j + 1 << "," << k + 1 << ","
	  << l + 1 << ")" << endl;
      }
    }
  }

  if (format == "JSON") {

  }

}


//Setter for 0F

void RGESolver::SetCoefficient(std::string name, double val) {
  *(Operators0F.at(name)) = val;
}

//Getter for 0F

double RGESolver::GetCoefficient(std::string name) {
  return * (Operators0F.at(name));
}



//Setter for 2F

void RGESolver::SetCoefficient(std::string name, double val,
			       int i, int j) {
  if ((i <= 2)&&(i >= 0)&&(j <= 2)&&(j >= 0)) {

    Setter2F.at(name)(i, j, val);

  } else {
    std::cout << "ERROR : INDICES OUT OF RANGE [0:2]"
	      << std::endl;
  }
}

//Getter for 2F

double RGESolver::GetCoefficient(std::string name,
				 int i, int j) {

  double val = 0.;
  if ((i <= 2)&&(i >= 0)&&(j <= 2)&&(j >= 0)) {

    val = Getter2F.at(name)(i, j);

  } else {
    std::cout << "ERROR : INDICES OUT OF RANGE [0:2]"
	      << std::endl;
  }
  return val;
}

//Setter for 4F

void RGESolver::SetCoefficient(std::string name, double val,
			       int i, int j, int k, int l) {
  if ((i <= 2)&&(i >= 0)&&(j <= 2)&&(j >= 0)
      &&(k <= 2)&&(k >= 0)&&(l <= 2)&&(l >= 0)) {


    Setter4F.at(name)(i, j, k, l, val);

  } else {
    std::cout << "ERROR : INDICES OUT OF RANGE [0:2]"
	      << std::endl;
  }
}
//Setter for 4F

double RGESolver::GetCoefficient(std::string name,
				 int i, int j, int k, int l) {
  double val = 0.;
  if ((i <= 2)&&(i >= 0)&&(j <= 2)&&(j >= 0)
      &&(k <= 2)&&(k >= 0)&&(l <= 2)&&(l >= 0)) {


    val = Getter4F.at(name)(i, j, k, l);

  } else {
    std::cout << "ERROR : INDICES OUT OF RANGE [0:2]"
	      << std::endl;
  }
  return val;
}



#include "BetaFunction.cpp"



