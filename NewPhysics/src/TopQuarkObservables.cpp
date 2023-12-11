#include "TopQuarkObservables.h"
#include "std_make_vector.h"


TopQuarkObservables::TopQuarkObservables(const NPSMEFTd6General& NP_i) : NP(NP_i) {};


//// ALL THE UNITS OF THE WILSON COEFFICIENTS ARE IN TeV^(-2). THIS FACTOR NEEDS TO BE CONSIDERED => NEEDS TO BE CHECKED CAREFULLY!!! 




////// ttbar differential cross section ///////////////////////////////////

sigma_tt_diff_LO::sigma_tt_diff_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tt_bin_250_400" << "SM_sigma_tt_bin_400_480" <<
            "SM_sigma_tt_bin_480_560" << "SM_sigma_tt_bin_560_640" << "SM_sigma_tt_bin_640_720" << "SM_sigma_tt_bin_720_800" <<
            "SM_sigma_tt_bin_800_900" << "SM_sigma_tt_bin_900_1000" << "SM_sigma_tt_bin_1000_1150" << "SM_sigma_tt_bin_1150_1300" <<
            "SM_sigma_tt_bin_1300_1500" << "SM_sigma_tt_bin_1500_1700" << "SM_sigma_tt_bin_1700_2000" << "SM_sigma_tt_bin_2000_2300" <<
            "SM_sigma_tt_bin_2300_3500");
            //"SM_sigma_tt_bin_2300_2600" << "SM_sigma_tt_bin_2600_3000" << "SM_sigma_tt_bin_3000_3500" << "SM_sigma_tt_bin_3500_4000");
    
};

double sigma_tt_diff_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    //double C_tG = ewgc("CuGR",2,2);

    
    
    
    if(b_min == 250 && b_max == 400){
        
        double SM_sigma_tt_bin_250_400 = SM.getOptionalParameter("SM_sigma_tt_bin_250_400");
        //double sigma_tt_bin_250_400_madgraph = 105600.0;//fb maybe over the bin width? Check!
        double sigma_tt_bin_250_400_madgraph = 78.586100;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return  (sigma_tt_bin_250_400_madgraph + 0.070920*ewgc("Cqd8R",2,2,2,2)+0.398652*ewgc("Cqd8R",2,2,0,0)+0.393474*ewgc("Cud8R",2,2,0,0)+0.768666*ewgc("Cqq3R",1,2,2,1)
                    +0.243408*ewgc("Cqq1R",1,2,2,1)+0.884946*ewgc("Cqu8R",0,0,2,2)+0.568908*ewgc("Cqu8R",2,2,0,0)-25.503700*ewgc("CuGR",2,2)+0.180606*ewgc("Cqd8R",2,2,1,1)
                    +2.084240*ewgc("Cqq1R",0,2,2,0)+4.638680*ewgc("Cqq3R",0,2,2,0)+0.122778*ewgc("Cqu8R",2,2,1,1)+0.070920*ewgc("Cud8R",2,2,2,2)+2.055090*ewgc("CuuR",0,2,2,0)
                    +0.185274*ewgc("Cud8R",2,2,1,1)+0.238692*ewgc("CuuR",1,2,2,1)-2.030250*ewgc("CG")+0.191568*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_250_400/sigma_tt_bin_250_400_madgraph);
            
        }
        
    } else if(b_min == 400 && b_max == 480){
        
        double SM_sigma_tt_bin_400_480 = SM.getOptionalParameter("SM_sigma_tt_bin_400_480");
        double sigma_tt_bin_400_480_madgraph = 108.694000;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_400_480_madgraph -0.030486*ewgc("Cqd8R", 2, 2, 2, 2) + 0.489354 * ewgc("Cqd8R", 2, 2, 0, 0) + 0.471660 * ewgc("Cud8R", 2, 2, 0, 0)
                    + 0.937392 * ewgc("Cqq3R", 1, 2, 2, 1) + 0.147936 * ewgc("Cqq1R", 1, 2, 2, 1) + 1.278050 * ewgc("Cqu8R", 0, 0, 2, 2) + 0.773802 * ewgc("Cqu8R", 2, 2, 0, 0)
                    -31.928800 * ewgc("CuGR", 2, 2) + 0.083934 * ewgc("Cqd8R", 2, 2, 1, 1) + 3.277280 * ewgc("Cqq1R", 0, 2, 2, 0) + 7.079870 * ewgc("Cqq3R", 0, 2, 2, 0) 
                    + 0.020028 * ewgc("Cqu8R", 2, 2, 1, 1)-0.030486 * ewgc("Cud8R", 2, 2, 2, 2) + 3.207280 * ewgc("CuuR", 0, 2, 2, 0) + 0.090372 * ewgc("Cud8R", 2, 2, 1, 1) 
                    + 0.141108 * ewgc("CuuR", 1, 2, 2, 1)-6.256640 * ewgc("CG") + 0.113838 * ewgc("Cqu8R", 1, 1, 2, 2))*(SM_sigma_tt_bin_400_480/sigma_tt_bin_400_480_madgraph);
            
        }
        
    } else if(b_min == 480 && b_max == 560){
        
        double SM_sigma_tt_bin_480_560 = SM.getOptionalParameter("SM_sigma_tt_bin_480_560");
        double sigma_tt_bin_480_560_madgraph = 64.345700;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_480_560_madgraph+0.010313*ewgc("Cqd8R",2,2,2,2)+0.376856*ewgc("Cqd8R",2,2,0,0)+0.377929*ewgc("Cud8R",2,2,0,0)
                    +0.645970*ewgc("Cqq3R",1,2,2,1)+0.137390*ewgc("Cqq1R",1,2,2,1)+0.932173*ewgc("Cqu8R",0,0,2,2)+0.578978*ewgc("Cqu8R",2,2,0,0)
                    -18.013300*ewgc("CuGR",2,2)+0.079471*ewgc("Cqd8R",2,2,1,1)+2.430400*ewgc("Cqq1R",0,2,2,0)+5.120730*ewgc("Cqq3R",0,2,2,0)
                    +0.037264*ewgc("Cqu8R",2,2,1,1)+0.010313*ewgc("Cud8R",2,2,2,2)+2.344900*ewgc("CuuR",0,2,2,0)+0.077563*ewgc("Cud8R",2,2,1,1)
                    +0.129405*ewgc("CuuR",1,2,2,1)-5.150190*ewgc("CG")+0.107117*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_480_560/sigma_tt_bin_480_560_madgraph);
            
        }
        
    } else if(b_min == 560 && b_max == 640){
        
        double SM_sigma_tt_bin_560_640 = SM.getOptionalParameter("SM_sigma_tt_bin_560_640");
        double sigma_tt_bin_560_640_madgraph = 35.478800;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_560_640_madgraph-0.013071*ewgc("Cqd8R",2,2,2,2)+0.251126*ewgc("Cqd8R",2,2,0,0)+0.249005*ewgc("Cud8R",2,2,0,0)
                    +0.384914*ewgc("Cqq3R",1,2,2,1)+0.074809*ewgc("Cqq1R",1,2,2,1)+0.665144*ewgc("Cqu8R",0,0,2,2)+0.396944*ewgc("Cqu8R",2,2,0,0)
                    -9.816180*ewgc("CuGR",2,2)+0.030503*ewgc("Cqd8R",2,2,1,1)+1.711280*ewgc("Cqq1R",0,2,2,0)+3.778140*ewgc("Cqq3R",0,2,2,0)
                    +0.008984*ewgc("Cqu8R",2,2,1,1)+-0.013071*ewgc("Cud8R",2,2,2,2)+1.668460*ewgc("CuuR",0,2,2,0)+0.030811*ewgc("Cud8R",2,2,1,1)
                    +0.070984*ewgc("CuuR",1,2,2,1)-3.307920*ewgc("CG")+0.051151*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_560_640/sigma_tt_bin_560_640_madgraph);
            
        }
        
    } else if(b_min == 640 && b_max == 720){
        
        double SM_sigma_tt_bin_640_720 = SM.getOptionalParameter("SM_sigma_tt_bin_640_720");
        double sigma_tt_bin_640_720_madgraph = 19.766700;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_640_720_madgraph+0.005380*ewgc("Cqd8R",2,2,2,2)+0.188524*ewgc("Cqd8R",2,2,0,0)+0.188168*ewgc("Cud8R",2,2,0,0)
                    +0.268617*ewgc("Cqq3R",1,2,2,1)+0.054757*ewgc("Cqq1R",1,2,2,1)+0.486146*ewgc("Cqu8R",0,0,2,2)+0.307089*ewgc("Cqu8R",2,2,0,0)
                    -5.471020*ewgc("CuGR",2,2)+0.035139*ewgc("Cqd8R",2,2,1,1)+1.271600*ewgc("Cqq1R",0,2,2,0)+2.593230*ewgc("Cqq3R",0,2,2,0)
                    +0.016561*ewgc("Cqu8R",2,2,1,1)+0.005380*ewgc("Cud8R",2,2,2,2)+1.245390*ewgc("CuuR",0,2,2,0)+0.034371*ewgc("Cud8R",2,2,1,1)
                    +0.053086*ewgc("CuuR",1,2,2,1)-2.134820*ewgc("CG")+0.045496*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_640_720/sigma_tt_bin_640_720_madgraph);
            
        }
        
    } else if(b_min == 720 && b_max == 800){
        
        double SM_sigma_tt_bin_720_800 = SM.getOptionalParameter("SM_sigma_tt_bin_720_800");
        double sigma_tt_bin_720_800_madgraph = 11.325800;//pb
        
        //PROBABLY THE OPERATORS PURELY FROM THE THIRD FAMILY NEED TO BE NEGLECTED SINCE THEIR CONTRIBUTION IS REALLY SMALL
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_720_800_madgraph-0.000354*ewgc("Cqd8R",2,2,2,2)+0.130586*ewgc("Cqd8R",2,2,0,0)+0.131137*ewgc("Cud8R",2,2,0,0)
                    +0.159377*ewgc("Cqq3R",1,2,2,1)+0.030718*ewgc("Cqq1R",1,2,2,1)+0.345576*ewgc("Cqu8R",0,0,2,2)+0.214910*ewgc("Cqu8R",2,2,0,0)
                    -3.155290*ewgc("CuGR",2,2)+0.018608*ewgc("Cqd8R",2,2,1,1)+0.914467*ewgc("Cqq1R",0,2,2,0)+1.958960*ewgc("Cqq3R",0,2,2,0)
                    +0.007071*ewgc("Cqu8R",2,2,1,1)-0.000354*ewgc("Cud8R",2,2,2,2)+0.910969*ewgc("CuuR",0,2,2,0)+0.018563*ewgc("Cud8R",2,2,1,1)
                    +0.031329*ewgc("CuuR",1,2,2,1)-1.307530*ewgc("CG")+0.023518*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_720_800/sigma_tt_bin_720_800_madgraph);
            
        }
        
    } else if(b_min == 800 && b_max == 900){
        
        double SM_sigma_tt_bin_800_900 = SM.getOptionalParameter("SM_sigma_tt_bin_800_900");
        double sigma_tt_bin_800_900_madgraph = 7.891270;//pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_800_900_madgraph+0.002275*ewgc("Cqd8R",2,2,2,2)+0.122276*ewgc("Cqd8R",2,2,0,0)+0.122208*ewgc("Cud8R",2,2,0,0)
                    +0.145699*ewgc("Cqq3R",1,2,2,1)+0.028882*ewgc("Cqq1R",1,2,2,1)+0.321551*ewgc("Cqu8R",0,0,2,2)+0.202060*ewgc("Cqu8R",2,2,0,0)
                    -2.209900*ewgc("CuGR",2,2)+0.016045*ewgc("Cqd8R",2,2,1,1)+0.851861*ewgc("Cqq1R",0,2,2,0)+1.762370*ewgc("Cqq3R",0,2,2,0)
                    +0.008177*ewgc("Cqu8R",2,2,1,1)+0.002275*ewgc("Cud8R",2,2,2,2)+0.826565*ewgc("CuuR",0,2,2,0)+0.015907*ewgc("Cud8R",2,2,1,1)
                    +0.027875*ewgc("CuuR",1,2,2,1)+-0.974747*ewgc("CG")+0.023889*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_800_900/sigma_tt_bin_800_900_madgraph);
            
        }
        
    } else if(b_min == 900 && b_max == 1000){
        
        double SM_sigma_tt_bin_900_1000 = SM.getOptionalParameter("SM_sigma_tt_bin_900_1000");
        double sigma_tt_bin_900_1000_madgraph = 4.270620;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_900_1000_madgraph+0.000384*ewgc("Cqd8R",2,2,2,2)+0.087383*ewgc("Cqd8R",2,2,0,0)+0.086814*ewgc("Cud8R",2,2,0,0)
                    +0.088968*ewgc("Cqq3R",1,2,2,1)+0.017917*ewgc("Cqq1R",1,2,2,1)+0.227313*ewgc("Cqu8R",0,0,2,2)+0.142724*ewgc("Cqu8R",2,2,0,0)
                    -1.201160*ewgc("CuGR",2,2)+0.009999*ewgc("Cqd8R",2,2,1,1)+0.611818*ewgc("Cqq1R",0,2,2,0)+1.261570*ewgc("Cqq3R",0,2,2,0)
                    +0.003874*ewgc("Cqu8R",2,2,1,1)+0.000384*ewgc("Cud8R",2,2,2,2)+0.595264*ewgc("CuuR",0,2,2,0)+0.009766*ewgc("Cud8R",2,2,1,1)
                    +0.017022*ewgc("CuuR",1,2,2,1)+-0.540272*ewgc("CG")+0.014666*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_900_1000/sigma_tt_bin_900_1000_madgraph);
            
        }
        
    } else if(b_min == 1000 && b_max == 1150){
        
        double SM_sigma_tt_bin_1000_1150 = SM.getOptionalParameter("SM_sigma_tt_bin_1000_1150");
        double sigma_tt_bin_1000_1150_madgraph = 3.191480;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_1000_1150_madgraph+0.000886*ewgc("Cqd8R",2,2,2,2)+0.085754*ewgc("Cqd8R",2,2,0,0)+0.085872*ewgc("Cud8R",2,2,0,0)
                    +0.080947*ewgc("Cqq3R",1,2,2,1)+0.015686*ewgc("Cqq1R",1,2,2,1)+0.229640*ewgc("Cqu8R",0,0,2,2)+0.143023*ewgc("Cqu8R",2,2,0,0)
                    -0.908071*ewgc("CuGR",2,2)+0.009051*ewgc("Cqd8R",2,2,1,1)+0.607012*ewgc("Cqq1R",0,2,2,0)+1.261550*ewgc("Cqq3R",0,2,2,0)
                    +0.004500*ewgc("Cqu8R",2,2,1,1)+0.000886*ewgc("Cud8R",2,2,2,2)+0.587639*ewgc("CuuR",0,2,2,0)+0.008831*ewgc("Cud8R",2,2,1,1)
                    +0.015380*ewgc("CuuR",1,2,2,1)+-0.391609*ewgc("CG")+0.013411*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_1000_1150/sigma_tt_bin_1000_1150_madgraph);
            
        }
        
    } else if(b_min == 1150 && b_max == 1300){
        
        double SM_sigma_tt_bin_1150_1300 = SM.getOptionalParameter("SM_sigma_tt_bin_1150_1300");
        double sigma_tt_bin_1150_1300_madgraph = 1.449930;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_1150_1300_madgraph-0.000060*ewgc("Cqd8R",2,2,2,2)+0.055487*ewgc("Cqd8R",2,2,0,0)+0.055200*ewgc("Cud8R",2,2,0,0)
                    +0.046290*ewgc("Cqq3R",1,2,2,1)+0.009274*ewgc("Cqq1R",1,2,2,1)+0.148907*ewgc("Cqu8R",0,0,2,2)+0.094147*ewgc("Cqu8R",2,2,0,0)
                    -0.413945*ewgc("CuGR",2,2)+0.005015*ewgc("Cqd8R",2,2,1,1)+0.403567*ewgc("Cqq1R",0,2,2,0)+0.790825*ewgc("Cqq3R",0,2,2,0)
                    +0.002144*ewgc("Cqu8R",2,2,1,1)-0.000060*ewgc("Cud8R",2,2,2,2)+0.388209*ewgc("CuuR",0,2,2,0)+0.005003*ewgc("Cud8R",2,2,1,1)
                    +0.008958*ewgc("CuuR",1,2,2,1)+-0.201398*ewgc("CG")+0.007457*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_1150_1300/sigma_tt_bin_1150_1300_madgraph);
            
        }
        
    } else if(b_min == 1300 && b_max == 1500){
        
        double SM_sigma_tt_bin_1300_1500 = SM.getOptionalParameter("SM_sigma_tt_bin_1300_1500");
        double sigma_tt_bin_1300_1500_madgraph = 0.844669;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_1300_1500_madgraph+0.000194*ewgc("Cqd8R",2,2,2,2)+0.044943*ewgc("Cqd8R",2,2,0,0)+0.045385*ewgc("Cud8R",2,2,0,0)
                    +0.033908*ewgc("Cqq3R",1,2,2,1)+0.007047*ewgc("Cqq1R",1,2,2,1)+0.121928*ewgc("Cqu8R",0,0,2,2)+0.077764*ewgc("Cqu8R",2,2,0,0)
                    -0.242809*ewgc("CuGR",2,2)+0.003933*ewgc("Cqd8R",2,2,1,1)+0.335057*ewgc("Cqq1R",0,2,2,0)+0.672494*ewgc("Cqq3R",0,2,2,0)
                    +0.001675*ewgc("Cqu8R",2,2,1,1)+0.000194*ewgc("Cud8R",2,2,2,2)+0.323460*ewgc("CuuR",0,2,2,0)+0.003993*ewgc("Cud8R",2,2,1,1)
                    +0.006851*ewgc("CuuR",1,2,2,1)-0.102284*ewgc("CG")+0.005307*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_1300_1500/sigma_tt_bin_1300_1500_madgraph);
            
        }
        
    } else if(b_min == 1500 && b_max == 1700){
        
        double SM_sigma_tt_bin_1500_1700 = SM.getOptionalParameter("SM_sigma_tt_bin_1500_1700");
        double sigma_tt_bin_1500_1700_madgraph = 0.348438;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_1500_1700_madgraph+0.000086*ewgc("Cqd8R",2,2,2,2)+0.027083*ewgc("Cqd8R",2,2,0,0)+0.026758*ewgc("Cud8R",2,2,0,0)
                    +0.017806*ewgc("Cqq3R",1,2,2,1)+0.003851*ewgc("Cqq1R",1,2,2,1)+0.073439*ewgc("Cqu8R",0,0,2,2)+0.046555*ewgc("Cqu8R",2,2,0,0)
                    -0.100864*ewgc("CuGR",2,2)+0.002036*ewgc("Cqd8R",2,2,1,1)+0.199995*ewgc("Cqq1R",0,2,2,0)+0.412358*ewgc("Cqq3R",0,2,2,0)
                    +0.001051*ewgc("Cqu8R",2,2,1,1)+0.000086*ewgc("Cud8R",2,2,2,2)+0.192436*ewgc("CuuR",0,2,2,0)+0.002041*ewgc("Cud8R",2,2,1,1)
                    +0.003724*ewgc("CuuR",1,2,2,1)+-0.044097*ewgc("CG")+0.002887*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_1500_1700/sigma_tt_bin_1500_1700_madgraph);
            
        }
        
    } else if(b_min == 1700 && b_max == 2000){
        
        double SM_sigma_tt_bin_1700_2000 = SM.getOptionalParameter("SM_sigma_tt_bin_1700_2000");
        double sigma_tt_bin_1700_2000_madgraph = 0.195635;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_1700_2000_madgraph-0.000003*ewgc("Cqd8R",2,2,2,2)+0.021610*ewgc("Cqd8R",2,2,0,0)+0.021827*ewgc("Cud8R",2,2,0,0)
                    +0.013315*ewgc("Cqq3R",1,2,2,1)+0.002866*ewgc("Cqq1R",1,2,2,1)+0.061115*ewgc("Cqu8R",0,0,2,2)+0.038769*ewgc("Cqu8R",2,2,0,0)
                    -0.057671*ewgc("CuGR",2,2)+0.001371*ewgc("Cqd8R",2,2,1,1)+0.165805*ewgc("Cqq1R",0,2,2,0)+0.335094*ewgc("Cqq3R",0,2,2,0)
                    +0.000659*ewgc("Cqu8R",2,2,1,1)+-0.000003*ewgc("Cud8R",2,2,2,2)+0.159666*ewgc("CuuR",0,2,2,0)+0.001359*ewgc("Cud8R",2,2,1,1)
                    +0.002774*ewgc("CuuR",1,2,2,1)+-0.032997*ewgc("CG")+0.002083*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_1700_2000/sigma_tt_bin_1700_2000_madgraph);
            
        }
        
    } else if(b_min == 2000 && b_max == 2300){
        
        double SM_sigma_tt_bin_2000_2300 = SM.getOptionalParameter("SM_sigma_tt_bin_2000_2300");
        double sigma_tt_bin_2000_2300_madgraph = 0.064115;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_2000_2300_madgraph+0.000015*ewgc("Cqd8R",2,2,2,2)+0.010755*ewgc("Cqd8R",2,2,0,0)+0.010716*ewgc("Cud8R",2,2,0,0)
                    +0.005763*ewgc("Cqq3R",1,2,2,1)+0.001393*ewgc("Cqq1R",1,2,2,1)+0.030715*ewgc("Cqu8R",0,0,2,2)+0.019954*ewgc("Cqu8R",2,2,0,0)
                    -0.018962*ewgc("CuGR",2,2)+0.000634*ewgc("Cqd8R",2,2,1,1)+0.088213*ewgc("Cqq1R",0,2,2,0)+0.168864*ewgc("Cqq3R",0,2,2,0)
                    +0.000325*ewgc("Cqu8R",2,2,1,1)+0.000015*ewgc("Cud8R",2,2,2,2)+0.085357*ewgc("CuuR",0,2,2,0)+0.000622*ewgc("Cud8R",2,2,1,1)
                    +0.001357*ewgc("CuuR",1,2,2,1)+-0.005505*ewgc("CG")+0.000949*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_2000_2300/sigma_tt_bin_2000_2300_madgraph);
            
        }
        
    } else if(b_min == 2300 && b_max == 3500){
        
        double SM_sigma_tt_bin_2300_3500 = SM.getOptionalParameter("SM_sigma_tt_bin_2300_3500");
        double sigma_tt_bin_2300_3500_madgraph = 0.036367;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_tt_bin_2300_3500_madgraph+0.000010*ewgc("Cqd8R",2,2,2,2)+0.011318*ewgc("Cqd8R",2,2,0,0)+0.011335*ewgc("Cud8R",2,2,0,0)
                    +0.005970*ewgc("Cqq3R",1,2,2,1)+0.001645*ewgc("Cqq1R",1,2,2,1)+0.034199*ewgc("Cqu8R",0,0,2,2)+0.022847*ewgc("Cqu8R",2,2,0,0)
                    -0.010997*ewgc("CuGR",2,2)+0.000580*ewgc("Cqd8R",2,2,1,1)+0.097656*ewgc("Cqq1R",0,2,2,0)+0.183216*ewgc("Cqq3R",0,2,2,0)
                    +0.000383*ewgc("Cqu8R",2,2,1,1)+0.000010*ewgc("Cud8R",2,2,2,2)+0.093843*ewgc("CuuR",0,2,2,0)+0.000574*ewgc("Cud8R",2,2,1,1)
                    +0.001574*ewgc("CuuR",1,2,2,1)+-0.004467*ewgc("CG")+0.000948*ewgc("Cqu8R",1,1,2,2))*(SM_sigma_tt_bin_2300_3500/sigma_tt_bin_2300_3500_madgraph);
            
        }
        
    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tt_diff_LO.\n");
    }

}



////// ttbar charge asymmetry differential ///////////////////////////////////

charge_asymmetry_tt_diff_mtt_LO::charge_asymmetry_tt_diff_mtt_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_charge_asymmetry_num_bin_mtt_0_500" << "SM_charge_asymmetry_deno_bin_mtt_0_500"
            << "SM_charge_asymmetry_num_bin_mtt_500_750" << "SM_charge_asymmetry_deno_bin_mtt_500_750" << "SM_charge_asymmetry_num_bin_mtt_750_1000" 
            << "SM_charge_asymmetry_deno_bin_mtt_750_1000" << "SM_charge_asymmetry_num_bin_mtt_1000_1500" << "SM_charge_asymmetry_deno_bin_mtt_1000_1500"
            << "SM_charge_asymmetry_num_bin_mtt_1500_3000" << "SM_charge_asymmetry_deno_bin_mtt_1500_3000");
            //"SM_sigma_tt_bin_2300_2600" << "SM_sigma_tt_bin_2600_3000" << "SM_sigma_tt_bin_3000_3500" << "SM_sigma_tt_bin_3500_4000");
    
};


double charge_asymmetry_tt_diff_mtt_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    
    //double C_tG = ewgc("CuGR",2,2);

    
    
    
    if(b_min == 0 && b_max == 500){
        
        double SM_charge_asymmetry_num_bin_mtt_0_500 = SM.getOptionalParameter("SM_charge_asymmetry_num_bin_mtt_0_500");
        double SM_charge_asymmetry_deno_bin_mtt_0_500 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_0_500");
        double SM_charge_asymmetry_bin_mtt_0_500 = SM_charge_asymmetry_num_bin_mtt_0_500/SM_charge_asymmetry_deno_bin_mtt_0_500;
        
        double SM_sigma_pos_bin_mtt_0_500 =0.5*(SM_charge_asymmetry_num_bin_mtt_0_500+SM_charge_asymmetry_deno_bin_mtt_0_500);
        double SM_sigma_neg_bin_mtt_0_500 =0.5*(SM_charge_asymmetry_deno_bin_mtt_0_500-SM_charge_asymmetry_num_bin_mtt_0_500);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_pos_bin_mtt_0_500_madgraph = 0.;
            double sigma_pos_bin_mtt_0_500_NP = 0.;
            double sigma_neg_bin_mtt_0_500_madgraph = 0.;
            double sigma_neg_bin_mtt_0_500_NP = 0.;
            
            double sigma_pos_bin_mtt_0_500_NP_Corrected = SM_sigma_pos_bin_mtt_0_500*sigma_pos_bin_mtt_0_500_NP/sigma_pos_bin_mtt_0_500_madgraph;
            double sigma_neg_bin_mtt_0_500_NP_Corrected = SM_sigma_neg_bin_mtt_0_500*sigma_neg_bin_mtt_0_500_NP/sigma_neg_bin_mtt_0_500_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected - sigma_neg_bin_mtt_0_500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected + sigma_neg_bin_mtt_0_500_NP_Corrected;;
            
            return  SM_charge_asymmetry_bin_mtt_0_500*(1+(NP_charge_asymmetry_num_bin_mtt_0_500-NP_charge_asymmetry_deno_bin_mtt_0_500)/SM_charge_asymmetry_deno_bin_mtt_0_500);            
        }
        
    } else if(b_min == 500 && b_max == 750){
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return 0.;
            
        }
        
    }
    
    else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tt_diff_LO.\n");
    }

}



sigma_tta_diff_LO_CMS_dilepton::sigma_tta_diff_LO_CMS_dilepton(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_35_CMS_dilepton" << "SM_sigma_tta_bin_35_50_CMS_dilepton"
            << "SM_sigma_tta_bin_50_70_CMS_dilepton" << "SM_sigma_tta_bin_70_130_CMS_dilepton" 
            << "SM_sigma_tta_bin_130_200_CMS_dilepton" << "SM_sigma_tta_bin_200_300_CMS_dilepton");

    
};


double sigma_tta_diff_LO_CMS_dilepton::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    //std::cout<<"\033[1;33m b_min \033[0m "<< b_min << std::endl;
    //std::cout<<"\033[1;33m b_max \033[0m "<< b_max << std::endl;
    
    
    if(b_min == 20 && b_max == 35){
        
        double SM_sigma_tta_bin_20_35 = SM.getOptionalParameter("SM_sigma_tta_bin_20_35_CMS_dilepton");
        double sigma_tta_bin_20_35_madgraph = 0.583675; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            //There is no dependence on ewgc("Cqd8R",2,2,2,2) neither on ewgc("Cud8R",2,2,2,2) with the precision considered
        
            return SM_sigma_tta_bin_20_35 +(-0.034106*ewgc("CG")+0.004794*ewgc("Cqd8R",2,2,0,0)
                    +0.000605*ewgc("Cqd8R",2,2,1,1)+0.082493*ewgc("Cqq1R",0,2,2,0)
                    +0.002529*ewgc("Cqq1R",1,2,2,1)+0.146282*ewgc("Cqq3R",0,2,2,0)+0.010782*ewgc("Cqq3R",1,2,2,1)
                    +0.036042*ewgc("Cqu8R",0,0,2,2)+0.001578*ewgc("Cqu8R",1,1,2,2)+0.031190*ewgc("Cqu8R",2,2,0,0)
                    +0.000906*ewgc("Cqu8R",2,2,1,1)+-0.007204*ewgc("CuBR",2,2)+0.008516*ewgc("Cud8R",2,2,0,0)
                    +0.001057*ewgc("Cud8R",2,2,1,1)+-0.154829*ewgc("CuGR",2,2)
                    +0.080748*ewgc("CuuR",0,2,2,0)+0.002498*ewgc("CuuR",1,2,2,1)+-0.005006*ewgc("CuWR",2,2))
                    *(SM_sigma_tta_bin_20_35/sigma_tta_bin_20_35_madgraph);
        }
    } else if(b_min == 35 && b_max == 50){
        
    
        double SM_sigma_tta_bin_35_50 = SM.getOptionalParameter("SM_sigma_tta_bin_35_50_CMS_dilepton");
        double sigma_tta_bin_35_50_madgraph = 0.305601; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return SM_sigma_tta_bin_35_50 +(-0.018555*ewgc("CG")+0.002547*ewgc("Cqd8R",2,2,0,0)
                    +0.000295*ewgc("Cqd8R",2,2,1,1)+0.040351*ewgc("Cqq1R",0,2,2,0)
                    +0.001256*ewgc("Cqq1R",1,2,2,1)+0.075386*ewgc("Cqq3R",0,2,2,0)+0.005401*ewgc("Cqq3R",1,2,2,1)
                    +0.018788*ewgc("Cqu8R",0,0,2,2)+0.000807*ewgc("Cqu8R",1,1,2,2)+0.016262*ewgc("Cqu8R",2,2,0,0)
                    +0.000470*ewgc("Cqu8R",2,2,1,1)+-0.006535*ewgc("CuBR",2,2)+0.004622*ewgc("Cud8R",2,2,0,0)
                    +0.000557*ewgc("Cud8R",2,2,1,1)-0.080445*ewgc("CuGR",2,2)
                    +0.039496*ewgc("CuuR",0,2,2,0)+0.001217*ewgc("CuuR",1,2,2,1)+-0.004069*ewgc("CuWR",2,2))
                    *(SM_sigma_tta_bin_35_50/sigma_tta_bin_35_50_madgraph);
            
        }
    } else if(b_min == 50 && b_max == 70){
        
    
        double SM_sigma_tta_bin_50_70 = SM.getOptionalParameter("SM_sigma_tta_bin_50_70_CMS_dilepton");
        double sigma_tta_bin_50_70_madgraph = 0.238500; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return SM_sigma_tta_bin_50_70 +(-0.014859*ewgc("CG")+0.001934*ewgc("Cqd8R",2,2,0,0)
                    +0.000192*ewgc("Cqd8R",2,2,1,1)+
                    +0.029858*ewgc("Cqq1R",0,2,2,0)+0.000825*ewgc("Cqq1R",1,2,2,1)
                    +0.057705*ewgc("Cqq3R",0,2,2,0)+0.004114*ewgc("Cqq3R",1,2,2,1)
                    +0.014644*ewgc("Cqu8R",0,0,2,2)+0.000602*ewgc("Cqu8R",1,1,2,2)
                    +0.012747*ewgc("Cqu8R",2,2,0,0)+-0.007537*ewgc("CuBR",2,2)
                    +0.003692*ewgc("Cud8R",2,2,0,0)+0.000409*ewgc("Cud8R",2,2,1,1)
                    -0.061847*ewgc("CuGR",2,2)
                    +0.029209*ewgc("CuuR",0,2,2,0)+0.000813*ewgc("CuuR",1,2,2,1)
                    +-0.004460*ewgc("CuWR",2,2))
                    *(SM_sigma_tta_bin_50_70/sigma_tta_bin_50_70_madgraph);
            
        }
    } else if(b_min == 70 && b_max == 130){
        
    
        double SM_sigma_tta_bin_70_130 = SM.getOptionalParameter("SM_sigma_tta_bin_70_130_CMS_dilepton");
        double sigma_tta_bin_70_100_madgraph = 0.198831; //pb
        double sigma_tta_bin_100_130_madgraph =0.110866; //pb
        double sigma_tta_bin_70_130_madgraph =sigma_tta_bin_70_100_madgraph+sigma_tta_bin_100_130_madgraph; //pb
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_tta_bin_70_100_NP = (-0.013285*ewgc("CG")+0.001630*ewgc("Cqd8R",2,2,0,0)
                    +0.000250*ewgc("Cqd8R",2,2,1,1)
                    +0.023767*ewgc("Cqq1R",0,2,2,0)+0.000676*ewgc("Cqq1R",1,2,2,1)
                    +0.047671*ewgc("Cqq3R",0,2,2,0)+0.003336*ewgc("Cqq3R",1,2,2,1)
                    +0.012379*ewgc("Cqu8R",0,0,2,2)+0.000535*ewgc("Cqu8R",1,1,2,2)
                    +0.010765*ewgc("Cqu8R",2,2,0,0)+0.000350*ewgc("Cqu8R",2,2,1,1)
                    -0.008505*ewgc("CuBR",2,2)+0.003205*ewgc("Cud8R",2,2,0,0)
                    +0.000431*ewgc("Cud8R",2,2,1,1)
                    -0.050838*ewgc("CuGR",2,2)+0.023185*ewgc("CuuR",0,2,2,0)
                    +0.000655*ewgc("CuuR",1,2,2,1)-0.004930*ewgc("CuWR",2,2));
             
            
            double sigma_tta_bin_100_130_NP = (-0.008092*ewgc("CG")+0.000916*ewgc("Cqd8R",2,2,0,0)+0.000088*ewgc("Cqd8R",2,2,1,1)
                    +0.012744*ewgc("Cqq1R",0,2,2,0)+0.000316*ewgc("Cqq1R",1,2,2,1)
                    +0.026925*ewgc("Cqq3R",0,2,2,0)+0.001792*ewgc("Cqq3R",1,2,2,1)+0.007060*ewgc("Cqu8R",0,0,2,2)
                    +0.000247*ewgc("Cqu8R",1,1,2,2)+0.006144*ewgc("Cqu8R",2,2,0,0)+0.000162*ewgc("Cqu8R",2,2,1,1)
                    +-0.005920*ewgc("CuBR",2,2)+0.001849*ewgc("Cud8R",2,2,0,0)+0.000200*ewgc("Cud8R",2,2,1,1)
                    +-0.027969*ewgc("CuGR",2,2)+0.012428*ewgc("CuuR",0,2,2,0)
                    +0.000308*ewgc("CuuR",1,2,2,1)+-0.003364*ewgc("CuWR",2,2));
            
            return SM_sigma_tta_bin_70_130 + (sigma_tta_bin_70_100_NP + sigma_tta_bin_100_130_NP)*(SM_sigma_tta_bin_70_130/sigma_tta_bin_70_130_madgraph);
            
        }
    } else if(b_min == 130 && b_max == 200){
        
    
        double SM_sigma_tta_bin_130_200 = SM.getOptionalParameter("SM_sigma_tta_bin_130_200_CMS_dilepton");
        double sigma_tta_bin_130_165_madgraph = 0.075543; //pb
        double sigma_tta_bin_165_200_madgraph = 0.044687; //pb
        double sigma_tta_bin_130_200_madgraph = sigma_tta_bin_130_165_madgraph+sigma_tta_bin_165_200_madgraph;
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_tta_bin_130_165_NP = (-0.006183*ewgc("CG")+0.000634*ewgc("Cqd8R",2,2,0,0)+0.000043*ewgc("Cqd8R",2,2,1,1)
                    +0.008637*ewgc("Cqq1R",0,2,2,0)+0.000173*ewgc("Cqq1R",1,2,2,1)
                    +0.018848*ewgc("Cqq3R",0,2,2,0)+0.001189*ewgc("Cqq3R",1,2,2,1)+0.005006*ewgc("Cqu8R",0,0,2,2)
                    +0.000145*ewgc("Cqu8R",1,1,2,2)+0.004346*ewgc("Cqu8R",2,2,0,0)+0.000081*ewgc("Cqu8R",2,2,1,1)
                    +-0.004429*ewgc("CuBR",2,2)+0.001330*ewgc("Cud8R",2,2,0,0)+0.000119*ewgc("Cud8R",2,2,1,1)
                    +-0.018971*ewgc("CuGR",2,2)+0.008402*ewgc("CuuR",0,2,2,0)
                    +0.000171*ewgc("CuuR",1,2,2,1)+-0.002499*ewgc("CuWR",2,2));
            
            double sigma_tta_bin_165_200_NP = (-0.004285*ewgc("CG")+0.000401*ewgc("Cqd8R",2,2,0,0)+0.000029*ewgc("Cqd8R",2,2,1,1)
                    +0.005252*ewgc("Cqq1R",0,2,2,0)+0.000104*ewgc("Cqq1R",1,2,2,1)
                    +0.011756*ewgc("Cqq3R",0,2,2,0)+0.000720*ewgc("Cqq3R",1,2,2,1)+0.003157*ewgc("Cqu8R",0,0,2,2)
                    +0.000089*ewgc("Cqu8R",1,1,2,2)+0.002752*ewgc("Cqu8R",2,2,0,0)+0.000051*ewgc("Cqu8R",2,2,1,1)
                    +-0.002716*ewgc("CuBR",2,2)+0.000855*ewgc("Cud8R",2,2,0,0)+0.000072*ewgc("Cud8R",2,2,1,1)
                    +-0.011315*ewgc("CuGR",2,2)+0.005106*ewgc("CuuR",0,2,2,0)
                    +0.000098*ewgc("CuuR",1,2,2,1)+-0.001515*ewgc("CuWR",2,2));
            
            return SM_sigma_tta_bin_130_200 + (sigma_tta_bin_130_165_NP+sigma_tta_bin_165_200_NP)*(SM_sigma_tta_bin_130_200/sigma_tta_bin_130_200_madgraph);
            
        }
    } else if(b_min == 200 && b_max == 300){
        
    
        double SM_sigma_tta_bin_200_300 = SM.getOptionalParameter("SM_sigma_tta_bin_200_300_CMS_dilepton");
        double sigma_tta_bin_200_250_madgraph = 0.035911; //pb
        double sigma_tta_bin_250_300_madgraph = 0.018900; //pb
        double sigma_tta_bin_200_300_madgraph = sigma_tta_bin_200_250_madgraph+sigma_tta_bin_250_300_madgraph; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            double sigma_tta_bin_200_250_NP = (-0.004232*ewgc("CG")+0.000362*ewgc("Cqd8R",2,2,0,0)+0.000028*ewgc("Cqd8R",2,2,1,1)
                    +0.004454*ewgc("Cqq1R",0,2,2,0)+0.000087*ewgc("Cqq1R",1,2,2,1)
                    +0.010395*ewgc("Cqq3R",0,2,2,0)+0.000590*ewgc("Cqq3R",1,2,2,1)+0.002802*ewgc("Cqu8R",0,0,2,2)
                    +0.000078*ewgc("Cqu8R",1,1,2,2)+0.002436*ewgc("Cqu8R",2,2,0,0)+0.000047*ewgc("Cqu8R",2,2,1,1)
                    -0.002135*ewgc("CuBR",2,2)+0.000773*ewgc("Cud8R",2,2,0,0)+0.000064*ewgc("Cud8R",2,2,1,1)
                    -0.009208*ewgc("CuGR",2,2)+0.004354*ewgc("CuuR",0,2,2,0)
                    +0.000086*ewgc("CuuR",1,2,2,1)+-0.001220*ewgc("CuWR",2,2));
            
            double sigma_tta_bin_250_300_NP = (-0.002904*ewgc("CG")+0.000216*ewgc("Cqd8R",2,2,0,0)+0.000014*ewgc("Cqd8R",2,2,1,1)
                    +0.002603*ewgc("Cqq1R",0,2,2,0)+0.000045*ewgc("Cqq1R",1,2,2,1)
                    +0.006131*ewgc("Cqq3R",0,2,2,0)+0.000331*ewgc("Cqq3R",1,2,2,1)+0.001689*ewgc("Cqu8R",0,0,2,2)
                    +0.000042*ewgc("Cqu8R",1,1,2,2)+0.001470*ewgc("Cqu8R",2,2,0,0)+0.000026*ewgc("Cqu8R",2,2,1,1)
                    -0.001064*ewgc("CuBR",2,2)+0.000469*ewgc("Cud8R",2,2,0,0)+0.000035*ewgc("Cud8R",2,2,1,1)
                    +-0.004994*ewgc("CuGR",2,2)+0.002538*ewgc("CuuR",0,2,2,0)
                    +0.000044*ewgc("CuuR",1,2,2,1)+-0.000597*ewgc("CuWR",2,2));
            
            return SM_sigma_tta_bin_200_300 + (sigma_tta_bin_200_250_NP + sigma_tta_bin_250_300_NP)*(SM_sigma_tta_bin_200_300/sigma_tta_bin_200_300_madgraph);
            
        }
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO.\n");
    }

}


//// ttz differential cross section //////


sigma_ttz_diff_LO_ATLAS_210312603::sigma_ttz_diff_LO_ATLAS_210312603(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttz_bin_0_40_ATLAS_210312603" << "SM_sigma_ttz_bin_40_70_ATLAS_210312603"
            << "SM_sigma_ttz_bin_70_110_ATLAS_210312603" << "SM_sigma_ttz_bin_110_160_ATLAS_210312603" << "SM_sigma_ttz_bin_160_220_ATLAS_210312603" 
            << "SM_sigma_ttz_bin_220_290_ATLAS_210312603" << "SM_sigma_ttz_bin_290_400_ATLAS_210312603");

    
};


double sigma_ttz_diff_LO_ATLAS_210312603::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
   
    if(b_min == 0 && b_max == 40){
        
        double SM_sigma_ttz_bin_0_40 = SM.getOptionalParameter("SM_sigma_ttz_bin_0_40_ATLAS_210312603");
        double sigma_ttz_bin_0_40_madgraph = 0.064177; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            return (sigma_ttz_bin_0_40_madgraph -0.010478*ewgc("CG")-0.005804*ewgc("CHq1R",2,2)
                    +0.005759*ewgc("CHq3R",2,2)+0.002967*ewgc("CHuR",2,2)+0.001126*ewgc("Cqd8R",2,2,0,0)
                    +0.000122*ewgc("Cqd8R",2,2,1,1)
                    +0.008742*ewgc("Cqq1R",0,2,2,0)+0.000211*ewgc("Cqq1R",1,2,2,1)
                    +0.045418*ewgc("Cqq3R",0,2,2,0)+0.004850*ewgc("Cqq3R",1,2,2,1)
                    +0.005860*ewgc("Cqu8R",0,0,2,2)+0.000532*ewgc("Cqu8R",1,1,2,2)
                    +0.000976*ewgc("Cqu8R",2,2,0,0)+0.000001*ewgc("Cqu8R",2,2,1,1)
                    -0.000170*ewgc("CuBR",2,2)+0.000666*ewgc("Cud8R",2,2,0,0)+0.000058*ewgc("Cud8R",2,2,1,1)
                    -0.018690*ewgc("CuGR",2,2)+0.002251*ewgc("CuuR",0,2,2,0)
                    +0.000036*ewgc("CuuR",1,2,2,1)-0.000098*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_0_40/sigma_ttz_bin_0_40_madgraph);
        }
        
    } else if(b_min == 40 && b_max == 70){
        
        double SM_sigma_ttz_bin_40_70 = SM.getOptionalParameter("SM_sigma_ttz_bin_40_70_ATLAS_210312603");
        double sigma_ttz_bin_40_70_madgraph = 0.096642; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            return (sigma_ttz_bin_40_70_madgraph -0.018049*ewgc("CG")-0.009212*ewgc("CHq1R",2,2)
                    +0.009177*ewgc("CHq3R",2,2)+0.005078*ewgc("CHuR",2,2)+0.001840*ewgc("Cqd8R",2,2,0,0)
                    +0.000234*ewgc("Cqd8R",2,2,1,1)
                    +0.013449*ewgc("Cqq1R",0,2,2,0)+0.000372*ewgc("Cqq1R",1,2,2,1)
                    +0.065355*ewgc("Cqq3R",0,2,2,0)+0.006633*ewgc("Cqq3R",1,2,2,1)
                    +0.008494*ewgc("Cqu8R",0,0,2,2)+0.000775*ewgc("Cqu8R",1,1,2,2)
                    +0.001722*ewgc("Cqu8R",2,2,0,0)+0.000055*ewgc("Cqu8R",2,2,1,1)
                    -0.000190*ewgc("CuBR",2,2)+0.001097*ewgc("Cud8R",2,2,0,0)
                    +0.000148*ewgc("Cud8R",2,2,1,1)+-0.029042*ewgc("CuGR",2,2)
                    +0.003988*ewgc("CuuR",0,2,2,0)+0.000126*ewgc("CuuR",1,2,2,1)
                    -0.000106*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_40_70/sigma_ttz_bin_40_70_madgraph);
        }
        
    } else if(b_min == 70 && b_max == 110){
        
        double SM_sigma_ttz_bin_70_110 = SM.getOptionalParameter("SM_sigma_ttz_bin_70_110_ATLAS_210312603");
        double sigma_ttz_bin_70_110_madgraph = 0.127419; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            return (sigma_ttz_bin_70_110_madgraph-0.028383*ewgc("CG")+-0.012875*ewgc("CHq1R",2,2)+0.012884*ewgc("CHq3R",2,2)
                    +0.007778*ewgc("CHuR",2,2)+0.002468*ewgc("Cqd8R",2,2,0,0)+0.000298*ewgc("Cqd8R",2,2,1,1)
                    +0.017859*ewgc("Cqq1R",0,2,2,0)+0.000464*ewgc("Cqq1R",1,2,2,1)
                    +0.080694*ewgc("Cqq3R",0,2,2,0)+0.007774*ewgc("Cqq3R",1,2,2,1)+0.010371*ewgc("Cqu8R",0,0,2,2)
                    +0.000893*ewgc("Cqu8R",1,1,2,2)+0.002610*ewgc("Cqu8R",2,2,0,0)+0.000063*ewgc("Cqu8R",2,2,1,1)
                    -0.000135*ewgc("CuBR",2,2)+0.001481*ewgc("Cud8R",2,2,0,0)+0.000183*ewgc("Cud8R",2,2,1,1)
                    -0.039573*ewgc("CuGR",2,2)+0.006317*ewgc("CuuR",0,2,2,0)+0.000171*ewgc("CuuR",1,2,2,1)
                    -0.000281*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_70_110/sigma_ttz_bin_70_110_madgraph);
        }
        
    } else if(b_min == 110 && b_max == 160){
        
        double SM_sigma_ttz_bin_110_160 = SM.getOptionalParameter("SM_sigma_ttz_bin_110_160_ATLAS_210312603");
        double sigma_ttz_bin_110_160_madgraph = 0.118917; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            return (sigma_ttz_bin_110_160_madgraph -0.032594*ewgc("CG")-0.012662*ewgc("CHq1R",2,2)
                    +0.012643*ewgc("CHq3R",2,2)+0.008236*ewgc("CHuR",2,2)+0.002540*ewgc("Cqd8R",2,2,0,0)
                    +0.000253*ewgc("Cqd8R",2,2,1,1)+0.017566*ewgc("Cqq1R",0,2,2,0)+0.000417*ewgc("Cqq1R",1,2,2,1)
                    +0.071863*ewgc("Cqq3R",0,2,2,0)+0.006279*ewgc("Cqq3R",1,2,2,1)+0.009224*ewgc("Cqu8R",0,0,2,2)
                    +0.000682*ewgc("Cqu8R",1,1,2,2)+0.002931*ewgc("Cqu8R",2,2,0,0)+0.000049*ewgc("Cqu8R",2,2,1,1)
                    +-0.000017*ewgc("CuBR",2,2)+0.001474*ewgc("Cud8R",2,2,0,0)+0.000150*ewgc("Cud8R",2,2,1,1)
                    -0.038866*ewgc("CuGR",2,2)+0.007314*ewgc("CuuR",0,2,2,0)
                    +0.000171*ewgc("CuuR",1,2,2,1)+-0.000599*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_110_160/sigma_ttz_bin_110_160_madgraph);
        }
        
    } else if(b_min == 160 && b_max == 220){
        
        double SM_sigma_ttz_bin_160_220 = SM.getOptionalParameter("SM_sigma_ttz_bin_160_220_ATLAS_210312603");
        double sigma_ttz_bin_160_220_madgraph = 0.086169; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            return (sigma_ttz_bin_160_220_madgraph -0.028727*ewgc("CG")+-0.009451*ewgc("CHq1R",2,2)
                    +0.009451*ewgc("CHq3R",2,2)+0.006505*ewgc("CHuR",2,2)+0.002109*ewgc("Cqd8R",2,2,0,0)
                    +0.000210*ewgc("Cqd8R",2,2,1,1)+0.014508*ewgc("Cqq1R",0,2,2,0)
                    +0.000323*ewgc("Cqq1R",1,2,2,1)+0.054025*ewgc("Cqq3R",0,2,2,0)+0.004204*ewgc("Cqq3R",1,2,2,1)
                    +0.006853*ewgc("Cqu8R",0,0,2,2)+0.000457*ewgc("Cqu8R",1,1,2,2)+0.002647*ewgc("Cqu8R",2,2,0,0)
                    +0.000061*ewgc("Cqu8R",2,2,1,1)+0.000070*ewgc("CuBR",2,2)+0.001274*ewgc("Cud8R",2,2,0,0)
                    +0.000130*ewgc("Cud8R",2,2,1,1)+-0.029532*ewgc("CuGR",2,2)
                    +0.006980*ewgc("CuuR",0,2,2,0)+0.000163*ewgc("CuuR",1,2,2,1)+-0.000537*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_160_220/sigma_ttz_bin_160_220_madgraph);
        }
        
    } else if(b_min == 220 && b_max == 290){
        
        double SM_sigma_ttz_bin_220_290 = SM.getOptionalParameter("SM_sigma_ttz_bin_220_290_ATLAS_210312603");
        double sigma_ttz_bin_220_290_madgraph = 0.051619; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            return (sigma_ttz_bin_220_290_madgraph -0.020633*ewgc("CG")+-0.005708*ewgc("CHq1R",2,2)
                    +0.005730*ewgc("CHq3R",2,2)+0.004111*ewgc("CHuR",2,2)+0.001568*ewgc("Cqd8R",2,2,0,0)
                    +0.000156*ewgc("Cqd8R",2,2,1,1)+0.010800*ewgc("Cqq1R",0,2,2,0)
                    +0.000238*ewgc("Cqq1R",1,2,2,1)+0.036604*ewgc("Cqq3R",0,2,2,0)+0.002537*ewgc("Cqq3R",1,2,2,1)
                    +0.004688*ewgc("Cqu8R",0,0,2,2)+0.000287*ewgc("Cqu8R",1,1,2,2)+0.002098*ewgc("Cqu8R",2,2,0,0)
                    +0.000060*ewgc("Cqu8R",2,2,1,1)+0.000066*ewgc("CuBR",2,2)+0.000986*ewgc("Cud8R",2,2,0,0)
                    +0.000107*ewgc("Cud8R",2,2,1,1)-0.018706*ewgc("CuGR",2,2)
                    +0.005776*ewgc("CuuR",0,2,2,0)+0.000142*ewgc("CuuR",1,2,2,1)-0.000356*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_220_290/sigma_ttz_bin_220_290_madgraph);
        }
        
    } else if(b_min == 290 && b_max == 400){
        
        double SM_sigma_ttz_bin_290_400 = SM.getOptionalParameter("SM_sigma_ttz_bin_290_400_ATLAS_210312603");
        double sigma_ttz_bin_290_400_madgraph = 0.032507; //pb
        
        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
                    
            return (sigma_ttz_bin_290_400_madgraph -0.015551*ewgc("CG")+-0.003576*ewgc("CHq1R",2,2)
                    +0.003588*ewgc("CHq3R",2,2)+0.002640*ewgc("CHuR",2,2)+0.001358*ewgc("Cqd8R",2,2,0,0)
                    +0.000110*ewgc("Cqd8R",2,2,1,1)+0.009422*ewgc("Cqq1R",0,2,2,0)
                    +0.000177*ewgc("Cqq1R",1,2,2,1)+0.029289*ewgc("Cqq3R",0,2,2,0)+0.001743*ewgc("Cqq3R",1,2,2,1)
                    +0.003812*ewgc("Cqu8R",0,0,2,2)+0.000192*ewgc("Cqu8R",1,1,2,2)+0.001922*ewgc("Cqu8R",2,2,0,0)
                    +0.000038*ewgc("Cqu8R",2,2,1,1)+0.000016*ewgc("CuBR",2,2)+0.000882*ewgc("Cud8R",2,2,0,0)
                    +0.000076*ewgc("Cud8R",2,2,1,1)-0.012764*ewgc("CuGR",2,2)
                    +0.005527*ewgc("CuuR",0,2,2,0)+0.000108*ewgc("CuuR",1,2,2,1)+-0.000264*ewgc("CuWR",2,2))
                    *(SM_sigma_ttz_bin_290_400/sigma_ttz_bin_290_400_madgraph);
        }
        
    }
}



//// s-channel 13 TeV ////

sigma_tb_13_LO::sigma_tb_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tb_13");
};

double sigma_tb_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tb_13_LO_madgraph = 8.005990; //pb?
    double sigma_tb_13_SM = SM.getOptionalParameter("SM_sigma_tb_13"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return sigma_tb_13_SM + (0.970643*ewgc("CHq3R",2,2)-3.718540*ewgc("CuWR",2,2) )*(sigma_tb_13_SM/sigma_tb_13_LO_madgraph);
        }

}




//// t-channel 13 TeV ////

sigma_tq_13_LO::sigma_tq_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tq_13");
};

double sigma_tq_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tq_13_LO_madgraph = 186.321000; //pb?
    double sigma_tq_13_SM = SM.getOptionalParameter("SM_sigma_tq_13"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return sigma_tq_13_SM +(22.491500*ewgc("CHq3R",2,2)-4.598280*ewgc("CuWR",2,2) )*(sigma_tq_13_SM/sigma_tq_13_LO_madgraph);
        }

}



//// taq ////

sigma_taq_LO_CMS::sigma_taq_LO_CMS(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_taq_CMS");
};

double sigma_taq_LO_CMS::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_taq_LO_madgraph = 2.176180; //pb?
    double sigma_taq_SM_CMS = SM.getOptionalParameter("SM_sigma_taq_CMS"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
           
            return sigma_taq_SM_CMS +(0.262660*ewgc("CHq3R",2,2)-0.003262*ewgc("CuBR",2,2)-0.030453*ewgc("CuWR",2,2))*(sigma_taq_SM_CMS/sigma_taq_LO_madgraph);
        }

}


sigma_taq_LO_ATLAS::sigma_taq_LO_ATLAS(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_taq_ATLAS");
};

double sigma_taq_LO_ATLAS::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_taq_LO_madgraph = 2.176180; //pb?
    double sigma_taq_SM_ATLAS = SM.getOptionalParameter("SM_sigma_taq_ATLAS"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return sigma_taq_SM_ATLAS +(0.262660*ewgc("CHq3R",2,2)-0.003262*ewgc("CuBR",2,2)-0.030453*ewgc("CuWR",2,2))
                    *(sigma_taq_SM_ATLAS/sigma_taq_LO_madgraph);
        }

}



//// tzq ////

sigma_tzq_LO::sigma_tzq_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tzq");
};

double sigma_tzq_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tzq_LO_madgraph = 0.725915; //pb?
    double sigma_tzq_SM = SM.getOptionalParameter("SM_sigma_tzq"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return sigma_tzq_SM +(0.015938*ewgc("CHq1R",2,2)+0.141467*ewgc("CHq3R",2,2)
                    +0.004390*ewgc("CHuR",2,2)+-0.002539*ewgc("CuBR",2,2)-0.001949*ewgc("CuWR",2,2))
                    *(sigma_tzq_SM/sigma_tzq_LO_madgraph);
        }

}



//// tw  ////

sigma_tw_13_LO::sigma_tw_13_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tw_13");
};

double sigma_tw_13_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_tw_13_LO_madgraph = 51.021400; //pb?
    double sigma_tw_13_SM = SM.getOptionalParameter("SM_sigma_tw_13"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return sigma_tw_13_SM + (6.152940*ewgc("CHq3R",2,2)+4.361620*ewgc("CuWR",2,2))
                    *(sigma_tw_13_SM/sigma_tw_13_LO_madgraph);
        }

}




//// ttw  ////

sigma_ttw_LO::sigma_ttw_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttw");
};

double sigma_ttw_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double sigma_ttw_LO_madgraph = 0.516319; //pb?
    double sigma_ttw_SM = SM.getOptionalParameter("SM_sigma_ttw"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return sigma_ttw_SM +(-0.000325*ewgc("CHq1R",2,2)+0.000604*ewgc("CHq3R",2,2)-0.000146*ewgc("CHuR",2,2)
                    +0.184757*ewgc("Cqq1R",0,2,2,0)+0.010580*ewgc("Cqq1R",1,2,2,1)+0.872452*ewgc("Cqq3R",0,2,2,0)
                    +0.054539*ewgc("Cqq3R",1,2,2,1)+0.133233*ewgc("Cqu8R",0,0,2,2)+0.008296*ewgc("Cqu8R",1,1,2,2)
                    +0.000098*ewgc("CuBR",2,2)-0.125985*ewgc("CuGR",2,2)-0.005842*ewgc("CuWR",2,2))
                    *(sigma_ttw_SM/sigma_ttw_LO_madgraph);
        }

}




/////// Helicities  /////////////////////////////////////////////////////////////////////////////////////



F0_LO::F0_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_F0");
};

double F0_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double F0_madgraph = 0.70381;
    double F0_SM = SM.getOptionalParameter("SM_F0"); 

    
    
    if(flag_Quadratic){
        return  0.;
        }
        else{
            return F0_SM + (-0.0569231*getSMEFTCoeffEW("CuWR",2,2) )*(F0_SM/F0_madgraph);
        }

}



FL_LO::FL_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "SM_FL");
};

double FL_LO::computeThValue()
{
    
    bool   flag_Quadratic= false;//Needs to be properly defined
    
    double FL_madgraph = 0.29619;
    double FL_SM = SM.getOptionalParameter("SM_FL");
    

    if(flag_Quadratic){
        return  0.;
        }
    else{
        return FL_SM + (0.0569227*(getSMEFTCoeffEW("CuWR",2,2)) )*(FL_SM/FL_madgraph);
        }

    
}




Test_direct::Test_direct(const StandardModel& SM_i) 
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{ 

    setParametersForObservable(make_vector<std::string>() << "Test_direct_in");
};

double Test_direct::computeThValue()
{
    double Test_direct_in = SM.getOptionalParameter("Test_direct_in");
    
    return Test_direct_in;

    
}



