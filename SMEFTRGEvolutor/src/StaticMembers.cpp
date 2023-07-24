//StaticMembers.cc

/*
 Non-integral static members can't be inizialized 
 directly in the declaration in the header file, 
 so are initialized here.
 */

const double RGESolver::TWO_THIRDS = (2. / 3.); 
const double RGESolver::FOUR_THIRDS = (4. / 3.);  
const double RGESolver::EIGHT_THIRDS = (8. / 3.);  
const double RGESolver::ONE_THIRD = (1. / 3.);  
const double RGESolver::ONE_SIXTH = (1. / 6.);  
const double RGESolver::TEN_THIRDS = (10. / 3.);  

const double RGESolver::NC = 3.; 
const double RGESolver::NC2 = 9.; 

const double RGESolver::b01 = (- 1. / 6. - 3. * 20. / 9.); 
const double RGESolver::b02 = (43. / 6. - 4. * 3. / 3.); 
const double RGESolver::b03 = (11. - 4. * 3. / 3.);

//Casimirs
const double RGESolver::cA2 = double(2.); 
const double RGESolver::cA3 = double(NC); 
const double RGESolver::cF2 = double(3. / 4.); 
const double RGESolver::cF3 = double(0.5 * (NC*NC - 1) / NC); 

//Hypercharges (and product of hypercharges)
//We use Q = T3 + Y (the other possibility 
//is Q = T3 + Y/2) as in https://arxiv.org/abs/1308.2627

const double RGESolver::Yh = (0.5); 
const double RGESolver::Yh2 = Yh*Yh;

const double RGESolver::Yq = (1. / 6.); 
const double RGESolver::Yq2 = Yq*Yq;
const double RGESolver::Yl = (- 0.5); 
const double RGESolver::Yl2 = Yl*Yl;

const double RGESolver::Yu = (2. / 3.); 
const double RGESolver::Yu2 = Yu*Yu;
const double RGESolver::Yd = (- 1. / 3.); 
const double RGESolver::Yd2 = Yd*Yd;
const double RGESolver::Ye = (- 1.); 
const double RGESolver::Ye2 = Ye*Ye;

const double RGESolver::YhYu = Yh*Yu;
const double RGESolver::YhYd = Yh*Yd;
const double RGESolver::YhYe = Yh*Ye;
const double RGESolver::YhYq = Yh*Yq;
const double RGESolver::YhYl = Yh*Yl;

const double RGESolver::YuYd = Yu*Yd;
const double RGESolver::YuYe = Yu*Ye;
const double RGESolver::YuYq = Yu*Yq;
const double RGESolver::YuYl = Yu*Yl;

const double RGESolver::YdYe = Yd*Ye;
const double RGESolver::YdYq = Yd*Yq;
const double RGESolver::YdYl = Yd*Yl;

const double RGESolver::YeYq = Ye*Yq;
const double RGESolver::YeYl = Ye*Yl;

const double RGESolver::YlYq = Yl*Yq;



const double RGESolver::delta[3][3] = {
    {1., 0., 0.},
    { 0., 1., 0.},
    { 0., 0., 1.}
};

const int RGESolver::WC2R_indices[DWC2R][2] = {
    {0, 0},
    {0, 1},
    {0, 2},
    {1, 1},
    {1, 2},
    {2, 2}
};



const int RGESolver::WC2I_indices[DWC2I][2] = {
    {0, 1},
    {0, 2},
    {1, 2}
};


const int RGESolver::WC6R_indices[DWC6R][4] = {
    {0, 0, 0, 0},
    {0, 0, 0, 1},
    {0, 0, 0, 2},
    {0, 0, 1, 1},
    {0, 0, 1, 2},
    {0, 0, 2, 2},
    {0, 1, 0, 1},
    {0, 1, 0, 2},
    {0, 1, 1, 0},
    {0, 1, 1, 1},
    {0, 1, 1, 2},
    {0, 1, 2, 0},
    {0, 1, 2, 1},
    {0, 1, 2, 2},
    {0, 2, 0, 2},
    {0, 2, 1, 1},
    {0, 2, 1, 2},
    {0, 2, 2, 0},
    {0, 2, 2, 1},
    {0, 2, 2, 2},
    {1, 1, 1, 1},
    {1, 1, 1, 2},
    {1, 1, 2, 2},
    {1, 2, 1, 2},
    {1, 2, 2, 1},
    {1, 2, 2, 2},
    {2, 2, 2, 2}
};
const int RGESolver::WC6I_indices[DWC6I][4] = {
    {0, 0, 0, 1},
    {0, 0, 0, 2},
    {0, 0, 1, 2},
    {0, 1, 0, 1},
    {0, 1, 0, 2},
    {0, 1, 1, 1},
    {0, 1, 1, 2},
    {0, 1, 2, 0},
    {0, 1, 2, 1},
    {0, 1, 2, 2},
    {0, 2, 0, 2},
    {0, 2, 1, 1},
    {0, 2, 1, 2},
    {0, 2, 2, 1},
    {0, 2, 2, 2},
    {1, 1, 1, 2},
    {1, 2, 1, 2},
    {1, 2, 2, 2}
};

const int RGESolver::WC7R_indices[DWC7R][4] = {
    {0, 0, 0, 0},
    {0, 0, 0, 1},
    {0, 0, 0, 2},
    {0, 0, 1, 1},
    {0, 0, 1, 2},
    {0, 0, 2, 2},
    {0, 1, 0, 0},
    {0, 1, 0, 1},
    {0, 1, 0, 2},
    {0, 1, 1, 0},
    {0, 1, 1, 1},
    {0, 1, 1, 2},
    {0, 1, 2, 0},
    {0, 1, 2, 1},
    {0, 1, 2, 2},
    {0, 2, 0, 0},
    {0, 2, 0, 1},
    {0, 2, 0, 2},
    {0, 2, 1, 0},
    {0, 2, 1, 1},
    {0, 2, 1, 2},
    {0, 2, 2, 0},
    {0, 2, 2, 1},
    {0, 2, 2, 2},
    {1, 1, 0, 0},
    {1, 1, 0, 1},
    {1, 1, 0, 2},
    {1, 1, 1, 1},
    {1, 1, 1, 2},
    {1, 1, 2, 2},
    {1, 2, 0, 0},
    {1, 2, 0, 1},
    {1, 2, 0, 2},
    {1, 2, 1, 0},
    {1, 2, 1, 1},
    {1, 2, 1, 2},
    {1, 2, 2, 0},
    {1, 2, 2, 1},
    {1, 2, 2, 2},
    {2, 2, 0, 0},
    {2, 2, 0, 1},
    {2, 2, 0, 2},
    {2, 2, 1, 1},
    {2, 2, 1, 2},
    {2, 2, 2, 2}
};

const int RGESolver::WC7I_indices[DWC7I][4] = {
    {0, 0, 0, 1},
    {0, 0, 0, 2},
    {0, 0, 1, 2},
    {0, 1, 0, 0},
    {0, 1, 0, 1},
    {0, 1, 0, 2},
    {0, 1, 1, 0},
    {0, 1, 1, 1},
    {0, 1, 1, 2},
    {0, 1, 2, 0},
    {0, 1, 2, 1},
    {0, 1, 2, 2},
    {0, 2, 0, 0},
    {0, 2, 0, 1},
    {0, 2, 0, 2},
    {0, 2, 1, 0},
    {0, 2, 1, 1},
    {0, 2, 1, 2},
    {0, 2, 2, 0},
    {0, 2, 2, 1},
    {0, 2, 2, 2},
    {1, 1, 0, 1},
    {1, 1, 0, 2},
    {1, 1, 1, 2},
    {1, 2, 0, 0},
    {1, 2, 0, 1},
    {1, 2, 0, 2},
    {1, 2, 1, 0},
    {1, 2, 1, 1},
    {1, 2, 1, 2},
    {1, 2, 2, 0},
    {1, 2, 2, 1},
    {1, 2, 2, 2},
    {2, 2, 0, 1},
    {2, 2, 0, 2},
    {2, 2, 1, 2}
};


const int RGESolver::WC8R_indices[DWC8R][4] = {
    {0, 0, 0, 0},
    {0, 0, 0, 1},
    {0, 0, 0, 2},
    {0, 0, 1, 1},
    {0, 0, 1, 2},
    {0, 0, 2, 2},
    {0, 1, 0, 1},
    {0, 1, 0, 2},
    {0, 1, 1, 1},
    {0, 1, 1, 2},
    {0, 1, 2, 1},
    {0, 1, 2, 2},
    {0, 2, 0, 2},
    {0, 2, 1, 2},
    {0, 2, 2, 2},
    {1, 1, 1, 1},
    {1, 1, 1, 2},
    {1, 1, 2, 2},
    {1, 2, 1, 2},
    {1, 2, 2, 2},
    {2, 2, 2, 2}
};
const int RGESolver::WC8I_indices[DWC8I][4] = {
    {0, 0, 0, 1},
    {0, 0, 0, 2},
    {0, 0, 1, 2},
    {0, 1, 0, 1},
    {0, 1, 0, 2},
    {0, 1, 1, 1},
    {0, 1, 1, 2},
    {0, 1, 2, 1},
    {0, 1, 2, 2},
    {0, 2, 0, 2},
    {0, 2, 1, 2},
    {0, 2, 2, 2},
    {1, 1, 1, 2},
    {1, 2, 1, 2},
    {1, 2, 2, 2}
};



