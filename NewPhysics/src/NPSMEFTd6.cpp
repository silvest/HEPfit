/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6.h"
#include <limits>
#include <gsl/gsl_sf.h>
#include <boost/bind.hpp>
#include "gslpp_function_adapter.h"

const std::string NPSMEFTd6::NPSMEFTd6Vars[NNPSMEFTd6Vars]
        = {"CG", "CW", "CHG", "CHW", "CHB", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
    "CHL1_11", "CHL1_12r", "CHL1_13r", "CHL1_22", "CHL1_23r", "CHL1_33",
    "CHL1_12i", "CHL1_13i", "CHL1_23i",
    "CHL3_11", "CHL3_12r", "CHL3_13r", "CHL3_22", "CHL3_23r", "CHL3_33",
    "CHL3_12i", "CHL3_13i", "CHL3_23i",
    "CHe_11", "CHe_12r", "CHe_13r", "CHe_22", "CHe_23r", "CHe_33",
    "CHe_12i", "CHe_13i", "CHe_23i",
    "CHQ1_11", "CHQ1_12r", "CHQ1_13r", "CHQ1_22", "CHQ1_23r", "CHQ1_33",
    "CHQ1_12i", "CHQ1_13i", "CHQ1_23i",
    "CHQ3_11", "CHQ3_12r", "CHQ3_13r", "CHQ3_22", "CHQ3_23r", "CHQ3_33",
    "CHQ3_12i", "CHQ3_13i", "CHQ3_23i",
    "CHu_11", "CHu_12r", "CHu_13r", "CHu_22", "CHu_23r", "CHu_33",
    "CHu_12i", "CHu_13i", "CHu_23i",
    "CHd_11", "CHd_12r", "CHd_13r", "CHd_22", "CHd_23r", "CHd_33",
    "CHd_12i", "CHd_13i", "CHd_23i",
    "CHud_11r", "CHud_12r", "CHud_13r", "CHud_22r", "CHud_23r", "CHud_33r",
    "CHud_11i", "CHud_12i", "CHud_13i", "CHud_22i", "CHud_23i", "CHud_33i",
    "CeH_11r", "CeH_12r", "CeH_13r", "CeH_22r", "CeH_23r", "CeH_33r",
    "CeH_11i", "CeH_12i", "CeH_13i", "CeH_22i", "CeH_23i", "CeH_33i",
    "CuH_11r", "CuH_12r", "CuH_13r", "CuH_22r", "CuH_23r", "CuH_33r",
    "CuH_11i", "CuH_12i", "CuH_13i", "CuH_22i", "CuH_23i", "CuH_33i",
    "CdH_11r", "CdH_12r", "CdH_13r", "CdH_22r", "CdH_23r", "CdH_33r",
    "CdH_11i", "CdH_12i", "CdH_13i", "CdH_22i", "CdH_23i", "CdH_33i",
    "CuG_11r", "CuG_12r", "CuG_13r", "CuG_22r", "CuG_23r", "CuG_33r",
    "CuG_11i", "CuG_12i", "CuG_13i", "CuG_22i", "CuG_23i", "CuG_33i",
    "CuW_11r", "CuW_12r", "CuW_13r", "CuW_22r", "CuW_23r", "CuW_33r",
    "CuW_11i", "CuW_12i", "CuW_13i", "CuW_22i", "CuW_23i", "CuW_33i",
    "CuB_11r", "CuB_12r", "CuB_13r", "CuB_22r", "CuB_23r", "CuB_33r",
    "CuB_11i", "CuB_12i", "CuB_13i", "CuB_22i", "CuB_23i", "CuB_33i",
    "CLL_1111","CLL_1221","CLL_1122",
    "CLL_1133","CLL_1331",
    "CLQ1_1111","CLQ1_1122","CLQ1_2211","CLQ1_1221","CLQ1_2112",
    "CLQ1_1133","CLQ1_3311","CLQ1_1331","CLQ1_3113",
    "CLQ1_1123","CLQ1_2223","CLQ1_3323",
    "CLQ1_1132","CLQ1_2232","CLQ1_3332",        
    "CLQ3_1111","CLQ3_1122","CLQ3_2211","CLQ3_1221","CLQ3_2112",
    "CLQ3_1133","CLQ3_3311","CLQ3_1331","CLQ3_3113",
    "CLQ3_1123","CLQ3_2223","CLQ3_3323",
    "CLQ3_1132","CLQ3_2232","CLQ3_3332",
    "Cee_1111","Cee_1122","Cee_1133",
    "Ceu_1111","Ceu_1122","Ceu_2211","Ceu_1133","Ceu_2233","Ceu_3311",
    "Ced_1111","Ced_1122","Ced_2211","Ced_1133","Ced_3311",
    "Ced_1123","Ced_2223","Ced_3323",
    "Ced_1132","Ced_2232","Ced_3332",
    "CLe_1111","CLe_1122","CLe_2211","CLe_1133","CLe_3311",
    "CLu_1111","CLu_1122","CLu_2211","CLu_1133","CLu_2233","CLu_3311",
    "CLd_1111","CLd_1122","CLd_2211","CLd_1133","CLd_3311",
    "CLd_1123","CLd_2223","CLd_3323",
    "CLd_1132","CLd_2232","CLd_3332",
    "CQe_1111","CQe_1122","CQe_2211","CQe_1133","CQe_3311", 
    "CQe_2311","CQe_2322","CQe_2333",
    "CQe_3211","CQe_3222","CQe_3233",
    "CLedQ_11","CLedQ_22","CpLedQ_11","CpLedQ_22",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

const std::string NPSMEFTd6::NPSMEFTd6VarsRot[NNPSMEFTd6Vars]
        = {"CG", "CW", "CHG", "CHWHB_gaga", "CHWHB_gagaorth", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
    "CHL1_11", "CHL1_12r", "CHL1_13r", "CHL1_22", "CHL1_23r", "CHL1_33",
    "CHL1_12i", "CHL1_13i", "CHL1_23i",
    "CHL3_11", "CHL3_12r", "CHL3_13r", "CHL3_22", "CHL3_23r", "CHL3_33",
    "CHL3_12i", "CHL3_13i", "CHL3_23i",
    "CHe_11", "CHe_12r", "CHe_13r", "CHe_22", "CHe_23r", "CHe_33",
    "CHe_12i", "CHe_13i", "CHe_23i",
    "CHQ1_11", "CHQ1_12r", "CHQ1_13r", "CHQ1_22", "CHQ1_23r", "CHQ1_33",
    "CHQ1_12i", "CHQ1_13i", "CHQ1_23i",
    "CHQ3_11", "CHQ3_12r", "CHQ3_13r", "CHQ3_22", "CHQ3_23r", "CHQ3_33",
    "CHQ3_12i", "CHQ3_13i", "CHQ3_23i",
    "CHu_11", "CHu_12r", "CHu_13r", "CHu_22", "CHu_23r", "CHu_33",
    "CHu_12i", "CHu_13i", "CHu_23i",
    "CHd_11", "CHd_12r", "CHd_13r", "CHd_22", "CHd_23r", "CHd_33",
    "CHd_12i", "CHd_13i", "CHd_23i",
    "CHud_11r", "CHud_12r", "CHud_13r", "CHud_22r", "CHud_23r", "CHud_33r",
    "CHud_11i", "CHud_12i", "CHud_13i", "CHud_22i", "CHud_23i", "CHud_33i",
    "CeH_11r", "CeH_12r", "CeH_13r", "CeH_22r", "CeH_23r", "CeH_33r",
    "CeH_11i", "CeH_12i", "CeH_13i", "CeH_22i", "CeH_23i", "CeH_33i",
    "CuH_11r", "CuH_12r", "CuH_13r", "CuH_22r", "CuH_23r", "CuH_33r",
    "CuH_11i", "CuH_12i", "CuH_13i", "CuH_22i", "CuH_23i", "CuH_33i",
    "CdH_11r", "CdH_12r", "CdH_13r", "CdH_22r", "CdH_23r", "CdH_33r",
    "CdH_11i", "CdH_12i", "CdH_13i", "CdH_22i", "CdH_23i", "CdH_33i",
    "CuG_11r", "CuG_12r", "CuG_13r", "CuG_22r", "CuG_23r", "CuG_33r",
    "CuG_11i", "CuG_12i", "CuG_13i", "CuG_22i", "CuG_23i", "CuG_33i",
    "CuW_11r", "CuW_12r", "CuW_13r", "CuW_22r", "CuW_23r", "CuW_33r",
    "CuW_11i", "CuW_12i", "CuW_13i", "CuW_22i", "CuW_23i", "CuW_33i",
    "CuB_11r", "CuB_12r", "CuB_13r", "CuB_22r", "CuB_23r", "CuB_33r",
    "CuB_11i", "CuB_12i", "CuB_13i", "CuB_22i", "CuB_23i", "CuB_33i",
    "CLL_1111","CLL_1221","CLL_1122",
    "CLL_1133","CLL_1331",
    "CLQ1_1111","CLQ1_1122","CLQ1_2211","CLQ1_1221","CLQ1_2112",
    "CLQ1_1133","CLQ1_3311","CLQ1_1331","CLQ1_3113",
    "CLQ1_1123","CLQ1_2223","CLQ1_3323",
    "CLQ1_1132","CLQ1_2232","CLQ1_3332",
    "CLQ3_1111","CLQ3_1122","CLQ3_2211","CLQ3_1221","CLQ3_2112",
    "CLQ3_1133","CLQ3_3311","CLQ3_1331","CLQ3_3113",
    "CLQ3_1123","CLQ3_2223","CLQ3_3323",
    "CLQ3_1132","CLQ3_2232","CLQ3_3332",
    "Cee_1111","Cee_1122","Cee_1133",
    "Ceu_1111","Ceu_1122","Ceu_2211","Ceu_1133","Ceu_2233","Ceu_3311",
    "Ced_1111","Ced_1122","Ced_2211","Ced_1133","Ced_3311",
    "Ced_1123","Ced_2223","Ced_3323",
    "Ced_1132","Ced_2232","Ced_3332",
    "CLe_1111","CLe_1122","CLe_2211","CLe_1133","CLe_3311",
    "CLu_1111","CLu_1122","CLu_2211","CLu_1133","CLu_2233","CLu_3311",
    "CLd_1111","CLd_1122","CLd_2211","CLd_1133","CLd_3311",
    "CLd_1123","CLd_2223","CLd_3323",
    "CLd_1132","CLd_2232","CLd_3332",
    "CQe_1111","CQe_1122","CQe_2211","CQe_1133","CQe_3311",
    "CQe_2311","CQe_2322","CQe_2333",
    "CQe_3211","CQe_3222","CQe_3233",
    "CLedQ_11","CLedQ_22","CpLedQ_11","CpLedQ_22",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

const std::string NPSMEFTd6::NPSMEFTd6Vars_LFU_QFU[NNPSMEFTd6Vars_LFU_QFU]
        = {"CG", "CW", "CHG", "CHW", "CHB", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
    "CHL1", "CHL3", "CHe", "CHQ1", "CHQ3", "CHu", "CHd", "CHud_r", "CHud_i",
    "CeH_11r", "CeH_22r", "CeH_33r", "CeH_11i", "CeH_22i", "CeH_33i", 
    "CuH_11r", "CuH_22r", "CuH_33r", "CuH_11i", "CuH_22i", "CuH_33i", 
    "CdH_11r", "CdH_22r", "CdH_33r", "CdH_11i", "CdH_22i", "CdH_33i",
    "CuG_r", "CuG_i", "CuW_r", "CuW_i", "CuB_r", "CuB_i",
    "CLL", "CLQ1", "CLQ3",
    "Cee", "Ceu", "Ced", "CLe", "CLu", "CLd", "CQe",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

const std::string NPSMEFTd6::NPSMEFTd6VarsRot_LFU_QFU[NNPSMEFTd6Vars_LFU_QFU]
        = {"CG", "CW", "CHG", "CHWHB_gaga", "CHWHB_gagaorth", "CDHB", "CDHW", "CHWB", "CHD", "CHbox", "CH",
    "CHL1", "CHL3", "CHe", "CHQ1", "CHQ3", "CHu", "CHd", "CHud_r", "CHud_i",
    "CeH_11r", "CeH_22r", "CeH_33r", "CeH_11i", "CeH_22i", "CeH_33i", 
    "CuH_11r", "CuH_22r", "CuH_33r", "CuH_11i", "CuH_22i", "CuH_33i", 
    "CdH_11r", "CdH_22r", "CdH_33r", "CdH_11i", "CdH_22i", "CdH_33i",
    "CuG_r", "CuG_i", "CuW_r", "CuW_i", "CuB_r", "CuB_i",
    "CLL", "CLQ1", "CLQ3",
    "Cee", "Ceu", "Ced", "CLe", "CLu", "CLd", "CQe",
    "Lambda_NP",
    "BrHinv","BrHexo",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eVBF_2_Hbox", "eVBF_2_HQ1_11", "eVBF_2_Hu_11", "eVBF_2_Hd_11", "eVBF_2_HQ3_11",
    "eVBF_2_HD", "eVBF_2_HB", "eVBF_2_HW", "eVBF_2_HWB", "eVBF_2_HG", "eVBF_2_DHB",
    "eVBF_2_DHW", "eVBF_2_DeltaGF",
    "eVBF_78_Hbox", "eVBF_78_HQ1_11", "eVBF_78_Hu_11", "eVBF_78_Hd_11", "eVBF_78_HQ3_11",
    "eVBF_78_HD", "eVBF_78_HB", "eVBF_78_HW", "eVBF_78_HWB", "eVBF_78_HG", "eVBF_78_DHB",
    "eVBF_78_DHW", "eVBF_78_DeltaGF",
    "eVBF_1314_Hbox", "eVBF_1314_HQ1_11", "eVBF_1314_Hu_11", "eVBF_1314_Hd_11", "eVBF_1314_HQ3_11",
    "eVBF_1314_HD", "eVBF_1314_HB", "eVBF_1314_HW", "eVBF_1314_HWB", "eVBF_1314_HG", "eVBF_1314_DHB",
    "eVBF_1314_DHW", "eVBF_1314_DeltaGF",
    "eWH_2_Hbox", "eWH_2_HQ3_11", "eWH_2_HD", "eWH_2_HW", "eWH_2_HWB", "eWH_2_DHW", "eWH_2_DeltaGF",
    "eWH_78_Hbox", "eWH_78_HQ3_11", "eWH_78_HD", "eWH_78_HW", "eWH_78_HWB", "eWH_78_DHW", "eWH_78_DeltaGF",
    "eWH_1314_Hbox", "eWH_1314_HQ3_11", "eWH_1314_HD", "eWH_1314_HW", "eWH_1314_HWB", "eWH_1314_DHW", "eWH_1314_DeltaGF",
    "eZH_2_Hbox", "eZH_2_HQ1_11", "eZH_2_Hu_11", "eZH_2_Hd_11", "eZH_2_HQ3_11", "eZH_2_HD", "eZH_2_HB", "eZH_2_HW", "eZH_2_HWB", "eZH_2_DHB", "eZH_2_DHW", "eZH_2_DeltaGF",
    "eZH_78_Hbox", "eZH_78_HQ1_11", "eZH_78_Hu_11", "eZH_78_Hd_11", "eZH_78_HQ3_11", "eZH_78_HD", "eZH_78_HB", "eZH_78_HW", "eZH_78_HWB", "eZH_78_DHB", "eZH_78_DHW", "eZH_78_DeltaGF",
    "eZH_1314_Hbox", "eZH_1314_HQ1_11", "eZH_1314_Hu_11", "eZH_1314_Hd_11", "eZH_1314_HQ3_11", "eZH_1314_HD", "eZH_1314_HB", "eZH_1314_HW", "eZH_1314_HWB", "eZH_1314_DHB", "eZH_1314_DHW", "eZH_1314_DeltaGF",
    "ettH_2_HG", "ettH_2_G", "ettH_2_uG_33r", "ettH_2_DeltagHt",
    "ettH_78_HG", "ettH_78_G", "ettH_78_uG_33r", "ettH_78_DeltagHt",
    "ettH_1314_HG", "ettH_1314_G", "ettH_1314_uG_33r", "ettH_1314_DeltagHt"};

NPSMEFTd6::NPSMEFTd6(const bool FlagLeptonUniversal_in, const bool FlagQuarkUniversal_in)
: NPbase(), NPSMEFTd6M(*this), FlagLeptonUniversal(FlagLeptonUniversal_in), FlagQuarkUniversal(FlagQuarkUniversal_in)
{
    if ((!FlagLeptonUniversal && FlagQuarkUniversal)
            || (FlagLeptonUniversal && !FlagQuarkUniversal))
        throw std::runtime_error("Invalid arguments for NPSMEFTd6::NPSMEFTd6()");

    FlagQuadraticTerms = false;
    FlagRotateCHWCHB = false;
    FlagPartialQFU = false;
    setModelLinearized();
    
    SMM.setObj((StandardModelMatching&) NPSMEFTd6M.getObj());

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CG", boost::cref(CG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CW", boost::cref(CW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHG", boost::cref(CHG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHW", boost::cref(CHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHB", boost::cref(CHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWHB_gaga", boost::cref(CHWHB_gaga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWHB_gagaorth", boost::cref(CHWHB_gagaorth)));        
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CDHB", boost::cref(CDHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CDHW", boost::cref(CDHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWB", boost::cref(CHWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHD", boost::cref(CHD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHbox", boost::cref(CHbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CH", boost::cref(CH)));
    if (FlagLeptonUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11r", boost::cref(CeH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22r", boost::cref(CeH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33r", boost::cref(CeH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11i", boost::cref(CeH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22i", boost::cref(CeH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33i", boost::cref(CeH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL", boost::cref(CLL_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee", boost::cref(Cee_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe", boost::cref(CLe_1111)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_11", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12r", boost::cref(CHL1_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13r", boost::cref(CHL1_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_22", boost::cref(CHL1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23r", boost::cref(CHL1_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_33", boost::cref(CHL1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12i", boost::cref(CHL1_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13i", boost::cref(CHL1_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23i", boost::cref(CHL1_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_11", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12r", boost::cref(CHL3_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13r", boost::cref(CHL3_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_22", boost::cref(CHL3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23r", boost::cref(CHL3_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_33", boost::cref(CHL3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12i", boost::cref(CHL3_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13i", boost::cref(CHL3_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23i", boost::cref(CHL3_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_11", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12r", boost::cref(CHe_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13r", boost::cref(CHe_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_22", boost::cref(CHe_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23r", boost::cref(CHe_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_33", boost::cref(CHe_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12i", boost::cref(CHe_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13i", boost::cref(CHe_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23i", boost::cref(CHe_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11r", boost::cref(CeH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_12r", boost::cref(CeH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_13r", boost::cref(CeH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22r", boost::cref(CeH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_23r", boost::cref(CeH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33r", boost::cref(CeH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11i", boost::cref(CeH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_12i", boost::cref(CeH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_13i", boost::cref(CeH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22i", boost::cref(CeH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_23i", boost::cref(CeH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33i", boost::cref(CeH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1111", boost::cref(CLL_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1221", boost::cref(CLL_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1122", boost::cref(CLL_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1331", boost::cref(CLL_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1133", boost::cref(CLL_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1111", boost::cref(Cee_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1122", boost::cref(Cee_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1133", boost::cref(Cee_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1111", boost::cref(CLe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1122", boost::cref(CLe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_2211", boost::cref(CLe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1133", boost::cref(CLe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_3311", boost::cref(CLe_3311)));       
    }
    if (FlagQuarkUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_r", boost::cref(CHud_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_i", boost::cref(CHud_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22r", boost::cref(CuH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33r", boost::cref(CuH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22i", boost::cref(CuH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33i", boost::cref(CuH_33i)));        
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22r", boost::cref(CdH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33r", boost::cref(CdH_33r)));        
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11i", boost::cref(CdH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22i", boost::cref(CdH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33i", boost::cref(CdH_33i)));        
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_r", boost::cref(CuG_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_i", boost::cref(CuG_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_r", boost::cref(CuW_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_i", boost::cref(CuW_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_r", boost::cref(CuB_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_i", boost::cref(CuB_11i)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_11", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12r", boost::cref(CHQ1_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13r", boost::cref(CHQ1_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_22", boost::cref(CHQ1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23r", boost::cref(CHQ1_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_33", boost::cref(CHQ1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12i", boost::cref(CHQ1_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13i", boost::cref(CHQ1_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23i", boost::cref(CHQ1_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_11", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12r", boost::cref(CHQ3_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13r", boost::cref(CHQ3_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_22", boost::cref(CHQ3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23r", boost::cref(CHQ3_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_33", boost::cref(CHQ3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12i", boost::cref(CHQ3_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13i", boost::cref(CHQ3_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23i", boost::cref(CHQ3_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_11", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12r", boost::cref(CHu_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13r", boost::cref(CHu_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_22", boost::cref(CHu_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23r", boost::cref(CHu_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_33", boost::cref(CHu_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12i", boost::cref(CHu_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13i", boost::cref(CHu_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23i", boost::cref(CHu_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_11", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12r", boost::cref(CHd_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13r", boost::cref(CHd_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_22", boost::cref(CHd_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23r", boost::cref(CHd_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_33", boost::cref(CHd_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12i", boost::cref(CHd_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13i", boost::cref(CHd_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23i", boost::cref(CHd_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_11r", boost::cref(CHud_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_12r", boost::cref(CHud_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_13r", boost::cref(CHud_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_22r", boost::cref(CHud_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_23r", boost::cref(CHud_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_33r", boost::cref(CHud_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_11i", boost::cref(CHud_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_12i", boost::cref(CHud_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_13i", boost::cref(CHud_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_22i", boost::cref(CHud_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_23i", boost::cref(CHud_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_33i", boost::cref(CHud_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_12r", boost::cref(CuH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_13r", boost::cref(CuH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22r", boost::cref(CuH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_23r", boost::cref(CuH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33r", boost::cref(CuH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_12i", boost::cref(CuH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_13i", boost::cref(CuH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22i", boost::cref(CuH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_23i", boost::cref(CuH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33i", boost::cref(CuH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_12r", boost::cref(CdH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_13r", boost::cref(CdH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22r", boost::cref(CdH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_23r", boost::cref(CdH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33r", boost::cref(CdH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11i", boost::cref(CdH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_12i", boost::cref(CdH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_13i", boost::cref(CdH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22i", boost::cref(CdH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_23i", boost::cref(CdH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33i", boost::cref(CdH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_11r", boost::cref(CuG_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_12r", boost::cref(CuG_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_13r", boost::cref(CuG_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_22r", boost::cref(CuG_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_23r", boost::cref(CuG_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_33r", boost::cref(CuG_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_11i", boost::cref(CuG_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_12i", boost::cref(CuG_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_13i", boost::cref(CuG_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_22i", boost::cref(CuG_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_23i", boost::cref(CuG_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_33i", boost::cref(CuG_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_11r", boost::cref(CuW_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_12r", boost::cref(CuW_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_13r", boost::cref(CuW_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_22r", boost::cref(CuW_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_23r", boost::cref(CuW_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_33r", boost::cref(CuW_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_11i", boost::cref(CuW_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_12i", boost::cref(CuW_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_13i", boost::cref(CuW_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_22i", boost::cref(CuW_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_23i", boost::cref(CuW_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_33i", boost::cref(CuW_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_11r", boost::cref(CuB_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_12r", boost::cref(CuB_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_13r", boost::cref(CuB_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_22r", boost::cref(CuB_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_23r", boost::cref(CuB_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_33r", boost::cref(CuB_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_11i", boost::cref(CuB_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_12i", boost::cref(CuB_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_13i", boost::cref(CuB_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_22i", boost::cref(CuB_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_23i", boost::cref(CuB_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_33i", boost::cref(CuB_33i)));
    }
    if(FlagLeptonUniversal && FlagQuarkUniversal){
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1", boost::cref(CLQ1_1111)));  
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3", boost::cref(CLQ3_1111))); 
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu", boost::cref(Ceu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced", boost::cref(Ced_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu", boost::cref(CLu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd", boost::cref(CLd_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe", boost::cref(CQe_1111)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1111", boost::cref(CLQ1_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1122", boost::cref(CLQ1_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2211", boost::cref(CLQ1_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1221", boost::cref(CLQ1_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2112", boost::cref(CLQ1_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1133", boost::cref(CLQ1_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3311", boost::cref(CLQ1_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1331", boost::cref(CLQ1_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3113", boost::cref(CLQ1_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1123", boost::cref(CLQ1_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2223", boost::cref(CLQ1_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3323", boost::cref(CLQ1_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1132", boost::cref(CLQ1_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2232", boost::cref(CLQ1_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3332", boost::cref(CLQ1_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1111", boost::cref(CLQ3_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1122", boost::cref(CLQ3_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2211", boost::cref(CLQ3_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1221", boost::cref(CLQ3_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2112", boost::cref(CLQ3_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1133", boost::cref(CLQ3_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3311", boost::cref(CLQ3_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1331", boost::cref(CLQ3_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3113", boost::cref(CLQ3_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1123", boost::cref(CLQ3_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2223", boost::cref(CLQ3_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3323", boost::cref(CLQ3_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1132", boost::cref(CLQ3_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2232", boost::cref(CLQ3_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3332", boost::cref(CLQ3_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1111", boost::cref(Ceu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1122", boost::cref(Ceu_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_2211", boost::cref(Ceu_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1133", boost::cref(Ceu_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_2233", boost::cref(Ceu_2233)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_3311", boost::cref(Ceu_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1111", boost::cref(Ced_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1122", boost::cref(Ced_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2211", boost::cref(Ced_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1133", boost::cref(Ced_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3311", boost::cref(Ced_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1123", boost::cref(Ced_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2223", boost::cref(Ced_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3323", boost::cref(Ced_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1132", boost::cref(Ced_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2232", boost::cref(Ced_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3332", boost::cref(Ced_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1111", boost::cref(CLu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1122", boost::cref(CLu_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_2211", boost::cref(CLu_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1133", boost::cref(CLu_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_2233", boost::cref(CLu_2233)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_3311", boost::cref(CLu_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1111", boost::cref(CLd_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1122", boost::cref(CLd_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2211", boost::cref(CLd_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1133", boost::cref(CLd_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3311", boost::cref(CLd_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1123", boost::cref(CLd_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2223", boost::cref(CLd_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3323", boost::cref(CLd_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1132", boost::cref(CLd_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2232", boost::cref(CLd_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3332", boost::cref(CLd_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1111", boost::cref(CQe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1122", boost::cref(CQe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2211", boost::cref(CQe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1133", boost::cref(CQe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3311", boost::cref(CQe_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2311", boost::cref(CQe_2311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2322", boost::cref(CQe_2322)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2333", boost::cref(CQe_2333)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3211", boost::cref(CQe_3211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3222", boost::cref(CQe_3222)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3233", boost::cref(CQe_3233)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLedQ_11", boost::cref(CLedQ_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLedQ_22", boost::cref(CLedQ_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CpLedQ_11", boost::cref(CpLedQ_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CpLedQ_22", boost::cref(CpLedQ_22)));
    }
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(Lambda_NP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHexo", boost::cref(BrHexo)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eggFint", boost::cref(eggFint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eggFpar", boost::cref(eggFpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettHint", boost::cref(ettHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettHpar", boost::cref(ettHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBFint", boost::cref(eVBFint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBFpar", boost::cref(eVBFpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWHint", boost::cref(eWHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWHpar", boost::cref(eWHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZHint", boost::cref(eZHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZHpar", boost::cref(eZHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeWBFint", boost::cref(eeeWBFint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeWBFpar", boost::cref(eeeWBFpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeZHint", boost::cref(eeeZHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeeZHpar", boost::cref(eeeZHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeettHint", boost::cref(eeettHint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eeettHpar", boost::cref(eeettHpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHggint", boost::cref(eHggint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHggpar", boost::cref(eHggpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHWWint", boost::cref(eHWWint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHWWpar", boost::cref(eHWWpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZZint", boost::cref(eHZZint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZZpar", boost::cref(eHZZpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZgaint", boost::cref(eHZgaint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHZgapar", boost::cref(eHZgapar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHgagaint", boost::cref(eHgagaint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHgagapar", boost::cref(eHgagapar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHmumuint", boost::cref(eHmumuint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHmumupar", boost::cref(eHmumupar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHtautauint", boost::cref(eHtautauint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHtautaupar", boost::cref(eHtautaupar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHccint", boost::cref(eHccint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHccpar", boost::cref(eHccpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHbbint", boost::cref(eHbbint)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eHbbpar", boost::cref(eHbbpar)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_Hbox", boost::cref(eVBF_2_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HQ1_11", boost::cref(eVBF_2_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_Hu_11", boost::cref(eVBF_2_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_Hd_11", boost::cref(eVBF_2_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HQ3_11", boost::cref(eVBF_2_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HD", boost::cref(eVBF_2_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HB", boost::cref(eVBF_2_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HW", boost::cref(eVBF_2_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HWB", boost::cref(eVBF_2_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_HG", boost::cref(eVBF_2_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_DHB", boost::cref(eVBF_2_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_DHW", boost::cref(eVBF_2_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_2_DeltaGF", boost::cref(eVBF_2_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_Hbox", boost::cref(eVBF_78_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HQ1_11", boost::cref(eVBF_78_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_Hu_11", boost::cref(eVBF_78_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_Hd_11", boost::cref(eVBF_78_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HQ3_11", boost::cref(eVBF_78_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HD", boost::cref(eVBF_78_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HB", boost::cref(eVBF_78_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HW", boost::cref(eVBF_78_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HWB", boost::cref(eVBF_78_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_HG", boost::cref(eVBF_78_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_DHB", boost::cref(eVBF_78_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_DHW", boost::cref(eVBF_78_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_78_DeltaGF", boost::cref(eVBF_78_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_Hbox", boost::cref(eVBF_1314_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HQ1_11", boost::cref(eVBF_1314_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_Hu_11", boost::cref(eVBF_1314_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_Hd_11", boost::cref(eVBF_1314_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HQ3_11", boost::cref(eVBF_1314_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HD", boost::cref(eVBF_1314_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HB", boost::cref(eVBF_1314_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HW", boost::cref(eVBF_1314_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HWB", boost::cref(eVBF_1314_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_HG", boost::cref(eVBF_1314_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_DHB", boost::cref(eVBF_1314_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_DHW", boost::cref(eVBF_1314_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF_1314_DeltaGF", boost::cref(eVBF_1314_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_Hbox", boost::cref(eWH_2_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HQ3_11", boost::cref(eWH_2_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HD", boost::cref(eWH_2_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HW", boost::cref(eWH_2_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_HWB", boost::cref(eWH_2_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_DHW", boost::cref(eWH_2_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_2_DeltaGF", boost::cref(eWH_2_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_Hbox", boost::cref(eWH_78_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HQ3_11", boost::cref(eWH_78_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HD", boost::cref(eWH_78_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HW", boost::cref(eWH_78_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_HWB", boost::cref(eWH_78_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_DHW", boost::cref(eWH_78_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_78_DeltaGF", boost::cref(eWH_78_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_Hbox", boost::cref(eWH_1314_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HQ3_11", boost::cref(eWH_1314_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HD", boost::cref(eWH_1314_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HW", boost::cref(eWH_1314_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_HWB", boost::cref(eWH_1314_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_DHW", boost::cref(eWH_1314_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH_1314_DeltaGF", boost::cref(eWH_1314_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_Hbox", boost::cref(eZH_2_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HQ1_11", boost::cref(eZH_2_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_Hu_11", boost::cref(eZH_2_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_Hd_11", boost::cref(eZH_2_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HQ3_11", boost::cref(eZH_2_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HD", boost::cref(eZH_2_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HB", boost::cref(eZH_2_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HW", boost::cref(eZH_2_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_HWB", boost::cref(eZH_2_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_DHB", boost::cref(eZH_2_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_DHW", boost::cref(eZH_2_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_2_DeltaGF", boost::cref(eZH_2_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_Hbox", boost::cref(eZH_78_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HQ1_11", boost::cref(eZH_78_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_Hu_11", boost::cref(eZH_78_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_Hd_11", boost::cref(eZH_78_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HQ3_11", boost::cref(eZH_78_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HD", boost::cref(eZH_78_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HB", boost::cref(eZH_78_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HW", boost::cref(eZH_78_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_HWB", boost::cref(eZH_78_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_DHB", boost::cref(eZH_78_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_DHW", boost::cref(eZH_78_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_78_DeltaGF", boost::cref(eZH_78_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_Hbox", boost::cref(eZH_1314_Hbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HQ1_11", boost::cref(eZH_1314_HQ1_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_Hu_11", boost::cref(eZH_1314_Hu_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_Hd_11", boost::cref(eZH_1314_Hd_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HQ3_11", boost::cref(eZH_1314_HQ3_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HD", boost::cref(eZH_1314_HD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HB", boost::cref(eZH_1314_HB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HW", boost::cref(eZH_1314_HW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_HWB", boost::cref(eZH_1314_HWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_DHB", boost::cref(eZH_1314_DHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_DHW", boost::cref(eZH_1314_DHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH_1314_DeltaGF", boost::cref(eZH_1314_DeltaGF)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_HG", boost::cref(ettH_2_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_G", boost::cref(ettH_2_G)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_uG_33r", boost::cref(ettH_2_uG_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_2_DeltagHt", boost::cref(ettH_2_DeltagHt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_HG", boost::cref(ettH_78_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_G", boost::cref(ettH_78_G)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_uG_33r", boost::cref(ettH_78_uG_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_78_DeltagHt", boost::cref(ettH_78_DeltagHt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_HG", boost::cref(ettH_1314_HG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_G", boost::cref(ettH_1314_G)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_uG_33r", boost::cref(ettH_1314_uG_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH_1314_DeltagHt", boost::cref(ettH_1314_DeltagHt)));
    
    if (FlagLeptonUniversal) {
        CeH_12r = 0.0;
        CeH_13r = 0.0;
        CeH_23r = 0.0;
        CeH_12i = 0.0;
        CeH_13i = 0.0;
        CeH_23i = 0.0;  
        
//  bsll/sbll entries only interesting (for the moment) if non-lepton universal. Set to 0 otherwise
        CLQ1_1123 = 0.0;
        CLQ1_2223 = 0.0;
        CLQ1_3323 = 0.0;
        CLQ1_1132 = 0.0;
        CLQ1_2232 = 0.0;
        CLQ1_3332 = 0.0;
        
        CLQ3_1123 = 0.0;
        CLQ3_2223 = 0.0;
        CLQ3_3323 = 0.0;
        CLQ3_1132 = 0.0;
        CLQ3_2232 = 0.0;
        CLQ3_3332 = 0.0;
        
        Ced_1123 = 0.0;
        Ced_2223 = 0.0;
        Ced_3323 = 0.0;
        Ced_1132 = 0.0;
        Ced_2232 = 0.0;
        Ced_3332 = 0.0;  
        
        CLd_1123 = 0.0;
        CLd_2223 = 0.0;
        CLd_3323 = 0.0;
        CLd_1132 = 0.0;
        CLd_2232 = 0.0;
        CLd_3332 = 0.0;
        
        CQe_2311 = 0.0;
        CQe_2322 = 0.0;
        CQe_2333 = 0.0;
        CQe_3211 = 0.0;
        CQe_3222 = 0.0;
        CQe_3233 = 0.0;
    }
    if (FlagQuarkUniversal) {
        CuH_12r = 0.0;
        CuH_13r = 0.0;
        CuH_23r = 0.0;
        CuH_12i = 0.0;
        CuH_13i = 0.0;
        CuH_23i = 0.0;

        CdH_12r = 0.0;
        CdH_13r = 0.0;
        CdH_23r = 0.0;
        CdH_12i = 0.0;
        CdH_13i = 0.0;
        CdH_23i = 0.0;        
    }
}

bool NPSMEFTd6::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

    LambdaNP2 = Lambda_NP * Lambda_NP;
    v2_over_LambdaNP2 = v() * v() / LambdaNP2;
    cW_tree = Mw_tree() / Mz;
    cW2_tree = cW_tree * cW_tree;
    sW2_tree = 1.0 - cW2_tree;
    sW_tree = sqrt(sW2_tree);
      
    if (FlagRotateCHWCHB) {
        CHW = sW2_tree * CHWHB_gaga - cW2_tree * CHWHB_gagaorth;
        CHB = cW2_tree * CHWHB_gaga + sW2_tree * CHWHB_gagaorth;
    } else {
        CHWHB_gaga = sW2_tree * CHW + cW2_tree * CHB;
        CHWHB_gagaorth = - cW2_tree * CHW + sW2_tree * CHB;
    }

    delta_ZZ = (cW2_tree * CHW + sW2_tree * CHB + sW_tree * cW_tree * CHWB) * v2_over_LambdaNP2;
    delta_AA = (sW2_tree * CHW + cW2_tree * CHB - sW_tree * cW_tree * CHWB) * v2_over_LambdaNP2;
    delta_AZ = 2.0 * sW_tree * cW_tree * (CHW - CHB) * v2_over_LambdaNP2
            - (cW2_tree - sW2_tree) * CHWB * v2_over_LambdaNP2;
    delta_h = (-CHD / 4.0 + CHbox) * v2_over_LambdaNP2;
    
    NPSMEFTd6M.getObj().updateNPSMEFTd6Parameters();

    return (true);
}

void NPSMEFTd6::setParameter(const std::string name, const double& value)
{
    if (name.compare("CG") == 0)
        CG = value;
    else if (name.compare("CW") == 0)
        CW = value;
    else if (name.compare("CHG") == 0)
        CHG = value;
    else if (name.compare("CHW") == 0)
        CHW = value;
    else if (name.compare("CHB") == 0)
        CHB = value;
    else if (name.compare("CHWHB_gaga") == 0)
        CHWHB_gaga = value;
    else if (name.compare("CHWHB_gagaorth") == 0)
        CHWHB_gagaorth = value;
    else if (name.compare("CDHB") == 0)
        CDHB = value;
    else if (name.compare("CDHW") == 0)
        CDHW = value;
    else if (name.compare("CHWB") == 0)
        CHWB = value;
    else if (name.compare("CHD") == 0)
        CHD = value;
    else if (name.compare("CHbox") == 0)
        CHbox = value;
    else if (name.compare("CH") == 0)
        CH = value;
    else if (name.compare("CHL1_11") == 0)
        CHL1_11 = value;
    else if (name.compare("CHL1_12r") == 0)
        CHL1_12r = value;
    else if (name.compare("CHL1_13r") == 0)
        CHL1_13r = value;
    else if (name.compare("CHL1_22") == 0)
        CHL1_22 = value;
    else if (name.compare("CHL1_23r") == 0)
        CHL1_23r = value;
    else if (name.compare("CHL1_33") == 0)
        CHL1_33 = value;
    else if (name.compare("CHL1_12i") == 0)
        CHL1_12i = value;
    else if (name.compare("CHL1_13i") == 0)
        CHL1_13i = value;
    else if (name.compare("CHL1_23i") == 0)
        CHL1_23i = value;
    else if (name.compare("CHL1") == 0) {
        CHL1_11 = value;
        CHL1_12r = 0.0;
        CHL1_13r = 0.0;
        CHL1_22 = value;
        CHL1_23r = 0.0;
        CHL1_33 = value;
        CHL1_12i = 0.0;
        CHL1_13i = 0.0;
        CHL1_23i = 0.0;
    } else if (name.compare("CHL3_11") == 0)
        CHL3_11 = value;
    else if (name.compare("CHL3_12r") == 0)
        CHL3_12r = value;
    else if (name.compare("CHL3_13r") == 0)
        CHL3_13r = value;
    else if (name.compare("CHL3_22") == 0)
        CHL3_22 = value;
    else if (name.compare("CHL3_23r") == 0)
        CHL3_23r = value;
    else if (name.compare("CHL3_33") == 0)
        CHL3_33 = value;
    else if (name.compare("CHL3_12i") == 0)
        CHL3_12i = value;
    else if (name.compare("CHL3_13i") == 0)
        CHL3_13i = value;
    else if (name.compare("CHL3_23i") == 0)
        CHL3_23i = value;
    else if (name.compare("CHL3") == 0) {
        CHL3_11 = value;
        CHL3_12r = 0.0;
        CHL3_13r = 0.0;
        CHL3_22 = value;
        CHL3_23r = 0.0;
        CHL3_33 = value;
        CHL3_12i = 0.0;
        CHL3_13i = 0.0;
        CHL3_23i = 0.0;
    } else if (name.compare("CHe_11") == 0)
        CHe_11 = value;
    else if (name.compare("CHe_12r") == 0)
        CHe_12r = value;
    else if (name.compare("CHe_13r") == 0)
        CHe_13r = value;
    else if (name.compare("CHe_22") == 0)
        CHe_22 = value;
    else if (name.compare("CHe_23r") == 0)
        CHe_23r = value;
    else if (name.compare("CHe_33") == 0)
        CHe_33 = value;
    else if (name.compare("CHe_12i") == 0)
        CHe_12i = value;
    else if (name.compare("CHe_13i") == 0)
        CHe_13i = value;
    else if (name.compare("CHe_23i") == 0)
        CHe_23i = value;
    else if (name.compare("CHe") == 0) {
        CHe_11 = value;
        CHe_12r = 0.0;
        CHe_13r = 0.0;
        CHe_22 = value;
        CHe_23r = 0.0;
        CHe_33 = value;
        CHe_12i = 0.0;
        CHe_13i = 0.0;
        CHe_23i = 0.0;
    } else if (name.compare("CHQ1_11") == 0) {
        CHQ1_11 = value;
        if (FlagPartialQFU){
            CHQ1_22 = value;
        }
    } else if (name.compare("CHQ1_12r") == 0)
        CHQ1_12r = value;
    else if (name.compare("CHQ1_13r") == 0)
        CHQ1_13r = value;
    else if (name.compare("CHQ1_22") == 0) {
        if (!FlagPartialQFU){
            CHQ1_22 = value;
        }
    } else if (name.compare("CHQ1_23r") == 0)
        CHQ1_23r = value;
    else if (name.compare("CHQ1_33") == 0)
        CHQ1_33 = value;
    else if (name.compare("CHQ1_12i") == 0)
        CHQ1_12i = value;
    else if (name.compare("CHQ1_13i") == 0)
        CHQ1_13i = value;
    else if (name.compare("CHQ1_23i") == 0)
        CHQ1_23i = value;
    else if (name.compare("CHQ1") == 0) {
        CHQ1_11 = value;
        CHQ1_12r = 0.0;
        CHQ1_13r = 0.0;
        CHQ1_22 = value;
        CHQ1_23r = 0.0;
        CHQ1_33 = value;
        CHQ1_12i = 0.0;
        CHQ1_13i = 0.0;
        CHQ1_23i = 0.0;
    } else if (name.compare("CHQ3_11") == 0){
        CHQ3_11 = value;
        if (FlagPartialQFU){
            CHQ3_22 = value;
        }
    } else if (name.compare("CHQ3_12r") == 0)
        CHQ3_12r = value;
    else if (name.compare("CHQ3_13r") == 0)
        CHQ3_13r = value;
    else if (name.compare("CHQ3_22") == 0){
        if (!FlagPartialQFU){
            CHQ3_22 = value;
        }
    } else if (name.compare("CHQ3_23r") == 0)
        CHQ3_23r = value;
    else if (name.compare("CHQ3_33") == 0)
        CHQ3_33 = value;
    else if (name.compare("CHQ3_12i") == 0)
        CHQ3_12i = value;
    else if (name.compare("CHQ3_13i") == 0)
        CHQ3_13i = value;
    else if (name.compare("CHQ3_23i") == 0)
        CHQ3_23i = value;
    else if (name.compare("CHQ3") == 0) {
        CHQ3_11 = value;
        CHQ3_12r = 0.0;
        CHQ3_13r = 0.0;
        CHQ3_22 = value;
        CHQ3_23r = 0.0;
        CHQ3_33 = value;
        CHQ3_12i = 0.0;
        CHQ3_13i = 0.0;
        CHQ3_23i = 0.0;
    } else if (name.compare("CHu_11") == 0){
        CHu_11 = value;
        if (FlagPartialQFU){
            CHu_22 = value;
        }
    } else if (name.compare("CHu_12r") == 0)
        CHu_12r = value;
    else if (name.compare("CHu_13r") == 0)
        CHu_13r = value;
    else if (name.compare("CHu_22") == 0){
        if (!FlagPartialQFU){
            CHu_22 = value;
        }
    } else if (name.compare("CHu_23r") == 0)
        CHu_23r = value;
    else if (name.compare("CHu_33") == 0)
        CHu_33 = value;
    else if (name.compare("CHu_12i") == 0)
        CHu_12i = value;
    else if (name.compare("CHu_13i") == 0)
        CHu_13i = value;
    else if (name.compare("CHu_23i") == 0)
        CHu_23i = value;
    else if (name.compare("CHu") == 0) {
        CHu_11 = value;
        CHu_12r = 0.0;
        CHu_13r = 0.0;
        CHu_22 = value;
        CHu_23r = 0.0;
        CHu_33 = value;
        CHu_12i = 0.0;
        CHu_13i = 0.0;
        CHu_23i = 0.0;
    } else if (name.compare("CHd_11") == 0){
        CHd_11 = value;
        if (FlagPartialQFU){
            CHd_22 = value;
        }
    } else if (name.compare("CHd_12r") == 0)
        CHd_12r = value;
    else if (name.compare("CHd_13r") == 0)
        CHd_13r = value;
    else if (name.compare("CHd_22") == 0){
        if (!FlagPartialQFU){
            CHd_22 = value;
        }
    } else if (name.compare("CHd_23r") == 0)
        CHd_23r = value;
    else if (name.compare("CHd_33") == 0)
        CHd_33 = value;
    else if (name.compare("CHd_12i") == 0)
        CHd_12i = value;
    else if (name.compare("CHd_13i") == 0)
        CHd_13i = value;
    else if (name.compare("CHd_23i") == 0)
        CHd_23i = value;
    else if (name.compare("CHd") == 0) {
        CHd_11 = value;
        CHd_12r = 0.0;
        CHd_13r = 0.0;
        CHd_22 = value;
        CHd_23r = 0.0;
        CHd_33 = value;
        CHd_12i = 0.0;
        CHd_13i = 0.0;
        CHd_23i = 0.0;
    } else if (name.compare("CHud_11r") == 0){
        CHud_11r = value;
        if (FlagPartialQFU){
            CHud_22r = value;
        }
    } else if (name.compare("CHud_12r") == 0)
        CHud_12r = value;
    else if (name.compare("CHud_13r") == 0)
        CHud_13r = value;
    else if (name.compare("CHud_22r") == 0){
        if (!FlagPartialQFU){
            CHud_22r = value;
        }
    } else if (name.compare("CHud_23r") == 0)
        CHud_23r = value;
    else if (name.compare("CHud_33r") == 0)
        CHud_33r = value;
    else if (name.compare("CHud_r") == 0) {
        CHud_11r = value;
        CHud_12r = 0.0;
        CHud_13r = 0.0;
        CHud_22r = value;
        CHud_23r = 0.0;
        CHud_33r = value;
    } else if (name.compare("CHud_11i") == 0){
        CHud_11i = value;
        if (FlagPartialQFU){
            CHud_22i = value;
        }
    } else if (name.compare("CHud_12i") == 0)
        CHud_12i = value;
    else if (name.compare("CHud_13i") == 0)
        CHud_13i = value;
    else if (name.compare("CHud_22i") == 0){
        if (!FlagPartialQFU){
            CHud_22i = value;
        }
    } else if (name.compare("CHud_23i") == 0)
        CHud_23i = value;
    else if (name.compare("CHud_33i") == 0)
        CHud_33i = value;
    else if (name.compare("CHud_i") == 0) {
        CHud_11i = value;
        CHud_12i = 0.0;
        CHud_13i = 0.0;
        CHud_22i = value;
        CHud_23i = 0.0;
        CHud_33i = value;
    } else if (name.compare("CeH_11r") == 0)
        CeH_11r = value;
    else if (name.compare("CeH_12r") == 0)
        CeH_12r = value;
    else if (name.compare("CeH_13r") == 0)
        CeH_13r = value;
    else if (name.compare("CeH_22r") == 0)
        CeH_22r = value;
    else if (name.compare("CeH_23r") == 0)
        CeH_23r = value;
    else if (name.compare("CeH_33r") == 0)
        CeH_33r = value;
    else if (name.compare("CeH_11i") == 0)
        CeH_11i = value;
    else if (name.compare("CeH_12i") == 0)
        CeH_12i = value;
    else if (name.compare("CeH_13i") == 0)
        CeH_13i = value;
    else if (name.compare("CeH_22i") == 0)
        CeH_22i = value;
    else if (name.compare("CeH_23i") == 0)
        CeH_23i = value;
    else if (name.compare("CeH_33i") == 0)
        CeH_33i = value;
    else if (name.compare("CuH_11r") == 0)
        CuH_11r = value;
    else if (name.compare("CuH_12r") == 0)
        CuH_12r = value;
    else if (name.compare("CuH_13r") == 0)
        CuH_13r = value;
    else if (name.compare("CuH_22r") == 0)
        CuH_22r = value;
    else if (name.compare("CuH_23r") == 0)
        CuH_23r = value;
    else if (name.compare("CuH_33r") == 0)
        CuH_33r = value;
    else if (name.compare("CuH_11i") == 0)
        CuH_11i = value;
    else if (name.compare("CuH_12i") == 0)
        CuH_12i = value;
    else if (name.compare("CuH_13i") == 0)
        CuH_13i = value;
    else if (name.compare("CuH_22i") == 0)
        CuH_22i = value;
    else if (name.compare("CuH_23i") == 0)
        CuH_23i = value;
    else if (name.compare("CuH_33i") == 0)
        CuH_33i = value;
    else if (name.compare("CdH_11r") == 0)
        CdH_11r = value;
    else if (name.compare("CdH_12r") == 0)
        CdH_12r = value;
    else if (name.compare("CdH_13r") == 0)
        CdH_13r = value;
    else if (name.compare("CdH_22r") == 0)
        CdH_22r = value;
    else if (name.compare("CdH_23r") == 0)
        CdH_23r = value;
    else if (name.compare("CdH_33r") == 0)
        CdH_33r = value;
    else if (name.compare("CdH_11i") == 0)
        CdH_11i = value;
    else if (name.compare("CdH_12i") == 0)
        CdH_12i = value;
    else if (name.compare("CdH_13i") == 0)
        CdH_13i = value;
    else if (name.compare("CdH_22i") == 0)
        CdH_22i = value;
    else if (name.compare("CdH_23i") == 0)
        CdH_23i = value;
    else if (name.compare("CdH_33i") == 0)
        CdH_33i = value;
    else if (name.compare("CuG_11r") == 0)
        CuG_11r = value;
    else if (name.compare("CuG_12r") == 0)
        CuG_12r = value;
    else if (name.compare("CuG_13r") == 0)
        CuG_13r = value;
    else if (name.compare("CuG_22r") == 0)
        CuG_22r = value;
    else if (name.compare("CuG_23r") == 0)
        CuG_23r = value;
    else if (name.compare("CuG_33r") == 0)
        CuG_33r = value;
    else if (name.compare("CuG_r") == 0) {
        CuG_11r = value;
        CuG_12r = 0.0;
        CuG_13r = 0.0;
        CuG_22r = value;
        CuG_23r = 0.0;
        CuG_33r = value;
    } else if (name.compare("CuG_11i") == 0)
        CuG_11i = value;
    else if (name.compare("CuG_12i") == 0)
        CuG_12i = value;
    else if (name.compare("CuG_13i") == 0)
        CuG_13i = value;
    else if (name.compare("CuG_22i") == 0)
        CuG_22i = value;
    else if (name.compare("CuG_23i") == 0)
        CuG_23i = value;
    else if (name.compare("CuG_33i") == 0)
        CuG_33i = value;
    else if (name.compare("CuG_i") == 0) {
        CuG_11i = value;
        CuG_12i = 0.0;
        CuG_13i = 0.0;
        CuG_22i = value;
        CuG_23i = 0.0;
        CuG_33i = value;
    } else if (name.compare("CuW_11r") == 0)
        CuW_11r = value;
    else if (name.compare("CuW_12r") == 0)
        CuW_12r = value;
    else if (name.compare("CuW_13r") == 0)
        CuW_13r = value;
    else if (name.compare("CuW_22r") == 0)
        CuW_22r = value;
    else if (name.compare("CuW_23r") == 0)
        CuW_23r = value;
    else if (name.compare("CuW_33r") == 0)
        CuW_33r = value;
    else if (name.compare("CuW_r") == 0) {
        CuW_11r = value;
        CuW_12r = 0.0;
        CuW_13r = 0.0;
        CuW_22r = value;
        CuW_23r = 0.0;
        CuW_33r = value;
    } else if (name.compare("CuW_11i") == 0)
        CuW_11i = value;
    else if (name.compare("CuW_12i") == 0)
        CuW_12i = value;
    else if (name.compare("CuW_13i") == 0)
        CuW_13i = value;
    else if (name.compare("CuW_22i") == 0)
        CuW_22i = value;
    else if (name.compare("CuW_23i") == 0)
        CuW_23i = value;
    else if (name.compare("CuW_33i") == 0)
        CuW_33i = value;
    else if (name.compare("CuW_i") == 0) {
        CuW_11i = value;
        CuW_12i = 0.0;
        CuW_13i = 0.0;
        CuW_22i = value;
        CuW_23i = 0.0;
        CuW_33i = value;
    } else if (name.compare("CuB_11r") == 0)
        CuB_11r = value;
    else if (name.compare("CuB_12r") == 0)
        CuB_12r = value;
    else if (name.compare("CuB_13r") == 0)
        CuB_13r = value;
    else if (name.compare("CuB_22r") == 0)
        CuB_22r = value;
    else if (name.compare("CuB_23r") == 0)
        CuB_23r = value;
    else if (name.compare("CuB_33r") == 0)
        CuB_33r = value;
    else if (name.compare("CuB_r") == 0) {
        CuB_11r = value;
        CuB_12r = 0.0;
        CuB_13r = 0.0;
        CuB_22r = value;
        CuB_23r = 0.0;
        CuB_33r = value;
    } else if (name.compare("CuB_11i") == 0)
        CuB_11i = value;
    else if (name.compare("CuB_12i") == 0)
        CuB_12i = value;
    else if (name.compare("CuB_13i") == 0)
        CuB_13i = value;
    else if (name.compare("CuB_22i") == 0)
        CuB_22i = value;
    else if (name.compare("CuB_23i") == 0)
        CuB_23i = value;
    else if (name.compare("CuB_33i") == 0)
        CuB_33i = value;
    else if (name.compare("CuB_i") == 0) {
        CuB_11i = value;
        CuB_12i = 0.0;
        CuB_13i = 0.0;
        CuB_22i = value;
        CuB_23i = 0.0;
        CuB_33i = value;
//  Several redundancies for the 4-fermionn operators below
    } else if (name.compare("CLL_1111") == 0) {
        CLL_1111 = value;
    } else if (name.compare("CLL_1122") == 0) {
        CLL_1122 = value;
        CLL_2211 = value;
    } else if (name.compare("CLL_1133") == 0) {
        CLL_1133 = value;
        CLL_3311 = value;
    } else if (name.compare("CLL_1221") == 0) {
        CLL_1221 = value;
        CLL_2112 = value;
    } else if (name.compare("CLL_1331") == 0) {
        CLL_1331 = value;
        CLL_3113 = value;
    } else if (name.compare("CLL") == 0) {
        CLL_1111 = value;
        CLL_1221 = value;
        CLL_2112 = value;
        CLL_2211 = value;
        CLL_1122 = value;
        CLL_3311 = value;
        CLL_1133 = value;
        CLL_1331 = value;
        CLL_3113 = value;
    } else if (name.compare("CLQ1_1111") == 0) {
        CLQ1_1111 = value;
    } else if (name.compare("CLQ1_1122") == 0) {
        CLQ1_1122 = value;
    } else if (name.compare("CLQ1_2211") == 0) {
        CLQ1_2211 = value;
    } else if (name.compare("CLQ1_2112") == 0) {
        CLQ1_2112 = value;
    } else if (name.compare("CLQ1_1221") == 0) {
        CLQ1_1221 = value;
    } else if (name.compare("CLQ1_1133") == 0) {
        CLQ1_1133 = value;
    } else if (name.compare("CLQ1_3311") == 0) {
        CLQ1_3311 = value;
    } else if (name.compare("CLQ1_3113") == 0) {
        CLQ1_3113 = value;
    } else if (name.compare("CLQ1_1331") == 0) {
        CLQ1_1331 = value;
    } else if (name.compare("CLQ1_1123") == 0) {
        CLQ1_1123 = value;
    } else if (name.compare("CLQ1_2223") == 0) {
        CLQ1_2223 = value;
    } else if (name.compare("CLQ1_3323") == 0) {
        CLQ1_3323 = value;
    } else if (name.compare("CLQ1_1132") == 0) {
        CLQ1_1132 = value;
    } else if (name.compare("CLQ1_2232") == 0) {
        CLQ1_2232 = value;
    } else if (name.compare("CLQ1_3332") == 0) {
        CLQ1_3332 = value;
    } else if (name.compare("CLQ1") == 0) {
        CLQ1_1111 = value;
        CLQ1_1122 = value;
        CLQ1_2211 = value;
        CLQ1_1221 = value;
        CLQ1_2112 = value;
        CLQ1_1133 = value;
        CLQ1_3311 = value;
        CLQ1_1331 = value;
        CLQ1_3113 = value;
    } else if (name.compare("CLQ3_1111") == 0) {
        CLQ3_1111 = value;
    } else if (name.compare("CLQ3_1122") == 0) {
        CLQ3_1122 = value;
    } else if (name.compare("CLQ3_2211") == 0) {
        CLQ3_2211 = value;
    } else if (name.compare("CLQ3_2112") == 0) {
        CLQ3_2112 = value;
    } else if (name.compare("CLQ3_1221") == 0) {
        CLQ3_1221 = value;
    } else if (name.compare("CLQ3_1133") == 0) {
        CLQ3_1133 = value;
    } else if (name.compare("CLQ3_3311") == 0) {
        CLQ3_3311 = value;
    } else if (name.compare("CLQ3_3113") == 0) {
        CLQ3_3113 = value;
    } else if (name.compare("CLQ3_1331") == 0) {
        CLQ3_1331 = value;
    } else if (name.compare("CLQ3_1123") == 0) {
        CLQ3_1123 = value;
    } else if (name.compare("CLQ3_2223") == 0) {
        CLQ3_2223 = value;
    } else if (name.compare("CLQ3_3323") == 0) {
        CLQ3_3323 = value;
    } else if (name.compare("CLQ3_1132") == 0) {
        CLQ3_1132 = value;
    } else if (name.compare("CLQ3_2232") == 0) {
        CLQ3_2232 = value;
    } else if (name.compare("CLQ3_3332") == 0) {
        CLQ3_3332 = value;
    } else if (name.compare("CLQ3") == 0) {
        CLQ3_1111 = value;
        CLQ3_1122 = value;
        CLQ3_2211 = value;
        CLQ3_1221 = value;
        CLQ3_2112 = value;
        CLQ3_1133 = value;
        CLQ3_3311 = value;
        CLQ3_1331 = value;
        CLQ3_3113 = value;
    } else if (name.compare("Cee") == 0) {
        Cee_1111 = value;
        Cee_1122 = value;
        Cee_2211 = value;
        Cee_1133 = value;
        Cee_3311 = value;
    } else if (name.compare("Cee_1111") == 0) {
        Cee_1111 = value;
    } else if (name.compare("Cee_1122") == 0) {
        Cee_1122 = value;
        Cee_2211 = value;
    } else if (name.compare("Cee_1133") == 0) {
        Cee_1133 = value;
        Cee_3311 = value;
    } else if (name.compare("Ceu") == 0) {
        Ceu_1111 = value;
        Ceu_1122 = value;
        Ceu_2211 = value;
        Ceu_1133 = value;
        Ceu_2233 = value;
        Ceu_3311 = value;
    } else if (name.compare("Ceu_1111") == 0) {
        Ceu_1111 = value;
    } else if (name.compare("Ceu_1122") == 0) {
        Ceu_1122 = value;
    } else if (name.compare("Ceu_2211") == 0) {
        Ceu_2211 = value;
    } else if (name.compare("Ceu_1133") == 0) {
        Ceu_1133 = value;
    } else if (name.compare("Ceu_2233") == 0) {
        Ceu_2233 = value;
    } else if (name.compare("Ceu_3311") == 0) {
        Ceu_3311 = value;
    } else if (name.compare("Ced") == 0) {
        Ced_1111 = value;
        Ced_1122 = value;
        Ced_2211 = value;
        Ced_1133 = value;
        Ced_3311 = value;
    } else if (name.compare("Ced_1111") == 0) {
        Ced_1111 = value;
    } else if (name.compare("Ced_1122") == 0) {
        Ced_1122 = value;
    } else if (name.compare("Ced_2211") == 0) {
        Ced_2211 = value;
    } else if (name.compare("Ced_1133") == 0) {
        Ced_1133 = value;
    } else if (name.compare("Ced_3311") == 0) {
        Ced_3311 = value;
    } else if (name.compare("Ced_1123") == 0) {
        Ced_1123 = value;
    } else if (name.compare("Ced_2223") == 0) {
        Ced_2223 = value;
    } else if (name.compare("Ced_3323") == 0) {
        Ced_3323 = value;
    } else if (name.compare("Ced_1132") == 0) {
        Ced_1132 = value;
    } else if (name.compare("Ced_2232") == 0) {
        Ced_2232 = value;
    } else if (name.compare("Ced_3332") == 0) {
        Ced_3332 = value;
    } else if (name.compare("CLe") == 0) {
        CLe_1111 = value;
        CLe_1122 = value;
        CLe_2211 = value;
        CLe_1133 = value;
        CLe_3311 = value;
    } else if (name.compare("CLe_1111") == 0) {
        CLe_1111 = value;
    } else if (name.compare("CLe_1122") == 0) {
        CLe_1122 = value;
    } else if (name.compare("CLe_2211") == 0) {
        CLe_2211 = value;
    } else if (name.compare("CLe_1133") == 0) {
        CLe_1133 = value;
    } else if (name.compare("CLe_3311") == 0) {
        CLe_3311 = value;
    } else if (name.compare("CLu") == 0) {
        CLu_1111 = value;
        CLu_1122 = value;
        CLu_2211 = value;
        CLu_1133 = value;
        CLu_2233 = value;
        CLu_3311 = value;
    } else if (name.compare("CLu_1111") == 0) {
        CLu_1111 = value;
    } else if (name.compare("CLu_1122") == 0) {
        CLu_1122 = value;
    } else if (name.compare("CLu_2211") == 0) {
        CLu_2211 = value;
    } else if (name.compare("CLu_1133") == 0) {
        CLu_1133 = value;
    } else if (name.compare("CLu_2233") == 0) {
        CLu_2233 = value;
    } else if (name.compare("CLu_3311") == 0) {
        CLu_3311 = value;
    } else if (name.compare("CLd") == 0) {
        CLd_1111 = value;
        CLd_1122 = value;
        CLd_2211 = value;
        CLd_1133 = value;
        CLd_3311 = value;
    } else if (name.compare("CLd_1111") == 0) {
        CLd_1111 = value;
    } else if (name.compare("CLd_1122") == 0) {
        CLd_1122 = value;
    } else if (name.compare("CLd_2211") == 0) {
        CLd_2211 = value;
    } else if (name.compare("CLd_1133") == 0) {
        CLd_1133 = value;
    } else if (name.compare("CLd_3311") == 0) {
        CLd_3311 = value;
    } else if (name.compare("CLd_1123") == 0) {
        CLd_1123 = value;
    } else if (name.compare("CLd_2223") == 0) {
        CLd_2223 = value;
    } else if (name.compare("CLd_3323") == 0) {
        CLd_3323 = value;
    } else if (name.compare("CLd_1132") == 0) {
        CLd_1132 = value;
    } else if (name.compare("CLd_2232") == 0) {
        CLd_2232 = value;
    } else if (name.compare("CLd_3332") == 0) {
        CLd_3332 = value;
    } else if (name.compare("CQe") == 0) {
        CQe_1111 = value;
        CQe_1122 = value;
        CQe_2211 = value;
        CQe_1133 = value;
        CQe_3311 = value;
    } else if (name.compare("CQe_1111") == 0) {
        CQe_1111 = value;
    } else if (name.compare("CQe_1122") == 0) {
        CQe_1122 = value;
    } else if (name.compare("CQe_2211") == 0) {
        CQe_2211 = value;
    } else if (name.compare("CQe_1133") == 0) {
        CQe_1133 = value;
    } else if (name.compare("CQe_3311") == 0) {
        CQe_3311 = value;
    } else if (name.compare("CQe_2311") == 0) {
        CQe_2311 = value;
    } else if (name.compare("CQe_2322") == 0) {
        CQe_2322 = value;
    } else if (name.compare("CQe_2333") == 0) {
        CQe_2333 = value;
    } else if (name.compare("CQe_3211") == 0) {
        CQe_3211 = value;
    } else if (name.compare("CQe_3222") == 0) {
        CQe_3222 = value;
    } else if (name.compare("CLedQ_11") == 0) {
        CLedQ_11 = value;
    } else if (name.compare("CLedQ_22") == 0) {
        CLedQ_22 = value; 
    } else if (name.compare("CpLedQ_11") == 0) {
        CpLedQ_11 = value;
    } else if (name.compare("CpLedQ_22") == 0) {
        CpLedQ_22 = value; 
    } else if (name.compare("CQe_3233") == 0) {
        CQe_3233 = value;
    } else if (name.compare("Lambda_NP") == 0) {
        Lambda_NP = value;        
    } else if (name.compare("BrHinv") == 0) {
        BrHinv = value;       
    } else if (name.compare("BrHexo") == 0) {
        BrHexo = value;   
    } else if (name.compare("eggFint") == 0) {
        eggFint = value;
    } else if (name.compare("eggFpar") == 0) {
        eggFpar = value;
    } else if (name.compare("ettHint") == 0) {
        ettHint = value;
    } else if (name.compare("ettHpar") == 0) {
        ettHpar = value;
    } else if (name.compare("eVBFint") == 0) {
        eVBFint = value;
    } else if (name.compare("eVBFpar") == 0) {
        eVBFpar = value;
    } else if (name.compare("eWHint") == 0) {
        eWHint = value;
    } else if (name.compare("eWHpar") == 0) {
        eWHpar = value;
    } else if (name.compare("eZHint") == 0) {
        eZHint = value;
    } else if (name.compare("eZHpar") == 0) {
        eZHpar = value;
    } else if (name.compare("eeeWBFint") == 0) {
        eeeWBFint = value;
    } else if (name.compare("eeeWBFpar") == 0) {
        eeeWBFpar = value;
    } else if (name.compare("eeeZHint") == 0) {
        eeeZHint = value;
    } else if (name.compare("eeeZHpar") == 0) {
        eeeZHpar = value;
    } else if (name.compare("eeettHint") == 0) {
        eeettHint = value;
    } else if (name.compare("eeettHpar") == 0) {
        eeettHpar = value;
    } else if (name.compare("eHggint") == 0) {
        eHggint = value;
    } else if (name.compare("eHggpar") == 0) {
        eHggpar = value;
    } else if (name.compare("eHWWint") == 0) {
        eHWWint = value;
    } else if (name.compare("eHWWpar") == 0) {
        eHWWpar = value;
    } else if (name.compare("eHZZint") == 0) {
        eHZZint = value;
    } else if (name.compare("eHZZpar") == 0) {
        eHZZpar = value;
    } else if (name.compare("eHZgaint") == 0) {
        eHZgaint = value;
    } else if (name.compare("eHZgapar") == 0) {
        eHZgapar = value;
    } else if (name.compare("eHgagaint") == 0) {
        eHgagaint = value;
    } else if (name.compare("eHgagapar") == 0) {
        eHgagapar = value;
    } else if (name.compare("eHmumuint") == 0) {
        eHmumuint = value;
    } else if (name.compare("eHmumupar") == 0) {
        eHmumupar = value;
    } else if (name.compare("eHtautauint") == 0) {
        eHtautauint = value;
    } else if (name.compare("eHtautaupar") == 0) {
        eHtautaupar = value;
    } else if (name.compare("eHccint") == 0) {
        eHccint = value;
    } else if (name.compare("eHccpar") == 0) {
        eHccpar = value;
    } else if (name.compare("eHbbint") == 0) {
        eHbbint = value;
    } else if (name.compare("eHbbpar") == 0) {
        eHbbpar = value;
    } else if (name.compare("eVBF_2_Hbox") == 0) {
         eVBF_2_Hbox = value;
    } else if (name.compare("eVBF_2_HQ1_11") == 0) {
         eVBF_2_HQ1_11 = value;
    } else if (name.compare("eVBF_2_Hu_11") == 0) {
         eVBF_2_Hu_11 = value;
    } else if (name.compare("eVBF_2_Hd_11") == 0) {
         eVBF_2_Hd_11 = value;
    } else if (name.compare("eVBF_2_HQ3_11") == 0) {
         eVBF_2_HQ3_11 = value;
    } else if (name.compare("eVBF_2_HD") == 0) {
         eVBF_2_HD = value;
    } else if (name.compare("eVBF_2_HB") == 0) {
         eVBF_2_HB = value;
    } else if (name.compare("eVBF_2_HW") == 0) {
         eVBF_2_HW = value;
    } else if (name.compare("eVBF_2_HWB") == 0) {
         eVBF_2_HWB = value;
    } else if (name.compare("eVBF_2_HG") == 0) {
         eVBF_2_HG = value;
    } else if (name.compare("eVBF_2_DHB") == 0) {
         eVBF_2_DHB = value;
    } else if (name.compare("eVBF_2_DHW") == 0) {
         eVBF_2_DHW = value;
    } else if (name.compare("eVBF_2_DeltaGF") == 0) {
         eVBF_2_DeltaGF = value;
    } else if (name.compare("eVBF_78_Hbox") == 0) {
         eVBF_78_Hbox = value;
    } else if (name.compare("eVBF_78_HQ1_11") == 0) {
         eVBF_78_HQ1_11 = value;
    } else if (name.compare("eVBF_78_Hu_11") == 0) {
         eVBF_78_Hu_11 = value;
    } else if (name.compare("eVBF_78_Hd_11") == 0) {
         eVBF_78_Hd_11 = value;
    } else if (name.compare("eVBF_78_HQ3_11") == 0) {
         eVBF_78_HQ3_11 = value;
    } else if (name.compare("eVBF_78_HD") == 0) {
         eVBF_78_HD = value;
    } else if (name.compare("eVBF_78_HB") == 0) {
         eVBF_78_HB = value;
    } else if (name.compare("eVBF_78_HW") == 0) {
         eVBF_78_HW = value;
    } else if (name.compare("eVBF_78_HWB") == 0) {
         eVBF_78_HWB = value;
    } else if (name.compare("eVBF_78_HG") == 0) {
         eVBF_78_HG = value;
    } else if (name.compare("eVBF_78_DHB") == 0) {
         eVBF_78_DHB = value;
    } else if (name.compare("eVBF_78_DHW") == 0) {
         eVBF_78_DHW = value;
    } else if (name.compare("eVBF_78_DeltaGF") == 0) {
         eVBF_78_DeltaGF = value;
    } else if (name.compare("eVBF_1314_Hbox") == 0) {
         eVBF_1314_Hbox = value;
    } else if (name.compare("eVBF_1314_HQ1_11") == 0) {
         eVBF_1314_HQ1_11 = value;
    } else if (name.compare("eVBF_1314_Hu_11") == 0) {
         eVBF_1314_Hu_11 = value;
    } else if (name.compare("eVBF_1314_Hd_11") == 0) {
         eVBF_1314_Hd_11 = value;
    } else if (name.compare("eVBF_1314_HQ3_11") == 0) {
         eVBF_1314_HQ3_11 = value;
    } else if (name.compare("eVBF_1314_HD") == 0) {
         eVBF_1314_HD = value;
    } else if (name.compare("eVBF_1314_HB") == 0) {
         eVBF_1314_HB = value;
    } else if (name.compare("eVBF_1314_HW") == 0) {
         eVBF_1314_HW = value;
    } else if (name.compare("eVBF_1314_HWB") == 0) {
         eVBF_1314_HWB = value;
    } else if (name.compare("eVBF_1314_HG") == 0) {
         eVBF_1314_HG = value;
    } else if (name.compare("eVBF_1314_DHB") == 0) {
         eVBF_1314_DHB = value;
    } else if (name.compare("eVBF_1314_DHW") == 0) {
         eVBF_1314_DHW = value;
    } else if (name.compare("eVBF_1314_DeltaGF") == 0) {
         eVBF_1314_DeltaGF = value;
    } else if (name.compare("eWH_2_Hbox") == 0) {
         eWH_2_Hbox = value;
    } else if (name.compare("eWH_2_HQ3_11") == 0) {
         eWH_2_HQ3_11 = value;
    } else if (name.compare("eWH_2_HD") == 0) {
         eWH_2_HD = value;
    } else if (name.compare("eWH_2_HW") == 0) {
         eWH_2_HW = value;
    } else if (name.compare("eWH_2_HWB") == 0) {
         eWH_2_HWB = value;
    } else if (name.compare("eWH_2_DHW") == 0) {
         eWH_2_DHW = value;
    } else if (name.compare("eWH_2_DeltaGF") == 0) {
         eWH_2_DeltaGF = value;
    } else if (name.compare("eWH_78_Hbox") == 0) {
         eWH_78_Hbox = value;
    } else if (name.compare("eWH_78_HQ3_11") == 0) {
         eWH_78_HQ3_11 = value;
    } else if (name.compare("eWH_78_HD") == 0) {
         eWH_78_HD = value;
    } else if (name.compare("eWH_78_HW") == 0) {
         eWH_78_HW = value;
    } else if (name.compare("eWH_78_HWB") == 0) {
         eWH_78_HWB = value;
    } else if (name.compare("eWH_78_DHW") == 0) {
         eWH_78_DHW = value;
    } else if (name.compare("eWH_78_DeltaGF") == 0) {
         eWH_78_DeltaGF = value;
    } else if (name.compare("eWH_1314_Hbox") == 0) {
         eWH_1314_Hbox = value;
    } else if (name.compare("eWH_1314_HQ3_11") == 0) {
         eWH_1314_HQ3_11 = value;
    } else if (name.compare("eWH_1314_HD") == 0) {
         eWH_1314_HD = value;
    } else if (name.compare("eWH_1314_HW") == 0) {
         eWH_1314_HW = value;
    } else if (name.compare("eWH_1314_HWB") == 0) {
         eWH_1314_HWB = value;
    } else if (name.compare("eWH_1314_DHW") == 0) {
         eWH_1314_DHW = value;
    } else if (name.compare("eWH_1314_DeltaGF") == 0) {
         eWH_1314_DeltaGF = value;
    } else if (name.compare("eZH_2_Hbox") == 0) {
         eZH_2_Hbox = value;
    } else if (name.compare("eZH_2_HQ1_11") == 0) {
         eZH_2_HQ1_11 = value;
    } else if (name.compare("eZH_2_Hu_11") == 0) {
         eZH_2_Hu_11 = value;
    } else if (name.compare("eZH_2_Hd_11") == 0) {
         eZH_2_Hd_11 = value;
    } else if (name.compare("eZH_2_HQ3_11") == 0) {
         eZH_2_HQ3_11 = value;
    } else if (name.compare("eZH_2_HD") == 0) {
         eZH_2_HD = value;
    } else if (name.compare("eZH_2_HB") == 0) {
         eZH_2_HB = value;
    } else if (name.compare("eZH_2_HW") == 0) {
         eZH_2_HW = value;
    } else if (name.compare("eZH_2_HWB") == 0) {
         eZH_2_HWB = value;
    } else if (name.compare("eZH_2_DHB") == 0) {
         eZH_2_DHB = value;
    } else if (name.compare("eZH_2_DHW") == 0) {
         eZH_2_DHW = value;
    } else if (name.compare("eZH_2_DeltaGF") == 0) {
         eZH_2_DeltaGF = value;
    } else if (name.compare("eZH_78_Hbox") == 0) {
         eZH_78_Hbox = value;
    } else if (name.compare("eZH_78_HQ1_11") == 0) {
         eZH_78_HQ1_11 = value;
    } else if (name.compare("eZH_78_Hu_11") == 0) {
         eZH_78_Hu_11 = value;
    } else if (name.compare("eZH_78_Hd_11") == 0) {
         eZH_78_Hd_11 = value;
    } else if (name.compare("eZH_78_HQ3_11") == 0) {
         eZH_78_HQ3_11 = value;
    } else if (name.compare("eZH_78_HD") == 0) {
         eZH_78_HD = value;
    } else if (name.compare("eZH_78_HB") == 0) {
         eZH_78_HB = value;
    } else if (name.compare("eZH_78_HW") == 0) {
         eZH_78_HW = value;
    } else if (name.compare("eZH_78_HWB") == 0) {
         eZH_78_HWB = value;
    } else if (name.compare("eZH_78_DHB") == 0) {
         eZH_78_DHB = value;
    } else if (name.compare("eZH_78_DHW") == 0) {
         eZH_78_DHW = value;
    } else if (name.compare("eZH_78_DeltaGF") == 0) {
         eZH_78_DeltaGF = value;
    } else if (name.compare("eZH_1314_Hbox") == 0) {
         eZH_1314_Hbox = value;
    } else if (name.compare("eZH_1314_HQ1_11") == 0) {
         eZH_1314_HQ1_11 = value;
    } else if (name.compare("eZH_1314_Hu_11") == 0) {
         eZH_1314_Hu_11 = value;
    } else if (name.compare("eZH_1314_Hd_11") == 0) {
         eZH_1314_Hd_11 = value;
    } else if (name.compare("eZH_1314_HQ3_11") == 0) {
         eZH_1314_HQ3_11 = value;
    } else if (name.compare("eZH_1314_HD") == 0) {
         eZH_1314_HD = value;
    } else if (name.compare("eZH_1314_HB") == 0) {
         eZH_1314_HB = value;
    } else if (name.compare("eZH_1314_HW") == 0) {
         eZH_1314_HW = value;
    } else if (name.compare("eZH_1314_HWB") == 0) {
         eZH_1314_HWB = value;
    } else if (name.compare("eZH_1314_DHB") == 0) {
         eZH_1314_DHB = value;
    } else if (name.compare("eZH_1314_DHW") == 0) {
         eZH_1314_DHW = value;
    } else if (name.compare("eZH_1314_DeltaGF") == 0) {
         eZH_1314_DeltaGF = value;
    } else if (name.compare("ettH_2_HG") == 0) {
         ettH_2_HG = value;
    } else if (name.compare("ettH_2_G") == 0) {
         ettH_2_G = value;
    } else if (name.compare("ettH_2_uG_33r") == 0) {
         ettH_2_uG_33r = value;
    } else if (name.compare("ettH_2_DeltagHt") == 0) {
         ettH_2_DeltagHt = value;
    } else if (name.compare("ettH_78_HG") == 0) {
         ettH_78_HG = value;
    } else if (name.compare("ettH_78_G") == 0) {
         ettH_78_G = value;
    } else if (name.compare("ettH_78_uG_33r") == 0) {
         ettH_78_uG_33r = value;
    } else if (name.compare("ettH_78_DeltagHt") == 0) {
         ettH_78_DeltagHt = value;
    } else if (name.compare("ettH_1314_HG") == 0) {
         ettH_1314_HG = value;
    } else if (name.compare("ettH_1314_G") == 0) {
         ettH_1314_G = value;
    } else if (name.compare("ettH_1314_uG_33r") == 0) {
         ettH_1314_uG_33r = value;
    } else if (name.compare("ettH_1314_DeltagHt") == 0) {
         ettH_1314_DeltagHt = value;
    } else
        NPbase::setParameter(name, value);
}

bool NPSMEFTd6::CheckParameters(const std::map<std::string, double>& DPars)
{
    if (FlagLeptonUniversal && FlagQuarkUniversal) {
        if (FlagRotateCHWCHB) {
              for (int i = 0; i < NNPSMEFTd6Vars_LFU_QFU; i++) {
                    if (DPars.find(NPSMEFTd6VarsRot_LFU_QFU[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6_LFU_QFU parameter "
                          << NPSMEFTd6VarsRot_LFU_QFU[i] << std::endl;
                          return false;
                    }
              }
        } else {
              for (int i = 0; i < NNPSMEFTd6Vars_LFU_QFU; i++) {
                    if (DPars.find(NPSMEFTd6Vars_LFU_QFU[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6_LFU_QFU parameter "
                          << NPSMEFTd6Vars_LFU_QFU[i] << std::endl;
                          return false;
                    }
              }
        }
          
        //} else if (FlagLeptonUniversal && !FlagQuarkUniversal) {
        //} else if (!FlagLeptonUniversal && FlagQuarkUniversal) {
    } else if (!FlagLeptonUniversal && !FlagQuarkUniversal) {
        if (FlagRotateCHWCHB) {
              for (int i = 0; i < NNPSMEFTd6Vars; i++) {
                    if (DPars.find(NPSMEFTd6VarsRot[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6 parameter"
                          << NPSMEFTd6VarsRot[i] << std::endl;
                          return false;
                    }
              }
        } else {
              for (int i = 0; i < NNPSMEFTd6Vars; i++) {
                    if (DPars.find(NPSMEFTd6Vars[i]) == DPars.end()) {
                          std::cout << "ERROR: Missing mandatory NPSMEFTd6 parameter"
                          << NPSMEFTd6Vars[i] << std::endl;
                          return false;
                    }
              }
        }
          
    } else
        throw std::runtime_error("Error in NPSMEFTd6::CheckParameters()");

    return (NPbase::CheckParameters(DPars));
}

bool NPSMEFTd6::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("QuadraticTerms") == 0) {
        FlagQuadraticTerms = value;
        if(value) setModelLinearized(false);
        res = true;
    } else if (name.compare("RotateCHWCHB") == 0) {
        FlagRotateCHWCHB = value;
        res = true;
    } else if (name.compare("PartialQFU") == 0) {        
        FlagPartialQFU = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}


////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::CHF1_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CHL1_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CHL1_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CHL1_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CHQ1_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CHQ1_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CHQ1_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHF1_diag(): wrong argument");
}

double NPSMEFTd6::CHF3_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CHL3_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CHL3_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CHL3_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CHQ3_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CHQ3_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CHQ3_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHF3_diag(): wrong argument");
}

double NPSMEFTd6::CHf_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return CHe_11;
    else if (f.is("MU"))
        return CHe_22;
    else if (f.is("TAU"))
        return CHe_33;
    else if (f.is("UP"))
        return CHu_11;
    else if (f.is("CHARM"))
        return CHu_22;
    else if (f.is("TOP"))
        return CHu_33;
    else if (f.is("DOWN"))
        return CHd_11;
    else if (f.is("STRANGE"))
        return CHd_22;
    else if (f.is("BOTTOM"))
        return CHd_33;
    else
        throw std::runtime_error("NPSMEFTd6::CHf_diag(): wrong argument");
}

gslpp::complex NPSMEFTd6::CHud_diag(const Particle u) const
{
    if (!u.is("QUARK") || u.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::CHud_diag(): wrong argument");

    if (u.is("UP"))
        return gslpp::complex(CHud_11r, CHud_11i, false);
    else if (u.is("CHARM"))
        return gslpp::complex(CHud_22r, CHud_22i, false);
    else if (u.is("TOP"))
        return gslpp::complex(CHud_22r, CHud_33i, false);
    else
        throw std::runtime_error("NPSMEFTd6::CHud_diag(): wrong argument");
}

gslpp::complex NPSMEFTd6::CfH_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return gslpp::complex(CeH_11r, CeH_11i, false);
    else if (f.is("MU"))
        return gslpp::complex(CeH_22r, CeH_22i, false);
    else if (f.is("TAU"))
        return gslpp::complex(CeH_33r, CeH_33i, false);
    else if (f.is("UP"))
        return gslpp::complex(CuH_11r, CuH_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuH_22r, CuH_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuH_33r, CuH_33i, false);
    else if (f.is("DOWN"))
        return gslpp::complex(CdH_11r, CdH_11i, false);
    else if (f.is("STRANGE"))
        return gslpp::complex(CdH_22r, CdH_22i, false);
    else if (f.is("BOTTOM"))
        return gslpp::complex(CdH_33r, CdH_33i, false);
    else
        throw std::runtime_error("NPSMEFTd6::CfH_diag(): wrong argument");
}
      
gslpp::complex NPSMEFTd6::CfG_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CuG_11r, CuG_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuG_22r, CuG_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuG_33r, CuG_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPSMEFTd6::CfG_diag(): wrong argument");
}
      
gslpp::complex NPSMEFTd6::CfW_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CuW_11r, CuW_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuW_22r, CuW_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuW_33r, CuW_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPSMEFTd6::CfW_diag(): wrong argument");
}
      
gslpp::complex NPSMEFTd6::CfB_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CuB_11r, CuB_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuB_22r, CuB_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuB_33r, CuB_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPSMEFTd6::CfB_diag(): wrong argument");
}


////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::DeltaGF() const
{
    return ((CHL3_11 + CHL3_22 - 0.5 * (CLL_1221 + CLL_2112)) * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueS() const
{
    return (4.0 * sW2_tree * cW_tree * CHWB / alphaMz() * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueT() const
{
    return (-CHD / 2.0 / alphaMz() * v2_over_LambdaNP2);
}

double NPSMEFTd6::obliqueU() const
{
    return 0.0;
}

double NPSMEFTd6::deltaMz2() const
{
    return (2.0*(cW2_tree*CHW+sW2_tree*CHB+cW_tree*sW_tree*CHWB)*v2_over_LambdaNP2);
}

double NPSMEFTd6::Mw() const
{
    return (trueSM.Mw() - Mw_tree() / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPSMEFTd6::GammaW(const Particle fi, const Particle fj) const
{
    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double Gamma_Wij = trueSM.GammaW(fi, fj);
    double GammaW_tree;
    double CHF3ij;
    
    if (fj.getIndex() - fi.getIndex() == 1)
        CHF3ij = CHF3_diag(fi);
    else
        CHF3ij = 0.;
    
    if (fi.is("QUARK")) {        
        GammaW_tree = Nc * G0;
    } else {
        GammaW_tree = G0;        
    }
    
    Gamma_Wij = Gamma_Wij - 3.0 * GammaW_tree / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * (1.0 + cW2_tree) / 3.0 * DeltaGF());

    return (Gamma_Wij
            + 2.0 * GammaW_tree * CHF3ij * v2_over_LambdaNP2);          
}

double NPSMEFTd6::GammaW() const
{
    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double GammaW_tree = (3.0 + 2.0 * Nc) * G0;

    return (trueSM.GammaW()
            - 3.0 * GammaW_tree / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * (1.0 + cW2_tree) / 3.0 * DeltaGF())
            + 2.0 * G0 * (CHL3_11 + CHL3_22 + CHL3_33 + Nc*(CHQ3_11 + CHQ3_22)) * v2_over_LambdaNP2);          
//            + 2.0 * GammaW_tree / 3.0 * (CHL3_11 + CHQ3_11 + CHQ3_22) * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaGV_f(const Particle p) const
{
    return (deltaGL_f(p) + deltaGR_f(p));
}

double NPSMEFTd6::deltaGA_f(const Particle p) const
{
    return (deltaGL_f(p) - deltaGR_f(p));
}

double NPSMEFTd6::deltaGL_f(const Particle p) const
{
    double I3p = p.getIsospin(), Qp = p.getCharge();
    double CHF1 = CHF1_diag(p);
    double CHF3 = CHF3_diag(p);
    double NPindirect;

    NPindirect = -I3p / 4.0 * (CHD * v2_over_LambdaNP2 + 2.0 * DeltaGF())
                - Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = -0.5 * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPSMEFTd6::deltaGR_f(const Particle p) const
{
    double Qp = p.getCharge();
    double CHf = CHf_diag(p);
    double NPindirect;

    NPindirect = -Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = -0.5 * CHf*v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPSMEFTd6::deltaGL_Wff(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGL_Wff(): Not implemented");

    double CHF3 = CHF3_diag(pbar);
    double NPindirect;

    NPindirect = -cW2_tree / 4.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree / cW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());

    double NPdirect = CHF3 * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

gslpp::complex NPSMEFTd6::deltaGR_Wff(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGR_Wff(): Not implemented");

    gslpp::complex CHud = CHud_diag(pbar);
    return (0.5 * CHud * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaG_hgg() const
{
    return (CHG * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG_hggRatio() const
{    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double aSPiv = AlsMz / 16.0 / M_PI / v();
    gslpp::complex gSM, dg;
    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = deltaG_hff(quarks[CHARM]) / (-m_c / v());
    double deltaloc = deltaG_hgg();
    
    gSM = aSPiv * (AH_f(tau_t) + AH_f(tau_b) + AH_f(tau_c));
    
    dg = deltaloc/gSM + (aSPiv/gSM) * (dKappa_t*AH_f(tau_t) + dKappa_b*AH_f(tau_b) + dKappa_c*AH_f(tau_c));

    return dg.real();
}

double NPSMEFTd6::deltaG1_hWW() const
{
    return (( 2.0 * CHW - sqrt( M_PI * ale ) * CDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG2_hWW() const
{
    return ( - sqrt( M_PI * ale ) * ( CDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG3_hWW() const
{
    double NPindirect;

    NPindirect = 2.0 * cW2_tree * Mz * Mz / v()
                * (delta_h - 1.0 / 2.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree * cW_tree * CHWB + cW2_tree * CHD) * v2_over_LambdaNP2 + DeltaGF()));

    return NPindirect;
}

double NPSMEFTd6::deltaG1_hZZ() const
{
    return ( (delta_ZZ - 0.5 * sqrt( M_PI * ale ) * (CDHB / cW_tree + CDHW / sW_tree) * v2_over_LambdaNP2 )/ v());
}

double NPSMEFTd6::deltaG2_hZZ() const
{
    return ( - sqrt( M_PI * ale ) * ( CDHB / cW_tree + CDHW / sW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG3_hZZ() const
{
    double NPindirect = Mz * Mz / v() * (-0.5 * CHD * v2_over_LambdaNP2 + delta_h - 0.5 * DeltaGF());
    double NPdirect = Mz * Mz / v() * CHD * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPSMEFTd6::deltaG1_hZA() const
{
    return ( (delta_AZ + 0.5 * sqrt( M_PI * ale ) * (CDHB / sW_tree - CDHW / cW_tree) * v2_over_LambdaNP2 )/ v());
}

double NPSMEFTd6::deltaG1_hZARatio() const
{    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double m_tau = leptons[TAU].getMass();
    double m_mu = leptons[MU].getMass();
    
    double M_w_2 = (trueSM.Mw())*(trueSM.Mw());
    
    double Qt = quarks[TOP].getCharge();
    double Qb = quarks[BOTTOM].getCharge();
    double Qc = quarks[CHARM].getCharge();
    double Qtau = leptons[TAU].getCharge();
    double Qmu = leptons[MU].getCharge();

    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double tau_tau = 4.0 * m_tau * m_tau / mHl / mHl;
    double tau_mu = 4.0 * m_mu * m_mu / mHl / mHl;
    double tau_W = 4.0 * M_w_2 / mHl / mHl;    
    
    double lambda_t = 4.0 * m_t * m_t / Mz / Mz;
    double lambda_b = 4.0 * m_b * m_b / Mz / Mz;
    double lambda_c = 4.0 * m_c * m_c / Mz / Mz;
    double lambda_tau = 4.0 * m_tau * m_tau / Mz / Mz;
    double lambda_mu = 4.0 * m_mu * m_mu / Mz / Mz;
    double lambda_W = 4.0 * M_w_2 / Mz / Mz; 
    double alpha2 = sqrt(2.0)*GF*M_w_2 / M_PI;
    double aPiv = sqrt(ale*alpha2) / 4.0 / M_PI / v();
    
//  mod. of Higgs couplings
    gslpp::complex gSM, dg;
    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = deltaG_hff(quarks[CHARM]) / (-m_c / v());
    gslpp::complex dKappa_tau = deltaG_hff(leptons[TAU]) / (-m_tau / v());
    gslpp::complex dKappa_mu = deltaG_hff(leptons[MU]) / (-m_mu / v());
    double dKappa_W = (0.5 * v() / M_w_2)*deltaG3_hWW();
    
//  mod of EW vector couplings vf =2 gvf    
    double vSMt = 2.0*(quarks[TOP].getIsospin()) - 4.0 * Qt * sW2_tree;
    double vSMb = 2.0*(quarks[BOTTOM].getIsospin()) - 4.0 * Qb * sW2_tree;
    double vSMc = 2.0*(quarks[CHARM].getIsospin()) - 4.0 * Qc * sW2_tree;
    double vSMtau = 2.0*(leptons[TAU].getIsospin()) - 4.0 * Qtau * sW2_tree;
    double vSMmu = 2.0*(leptons[MU].getIsospin()) - 4.0 * Qmu * sW2_tree;
    
    double dvSMt = 2.0*deltaGV_f(quarks[TOP]);
    double dvSMb = 2.0*deltaGV_f(quarks[BOTTOM]);
    double dvSMc = 2.0*deltaGV_f(quarks[CHARM]);
    double dvSMtau = 2.0*deltaGV_f(leptons[TAU]);
    double dvSMmu = 2.0*deltaGV_f(leptons[MU]);
    
    double deltaloc = deltaG1_hZA();
    
    gSM = -aPiv * ((3.0*vSMt*Qt*AHZga_f(tau_t,lambda_t) + 
            3.0*vSMb*Qb*AHZga_f(tau_b,lambda_b) + 
            3.0*vSMc*Qc*AHZga_f(tau_c,lambda_c) +
            vSMtau*Qtau*AHZga_f(tau_tau,lambda_tau)+
            vSMmu*Qmu*AHZga_f(tau_mu,lambda_mu))/cW_tree + 
            AHZga_W(tau_W,lambda_W));
    
    dg = deltaloc/gSM - (aPiv/gSM) * (
            (3.0*vSMt*dKappa_t*Qt*AHZga_f(tau_t,lambda_t) + 
            3.0*vSMb*dKappa_b*Qb*AHZga_f(tau_b,lambda_b) + 
            3.0*vSMc*dKappa_c*Qc*AHZga_f(tau_c,lambda_c)+ 
            dKappa_tau*vSMtau*Qtau*AHZga_f(tau_tau,lambda_tau)+ 
            dKappa_mu*vSMmu*Qmu*AHZga_f(tau_mu,lambda_mu))/cW_tree + 
            dKappa_W*AHZga_W(tau_W,lambda_W) + 
            (3.0*dvSMt*Qt*AHZga_f(tau_t,lambda_t) + 
            3.0*dvSMb*Qb*AHZga_f(tau_b,lambda_b) + 
            3.0*dvSMc*Qc*AHZga_f(tau_c,lambda_c)+ 
            dvSMtau*Qtau*AHZga_f(tau_tau,lambda_tau)+ 
            dvSMmu*Qmu*AHZga_f(tau_mu,lambda_mu))/cW_tree
            );

    return dg.real();
}

double NPSMEFTd6::deltaG2_hZA() const
{
    return ( sqrt( M_PI * ale ) * ( CDHB / sW_tree - CDHW / cW_tree ) * v2_over_LambdaNP2 / v());
}

double NPSMEFTd6::deltaG_hAA() const
{
    return (delta_AA / v());
}

double NPSMEFTd6::deltaG_hAARatio() const
{    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double m_tau = leptons[TAU].getMass();
    double m_mu = leptons[MU].getMass();
    
    double M_w_2 = (trueSM.Mw())*(trueSM.Mw());
    
    double Qt = quarks[TOP].getCharge();
    double Qb = quarks[BOTTOM].getCharge();
    double Qc = quarks[CHARM].getCharge();
    double Qtau = leptons[TAU].getCharge();
    double Qmu = leptons[MU].getCharge();

    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double tau_tau = 4.0 * m_tau * m_tau / mHl / mHl;
    double tau_mu = 4.0 * m_mu * m_mu / mHl / mHl;
    double tau_W = 4.0 * M_w_2 / mHl / mHl;    

    double aPiv = ale / 8.0 / M_PI / v();
    gslpp::complex gSM, dg;
    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = deltaG_hff(quarks[CHARM]) / (-m_c / v());
    gslpp::complex dKappa_tau = deltaG_hff(leptons[TAU]) / (-m_tau / v());
    gslpp::complex dKappa_mu = deltaG_hff(leptons[MU]) / (-m_mu / v());
    double dKappa_W = (0.5 * v() / M_w_2)*deltaG3_hWW();
    
    double deltaloc = deltaG_hAA();
    
    gSM = aPiv * (3.0*Qt*Qt*AH_f(tau_t) + 
            3.0*Qb*Qb*AH_f(tau_b) + 
            3.0*Qc*Qc*AH_f(tau_c) + 
            Qtau*Qtau*AH_f(tau_tau) + 
            Qmu*Qmu*AH_f(tau_mu) + 
            AH_W(tau_W));
    
    dg = deltaloc/gSM + (aPiv/gSM) * (
            3.0*Qt*Qt*dKappa_t*AH_f(tau_t) + 
            3.0*Qb*Qb*dKappa_b*AH_f(tau_b) + 
            3.0*Qc*Qc*dKappa_c*AH_f(tau_c) +
            dKappa_tau*Qtau*Qtau*AH_f(tau_tau) +
            dKappa_mu*Qmu*Qmu*AH_f(tau_mu) +
            dKappa_W*AH_W(tau_W)
            );

    return dg.real();
}

gslpp::complex NPSMEFTd6::deltaG_hff(const Particle p) const
{
    /* The effects of the RG running are neglected. */
    double mf;
    if (p.is("TOP"))
        //mf = p.getMass(); // m_t(m_t)
        mf = mtpole; // pole mass
    else
        mf = p.getMass();
    gslpp::complex CfH = CfH_diag(p);
    return (-mf / v() * (delta_h - 0.5 * DeltaGF())
            + CfH * v2_over_LambdaNP2 / sqrt(2.0));
}

double NPSMEFTd6::deltaG_hhhRatio() const
{    
    double dg;
    
    dg = -0.5 * DeltaGF() + 3.0 * delta_h - 2.0 * CH * v2_over_LambdaNP2 * v()*v()/mHl/mHl;

    return dg;
}

gslpp::complex NPSMEFTd6::deltaGL_Wffh(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGL_Wffh(): Not implemented");

    double CHF3 = CHF3_diag(pbar);
    return (2.0 * sqrt(2.0) * Mz * cW_tree / v() / v() * CHF3 * v2_over_LambdaNP2);
}

gslpp::complex NPSMEFTd6::deltaGR_Wffh(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPSMEFTd6::deltaGR_Wffh(): Not implemented");

    gslpp::complex CHud = CHud_diag(pbar);
    return (sqrt(2.0) * Mz * cW_tree / v() / v() * CHud * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaGL_Zffh(const Particle p) const
{
    double I3p = p.getIsospin();
    double CHF1 = CHF1_diag(p);
    double CHF3 = CHF3_diag(p);
    return (-2.0 * Mz / v() / v() * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2);
}

double NPSMEFTd6::deltaGR_Zffh(const Particle p) const
{
    double CHf = CHf_diag(p);
    return (-2.0 * Mz / v() / v() * CHf * v2_over_LambdaNP2);
}
      
gslpp::complex NPSMEFTd6::deltaG_hGff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_hZff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_hAff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_Gff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_Zff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPSMEFTd6::deltaG_Aff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
double NPSMEFTd6::deltag3G() const
{
    /* Set to 0. for the moment */

    return 0.;
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPSMEFTd6::f_triangle(const double tau) const
{
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = asin(1.0 / sqrt(tau));
        return (tmp * tmp);
    } else {
        tmp = log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i();
        return (-0.25 * tmp * tmp);
    }
}

gslpp::complex NPSMEFTd6::g_triangle(const double tau) const
{
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = sqrt(tau -1.0) * asin(1.0 / sqrt(tau));
        return tmp;
    } else {
        tmp = sqrt(1.0 - tau) * ( log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i() );
        return 0.5 * tmp;
    }
}

gslpp::complex NPSMEFTd6::I_triangle_1(const double tau, const double lambda) const
{
    gslpp::complex tmp;
    
    tmp = ( tau*lambda * (f_triangle(tau)- f_triangle(lambda)) + 2.0 * tau * (g_triangle(tau)- g_triangle(lambda)) ) / (tau-lambda);
    
    tmp = tau*lambda * ( 1.0 + tmp ) / (2.0*(tau-lambda));
    
    return tmp;
}

gslpp::complex NPSMEFTd6::I_triangle_2(const double tau, const double lambda) const
{
    gslpp::complex tmp;
    
    tmp = - 0.5 * tau*lambda * (f_triangle(tau)- f_triangle(lambda)) / (tau-lambda);
    
    return tmp;
}

gslpp::complex NPSMEFTd6::AH_f(const double tau) const
{
    return (2.0 * tau * (1.0 + (1.0 - tau) * f_triangle(tau)));
}

gslpp::complex NPSMEFTd6::AH_W(const double tau) const
{
    return -( 2.0 + 3.0 * tau + 3.0 * tau * (2.0 - tau) * f_triangle(tau) );
}

gslpp::complex NPSMEFTd6::AHZga_f(const double tau, const double lambda) const
{
    return I_triangle_1(tau,lambda) - I_triangle_2(tau,lambda);
}

gslpp::complex NPSMEFTd6::AHZga_W(const double tau, const double lambda) const
{
    gslpp::complex tmp;
    
    double tan2w = trueSM.sW2() / trueSM.cW2();
    
    tmp = 4.0 * (3.0 - tan2w ) * I_triangle_2(tau,lambda);
    
    tmp = tmp + ((1.0 +2.0 / tau)* tan2w - (5.0 + 2.0/tau)) * I_triangle_1(tau,lambda);

    return sqrt(trueSM.cW2()) * tmp;
}

double NPSMEFTd6::muggH(const double sqrt_s) const
{
    double m_t = mtpole;
    //doulbe m_t = quarks[TOP].getMass();
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();

    /* L_eff_SM = (G_eff_t_SM + G_eff_b_SM)*hGG */
    gslpp::complex G_eff_t_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_t * m_t / mHl / mHl);
    gslpp::complex G_eff_b_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_b * m_b / mHl / mHl);
    gslpp::complex G_eff_c_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_c * m_c / mHl / mHl);
    gslpp::complex G_eff_SM = G_eff_t_SM + G_eff_b_SM + G_eff_c_SM;

    //double sigma_tt_SM = trueSM.computeSigmaggH_tt(sqrt_s);
    //double sigma_bb_SM = trueSM.computeSigmaggH_bb(sqrt_s);
    //double sigma_tb_SM = trueSM.computeSigmaggH_tb(sqrt_s);
    //gslpp::complex tmp = (2.0 * dKappa_t * sigma_tt_SM
    //        + 2.0 * dKappa_b * sigma_bb_SM
    //        + (dKappa_t + dKappa_b) * sigma_tb_SM)
    //        / (sigma_tt_SM + sigma_bb_SM + sigma_tb_SM);
    
    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-m_t / v());
    gslpp::complex dKappa_b = deltaG_hff(quarks[BOTTOM]) / (-m_b / v());
    gslpp::complex dKappa_c = deltaG_hff(quarks[CHARM]) / (-m_c / v());

    gslpp::complex tmpHG = CHG / v() * v2_over_LambdaNP2 / G_eff_SM;
    gslpp::complex tmpt = G_eff_t_SM * dKappa_t / G_eff_SM;
    gslpp::complex tmpb = G_eff_b_SM * dKappa_b / G_eff_SM;
    gslpp::complex tmpc = G_eff_c_SM * dKappa_c / G_eff_SM;
    
    double mu = (1.0 + 2.0 * ( tmpt.real() + tmpb.real() + tmpc.real() + tmpHG.real() ) );
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        gslpp::complex tmp2 = tmpt +tmpb +tmpc + tmpHG;
        
        mu += tmp2.abs2();
        
    }
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eggFint + eggFpar;
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muggHH(const double sqrt_s) const
{
    double mu = 1.0;
    double A1HH = 0.0, A2HH = 0.0, A3HH = 0.0, A4HH = 0.0, A5HH = 0.0;
    double A6HH = 0.0, A7HH = 0.0, A8HH = 0.0, A9HH = 0.0, A10HH = 0.0;
    double A11HH = 0.0, A12HH = 0.0, A13HH = 0.0, A14HH = 0.0, A15HH = 0.0;    
    double ct,c2t,c3,cg,c2g;
    
    if (sqrt_s == 14.0) {

        A1HH = 2.13;
        A2HH = 10.1;
        A3HH = 0.300;
        A4HH = 21.8;
        A5HH = 188;
        A6HH = -8.62;
        A7HH = -1.43;
        A8HH = 2.93;
        A9HH = 21.0;
        A10HH = 59.8;
        A11HH = -9.93;
        A12HH = -23.1;
        A13HH = 4.87;
        A14HH = 10.5;
        A15HH = 96.6;
        
    } else if (sqrt_s == 100.0) {

        A1HH = 1.95;
        A2HH = 11.2;
        A3HH = 0.229;
        A4HH = 16.0;
        A5HH = 386;
        A6HH = -8.32;
        A7HH = -1.18;
        A8HH = 2.55;
        A9HH = 16.9;
        A10HH = 52.4;
        A11HH = -7.49;
        A12HH = -17.3;
        A13HH = 3.55;
        A14HH = 8.46;
        A15HH = 87.5;
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muggHH()");
    
    ct= 1.0 - 0.5 * DeltaGF() + delta_h - v() * CuH_33r * v2_over_LambdaNP2 / sqrt(2.0)/ mtpole;
    c2t = delta_h - 3.0 *v() * CuH_33r * v2_over_LambdaNP2 / 2.0 /sqrt(2.0)/ mtpole;
    c3 = 1.0 + deltaG_hhhRatio();
    cg = M_PI * CHG * v2_over_LambdaNP2 / AlsMz; 
    c2g = cg;

// In the SM the Eq. returns 0.999. Fix that small offset by adding 0.0010    
    mu = 0.0010 + A1HH*ct*ct*ct*ct +
            A2HH*c2t*c2t +
            A3HH*ct*ct*c3*c3 +
            A4HH*cg*cg*c3*c3 +
            A5HH*c2g*c2g +
            A6HH*c2t*ct*ct +
            A7HH*ct*ct*ct*c3 +
            A8HH*c2t*ct*c3 +
            A9HH*c2t*cg*c3 +
            A10HH*c2t*c2g +
            A11HH*ct*ct*cg*c3 +
            A12HH*ct*ct*c2g +
            A13HH*ct*c3*c3*cg +
            A14HH*ct*c3*c2g +
            A15HH*cg*c3*c2g;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muVBF(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {

        mu += 
                +106687. * (1. + eVBF_2_Hbox ) * CHbox / LambdaNP2
                -24330.2 * (1. + eVBF_2_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +762.361 * (1. + eVBF_2_Hu_11 ) * CHu_11 / LambdaNP2
                +2984.1 * (1. + eVBF_2_Hd_11 ) * CHd_11 / LambdaNP2
                -165737. * (1. + eVBF_2_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -148464. * (1. + eVBF_2_HD ) * CHD / LambdaNP2
                -1546.5 * (1. + eVBF_2_HB ) * CHB / LambdaNP2
                -59464.2 * (1. + eVBF_2_HW ) * CHW / LambdaNP2
                -280832. * (1. + eVBF_2_HWB ) * CHWB / LambdaNP2
                +4339795. * (1. + eVBF_2_HG ) * CHG / LambdaNP2
                -631.671 * (1. + eVBF_2_DHB ) * CDHB / LambdaNP2
                -28599.4 * (1. + eVBF_2_DHW ) * CDHW / LambdaNP2
                -3.907 * (1. + eVBF_2_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  
                +3168590721.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +12137442823.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +577651928.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +509910695.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                -11415160096.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -8881398388.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +279459813.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -3521455021.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17100196036.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                -77226530168.  * CHbox / LambdaNP2 * CHG / LambdaNP2
                -5467218.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                -1847702026.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -235643.  * CHbox / LambdaNP2 * DeltaGF()
                +22132433021.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                -1184927031.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +206490961.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -31657590939.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +18605968199.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +33404487040.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -6092354607.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +24543672403.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +825177521237.  * CHQ1_11 / LambdaNP2 * CHG / LambdaNP2
                +2355456328.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                +6108952298.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                -208291.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +17837072533.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -192115008.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +6612938358.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +3073404487.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +20967109562.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +3742104117.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +15950773252.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +692247876280.  * CHu_11 / LambdaNP2 * CHG / LambdaNP2
                +2434306251.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +678980614.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                -24334.7  * CHu_11 / LambdaNP2 * DeltaGF()
                +5570028316.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -451317796.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +1123284687.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                +2665650185.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +1559137443.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                +698758440.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                +88091483337.  * CHd_11 / LambdaNP2 * CHG / LambdaNP2
                -375321281.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +386887388.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -8195.38  * CHd_11 / LambdaNP2 * DeltaGF()
                +99945545633.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +42592027881.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +42535395339.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +69228926160.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +66412546286.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +1289729906339.  * CHQ3_11 / LambdaNP2 * CHG / LambdaNP2
                -5335678501.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +4564626443.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                +443963.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +7164016554.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +514049227.  * CHD / LambdaNP2 * CHB / LambdaNP2
                +5632759747.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +27425397517.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -48549335657.  * CHD / LambdaNP2 * CHG / LambdaNP2
                -83881507.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +2653495970.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +348173.  * CHD / LambdaNP2 * DeltaGF()
                +22766281856.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +13890219996.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -24508821607.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -11875408408.  * CHB / LambdaNP2 * CHG / LambdaNP2
                +183816162.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +195208016.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +2683.89  * CHB / LambdaNP2 * DeltaGF()
                +49601394032.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -9695055543.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +100692659551.  * CHW / LambdaNP2 * CHG / LambdaNP2
                +1407950338.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -2954432585.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +142986.  * CHW / LambdaNP2 * DeltaGF()
                +41415813548.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -82801132651.  * CHWB / LambdaNP2 * CHG / LambdaNP2
                -2377717273.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +6063210630.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +673549.  * CHWB / LambdaNP2 * DeltaGF()
                +30718797647571.  * CHG / LambdaNP2 * CHG / LambdaNP2
                -35170137225.  * CHG / LambdaNP2 * CDHB / LambdaNP2
                -84278762797.  * CHG / LambdaNP2 * CDHW / LambdaNP2
                -441038.  * CHG / LambdaNP2 * DeltaGF()
                +319342191.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -338118057.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -3017.67  * CDHB / LambdaNP2 * DeltaGF()
                +1508124592.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +71253.5  * CDHW / LambdaNP2 * DeltaGF()
                +4.578  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 7.0) {

        mu += 
                +108269. * (1. + eVBF_78_Hbox ) * CHbox / LambdaNP2
                +9205. * (1. + eVBF_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -23061.1 * (1. + eVBF_78_Hu_11 ) * CHu_11 / LambdaNP2
                +7456.81 * (1. + eVBF_78_Hd_11 ) * CHd_11 / LambdaNP2
                -352628. * (1. + eVBF_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -146546. * (1. + eVBF_78_HD ) * CHD / LambdaNP2
                -2007.19 * (1. + eVBF_78_HB ) * CHB / LambdaNP2
                -77179.6 * (1. + eVBF_78_HW ) * CHW / LambdaNP2
                -277143. * (1. + eVBF_78_HWB ) * CHWB / LambdaNP2
                +3725307. * (1. + eVBF_78_HG ) * CHG / LambdaNP2
                -582.732 * (1. + eVBF_78_DHB ) * CDHB / LambdaNP2
                -47999.8 * (1. + eVBF_78_DHW ) * CDHW / LambdaNP2
                -3.966 * (1. + eVBF_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  
                +3275491815.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +2790313704.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                -1413987719.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +592120172.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                -22300603235.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -8894644264.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -18668998.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -4939514606.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -16910334855.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +2120172230.  * CHbox / LambdaNP2 * CHG / LambdaNP2
                -23794878.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                -2985852569.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -243039.  * CHbox / LambdaNP2 * DeltaGF()
                +37213894914.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                -485609764.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +178056914.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -22127618246.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +2125783720.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +6037748066.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -5580195755.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +3323729051.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +230411041687.  * CHQ1_11 / LambdaNP2 * CHG / LambdaNP2
                +2147204507.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -1557404470.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                -82041.4  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +23061068126.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -427336592.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +4640271078.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -71330679.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +5172283553.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +174819515.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +2094596782.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +186290696796.  * CHu_11 / LambdaNP2 * CHG / LambdaNP2
                +4786763357.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +265250844.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +4450.06  * CHu_11 / LambdaNP2 * DeltaGF()
                +13008946011.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1116902457.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +553487218.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                +1064024950.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +971219528.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                +428415725.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                +140227912849.  * CHd_11 / LambdaNP2 * CHG / LambdaNP2
                -1358897990.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +130359243.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -5268.88  * CHd_11 / LambdaNP2 * DeltaGF()
                +147247671771.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +37023859626.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +14989154715.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +78545760627.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +72184272718.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +505627677598.  * CHQ3_11 / LambdaNP2 * CHG / LambdaNP2
                -1188395006.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +32408300690.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                +841738.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +7080298272.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +108884500.  * CHD / LambdaNP2 * CHB / LambdaNP2
                +6890155071.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +27101233449.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +55509404642.  * CHD / LambdaNP2 * CHG / LambdaNP2
                -233578296.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +4099301801.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +347963.  * CHD / LambdaNP2 * DeltaGF()
                +22370423128.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +11801396398.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -21667907669.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +52089740684.  * CHB / LambdaNP2 * CHG / LambdaNP2
                -210700681.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +138992306.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +1726.63  * CHB / LambdaNP2 * DeltaGF()
                +74291819094.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -9475865195.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +132901680210.  * CHW / LambdaNP2 * CHG / LambdaNP2
                +2188211553.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -9101621937.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +193133.  * CHW / LambdaNP2 * DeltaGF()
                +40369279247.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +45961345463.  * CHWB / LambdaNP2 * CHG / LambdaNP2
                -2966805875.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +9328887306.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +674701.  * CHWB / LambdaNP2 * DeltaGF()
                +31292693191751.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +43681947187.  * CHG / LambdaNP2 * CDHB / LambdaNP2
                +48239178996.  * CHG / LambdaNP2 * CDHW / LambdaNP2
                -1103.62  * CHG / LambdaNP2 * DeltaGF()
                +562497977.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -400628055.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -7809.86  * CDHB / LambdaNP2 * DeltaGF()
                +3559195830.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +113979.  * CDHW / LambdaNP2 * DeltaGF()
                +4.68  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 8.0) {

        mu += 
                +108661. * (1. + eVBF_78_Hbox ) * CHbox / LambdaNP2
                +9068.97 * (1. + eVBF_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -23462.9 * (1. + eVBF_78_Hu_11 ) * CHu_11 / LambdaNP2
                +8153.76 * (1. + eVBF_78_Hd_11 ) * CHd_11 / LambdaNP2
                -363308. * (1. + eVBF_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -146851. * (1. + eVBF_78_HD ) * CHD / LambdaNP2
                -1996.84 * (1. + eVBF_78_HB ) * CHB / LambdaNP2
                -75547.1 * (1. + eVBF_78_HW ) * CHW / LambdaNP2
                -277228. * (1. + eVBF_78_HWB ) * CHWB / LambdaNP2
                +3724686. * (1. + eVBF_78_HG ) * CHG / LambdaNP2
                -998.419 * (1. + eVBF_78_DHB ) * CDHB / LambdaNP2
                -50170.6 * (1. + eVBF_78_DHW ) * CDHW / LambdaNP2
                -3.972 * (1. + eVBF_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

        mu +=  
                +3271653216.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +2671436892.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                -1436225976.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +614693402.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                -22699226225.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -8888093851.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -14310675.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -4815042849.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -16945835760.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                -33224228305.  * CHbox / LambdaNP2 * CHG / LambdaNP2
                -27855895.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                -3098094683.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -242241.  * CHbox / LambdaNP2 * DeltaGF()
                +40982610866.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                -518345952.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +145602796.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -23635077793.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +1957733589.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +6082868791.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -5419752059.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +3115067809.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +223598469091.  * CHQ1_11 / LambdaNP2 * CHG / LambdaNP2
                +2359197937.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -1900324486.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                -80327.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +25420584075.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -415175971.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +4947999002.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -84865629.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +5136034612.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +149762875.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +2113320576.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +154294866461.  * CHu_11 / LambdaNP2 * CHG / LambdaNP2
                +5366935685.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +290373575.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +3787.84  * CHu_11 / LambdaNP2 * DeltaGF()
                +14545303270.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1108245278.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +538314336.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                +1062484400.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +891089109.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                +395207588.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                +116081204759.  * CHd_11 / LambdaNP2 * CHG / LambdaNP2
                -1526316665.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +128962476.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -5764.11  * CHd_11 / LambdaNP2 * DeltaGF()
                +160812047591.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +37611282137.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +15110242117.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +79784507863.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +74765787503.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +508960812048.  * CHQ3_11 / LambdaNP2 * CHG / LambdaNP2
                -1091172311.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +36909892670.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                +863066.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +7076295865.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +99009901.  * CHD / LambdaNP2 * CHB / LambdaNP2
                +6679424245.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +27018054747.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +7067975705.  * CHD / LambdaNP2 * CHG / LambdaNP2
                -228371745.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +4254929695.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +347769.  * CHD / LambdaNP2 * DeltaGF()
                +22022630835.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +11426907397.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -21217239371.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -14534487062.  * CHB / LambdaNP2 * CHG / LambdaNP2
                -213395457.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +125634412.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +1193.99  * CHB / LambdaNP2 * DeltaGF()
                +75710957650.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -9697978201.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +102340460937.  * CHW / LambdaNP2 * CHG / LambdaNP2
                +2186970630.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -9435061153.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +189200.  * CHW / LambdaNP2 * DeltaGF()
                +40426824195.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +23638405857.  * CHWB / LambdaNP2 * CHG / LambdaNP2
                -3057226059.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +9564855645.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +675403.  * CHWB / LambdaNP2 * DeltaGF()
                +31218071386971.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +26015908145.  * CHG / LambdaNP2 * CDHB / LambdaNP2
                +37667027207.  * CHG / LambdaNP2 * CDHW / LambdaNP2
                -73327.7  * CHG / LambdaNP2 * DeltaGF()
                +608603045.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -411165654.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -7884.75  * CDHB / LambdaNP2 * DeltaGF()
                +3933771528.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +118356.  * CDHW / LambdaNP2 * DeltaGF()
                +4.676  * DeltaGF() * DeltaGF()
                ;

        }
    } else if (sqrt_s == 13.0) {

        mu += 
                +107703. * (1. + eVBF_1314_Hbox ) * CHbox / LambdaNP2
                +6661.29 * (1. + eVBF_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -25488.2 * (1. + eVBF_1314_Hu_11 ) * CHu_11 / LambdaNP2
                +7958.49 * (1. + eVBF_1314_Hd_11 ) * CHd_11 / LambdaNP2
                -403218. * (1. + eVBF_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -146724. * (1. + eVBF_1314_HD ) * CHD / LambdaNP2
                -2068.51 * (1. + eVBF_1314_HB ) * CHB / LambdaNP2
                -73169. * (1. + eVBF_1314_HW ) * CHW / LambdaNP2
                -277706. * (1. + eVBF_1314_HWB ) * CHWB / LambdaNP2
                +3717947. * (1. + eVBF_1314_HG ) * CHG / LambdaNP2
                -1332.26 * (1. + eVBF_1314_DHB ) * CDHB / LambdaNP2
                -57988.3 * (1. + eVBF_1314_DHW ) * CDHW / LambdaNP2
                -3.974 * (1. + eVBF_1314_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3271430074.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +2236405708.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                -1592083582.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +671352943.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                -25330820741.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -8873926305.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -21771903.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -4537075343.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -16875188444.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                -4323212846.  * CHbox / LambdaNP2 * CHG / LambdaNP2
                -50906286.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                -3496862181.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -241626.  * CHbox / LambdaNP2 * DeltaGF()
                +58840234197.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                -546927041.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +100269958.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -30771657960.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +1454966168.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +5895593030.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -4729516531.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +2594397504.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +257395785857.  * CHQ1_11 / LambdaNP2 * CHG / LambdaNP2
                +3340777618.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -3383234583.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                -69529.5  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +36500368124.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -425621428.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +5549907093.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -227886267.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +4945482593.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +94660450.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +2247309189.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +181569259896.  * CHu_11 / LambdaNP2 * CHG / LambdaNP2
                +7930056446.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +312028889.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +5476.54  * CHu_11 / LambdaNP2 * DeltaGF()
                +21853942432.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1456368545.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +517126530.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                +1164674123.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +759387161.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                +324650282.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                +81220769204.  * CHd_11 / LambdaNP2 * CHG / LambdaNP2
                -2369701644.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +105178277.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -5198.95  * CHd_11 / LambdaNP2 * DeltaGF()
                +224068295761.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +39746870946.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +15515899450.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +84735126039.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +82070609683.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +493948743120.  * CHQ3_11 / LambdaNP2 * CHG / LambdaNP2
                -911930723.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +57027661887.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                +953317.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +7057462399.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +70118852.  * CHD / LambdaNP2 * CHB / LambdaNP2
                +6202012411.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +27034323178.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +8368684921.  * CHD / LambdaNP2 * CHG / LambdaNP2
                -246853417.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +4789818743.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +347140.  * CHD / LambdaNP2 * DeltaGF()
                +20921361708.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +10134277601.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -19711811521.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +26759457280.  * CHB / LambdaNP2 * CHG / LambdaNP2
                -259825404.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +187567928.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +1335.88  * CHB / LambdaNP2 * DeltaGF()
                +81758580794.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -9637836132.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +111418854959.  * CHW / LambdaNP2 * CHG / LambdaNP2
                +2183816569.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -10812326894.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +179893.  * CHW / LambdaNP2 * DeltaGF()
                +40451565403.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +65778494548.  * CHWB / LambdaNP2 * CHG / LambdaNP2
                -3148722084.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +10696630789.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +674511.  * CHWB / LambdaNP2 * DeltaGF()
                +30975002629457.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +17413631105.  * CHG / LambdaNP2 * CDHB / LambdaNP2
                +17140553238.  * CHG / LambdaNP2 * CDHW / LambdaNP2
                -133316.  * CHG / LambdaNP2 * DeltaGF()
                +817270273.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -445289766.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -7314.97  * CDHB / LambdaNP2 * DeltaGF()
                +5648073485.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +132993.  * CDHW / LambdaNP2 * DeltaGF()
                +4.669  * DeltaGF() * DeltaGF()
                ;
        }
        
    } else if (sqrt_s == 14.0) {

        mu += 
                +107633. * (1. + eVBF_1314_Hbox ) * CHbox / LambdaNP2
                +6335.01 * (1. + eVBF_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                -26240.6 * (1. + eVBF_1314_Hu_11 ) * CHu_11 / LambdaNP2
                +8229.3 * (1. + eVBF_1314_Hd_11 ) * CHd_11 / LambdaNP2
                -409260. * (1. + eVBF_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -146699. * (1. + eVBF_1314_HD ) * CHD / LambdaNP2
                -1925.35 * (1. + eVBF_1314_HB ) * CHB / LambdaNP2
                -73070. * (1. + eVBF_1314_HW ) * CHW / LambdaNP2
                -277778. * (1. + eVBF_1314_HWB ) * CHWB / LambdaNP2
                +3719893. * (1. + eVBF_1314_HG ) * CHG / LambdaNP2
                -1552.7 * (1. + eVBF_1314_DHB ) * CDHB / LambdaNP2
                -58350.4 * (1. + eVBF_1314_DHW ) * CDHW / LambdaNP2
                -3.981 * (1. + eVBF_1314_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3269983231.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +2235885970.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                -1627228122.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +695298429.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                -25635053723.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -8867803242.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -22979939.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -4493509720.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -16879386374.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +2630271412.  * CHbox / LambdaNP2 * CHG / LambdaNP2
                -52356996.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                -3551021676.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -241134.  * CHbox / LambdaNP2 * DeltaGF()
                +62257002671.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                -537233712.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +150611763.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -31799267126.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +1427861624.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +5934414012.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -4630147196.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +2465685361.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +223439537917.  * CHQ1_11 / LambdaNP2 * CHG / LambdaNP2
                +3587168499.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -3635799019.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                -69689.5  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +38618719334.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -442519098.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +5561766350.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -245947457.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +4897211353.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +9316191.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +2232780573.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +157021303025.  * CHu_11 / LambdaNP2 * CHG / LambdaNP2
                +8440904292.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +308055400.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +6633.43  * CHu_11 / LambdaNP2 * DeltaGF()
                +23257872182.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1762312900.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +524812124.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                +1159865847.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +728215639.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                +299670828.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                +110274206571.  * CHd_11 / LambdaNP2 * CHG / LambdaNP2
                -2519284516.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +108067822.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -5019.89  * CHd_11 / LambdaNP2 * DeltaGF()
                +235907707596.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +40130426682.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +15849947208.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +85649959630.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +83119060928.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +490755232594.  * CHQ3_11 / LambdaNP2 * CHG / LambdaNP2
                -872181852.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +60673871188.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                +966714.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +7052978076.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +46891497.  * CHD / LambdaNP2 * CHB / LambdaNP2
                +6219489473.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +26985280417.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +4747531209.  * CHD / LambdaNP2 * CHG / LambdaNP2
                -245202161.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +4871126017.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +347254.  * CHD / LambdaNP2 * DeltaGF()
                +20787528725.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +9955903360.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -19486367306.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -1534066207.  * CHB / LambdaNP2 * CHG / LambdaNP2
                -270355879.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +153406621.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +1511.09  * CHB / LambdaNP2 * DeltaGF()
                +82889261537.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -9629836656.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +91012980560.  * CHW / LambdaNP2 * CHG / LambdaNP2
                +2168001988.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -11004906528.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +180281.  * CHW / LambdaNP2 * DeltaGF()
                +40485063040.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +52676852369.  * CHWB / LambdaNP2 * CHG / LambdaNP2
                -3185702751.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +10816719458.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +674827.  * CHWB / LambdaNP2 * DeltaGF()
                +30979442270666.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +9391155829.  * CHG / LambdaNP2 * CDHB / LambdaNP2
                -13734550649.  * CHG / LambdaNP2 * CDHW / LambdaNP2
                -3611.25  * CHG / LambdaNP2 * DeltaGF()
                +856654866.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -447363518.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -7122.1  * CDHB / LambdaNP2 * DeltaGF()
                +5971057698.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +135609.  * CDHW / LambdaNP2 * DeltaGF()
                +4.672  * DeltaGF() * DeltaGF()
                ;

        }

    } else if (sqrt_s == 100.0) {

        mu += 
                +104969. * CHbox / LambdaNP2
                -3525.87 * CHQ1_11 / LambdaNP2
                -35635.5 * CHu_11 / LambdaNP2
                +13726.7 * CHd_11 / LambdaNP2
                -579211. * CHQ3_11 / LambdaNP2
                -143430. * CHD / LambdaNP2
                -2260.86 * CHB / LambdaNP2
                -80475.9 * CHW / LambdaNP2
                -272191. * CHWB / LambdaNP2
                +4215589. * CHG / LambdaNP2
                -2960.65 * CDHB / LambdaNP2
                -84270.9 * CDHW / LambdaNP2
                -3.88 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3208268289.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +640577058.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                -2333530710.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +1340367121.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                -32249555902.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -8632529472.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -41179954.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -4804327932.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -16541960489.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +4866232438.  * CHbox / LambdaNP2 * CHG / LambdaNP2
                -163912365.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                -5211551919.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -236432.  * CHbox / LambdaNP2 * DeltaGF()
                +289982236098.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                +153415514.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +1208483609.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -122301771007.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +991548689.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +8047047424.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -4352155892.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +2239328202.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +528233837541.  * CHQ1_11 / LambdaNP2 * CHG / LambdaNP2
                +14150831674.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -19914948592.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                -55517.4  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +175789955321.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -172256016.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +14552941810.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -1258545513.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +5745814717.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                -1157345104.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +3924207353.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +322137589492.  * CHu_11 / LambdaNP2 * CHG / LambdaNP2
                +38080960327.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +244388222.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +8990.23  * CHu_11 / LambdaNP2 * DeltaGF()
                +114149216773.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -648651558.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -133498412.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                +1350594822.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +40372504.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -409108037.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                +222288313506.  * CHd_11 / LambdaNP2 * CHG / LambdaNP2
                -14241804382.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                -1486246434.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                +3884.67  * CHd_11 / LambdaNP2 * DeltaGF()
                +972762017549.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +53556548420.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +13995262960.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +134636916617.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +95553641600.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +574850621737.  * CHQ3_11 / LambdaNP2 * CHG / LambdaNP2
                -3828389945.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +289903644291.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                +1277700.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6840717016.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -10227701.  * CHD / LambdaNP2 * CHB / LambdaNP2
                +6619475696.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +26119394951.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -12223717500.  * CHD / LambdaNP2 * CHG / LambdaNP2
                -304946977.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +6801959412.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +337114.  * CHD / LambdaNP2 * DeltaGF()
                +20778920170.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +7600258384.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -16660924799.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -12704419443.  * CHB / LambdaNP2 * CHG / LambdaNP2
                -681218711.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +216934920.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +75.474  * CHB / LambdaNP2 * DeltaGF()
                +139277063035.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -9293750336.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +6063950046.  * CHW / LambdaNP2 * CHG / LambdaNP2
                +2547774129.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -21699413253.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +189816.  * CHW / LambdaNP2 * DeltaGF()
                +42708725844.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -15793723421.  * CHWB / LambdaNP2 * CHG / LambdaNP2
                -4234806481.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +14067933466.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +657734.  * CHWB / LambdaNP2 * DeltaGF()
                +35035258653173.  * CHG / LambdaNP2 * CHG / LambdaNP2
                -22157506594.  * CHG / LambdaNP2 * CDHB / LambdaNP2
                -4062012166.  * CHG / LambdaNP2 * CDHW / LambdaNP2
                +320518.  * CHG / LambdaNP2 * DeltaGF()
                +3317004899.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -873660979.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -5611.68  * CDHB / LambdaNP2 * DeltaGF()
                +25594014104.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +195108.  * CDHW / LambdaNP2 * DeltaGF()
                +4.555  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muVBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eVBFint + eVBFpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeWBF(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 0.24) {

        mu += 
                +121115. * CHbox / LambdaNP2
                -137868. * CHL3_11 / LambdaNP2
                -203737. * CHD / LambdaNP2
                -24698.9 * CHW / LambdaNP2
                -379818. * CHWB / LambdaNP2
                -18173.1 * CDHW / LambdaNP2
                -286375. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//
//        }
          
    } else if (sqrt_s == 0.25) {

        mu += 
                +121115. * CHbox / LambdaNP2
                -137868. * CHL3_11 / LambdaNP2
                -203737. * CHD / LambdaNP2
                -24698.9 * CHW / LambdaNP2
                -379818. * CHWB / LambdaNP2
                -18173.1 * CDHW / LambdaNP2
                -286375. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//
//        }
        
    } else if (sqrt_s == 0.35) {

        mu += 
                +121119. * CHbox / LambdaNP2
                -218560. * CHL3_11 / LambdaNP2
                -203738. * CHD / LambdaNP2
                -39726.3 * CHW / LambdaNP2
                -379800. * CHWB / LambdaNP2
                -28867.1 * CDHW / LambdaNP2
                -286322. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//
//        }
//        
//    } else if (sqrt_s == 0.5) {
//        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//        }
//        
//    } else if (sqrt_s == 1.0) {
//        
//        if (FlagQuadraticTerms) {
//            //Add contributions that are quadratic in the effective coefficients
//        }
//        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeWBF()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeWBFint + eeeWBFpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}


double NPSMEFTd6::muepWBF(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 3.5) {

        mu += 
                +122379. * CHbox / LambdaNP2
                -218680. * CHL3_11 / LambdaNP2
                -217481. * CHQ3_11 / LambdaNP2
                -202145. * CHD / LambdaNP2
                -65129.9 * CHW / LambdaNP2
                -377695. * CHWB / LambdaNP2
                -58836.6 * CDHW / LambdaNP2
                -283479. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
          
    } else if (sqrt_s == 5.0) {

        mu += 
                +122344. * CHbox / LambdaNP2
                -235801. * CHL3_11 / LambdaNP2
                -235975. * CHQ3_11 / LambdaNP2
                -201446. * CHD / LambdaNP2
                -59025.7 * CHW / LambdaNP2
                -376522. * CHWB / LambdaNP2
                -65505.7 * CDHW / LambdaNP2
                -284078. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muepWBF()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muepZBF(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 3.5) {

        mu +=              
                +121495. * CHbox / LambdaNP2
                -207460. * CHL1_11 / LambdaNP2
                +18181.8 * CHQ1_11 / LambdaNP2
                +177060. * CHe_11 / LambdaNP2
                -54375.5 * CHu_11 / LambdaNP2
                +19796.1 * CHd_11 / LambdaNP2
                -204315. * CHL3_11 / LambdaNP2
                -227935. * CHQ3_11 / LambdaNP2
                -215293. * CHD / LambdaNP2
                -9090.91 * CHB / LambdaNP2
                -80034. * CHW / LambdaNP2
                -103483. * CHWB / LambdaNP2
                -11130. * CDHB / LambdaNP2
                -62531.9 * CDHW / LambdaNP2
                -311810. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
          
    } else if (sqrt_s == 5.0) {

        mu += 
                +120062. * CHbox / LambdaNP2
                -225482. * CHL1_11 / LambdaNP2
                +13203.2 * CHQ1_11 / LambdaNP2
                +192294. * CHe_11 / LambdaNP2
                -58926.3 * CHu_11 / LambdaNP2
                +20703.8 * CHd_11 / LambdaNP2
                -225821. * CHL3_11 / LambdaNP2
                -247881. * CHQ3_11 / LambdaNP2
                -215515. * CHD / LambdaNP2
                -4931.93 * CHB / LambdaNP2
                -72694.6 * CHW / LambdaNP2
                -104495. * CHWB / LambdaNP2
                -13100.4 * CDHB / LambdaNP2
                -70793.7 * CDHW / LambdaNP2
                -309838. * DeltaGF() / v() / v()
                ;
        
//        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients

//        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muepZBF()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muWH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {

        mu += 
                +121384. * (1. + eWH_2_Hbox ) * CHbox / LambdaNP2
                +1566381. * (1. + eWH_2_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -161022. * (1. + eWH_2_HD ) * CHD / LambdaNP2
                +861269. * (1. + eWH_2_HW ) * CHW / LambdaNP2
                -286326. * (1. + eWH_2_HWB ) * CHWB / LambdaNP2
                +134952. * (1. + eWH_2_DHW ) * CDHW / LambdaNP2
                -3.309 * (1. + eWH_2_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3674901007.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +94906032040.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -9760577677.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +52323973387.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17369701112.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +8159031979.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -200843.  * CHbox / LambdaNP2 * DeltaGF()
                +861942473518.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -126262085184.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +741985085727.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                -224539205242.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +175132158148.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2594431.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6485598279.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -69298406062.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +23066571526.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -10845008141.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +266944.  * CHD / LambdaNP2 * DeltaGF()
                +262768587566.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -123224658801.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +43858414906.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1425716.  * CHW / LambdaNP2 * DeltaGF()
                +20543305662.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -19302124580.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +475103.  * CHWB / LambdaNP2 * DeltaGF()
                +11346505598.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -222852.  * CDHW / LambdaNP2 * DeltaGF()
                +2.747  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 7.0) {

        mu += 
                +121106. * (1. + eWH_78_Hbox ) * CHbox / LambdaNP2
                +1791439. * (1. + eWH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159841. * (1. + eWH_78_HD ) * CHD / LambdaNP2
                +879865. * (1. + eWH_78_HW ) * CHW / LambdaNP2
                -283882. * (1. + eWH_78_HWB ) * CHWB / LambdaNP2
                +168277. * (1. + eWH_78_DHW ) * CDHW / LambdaNP2
                -3.267 * (1. + eWH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3674771194.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +108047517796.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -9689400943.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +53388185264.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17218267542.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +10208468152.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -198475.  * CHbox / LambdaNP2 * DeltaGF()
                +1662475732643.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -143111075159.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +902306554498.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                -254041092724.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +411946473144.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2937063.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6404964408.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -70239669039.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +22747064805.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -13449662568.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +262691.  * CHD / LambdaNP2 * DeltaGF()
                +316099657946.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -124897152630.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +52604696311.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1441125.  * CHW / LambdaNP2 * DeltaGF()
                +20242442452.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -23877923639.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +467191.  * CHWB / LambdaNP2 * DeltaGF()
                +29996764352.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -275407.  * CDHW / LambdaNP2 * DeltaGF()
                +2.696  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 8.0) {

        mu += 
                +121149. * (1. + eWH_78_Hbox ) * CHbox / LambdaNP2
                +1829711. * (1. + eWH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159724. * (1. + eWH_78_HD ) * CHD / LambdaNP2
                +883387. * (1. + eWH_78_HW ) * CHW / LambdaNP2
                -283633. * (1. + eWH_78_HWB ) * CHWB / LambdaNP2
                +173880. * (1. + eWH_78_DHW ) * CDHW / LambdaNP2
                -3.276 * (1. + eWH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3676053676.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +110461160461.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -9688527689.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +53761103761.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17210546211.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +10544320544.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -198198.  * CHbox / LambdaNP2 * DeltaGF()
                +1856000756001.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -146584010584.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +928671328671.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                -258635418635.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +471054621055.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2988095.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6402570403.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -70511434511.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +22763938764.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -13917217917.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +262632.  * CHD / LambdaNP2 * DeltaGF()
                +324948024948.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -125239085239.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +54290304290.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1449583.  * CHW / LambdaNP2 * DeltaGF()
                +20230580231.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -24690984691.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +466976.  * CHWB / LambdaNP2 * DeltaGF()
                +34690984691.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -284586.  * CDHW / LambdaNP2 * DeltaGF()
                +2.696  * DeltaGF() * DeltaGF()
                ;

        }

    } else if (sqrt_s == 13.0) {

        mu += 
                +121024. * (1. + eWH_1314_Hbox ) * CHbox / LambdaNP2
                +1954150. * (1. + eWH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159652. * (1. + eWH_1314_HD ) * CHD / LambdaNP2
                +892062. * (1. + eWH_1314_HW ) * CHW / LambdaNP2
                -283618. * (1. + eWH_1314_HWB ) * CHWB / LambdaNP2
                +192569. * (1. + eWH_1314_DHW ) * CDHW / LambdaNP2
                -3.263 * (1. + eWH_1314_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3674754036.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +121552106218.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -9677036966.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +53914318464.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17189607412.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +11654312733.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -198130.  * CHbox / LambdaNP2 * DeltaGF()
                +2780876874582.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -156644665202.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +1021329639889.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                -272177858439.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +761681153883.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -3262360.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6391536919.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -71182825485.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +22706753272.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -15419142229.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +262029.  * CHD / LambdaNP2 * DeltaGF()
                +354733021301.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -126471487248.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +59063902952.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1460030.  * CHW / LambdaNP2 * DeltaGF()
                +20195816219.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -27380838667.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +465963.  * CHWB / LambdaNP2 * DeltaGF()
                +57785843920.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -314175.  * CDHW / LambdaNP2 * DeltaGF()
                +2.689  * DeltaGF() * DeltaGF()
                ;

        }
      
    } else if (sqrt_s == 14.0) {

        mu += 
                +121073. * (1. + eWH_1314_Hbox ) * CHbox / LambdaNP2
                +1971095. * (1. + eWH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -159723. * (1. + eWH_1314_HD ) * CHD / LambdaNP2
                +893033. * (1. + eWH_1314_HW ) * CHW / LambdaNP2
                -283739. * (1. + eWH_1314_HWB ) * CHWB / LambdaNP2
                +195240. * (1. + eWH_1314_DHW ) * CDHW / LambdaNP2
                -3.265 * (1. + eWH_1314_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3674340113.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +121639376893.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -9675093033.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +54136131545.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17184162700.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +11861791432.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -198029.  * CHbox / LambdaNP2 * DeltaGF()
                +2961055819991.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -158964517525.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +1033215058416.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                -273497187365.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +817931631328.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -3305902.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6389874513.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -71141064474.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +22710947642.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -15578970143.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +261976.  * CHD / LambdaNP2 * DeltaGF()
                +358823020338.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -126591086110.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +59740372133.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1461549.  * CHW / LambdaNP2 * DeltaGF()
                +20186932064.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -27704024232.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +465787.  * CHWB / LambdaNP2 * DeltaGF()
                +62293379489.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -319678.  * CDHW / LambdaNP2 * DeltaGF()
                +2.688  * DeltaGF() * DeltaGF()
                ;

        }
          
    } else if (sqrt_s == 100.0) {

        mu += 
                +121131. * CHbox / LambdaNP2
                +2294774. * CHQ3_11 / LambdaNP2
                -159215. * CHD / LambdaNP2
                +907558. * CHW / LambdaNP2
                -282770. * CHWB / LambdaNP2
                +245365. * CDHW / LambdaNP2
                -3.257 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3673685094.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +77695327573.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -9650801107.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +55062998071.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -17130442077.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +14865028102.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -197073.  * CHbox / LambdaNP2 * DeltaGF()
                +16975085982720.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -343585605234.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +1245348544585.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                -332991359785.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                +5309026088415.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -3380187.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +6365573358.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -72517741800.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +22625283114.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -19543997987.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +260717.  * CHD / LambdaNP2 * DeltaGF()
                +436557335794.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -128218270279.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +71697005285.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1475400.  * CHW / LambdaNP2 * DeltaGF()
                +20098146129.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -34853619663.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +463420.  * CHWB / LambdaNP2 * DeltaGF()
                +427908732489.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -396814.  * CDHW / LambdaNP2 * DeltaGF()
                +2.673  * DeltaGF() * DeltaGF()
                ;

        }
          
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muWH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eWHint + eWHpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muZH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {

        mu += 
                +121306. * (1. + eZH_2_Hbox ) * CHbox / LambdaNP2
                -809998. * (1. + eZH_2_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +529178. * (1. + eZH_2_Hu_11 ) * CHu_11 / LambdaNP2
                -69391.9 * (1. + eZH_2_Hd_11 ) * CHd_11 / LambdaNP2
                +1568049. * (1. + eZH_2_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -16522.6 * (1. + eZH_2_HD ) * CHD / LambdaNP2
                +79360.5 * (1. + eZH_2_HB ) * CHB / LambdaNP2
                +712302. * (1. + eZH_2_HW ) * CHW / LambdaNP2
                +189306. * (1. + eZH_2_HWB ) * CHWB / LambdaNP2
                +10001. * (1. + eZH_2_DHB ) * CDHB / LambdaNP2
                +131065. * (1. + eZH_2_DHW ) * CDHW / LambdaNP2
                -2.543 * (1. + eZH_2_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3673819465.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -49374049260.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +32203450173.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -4143120691.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                +95882771790.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -993704890.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +4809204777.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +43175389196.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +11494319837.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +607664174.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +7959672460.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -154064.  * CHbox / LambdaNP2 * DeltaGF()
                +1002718710554.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                +356021620.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +388387222.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1276790626922.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +20648363919.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +61166132634.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -426368255818.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +33312295692.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +51511797262.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -135829044891.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                +1461046.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +820726931417.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +80914005.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +809140046.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +22978687251.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +220987474512.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +15454574878.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +313323299997.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +135872738454.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                -10266692559.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +223803.  * CHu_11 / LambdaNP2 * DeltaGF()
                +182185972748.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +534032430.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -2929977020.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                -27224325986.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                -2039032916.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -39558856847.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                -14255105674.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +1239278894.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -30852.3  * CHd_11 / LambdaNP2 * DeltaGF()
                +1002427420138.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -25391704696.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                -27094863579.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +715344531832.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +48783053371.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                -43465061333.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +204754183254.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2437802.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +486374082.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +4648428650.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -11296485096.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +5826536557.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +2383483833.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2437696216.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +48519.4  * CHD / LambdaNP2 * DeltaGF()
                +26334271936.  * CHB / LambdaNP2 * CHB / LambdaNP2
                -11478460692.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +67101984011.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +16000258925.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -4774249927.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +74640.3  * CHB / LambdaNP2 * DeltaGF()
                +198757160889.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +45449396382.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -17394569052.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +47230799107.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1080173.  * CHW / LambdaNP2 * DeltaGF()
                +50595527074.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +21248341263.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -2937502023.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +2167.51  * CHWB / LambdaNP2 * DeltaGF()
                +7231770075.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -4569861152.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +68714.3  * CDHB / LambdaNP2 * DeltaGF()
                +12710295498.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -211642.  * CDHW / LambdaNP2 * DeltaGF()
                +2.071  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 7.0) {

        mu += 
                +121146. * (1. + eZH_78_Hbox ) * CHbox / LambdaNP2
                -181981. * (1. + eZH_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +421787. * (1. + eZH_78_Hu_11 ) * CHu_11 / LambdaNP2
                -139069. * (1. + eZH_78_Hd_11 ) * CHd_11 / LambdaNP2
                +1712878. * (1. + eZH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15306.8 * (1. + eZH_78_HD ) * CHD / LambdaNP2
                +87131. * (1. + eZH_78_HB ) * CHB / LambdaNP2
                +717936. * (1. + eZH_78_HW ) * CHW / LambdaNP2
                +203131. * (1. + eZH_78_HWB ) * CHWB / LambdaNP2
                +17618.1 * (1. + eZH_78_DHB ) * CDHB / LambdaNP2
                +153199. * (1. + eZH_78_DHW ) * CDHW / LambdaNP2
                -2.507 * (1. + eZH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675199512.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -11699271728.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +25499280450.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -8428590118.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                +104475339061.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -931677642.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +5253325193.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +43448781126.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +12298809472.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1071213641.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +9259256029.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -151897.  * CHbox / LambdaNP2 * DeltaGF()
                +1604029479744.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                -43609088.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -2180454407.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -740395098338.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +10688356374.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +72296018490.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -161545506083.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +80410797610.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +77094326458.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -96755483843.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                +605446.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +986917273560.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -523309058.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -523309058.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +18145510444.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +184153329554.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +10492346605.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +257768959051.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +183674937857.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                -8534298548.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +177444.  * CHu_11 / LambdaNP2 * DeltaGF()
                +616632506214.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -2485718024.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -6101142558.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                -60973354847.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                -3549779774.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -85085691858.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                -57396101348.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +2804064367.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -52741.3  * CHd_11 / LambdaNP2 * DeltaGF()
                +1602982861628.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -24765832279.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                -2668004012.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +793101042257.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +100383759976.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                -37684793511.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +383816667394.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2500622.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +388352013.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +3566992281.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -9834080502.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +4256015874.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +2062042650.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2384339976.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +41041.6  * CHD / LambdaNP2 * DeltaGF()
                +24734202608.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +1810649339.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +62584274563.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +14856308055.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -3134621255.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +30909.6  * CHB / LambdaNP2 * DeltaGF()
                +218978675156.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +75160263399.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -15404910383.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +51868649427.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1038440.  * CHW / LambdaNP2 * DeltaGF()
                +49020975971.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +19207622869.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +43609088.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -62452.3  * CHWB / LambdaNP2 * DeltaGF()
                +11999040600.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -3774366578.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +50791.9  * CDHB / LambdaNP2 * DeltaGF()
                +26937333740.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -231567.  * CDHW / LambdaNP2 * DeltaGF()
                +1.93  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 8.0) {

        mu += 
                +121330. * (1. + eZH_78_Hbox ) * CHbox / LambdaNP2
                -177103. * (1. + eZH_78_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +428308. * (1. + eZH_78_Hu_11 ) * CHu_11 / LambdaNP2
                -143029. * (1. + eZH_78_Hd_11 ) * CHd_11 / LambdaNP2
                +1747483. * (1. + eZH_78_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15210.6 * (1. + eZH_78_HD ) * CHD / LambdaNP2
                +87647.1 * (1. + eZH_78_HB ) * CHB / LambdaNP2
                +720891. * (1. + eZH_78_HW ) * CHW / LambdaNP2
                +204226. * (1. + eZH_78_HWB ) * CHWB / LambdaNP2
                +18685.3 * (1. + eZH_78_DHB ) * CDHB / LambdaNP2
                +158453. * (1. + eZH_78_DHW ) * CDHW / LambdaNP2
                -2.505 * (1. + eZH_78_DeltaGF ) * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675365196.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -10163806552.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +25327613105.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -8816479932.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                +106947241526.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -930860871.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +5288966104.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +43658346334.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +12375549567.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1123244930.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +9578074032.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -151703.  * CHbox / LambdaNP2 * DeltaGF()
                +1778081123245.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                +35455964.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -1737342221.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -811551552971.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +13086264360.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +74950716211.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -163111615374.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +83729258261.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +87074528436.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -106665721174.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                +591196.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +1092610977166.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -1347326620.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +106367891.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +17553715785.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +188551623883.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +10587150759.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +263065522621.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +207447525174.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                -8665437527.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +187653.  * CHu_11 / LambdaNP2 * DeltaGF()
                +686179265352.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1914622039.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -6343603744.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                -62547510991.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                -3417954900.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -88267621614.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                -64889731953.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +3141398383.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -62660.3  * CHd_11 / LambdaNP2 * DeltaGF()
                +1777691107644.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -24993440647.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                -3548787406.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +816784853212.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +105392852078.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                -39077790384.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +435314139838.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2557608.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +387001844.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +3579102255.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -9832470572.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +4229364629.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +2119734789.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2445575096.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +40930.4  * CHD / LambdaNP2 * DeltaGF()
                +25425117005.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +1812154304.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +63913274713.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +15354559637.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -3222592540.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +30452.9  * CHB / LambdaNP2 * DeltaGF()
                +225159551837.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +78620053893.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -15836406184.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +53421500496.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1041807.  * CHW / LambdaNP2 * DeltaGF()
                +50010636789.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +19701106226.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +42547156.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -63654.5  * CHWB / LambdaNP2 * DeltaGF()
                +13726776344.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -4143029358.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +51624.3  * CDHB / LambdaNP2 * DeltaGF()
                +31073606581.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -239342.  * CDHW / LambdaNP2 * DeltaGF()
                +1.928  * DeltaGF() * DeltaGF()
                ;

        }

    } else if (sqrt_s == 13.0) {

        mu += 
                +121156. * (1. + eZH_1314_Hbox ) * CHbox / LambdaNP2
                -152858. * (1. + eZH_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +448053. * (1. + eZH_1314_Hu_11 ) * CHu_11 / LambdaNP2
                -155982. * (1. + eZH_1314_Hd_11 ) * CHd_11 / LambdaNP2
                +1861934. * (1. + eZH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15515.2 * (1. + eZH_1314_HD ) * CHD / LambdaNP2
                +88537.4 * (1. + eZH_1314_HB ) * CHB / LambdaNP2
                +727998. * (1. + eZH_1314_HW ) * CHW / LambdaNP2
                +207593. * (1. + eZH_1314_HWB ) * CHWB / LambdaNP2
                +21551.9 * (1. + eZH_1314_DHB ) * CDHB / LambdaNP2
                +174974. * (1. + eZH_1314_DHW ) * CDHW / LambdaNP2
                -2.503 * (1. + eZH_1314_DeltaGF ) * DeltaGF()
                ;

        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3674586084.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -7734317083.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +27179369506.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -9675927560.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                +111409962227.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -926695379.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +5362728139.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +43894870618.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +12626116426.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1301232040.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +10649203940.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -151696.  * CHbox / LambdaNP2 * DeltaGF()
                +2613937233029.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                +3600804886.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +17651004.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1148197832457.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +15151604476.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +85083136231.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -160677092527.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +91877007802.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +135381438204.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -161000105906.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                +570543.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +1594697638296.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -794295195.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1871006460.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +19582006566.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +206416140078.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +9725703393.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +282225438628.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +325023828856.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                -8901401490.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +196703.  * CHu_11 / LambdaNP2 * DeltaGF()
                +1020210399972.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -2453489604.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -6894499947.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                -71287111237.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                -3706710912.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -96935785646.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                -102783563385.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +3348395524.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -58929.3  * CHd_11 / LambdaNP2 * DeltaGF()
                +2617290923854.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -30387986726.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                -6349066262.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +892911356656.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +117823984185.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                -55072898648.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +691190383733.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2687456.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +381279345.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +3539008720.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -10018727716.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +4179740177.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +2238129700.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2645903202.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +40473.4  * CHD / LambdaNP2 * DeltaGF()
                +27706781516.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +1664489710.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +67984608324.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +16902601758.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -3616690790.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +29057.  * CHB / LambdaNP2 * DeltaGF()
                +245931443499.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +89529424224.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -17741024464.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +58596039115.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1050013.  * CHW / LambdaNP2 * DeltaGF()
                +53214247891.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +21308292442.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +93550323.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -69672.8  * CHWB / LambdaNP2 * DeltaGF()
                +22242030572.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -5870724044.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +53572.1  * CDHB / LambdaNP2 * DeltaGF()
                +51563878985.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -264803.  * CDHW / LambdaNP2 * DeltaGF()
                +1.919  * DeltaGF() * DeltaGF()
                ;

        }
        
    } else if (sqrt_s == 14.0) {

        mu += 
                +120925. * (1. + eZH_1314_Hbox ) * CHbox / LambdaNP2
                -148239. * (1. + eZH_1314_HQ1_11 ) * CHQ1_11 / LambdaNP2
                +450417. * (1. + eZH_1314_Hu_11 ) * CHu_11 / LambdaNP2
                -157929. * (1. + eZH_1314_Hd_11 ) * CHd_11 / LambdaNP2
                +1877941. * (1. + eZH_1314_HQ3_11 ) * CHQ3_11 / LambdaNP2
                -15564.6 * (1. + eZH_1314_HD ) * CHD / LambdaNP2
                +88598.7 * (1. + eZH_1314_HB ) * CHB / LambdaNP2
                +728872. * (1. + eZH_1314_HW ) * CHW / LambdaNP2
                +207768. * (1. + eZH_1314_HWB ) * CHWB / LambdaNP2
                +21918.2 * (1. + eZH_1314_DHB ) * CDHB / LambdaNP2
                +177277. * (1. + eZH_1314_DHW ) * CDHW / LambdaNP2
                -2.5 * (1. + eZH_1314_DeltaGF ) * DeltaGF()
                ;

        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3673414003.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -7823983909.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +28413684110.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -9276683375.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                +110786532997.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -925593053.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +5366846525.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +44026212445.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +12652693081.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1313655375.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +10760991028.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -151747.  * CHbox / LambdaNP2 * DeltaGF()
                +2776252354650.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                +2394559561.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -1197279780.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -1216116982216.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +14306359950.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +85966284601.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -159940614923.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +94481338399.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +145297085023.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -171543054181.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                +566538.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +1691357236359.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -878005172.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -3352383385.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +20212940200.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +208887008716.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +10073113885.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +286365377862.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +347717186552.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                -9351553271.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +182089.  * CHu_11 / LambdaNP2 * DeltaGF()
                +1084496025031.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -4390025861.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -6462453306.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                -72298138629.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                -3224673542.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -98057214010.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                -109276523738.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +4680565755.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -67802.6  * CHd_11 / LambdaNP2 * DeltaGF()
                +2779924012643.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -32627007439.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                -7261900961.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +903563104626.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +118905845918.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                -57601928419.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +742540148782.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -2706544.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +380272022.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +3534034673.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -9878691613.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +4163005651.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +2261725360.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2718958526.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +40396.4  * CHD / LambdaNP2 * DeltaGF()
                +28014750487.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +1821461639.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +68525909134.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +17121100859.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -3676447112.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +28836.3  * CHB / LambdaNP2 * DeltaGF()
                +248698955972.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +91320519779.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -17852239711.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +59404233581.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1052360.  * CHW / LambdaNP2 * DeltaGF()
                +53662079755.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +21522301331.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +71836787.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -70514.8  * CHWB / LambdaNP2 * DeltaGF()
                +23918457265.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -6240222215.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +53983.5  * CDHB / LambdaNP2 * DeltaGF()
                +55646371444.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -267479.  * CDHW / LambdaNP2 * DeltaGF()
                +1.918  * DeltaGF() * DeltaGF()
                ;

        }
          
    } else if (sqrt_s == 100.0) {
 
        mu += 
                +121366. * CHbox / LambdaNP2
                -377.758 * CHQ1_11 / LambdaNP2
                +484028. * CHu_11 / LambdaNP2
                -197401. * CHd_11 / LambdaNP2
                +2175431. * CHQ3_11 / LambdaNP2
                -15095.2 * CHD / LambdaNP2
                +91689.3 * CHB / LambdaNP2
                +741765. * CHW / LambdaNP2
                +215382. * CHWB / LambdaNP2
                +31489.9 * CDHB / LambdaNP2
                +223769. * CDHW / LambdaNP2
                -2.503 * DeltaGF()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3673315201.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -24827742520.  * CHbox / LambdaNP2 * CHQ1_11 / LambdaNP2
                +5997280145.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -15610456331.  * CHbox / LambdaNP2 * CHd_11 / LambdaNP2
                +142896645512.  * CHbox / LambdaNP2 * CHQ3_11 / LambdaNP2
                -912541553.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +5543970988.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +45117860381.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +13056814748.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1872166818.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +13612873980.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -151281.  * CHbox / LambdaNP2 * DeltaGF()
                +15217588395286.  * CHQ1_11 / LambdaNP2 * CHQ1_11 / LambdaNP2
                +6648534300.  * CHQ1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -1813236627.  * CHQ1_11 / LambdaNP2 * CHd_11 / LambdaNP2
                -6146872166818.  * CHQ1_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +16261680266.  * CHQ1_11 / LambdaNP2 * CHD / LambdaNP2
                +92272589906.  * CHQ1_11 / LambdaNP2 * CHB / LambdaNP2
                -121607736476.  * CHQ1_11 / LambdaNP2 * CHW / LambdaNP2
                +125859776368.  * CHQ1_11 / LambdaNP2 * CHWB / LambdaNP2
                +851722574796.  * CHQ1_11 / LambdaNP2 * CDHB / LambdaNP2
                -974191598670.  * CHQ1_11 / LambdaNP2 * CDHW / LambdaNP2
                +986792.  * CHQ1_11 / LambdaNP2 * DeltaGF()
                +9145663342400.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                -7252946510.  * CHu_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +3626473255.  * CHu_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                +23967935932.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +236878210940.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                +2447869447.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +316249622242.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                +2092127530976.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                -9398609852.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                +261492.  * CHu_11 / LambdaNP2 * DeltaGF()
                +6071622846782.  * CHd_11 / LambdaNP2 * CHd_11 / LambdaNP2
                +10275007555.  * CHd_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -2021789060.  * CHd_11 / LambdaNP2 * CHD / LambdaNP2
                -99628286491.  * CHd_11 / LambdaNP2 * CHB / LambdaNP2
                +1087941976.  * CHd_11 / LambdaNP2 * CHW / LambdaNP2
                -125273496525.  * CHd_11 / LambdaNP2 * CHWB / LambdaNP2
                -699652462980.  * CHd_11 / LambdaNP2 * CDHB / LambdaNP2
                +5258386219.  * CHd_11 / LambdaNP2 * CDHW / LambdaNP2
                -197117.  * CHd_11 / LambdaNP2 * DeltaGF()
                +15207011181626.  * CHQ3_11 / LambdaNP2 * CHQ3_11 / LambdaNP2
                -15469960713.  * CHQ3_11 / LambdaNP2 * CHD / LambdaNP2
                +12187972197.  * CHQ3_11 / LambdaNP2 * CHB / LambdaNP2
                +1121970383802.  * CHQ3_11 / LambdaNP2 * CHW / LambdaNP2
                +181767905712.  * CHQ3_11 / LambdaNP2 * CHWB / LambdaNP2
                -195421577516.  * CHQ3_11 / LambdaNP2 * CDHB / LambdaNP2
                +4751707464491.  * CHQ3_11 / LambdaNP2 * CDHW / LambdaNP2
                -3325128.  * CHQ3_11 / LambdaNP2 * DeltaGF()
                +359655485.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +3343880326.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -9652493200.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +3828921124.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +2526412814.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -3321275310.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +38800.3  * CHD / LambdaNP2 * DeltaGF()
                +33142943488.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +4239951647.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +77742520399.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +20695074041.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -4372922333.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +20547.7  * CHB / LambdaNP2 * DeltaGF()
                +301420368691.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +123306134784.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -22121486854.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +72574796011.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                -1059252.  * CHW / LambdaNP2 * DeltaGF()
                +61508008462.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +25149592022.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                +791779994.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -87624.2  * CHWB / LambdaNP2 * DeltaGF()
                +156134783923.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -30417044424.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +50880.7  * CDHB / LambdaNP2 * DeltaGF()
                +379752190994.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -325976.  * CDHW / LambdaNP2 * DeltaGF()
                +1.888  * DeltaGF() * DeltaGF()
                ;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muZH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eZHint + eZHpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::mueeZH(const double sqrt_s) const
{
    double mu = 1.0;

    if (sqrt_s == 0.24) {

        mu += 
                +121274. * CHbox / LambdaNP2
                +898665. * CHL1_11 / LambdaNP2
                -771954. * CHe_11 / LambdaNP2
                +898707. * CHL3_11 / LambdaNP2
                -48862.4 * CHD / LambdaNP2
                +122450. * CHB / LambdaNP2
                +540173. * CHW / LambdaNP2
                +231005. * CHWB / LambdaNP2
                +17546.8 * CDHB / LambdaNP2
                +53438. * CDHW / LambdaNP2
                -219041. * DeltaGF() / v() / v()
                ;
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675426077.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +54485097809.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                -46770044497.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                +27925614978.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -2963689027.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +7410964655.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +32778272185.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +13976324406.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1060280413.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +3239190664.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +13278230207.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +350692637058.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                -16791201.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +458706237931.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -33418100915.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                -99361934346.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +397048946352.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -140470993200.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -37751909999.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +49003442196.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +121320627991.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +350692637058.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +115578037109.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +5508101755.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                -227197548485.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -28146251364.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -373692385190.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -60063218873.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                +8570229200.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -57791117454.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +159537402401.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -20950633868.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -102497691210.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                +255108723029.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                -153509361095.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -34641339938.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +33471580892.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                +69830408866.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +817983377.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +2951641340.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -19170766518.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +3694652002.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +1330198976.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2270329947.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -6235664512.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                +44610024347.  * CHB / LambdaNP2 * CHB / LambdaNP2
                -47892704223.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +142103937537.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +24488875829.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -9691041894.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +1511208127.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +114906389052.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -64104609185.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -18962555621.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +26755100327.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +71018386366.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +114192762992.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +39197968265.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -14538661741.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +6502392746.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                +3623793132.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -3355511712.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -1596675342.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +1827722274.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +7773906473.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +12874653682.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                ;
            }

    } else if (sqrt_s == 0.25) {

        mu += 
                +121285. * CHbox / LambdaNP2
                +975147. * CHL1_11 / LambdaNP2
                -837552. * CHe_11 / LambdaNP2
                +975232. * CHL3_11 / LambdaNP2
                -48866.5 * CHD / LambdaNP2
                +128421. * CHB / LambdaNP2
                +568430. * CHW / LambdaNP2
                +244081. * CHWB / LambdaNP2
                +19815.3 * CDHB / LambdaNP2
                +60117.5 * CDHW / LambdaNP2
                -219060. * DeltaGF() / v() / v()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675482788.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +59120486986.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                -50749370277.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                +32567170445.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -2963769941.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +7764483627.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +34473131822.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +14796389589.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +1194374475.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +3644584383.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +13278337532.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +412896725441.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                +25188917.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +562497900924.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -36263014274.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                -115906801008.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +455554156171.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -153052057095.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -45923173804.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +59726448363.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +131659109992.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +412896725441.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +125428211587.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +5978799328.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                -261733837112.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -29710327456.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -419063811923.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -73135600336.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                +10272208228.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -62706968934.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +201095717884.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -23790302267.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -118803526448.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                +305936188077.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                -167963895886.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -42451301427.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +42266330814.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                +80151973132.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +818010076.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +3191225861.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -20217674223.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +3648824517.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +1490134341.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -2550545760.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -6236146096.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                +51267842150.  * CHB / LambdaNP2 * CHB / LambdaNP2
                -56297229219.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +156855583543.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +28942485307.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -11600167926.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +1356842989.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +130331654072.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -68648194794.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -22669605374.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +31220990764.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +74887489505.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +121742233417.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +44950041981.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -16650545760.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +7423173804.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                +4592359362.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -4199454240.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -1781276238.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +2355835432.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +8736943745.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +12874895046.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                ;
        }
        
    } else if (sqrt_s == 0.35) {

        mu += 
                +121274. * CHbox / LambdaNP2
                +1911269. * CHL1_11 / LambdaNP2
                -1641658. * CHe_11 / LambdaNP2
                +1911346. * CHL3_11 / LambdaNP2
                -48852.8 * CHD / LambdaNP2
                +173326. * CHB / LambdaNP2
                +785469. * CHW / LambdaNP2
                +344545. * CHWB / LambdaNP2
                +59232.1 * CDHB / LambdaNP2
                +167941. * CDHW / LambdaNP2
                -219057. * DeltaGF() / v() / v()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675511160.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +115959263306.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                -99430154518.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                +89347744654.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -2963789605.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +10512096145.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +47526299360.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +20891368815.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +3583736538.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +10181208054.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +13278601530.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +1586155767130.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                +156079288.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +2656391446855.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -70975651631.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                -351880755424.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +1268113001405.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -311112845325.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -215339472452.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +308906664586.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +258078663961.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +1586155767130.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +245824879039.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +11824410801.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                -745122522241.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -41626346184.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -1028726393008.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -365799906352.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                +35689870454.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -122832839082.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +1121585765569.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -58567348213.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -352895270798.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                +1058498517247.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                -340533791166.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -208784142344.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +262082878102.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                +206650538474.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +818011550.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +5175433120.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -28678164508.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +3162010301.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +3480412049.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -6795848291.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -6236147963.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                +144295302013.  * CHB / LambdaNP2 * CHB / LambdaNP2
                -176767597940.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +331738723271.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +94267207741.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -40281723115.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +251287654.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +342000936476.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -104378024036.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -78660839707.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +93899641018.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +104473232402.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +208178554706.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +126224442017.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -47734509131.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +14618386140.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                +30151396910.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -24436553769.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -3377555798.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +17923365070.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +23790385516.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +12874980490.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                ;
        }
        
    } else if (sqrt_s == 0.5) {

        mu += 
                +121276. * CHbox / LambdaNP2
                +3900301. * CHL1_11 / LambdaNP2
                -3349835. * CHe_11 / LambdaNP2
                +3900547. * CHL3_11 / LambdaNP2
                -48890.4 * CHD / LambdaNP2
                +209412. * CHB / LambdaNP2
                +959487. * CHW / LambdaNP2
                +425160. * CHWB / LambdaNP2
                +169814. * CDHB / LambdaNP2
                +455465. * CDHW / LambdaNP2
                -219056. * DeltaGF() / v() / v()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675401700.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +236728674258.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                -202545537741.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                +210329983633.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -2963728199.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +12656588233.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +57974340473.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +25785537037.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +10333503458.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +27596487214.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +13278717376.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +6606008341986.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                +1231938896.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +12159236901850.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -144602840499.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                -927737984196.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +3207845690854.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -656781823622.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -1017123950652.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +1626755952905.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +526714594956.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +6606008341986.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +501751113145.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +24348779501.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                -1905897467486.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -50808679889.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -2442037274951.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -1855018391088.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                +98799739533.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -250286865771.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +5604266028405.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -132107460270.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -927034019113.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                +2950194470354.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                -697083824642.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -1005860509319.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +1501274176801.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                +475677126415.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +817972229.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +6749652417.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -35505851710.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +2807447951.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +7823199169.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -17672656236.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -6236215484.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                +367733760406.  * CHB / LambdaNP2 * CHB / LambdaNP2
                -471463015434.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +689815385157.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +256788863272.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -113234543567.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                -819239366.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +841713451013.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -133330986783.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -220235476320.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +246350820999.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +128499146442.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +376955702997.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +323683145316.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -125729923796.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +20282114007.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                +178701536404.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -136074690695.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -5307016772.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +113797715633.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +62975836398.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +12874641418.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                ;
        }
        
    } else if (sqrt_s == 1.0) {

        mu += 
                +121352. * CHbox / LambdaNP2
                +15601820. * CHL1_11 / LambdaNP2
                -13395670. * CHe_11 / LambdaNP2
                +15603075. * CHL3_11 / LambdaNP2
                -48870.4 * CHD / LambdaNP2
                +244038. * CHB / LambdaNP2
                +1129275. * CHW / LambdaNP2
                +503216. * CHWB / LambdaNP2
                +900690. * CDHB / LambdaNP2
                +2321462. * CDHW / LambdaNP2
                -219093. * DeltaGF() / v() / v()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675478506.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +945493018513.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                -795956620646.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                +929804283025.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -2963837465.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +14994116724.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +67943598996.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                +30682852212.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +55235723251.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +139954894886.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +13279337935.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +105702855349859.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                +23533103232.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +207201129588955.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -573456855977.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                -4551615939755.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +15263884530907.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -2711405710700.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -18089739566991.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +31865155318481.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +2105104330091.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +105695010982115.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +2016002510198.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +109003137747.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                -9117037966740.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -56165673047.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -11104879196737.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -35275964857233.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                +503373078130.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -985576561029.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +101537496077816.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -557768120489.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -4520238468779.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                +14965798556636.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                -2766316284907.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -18050517728271.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +31229761531221.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                +2065882491371.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +818010668.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +8124568560.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -41844053969.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +2711954816.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +34010982115.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -87341386884.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -6236350800.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                +1759021022906.  * CHB / LambdaNP2 * CHB / LambdaNP2
                -2322011295890.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +2752039535613.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                +1284358330719.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -578365233762.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                -2285064324.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +3931283338563.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -161829306558.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -1121901474741.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +1197128961406.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +152248980232.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +1322952620019.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +1560401631629.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -620567932225.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +25876215877.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                +3789457169752.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -2790084719172.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -11698305617.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +2501019767807.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +314784279887.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +12874960778.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                ;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueeZH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeeZHint + eeeZHpar;
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double mu = ((sigmaWH + sigmaZH) / (sigmaWH_SM + sigmaZH_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muVBFpVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaVBF_SM = computeSigmaVBF(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;
    double mu = ((sigmaWH + sigmaZH + sigmaVBF) / (sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muttH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {

        mu += 
                +418178. * (1. + ettH_2_HG ) * CHG / LambdaNP2
                -4002.18 * (1. + ettH_2_G ) * CG / LambdaNP2
                +566474. * (1. + ettH_2_uG_33r ) * CuG_33r / LambdaNP2
                -2.787 * (1. + ettH_2_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +53209042933.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +4957481622.  * CHG / LambdaNP2 * CG / LambdaNP2
                +131440554928.  * CHG / LambdaNP2 * CuG_33r / LambdaNP2
                -585682.  * CHG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +22878805337.  * CG / LambdaNP2 * CG / LambdaNP2
                +15950038025.  * CG / LambdaNP2 * CuG_33r / LambdaNP2
                +5886.55  * CG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +95095983223.  * CuG_33r / LambdaNP2 * CuG_33r / LambdaNP2
                -794534.  * CuG_33r / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1.97  * deltaG_hff(quarks[TOP]).real() * deltaG_hff(quarks[TOP]).real()
                ;

        }
        
    } else if (sqrt_s == 7.0) {

        mu += 
                +526637. * (1. + ettH_78_HG ) * CHG / LambdaNP2
                -84978.6 * (1. + ettH_78_G ) * CG / LambdaNP2
                +814927. * (1. + ettH_78_uG_33r ) * CuG_33r / LambdaNP2
                -2.772 * (1. + ettH_78_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +402731665820.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +165932157331.  * CHG / LambdaNP2 * CG / LambdaNP2
                +583503668608.  * CHG / LambdaNP2 * CuG_33r / LambdaNP2
                -747043.  * CHG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1450193177925.  * CG / LambdaNP2 * CG / LambdaNP2
                +817751101764.  * CG / LambdaNP2 * CuG_33r / LambdaNP2
                +109850.  * CG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +648865153538.  * CuG_33r / LambdaNP2 * CuG_33r / LambdaNP2
                -1154681.  * CuG_33r / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1.938  * deltaG_hff(quarks[TOP]).real() * deltaG_hff(quarks[TOP]).real()
                ;

        }
        
    } else if (sqrt_s == 8.0) {

        mu += 
                +532535. * (1. + ettH_78_HG ) * CHG / LambdaNP2
                -86299.4 * (1. + ettH_78_G ) * CG / LambdaNP2
                +830579. * (1. + ettH_78_uG_33r ) * CuG_33r / LambdaNP2
                -2.766 * (1. + ettH_78_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +467736526946.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +177213572854.  * CHG / LambdaNP2 * CG / LambdaNP2
                +641301397206.  * CHG / LambdaNP2 * CuG_33r / LambdaNP2
                -766329.  * CHG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1666746506986.  * CG / LambdaNP2 * CG / LambdaNP2
                +913029940120.  * CG / LambdaNP2 * CuG_33r / LambdaNP2
                +112569.  * CG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +722898203593.  * CuG_33r / LambdaNP2 * CuG_33r / LambdaNP2
                -1167504.  * CuG_33r / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1.936  * deltaG_hff(quarks[TOP]).real() * deltaG_hff(quarks[TOP]).real()
                ;

        }
        
    } else if (sqrt_s == 13.0) {

        mu += 
                +544404. * (1. + ettH_1314_HG ) * CHG / LambdaNP2
                -85840.1 * (1. + ettH_1314_G ) * CG / LambdaNP2
                +872840. * (1. + ettH_1314_uG_33r ) * CuG_33r / LambdaNP2
                -2.763 * (1. + ettH_1314_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +778438030560.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +211941143181.  * CHG / LambdaNP2 * CG / LambdaNP2
                +883418222977.  * CHG / LambdaNP2 * CuG_33r / LambdaNP2
                -778316.  * CHG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +2648774541814.  * CG / LambdaNP2 * CG / LambdaNP2
                +1310979060555.  * CG / LambdaNP2 * CuG_33r / LambdaNP2
                +54406.2  * CG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1044599712681.  * CuG_33r / LambdaNP2 * CuG_33r / LambdaNP2
                -1237329.  * CuG_33r / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1.929  * deltaG_hff(quarks[TOP]).real() * deltaG_hff(quarks[TOP]).real()
                ;

        }
        
    } else if (sqrt_s == 14.0) {

        mu += 
                +544972. * (1. + ettH_1314_HG ) * CHG / LambdaNP2
                -85093.7 * (1. + ettH_1314_G ) * CG / LambdaNP2
                +877274. * (1. + ettH_1314_uG_33r ) * CuG_33r / LambdaNP2
                -2.76 * (1. + ettH_1314_DeltagHt ) * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +837632364920.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +212960448909.  * CHG / LambdaNP2 * CG / LambdaNP2
                +923486288352.  * CHG / LambdaNP2 * CuG_33r / LambdaNP2
                -778881.  * CHG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +2829939361028.  * CG / LambdaNP2 * CG / LambdaNP2
                +1378821612816.  * CG / LambdaNP2 * CuG_33r / LambdaNP2
                +73649.5  * CG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1100678794461.  * CuG_33r / LambdaNP2 * CuG_33r / LambdaNP2
                -1260250.  * CuG_33r / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1.926  * deltaG_hff(quarks[TOP]).real() * deltaG_hff(quarks[TOP]).real()
                ;

        }
        
    } else if (sqrt_s == 100.0) {

        mu += 
                +536172. * CHG / LambdaNP2
                -33677.8 * CG / LambdaNP2
                +923472. * CuG_33r / LambdaNP2
                -2.766 * deltaG_hff(quarks[TOP]).real()
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +4963225970074.  * CHG / LambdaNP2 * CHG / LambdaNP2
                +174676642151.  * CHG / LambdaNP2 * CG / LambdaNP2
                +2558965254882.  * CHG / LambdaNP2 * CuG_33r / LambdaNP2
                -1088040.  * CHG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +15399759066701.  * CG / LambdaNP2 * CG / LambdaNP2
                +4186533096627.  * CG / LambdaNP2 * CuG_33r / LambdaNP2
                +50819.8  * CG / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +4060677149379.  * CuG_33r / LambdaNP2 * CuG_33r / LambdaNP2
                -1102830.  * CuG_33r / LambdaNP2 * deltaG_hff(quarks[TOP]).real()
                +1.92  * deltaG_hff(quarks[TOP]).real() * deltaG_hff(quarks[TOP]).real()
                ;

        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::muttH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += ettHint + ettHpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = computeSigmaggH(sqrt_s);
    double sigmattH_SM = computeSigmattH(sqrt_s);
    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    double mu = ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
} 

double NPSMEFTd6::mueettH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 0.5) {

        mu += 
                +121314. * CHbox / LambdaNP2
                +83662.1 * CHL1_11 / LambdaNP2
                +41909.2 * CHe_11 / LambdaNP2
                -31246.3 * CHu_11 / LambdaNP2
                +84209. * CHL3_11 / LambdaNP2
                -121041. * CuH_33r / LambdaNP2
                -52181.4 * CHD / LambdaNP2
                +138812. * CHB / LambdaNP2
                +130180. * CHW / LambdaNP2
                -251963. * CHWB / LambdaNP2
                -1952.9 * CDHB / LambdaNP2
                +1991.95 * CDHW / LambdaNP2
                -105027. * DeltaGF() / v() / v()
                +1384291. * CuW_33r / LambdaNP2
                +1696715. * CuB_33r / LambdaNP2
                -1.188 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3674999024.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +5071632231.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                +2536773034.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                -1898371285.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                -7627270242.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -7335014256.  * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
                -3158614225.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +8434128813.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +7919735968.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -15270124595.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                -86810139.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +118325977.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +6352732102.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +83853806195.  * CHbox / LambdaNP2 * CuW_33r / LambdaNP2
                +101781392806.  * CHbox / LambdaNP2 * CuB_33r / LambdaNP2
                -71233.1  * CHbox / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +10202320041.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                +8983322.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                -3612975042.  * CHL1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +5157989298.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -4768152170.  * CHL1_11 / LambdaNP2 * CuH_33r / LambdaNP2
                -3647502246.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                +2985978206.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +7840487443.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -13708159200.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -1263156661.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +2429422333.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +7628012342.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +76271921259.  * CHL1_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +38541967738.  * CHL1_11 / LambdaNP2 * CuB_33r / LambdaNP2
                -147461.  * CHL1_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +10199976565.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +3111236964.  * CHe_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +3860094520.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -2820724134.  * CHe_11 / LambdaNP2 * CuH_33r / LambdaNP2
                +812912549.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                +6149279381.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -124985353.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -4091317424.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -2862969183.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                +68761473.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -1918134594.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                -2193102371.  * CHe_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +82367300707.  * CHe_11 / LambdaNP2 * CuB_33r / LambdaNP2
                +126362.  * CHe_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3588759130.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +10190876069.  * CHu_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +1833300785.  * CHu_11 / LambdaNP2 * CuH_33r / LambdaNP2
                +3421864625.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +2496465258.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                -6649220794.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +10326016482.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                -213748389.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +89071593.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                -6902862946.  * CHu_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                -77173885873.  * CHu_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +30586923407.  * CHu_11 / LambdaNP2 * CuB_33r / LambdaNP2
                +219454.  * CHu_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +15778619693.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +7760457759.  * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
                +6652072023.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -3087528805.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                -14483458970.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                +22703980002.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -1509612155.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +2239600828.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                -13192594618.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                -158652501660.  * CHL3_11 / LambdaNP2 * CuW_33r / LambdaNP2
                -32807874077.  * CHL3_11 / LambdaNP2 * CuB_33r / LambdaNP2
                +348361.  * CHL3_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3664765848.  * CuH_33r / LambdaNP2 * CuH_33r / LambdaNP2
                +3120884271.  * CuH_33r / LambdaNP2 * CHD / LambdaNP2
                -8478772019.  * CuH_33r / LambdaNP2 * CHB / LambdaNP2
                -7829043471.  * CuH_33r / LambdaNP2 * CHW / LambdaNP2
                +15184587744.  * CuH_33r / LambdaNP2 * CHWB / LambdaNP2
                +106643753.  * CuH_33r / LambdaNP2 * CDHB / LambdaNP2
                -72003281.  * CuH_33r / LambdaNP2 * CDHW / LambdaNP2
                -6259344608.  * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
                -83225754794.  * CuH_33r / LambdaNP2 * CuW_33r / LambdaNP2
                -103535874702.  * CuH_33r / LambdaNP2 * CuB_33r / LambdaNP2
                +68393.2  * CuH_33r / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +1273952271.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -1526266453.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -5532632895.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +9041635746.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -53806195.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -28102176.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -5146193805.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                -58728156857.  * CHD / LambdaNP2 * CuW_33r / LambdaNP2
                -18146974964.  * CHD / LambdaNP2 * CuB_33r / LambdaNP2
                +122539.  * CHD / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +6790219896.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +5158379877.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -12620786627.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -285146272.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +144143265.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +3032066555.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +56093035972.  * CHB / LambdaNP2 * CuW_33r / LambdaNP2
                +163892903175.  * CHB / LambdaNP2 * CuB_33r / LambdaNP2
                +77798.1  * CHB / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +6275436472.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -21292426669.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                +52314182.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +17986174.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +11183845643.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +129255555990.  * CHW / LambdaNP2 * CuW_33r / LambdaNP2
                +64145998516.  * CHW / LambdaNP2 * CuB_33r / LambdaNP2
                -240408.  * CHW / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +19678553295.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +232371207.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -223001211.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -18211147131.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                -227701050658.  * CHWB / LambdaNP2 * CuW_33r / LambdaNP2
                -160599539117.  * CHWB / LambdaNP2 * CuB_33r / LambdaNP2
                +332731.  * CHWB / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +287489747.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -202862946.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +124180760.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +259321173.  * CDHB / LambdaNP2 * CuW_33r / LambdaNP2
                -3607413194.  * CDHB / LambdaNP2 * CuB_33r / LambdaNP2
                -7090.2  * CDHB / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +194488927.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +103913604.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                -233546850.  * CDHW / LambdaNP2 * CuW_33r / LambdaNP2
                +2227102293.  * CDHW / LambdaNP2 * CuB_33r / LambdaNP2
                +2396.  * CDHW / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +5202515330.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                +117790883881.  * CLL_1221 / LambdaNP2 * CuW_33r / LambdaNP2
                +36120767098.  * CLL_1221 / LambdaNP2 * CuB_33r / LambdaNP2
                -247816.  * CLL_1221 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +711283833926.  * CuW_33r / LambdaNP2 * CuW_33r / LambdaNP2
                +685271257275.  * CuW_33r / LambdaNP2 * CuB_33r / LambdaNP2
                -2549838.  * CuW_33r / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +1008280279655.  * CuB_33r / LambdaNP2 * CuB_33r / LambdaNP2
                +1016780.  * CuB_33r / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3.906  * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        }
        
    } else if (sqrt_s == 1.0) {

        mu += 
                +121235. * CHbox / LambdaNP2
                +880833. * CHL1_11 / LambdaNP2
                -537246. * CHe_11 / LambdaNP2
                -8082.36 * CHu_11 / LambdaNP2
                +880977. * CHL3_11 / LambdaNP2
                -115462. * CuH_33r / LambdaNP2
                -59799.9 * CHD / LambdaNP2
                +351246. * CHB / LambdaNP2
                +360483. * CHW / LambdaNP2
                -396132. * CHWB / LambdaNP2
                +36899.8 * CDHB / LambdaNP2
                +112864. * CDHW / LambdaNP2
                -125950. * DeltaGF() / v() / v()
                +2759405. * CuW_33r / LambdaNP2
                +3461224. * CuB_33r / LambdaNP2
                -2.562 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        mu +=  
                +3675406524.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +50495670163.  * CHbox / LambdaNP2 * CHL1_11 / LambdaNP2
                -35860531127.  * CHbox / LambdaNP2 * CHe_11 / LambdaNP2
                -481429809.  * CHbox / LambdaNP2 * CHu_11 / LambdaNP2
                +39430530165.  * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
                -7015779852.  * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
                -3626335033.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                +21250024055.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                +21793659194.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -23991965746.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +2184306745.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +6855720196.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                +7599008948.  * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
                +167689935534.  * CHbox / LambdaNP2 * CuW_33r / LambdaNP2
                +208967718657.  * CHbox / LambdaNP2 * CuB_33r / LambdaNP2
                -156041.  * CHbox / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3930866929664.  * CHL1_11 / LambdaNP2 * CHL1_11 / LambdaNP2
                -3800635043.  * CHL1_11 / LambdaNP2 * CHe_11 / LambdaNP2
                +5985759646.  * CHL1_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +7591311459636.  * CHL1_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -19688973347.  * CHL1_11 / LambdaNP2 * CuH_33r / LambdaNP2
                -51883142500.  * CHL1_11 / LambdaNP2 * CHD / LambdaNP2
                +4738766477.  * CHL1_11 / LambdaNP2 * CHB / LambdaNP2
                +97233715001.  * CHL1_11 / LambdaNP2 * CHW / LambdaNP2
                -177229866256.  * CHL1_11 / LambdaNP2 * CHWB / LambdaNP2
                -685268931011.  * CHL1_11 / LambdaNP2 * CDHB / LambdaNP2
                +1256788222842.  * CHL1_11 / LambdaNP2 * CDHW / LambdaNP2
                +130552294814.  * CHL1_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                +356153180025.  * CHL1_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +189213893967.  * CHL1_11 / LambdaNP2 * CuB_33r / LambdaNP2
                -1218291.  * CHL1_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3930145290099.  * CHe_11 / LambdaNP2 * CHe_11 / LambdaNP2
                -6618878091.  * CHe_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +137977484846.  * CHe_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +2441306649.  * CHe_11 / LambdaNP2 * CuH_33r / LambdaNP2
                +17442365053.  * CHe_11 / LambdaNP2 * CHD / LambdaNP2
                +10319445781.  * CHe_11 / LambdaNP2 * CHB / LambdaNP2
                -3555277591.  * CHe_11 / LambdaNP2 * CHW / LambdaNP2
                -183387857212.  * CHe_11 / LambdaNP2 * CHWB / LambdaNP2
                -1377850476282.  * CHe_11 / LambdaNP2 * CDHB / LambdaNP2
                -5407485808.  * CHe_11 / LambdaNP2 * CDHW / LambdaNP2
                -66551525065.  * CHe_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                -7312614260.  * CHe_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +387664774367.  * CHe_11 / LambdaNP2 * CuB_33r / LambdaNP2
                +906866.  * CHe_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +10996824786.  * CHu_11 / LambdaNP2 * CHu_11 / LambdaNP2
                +19312036948.  * CHu_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                +544356779.  * CHu_11 / LambdaNP2 * CuH_33r / LambdaNP2
                +2173915135.  * CHu_11 / LambdaNP2 * CHD / LambdaNP2
                +12422784567.  * CHu_11 / LambdaNP2 * CHB / LambdaNP2
                -13185798133.  * CHu_11 / LambdaNP2 * CHW / LambdaNP2
                +23536033869.  * CHu_11 / LambdaNP2 * CHWB / LambdaNP2
                -840950640.  * CHu_11 / LambdaNP2 * CDHB / LambdaNP2
                +1583758299.  * CHu_11 / LambdaNP2 * CDHW / LambdaNP2
                -4291831040.  * CHu_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                -230182815356.  * CHu_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +91475993457.  * CHu_11 / LambdaNP2 * CuB_33r / LambdaNP2
                +243084.  * CHu_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3688059270663.  * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
                -4919416915.  * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
                -34756230155.  * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
                -7962089868.  * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
                +31949389012.  * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
                -127725392091.  * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
                -683681323968.  * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
                +1223929567978.  * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
                +107507938035.  * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
                -98624073896.  * CHL3_11 / LambdaNP2 * CuW_33r / LambdaNP2
                +35119792168.  * CHL3_11 / LambdaNP2 * CuB_33r / LambdaNP2
                -162855.  * CHL3_11 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3476137785.  * CuH_33r / LambdaNP2 * CuH_33r / LambdaNP2
                +3347060522.  * CuH_33r / LambdaNP2 * CHD / LambdaNP2
                -21416097373.  * CuH_33r / LambdaNP2 * CHB / LambdaNP2
                -21517126912.  * CuH_33r / LambdaNP2 * CHW / LambdaNP2
                +24095304532.  * CuH_33r / LambdaNP2 * CHWB / LambdaNP2
                -21889733.  * CuH_33r / LambdaNP2 * CDHB / LambdaNP2
                -1181323968.  * CuH_33r / LambdaNP2 * CDHW / LambdaNP2
                -6791879149.  * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
                -165123400366.  * CuH_33r / LambdaNP2 * CuW_33r / LambdaNP2
                -210346146445.  * CuH_33r / LambdaNP2 * CuB_33r / LambdaNP2
                +148578.  * CuH_33r / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +1560810161.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -3975897239.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -17148224767.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +11933849707.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +531944578.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                -6131194073.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -6525353603.  * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
                -114377128837.  * CHD / LambdaNP2 * CuW_33r / LambdaNP2
                -40529346676.  * CHD / LambdaNP2 * CuB_33r / LambdaNP2
                +228636.  * CHD / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +66366785336.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +30890984316.  * CHB / LambdaNP2 * CHW / LambdaNP2
                -65241027615.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -2078321947.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                +1005484461.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +7883190609.  * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +339964399115.  * CHB / LambdaNP2 * CuW_33r / LambdaNP2
                +956437024921.  * CHB / LambdaNP2 * CuB_33r / LambdaNP2
                +215198.  * CHB / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +70234773405.  * CHW / LambdaNP2 * CHW / LambdaNP2
                -105796690080.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -1361493313.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                +1871451939.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +34554988935.  * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +715293947850.  * CHW / LambdaNP2 * CuW_33r / LambdaNP2
                +379106129125.  * CHW / LambdaNP2 * CuB_33r / LambdaNP2
                -1332008.  * CHW / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +78702010969.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                +34874434716.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -21124795535.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                -23936303281.  * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
                -776768016934.  * CHWB / LambdaNP2 * CuW_33r / LambdaNP2
                -782829789281.  * CHWB / LambdaNP2 * CuB_33r / LambdaNP2
                +389940.  * CHWB / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +151399980756.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                -109573751564.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                +1239295680.  * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
                +5725007216.  * CDHB / LambdaNP2 * CuW_33r / LambdaNP2
                -34975464255.  * CDHB / LambdaNP2 * CuB_33r / LambdaNP2
                -51914.8  * CDHB / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +101481766574.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                +17966900799.  * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
                +1953237756.  * CDHW / LambdaNP2 * CuW_33r / LambdaNP2
                +16386029058.  * CDHW / LambdaNP2 * CuB_33r / LambdaNP2
                -26362.1  * CDHW / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +6944096988.  * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2
                +226482247667.  * CLL_1221 / LambdaNP2 * CuW_33r / LambdaNP2
                +77247185606.  * CLL_1221 / LambdaNP2 * CuB_33r / LambdaNP2
                -464787.  * CLL_1221 / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +3677956316752.  * CuW_33r / LambdaNP2 * CuW_33r / LambdaNP2
                +3862214952372.  * CuW_33r / LambdaNP2 * CuB_33r / LambdaNP2
                -7581065.  * CuW_33r / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +5414702203406.  * CuB_33r / LambdaNP2 * CuB_33r / LambdaNP2
                +3001867.  * CuB_33r / LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                +11.97  * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2 * 0.5 * (CHQ1_33 - CHQ3_33) * v2_over_LambdaNP2
                ;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPSMEFTd6::mueettH()");
      
    //Add intrinsic and parametric relative theory errors (free par). (Assume they are constant in energy.)
    mu += eeettHint + eeettHpar;

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPSMEFTd6::BrHggRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHggRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHggRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHggRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHWWRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHWWRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHWWRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHWWRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHZZRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZZRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZZRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZZRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHZgaRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZgaRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHZgaRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZgaRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHgagaRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHgagaRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHgagaRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHgagaRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}
      
double NPSMEFTd6::BrHmumuRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHmumuRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHmumuRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHmumuRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHtautauRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHtautauRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHtautauRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHtautauRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHccRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHccRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHccRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHccRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::BrHbbRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHbbRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        Br += - deltaGammaHbbRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHbbRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPSMEFTd6::computeGammaTotalRatio() const
{
    double width = 1.0;

    width += deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaTotalRatio2();
        }
    
    if (width < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return width;

}

double NPSMEFTd6::deltaGammaTotalRatio1() const
{
    return (trueSM.computeBrHtogg() * deltaGammaHggRatio1()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio1()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio1()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio1()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio1()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio1()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio1()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio1()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio1()
            + BrHinv + BrHexo);
}

double NPSMEFTd6::deltaGammaTotalRatio2() const
{
    double delta2SM;
    
    delta2SM = trueSM.computeBrHtogg() * deltaGammaHggRatio2()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio2()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio2()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio2()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio2()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio2()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio2()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio2()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio2();

    return (delta2SM + (BrHinv + BrHexo)*(BrHinv + BrHexo + delta2SM));
}

double NPSMEFTd6::GammaHggRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHggint + eHggpar;

    width += deltaGammaHggRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHggRatio2();            
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHggRatio1() const
{
    return ( +121249. * CHbox / LambdaNP2
            +173400. * CuH_22r / LambdaNP2
            -128860. * CuH_33r / LambdaNP2
            +248587. * CdH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            +37390592. * CHG / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
}

double NPSMEFTd6::deltaGammaHggRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3675338322. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -7350676644. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            +10512328996. * CHbox / LambdaNP2 * CuH_22r / LambdaNP2
            -7812106277. * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
            +15070455306. * CHbox / LambdaNP2 * CdH_33r / LambdaNP2
            -1837669161. * CHbox / LambdaNP2 * CHD / LambdaNP2
            +2266788502923. * CHbox / LambdaNP2 * CHG / LambdaNP2
            +3675338322. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +3675338322. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            -10512328996. * CHL3_11 / LambdaNP2 * CuH_22r / LambdaNP2
            +7812106277. * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
            -15070455306. * CHL3_11 / LambdaNP2 * CdH_33r / LambdaNP2
            +1837669161. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            -2266788502922. * CHL3_11 / LambdaNP2 * CHG / LambdaNP2
            -3675338322. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +14576143481. * CuH_22r / LambdaNP2 * CuH_22r / LambdaNP2
            -12069614762. * CuH_22r / LambdaNP2 * CuH_33r / LambdaNP2
            +54515542771. * CuH_22r / LambdaNP2 * CdH_33r / LambdaNP2
            -2628082249. * CuH_22r / LambdaNP2 * CHD / LambdaNP2
            +3502162286037. * CuH_22r / LambdaNP2 * CHG / LambdaNP2
            +5256164498. * CuH_22r / LambdaNP2 * CLL_1221 / LambdaNP2
            +4179770404. * CuH_33r / LambdaNP2 * CuH_33r / LambdaNP2
            -18111661684. * CuH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +1953026569. * CuH_33r / LambdaNP2 * CHD / LambdaNP2
            -2425634050957. * CuH_33r / LambdaNP2 * CHG / LambdaNP2
            -3906053138. * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +53929075760. * CdH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            -3767613827. * CdH_33r / LambdaNP2 * CHD / LambdaNP2
            +5255344080213. * CdH_33r / LambdaNP2 * CHG / LambdaNP2
            +7535227653. * CdH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +229708645. * CHD / LambdaNP2 * CHD / LambdaNP2
            -566697125731. * CHD / LambdaNP2 * CHG / LambdaNP2
            -918834581. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +228376225497250. * CHG / LambdaNP2 * CHG / LambdaNP2
            +1133394251461. * CHG / LambdaNP2 * CLL_1221 / LambdaNP2
            +918834581. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );
    
}

double NPSMEFTd6::GammaHWWRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHWWint + eHWWpar;

    width += deltaGammaHWWRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHWWRatio2();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHWWRatio1() const
{

    return ( +121100. * CHbox / LambdaNP2
                -117075. * CHD / LambdaNP2
                -91657.4 * CHW / LambdaNP2
                -190618. * CHWB / LambdaNP2
                +38345.4 * CDHW / LambdaNP2
                -1.849 * DeltaGF() );
    
}

double NPSMEFTd6::deltaGammaHWWRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3675070317.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                -7091873661.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -5530926638.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -11515676206.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +2241942559.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -112676.  * CHbox / LambdaNP2 * DeltaGF()
                +3421412536.  * CHD / LambdaNP2 * CHD / LambdaNP2
                +5327587542.  * CHD / LambdaNP2 * CHW / LambdaNP2
                +11099269734.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                -2166464532.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                +108782.  * CHD / LambdaNP2 * DeltaGF()
                +2929837690.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +8664103424.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -2001161200.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +84850.3  * CHW / LambdaNP2 * DeltaGF()
                +9020462932.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -3515805228.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +176601.  * CHWB / LambdaNP2 * DeltaGF()
                +372822749.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -34333.3  * CDHW / LambdaNP2 * DeltaGF()
                +0.865  * DeltaGF() * DeltaGF() );
    
}

double NPSMEFTd6::GammaHZZRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHZZint + eHZZpar;

    width += deltaGammaHZZRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZZRatio2();            
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZZRatio1() const
{

    return ( +120032. * CHbox / LambdaNP2
                +29732.2 * CHD / LambdaNP2
                -13963.3 * CHB / LambdaNP2
                -46223.3 * CHW / LambdaNP2
                -25579.3 * CHWB / LambdaNP2
                +15260.6 * CDHB / LambdaNP2
                +28274.3 * CDHW / LambdaNP2
                -0.992 * DeltaGF() );
    
}

double NPSMEFTd6::deltaGammaHZZRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3675794964.  * CHbox / LambdaNP2 * CHbox / LambdaNP2
                +1836272445.  * CHbox / LambdaNP2 * CHD / LambdaNP2
                -837280012.  * CHbox / LambdaNP2 * CHB / LambdaNP2
                -2787840082.  * CHbox / LambdaNP2 * CHW / LambdaNP2
                -1528719029.  * CHbox / LambdaNP2 * CHWB / LambdaNP2
                +948225722.  * CHbox / LambdaNP2 * CDHB / LambdaNP2
                +1730083258.  * CHbox / LambdaNP2 * CDHW / LambdaNP2
                -60644.3  * CHbox / LambdaNP2 * DeltaGF()
                +229531548.  * CHD / LambdaNP2 * CHD / LambdaNP2
                -209044482.  * CHD / LambdaNP2 * CHB / LambdaNP2
                -696216048.  * CHD / LambdaNP2 * CHW / LambdaNP2
                -381621025.  * CHD / LambdaNP2 * CHWB / LambdaNP2
                +237582729.  * CHD / LambdaNP2 * CDHB / LambdaNP2
                +433363427.  * CHD / LambdaNP2 * CDHW / LambdaNP2
                -15151.7  * CHD / LambdaNP2 * DeltaGF()
                +67092431.  * CHB / LambdaNP2 * CHB / LambdaNP2
                +447244235.  * CHB / LambdaNP2 * CHW / LambdaNP2
                +245015102.  * CHB / LambdaNP2 * CHWB / LambdaNP2
                -124762787.  * CHB / LambdaNP2 * CDHB / LambdaNP2
                -227788589.  * CHB / LambdaNP2 * CDHW / LambdaNP2
                +6903.17  * CHB / LambdaNP2 * DeltaGF()
                +744004191.  * CHW / LambdaNP2 * CHW / LambdaNP2
                +815851361.  * CHW / LambdaNP2 * CHWB / LambdaNP2
                -415479542.  * CHW / LambdaNP2 * CDHB / LambdaNP2
                -758148037.  * CHW / LambdaNP2 * CDHW / LambdaNP2
                +22981.7  * CHW / LambdaNP2 * DeltaGF()
                +223472766.  * CHWB / LambdaNP2 * CHWB / LambdaNP2
                -227707672.  * CHWB / LambdaNP2 * CDHB / LambdaNP2
                -415590330.  * CHWB / LambdaNP2 * CDHW / LambdaNP2
                +12595.7  * CHWB / LambdaNP2 * DeltaGF()
                +65918124.  * CDHB / LambdaNP2 * CDHB / LambdaNP2
                +240411498.  * CDHB / LambdaNP2 * CDHW / LambdaNP2
                -7829.46  * CDHB / LambdaNP2 * DeltaGF()
                +219560638.  * CDHW / LambdaNP2 * CDHW / LambdaNP2
                -14272.5  * CDHW / LambdaNP2 * DeltaGF()
                +0.25  * DeltaGF() * DeltaGF() );
    
}

double NPSMEFTd6::GammaHZgaRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHZgaint + eHZgapar;

    width += deltaGammaHZgaRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHZgaRatio2();            
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHZgaRatio1() const
{

    return ( +119538. * CHbox / LambdaNP2
            -321.71 * CeH_33r / LambdaNP2
            -2910.68 * CuH_22r / LambdaNP2
            +6522.67 * CuH_33r / LambdaNP2
            -6118.03 * CdH_33r / LambdaNP2
            -120435. * CHD / LambdaNP2
            +14997549. * CHB / LambdaNP2
            -14997549. * CHW / LambdaNP2
            +9237644. * CHWB / LambdaNP2
            -2872389. * CDHB / LambdaNP2
            +1586252. * CDHW / LambdaNP2
            -115000. * DeltaGF() / v() / v() );
    
}

double NPSMEFTd6::deltaGammaHZgaRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3572354374. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -13746871780. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            -19228428. * CHbox / LambdaNP2 * CeH_33r / LambdaNP2
            -173969500. * CHbox / LambdaNP2 * CuH_22r / LambdaNP2
            +389854519. * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
            -365671769. * CHbox / LambdaNP2 * CdH_33r / LambdaNP2
            -7198310278. * CHbox / LambdaNP2 * CHD / LambdaNP2
            +896391115001. * CHbox / LambdaNP2 * CHB / LambdaNP2
            -896391115001. * CHbox / LambdaNP2 * CHW / LambdaNP2
            +552126360087. * CHbox / LambdaNP2 * CHWB / LambdaNP2
            -171680296875. * CHbox / LambdaNP2 * CDHB / LambdaNP2
            +94808952497. * CHbox / LambdaNP2 * CDHW / LambdaNP2
            +6873435890. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +13224927460. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            +36990810. * CHL3_11 / LambdaNP2 * CeH_33r / LambdaNP2
            +334679379. * CHL3_11 / LambdaNP2 * CuH_22r / LambdaNP2
            -750104982. * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
            +703428859. * CHL3_11 / LambdaNP2 * CdH_33r / LambdaNP2
            +13850007508. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            -1724713728406. * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
            +1724713728406. * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
            -1062326363035. * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
            +330323850787. * CHL3_11 / LambdaNP2 * CDHB / LambdaNP2
            -182418476948. * CHL3_11 / LambdaNP2 * CDHW / LambdaNP2
            -13224927460. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +30650.9 * CeH_33r / LambdaNP2 * CeH_33r / LambdaNP2
            +547410. * CeH_33r / LambdaNP2 * CuH_22r / LambdaNP2
            -1048860. * CeH_33r / LambdaNP2 * CuH_33r / LambdaNP2
            +1221618. * CeH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +19367846. * CeH_33r / LambdaNP2 * CHD / LambdaNP2
            -2411640054. * CeH_33r / LambdaNP2 * CHB / LambdaNP2
            +2411640054. * CeH_33r / LambdaNP2 * CHW / LambdaNP2
            -1485434229. * CeH_33r / LambdaNP2 * CHWB / LambdaNP2
            +461886640. * CeH_33r / LambdaNP2 * CDHB / LambdaNP2
            -255072885. * CeH_33r / LambdaNP2 * CDHW / LambdaNP2
            -18495405. * CeH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +2446415. * CuH_22r / LambdaNP2 * CuH_22r / LambdaNP2
            -9489840. * CuH_22r / LambdaNP2 * CuH_33r / LambdaNP2
            +10873108. * CuH_22r / LambdaNP2 * CdH_33r / LambdaNP2
            +175234540. * CuH_22r / LambdaNP2 * CHD / LambdaNP2
            -21819956015. * CuH_22r / LambdaNP2 * CHB / LambdaNP2
            +21819956015. * CuH_22r / LambdaNP2 * CHW / LambdaNP2
            -13439862009. * CuH_22r / LambdaNP2 * CHWB / LambdaNP2
            +4179042456. * CuH_22r / LambdaNP2 * CDHB / LambdaNP2
            -2307839891. * CuH_22r / LambdaNP2 * CDHW / LambdaNP2
            -167339689. * CuH_22r / LambdaNP2 * CLL_1221 / LambdaNP2
            +10636307. * CuH_33r / LambdaNP2 * CuH_33r / LambdaNP2
            -19944392. * CuH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            -392779489. * CuH_33r / LambdaNP2 * CHD / LambdaNP2
            +48912047811. * CuH_33r / LambdaNP2 * CHB / LambdaNP2
            -48912047811. * CuH_33r / LambdaNP2 * CHW / LambdaNP2
            +30127062251. * CuH_33r / LambdaNP2 * CHWB / LambdaNP2
            -9367824768. * CuH_33r / LambdaNP2 * CDHB / LambdaNP2
            +5173299846. * CuH_33r / LambdaNP2 * CDHW / LambdaNP2
            +375052491. * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +12309844. * CdH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +368294826. * CdH_33r / LambdaNP2 * CHD / LambdaNP2
            -45858069994. * CdH_33r / LambdaNP2 * CHB / LambdaNP2
            +45858069994. * CdH_33r / LambdaNP2 * CHW / LambdaNP2
            -28245984194. * CdH_33r / LambdaNP2 * CHWB / LambdaNP2
            +8782915113. * CdH_33r / LambdaNP2 * CDHB / LambdaNP2
            -4850288570. * CdH_33r / LambdaNP2 * CDHW / LambdaNP2
            -351714430. * CdH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +3626158197. * CHD / LambdaNP2 * CHD / LambdaNP2
            -903116488795. * CHD / LambdaNP2 * CHB / LambdaNP2
            +903116488795. * CHD / LambdaNP2 * CHW / LambdaNP2
            -556268810955. * CHD / LambdaNP2 * CHWB / LambdaNP2
            +172968366503. * CHD / LambdaNP2 * CDHB / LambdaNP2
            -95520277759. * CHD / LambdaNP2 * CDHW / LambdaNP2
            -6925003754. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +56231649181760. * CHB / LambdaNP2 * CHB / LambdaNP2
            -112463298363521. * CHB / LambdaNP2 * CHW / LambdaNP2
            +69271047570247. * CHB / LambdaNP2 * CHWB / LambdaNP2
            -21539406323346. * CHB / LambdaNP2 * CDHB / LambdaNP2
            +11894950021056. * CHB / LambdaNP2 * CDHW / LambdaNP2
            +862356864203. * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
            +56231649181760. * CHW / LambdaNP2 * CHW / LambdaNP2
            -69271047570247. * CHW / LambdaNP2 * CHWB / LambdaNP2
            +21539406323346. * CHW / LambdaNP2 * CDHB / LambdaNP2
            -11894950021056. * CHW / LambdaNP2 * CDHW / LambdaNP2
            -862356864203. * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
            +21333528810302. * CHWB / LambdaNP2 * CHWB / LambdaNP2
            -13267059225282. * CHWB / LambdaNP2 * CDHB / LambdaNP2
            +7326618201175. * CHWB / LambdaNP2 * CDHW / LambdaNP2
            +531163181518. * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
            +2062655246259. * CDHB / LambdaNP2 * CDHB / LambdaNP2
            -2278166881352. * CDHB / LambdaNP2 * CDHW / LambdaNP2
            -165161925393. * CDHB / LambdaNP2 * CLL_1221 / LambdaNP2
            +629048934463. * CDHW / LambdaNP2 * CDHW / LambdaNP2
            +91209238474. * CDHW / LambdaNP2 * CLL_1221 / LambdaNP2
            +3306231865. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );
    
}

double NPSMEFTd6::GammaHgagaRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHgagaint + eHgagapar;

    width += deltaGammaHgagaRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHgagaRatio2();            
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHgagaRatio1() const
{
    return ( +119212. * CHbox / LambdaNP2
            -42570.4 * CeH_33r / LambdaNP2
            -48874.1 * CuH_22r / LambdaNP2
            +31992.9 * CuH_33r / LambdaNP2
            -18430.1 * CdH_33r / LambdaNP2
            -137624. * CHD / LambdaNP2
            -48142626. * CHB / LambdaNP2
            -14682079. * CHW / LambdaNP2
            +26348174. * CHWB / LambdaNP2
            -125370. * DeltaGF() / v() / v() );
    
}

double NPSMEFTd6::deltaGammaHgagaRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3552878939. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -14945652252. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            -2537854461. * CHbox / LambdaNP2 * CeH_33r / LambdaNP2
            -2913607301. * CHbox / LambdaNP2 * CuH_22r / LambdaNP2
            +1906965841. * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
            -1098797271. * CHbox / LambdaNP2 * CdH_33r / LambdaNP2
            -8203204692. * CHbox / LambdaNP2 * CHD / LambdaNP2
            -2869585490784. * CHbox / LambdaNP2 * CHB / LambdaNP2
            -875138828266. * CHbox / LambdaNP2 * CHW / LambdaNP2
            +1570507157057. * CHbox / LambdaNP2 * CHWB / LambdaNP2
            +7472826126. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +15718197048. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            +5311406908. * CHL3_11 / LambdaNP2 * CeH_33r / LambdaNP2
            +6100926530. * CHL3_11 / LambdaNP2 * CuH_22r / LambdaNP2
            -4011186872. * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
            +2294519933. * CHL3_11 / LambdaNP2 * CdH_33r / LambdaNP2
            +17254723530. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            +6035998863335. * CHL3_11 / LambdaNP2 * CHB / LambdaNP2
            +1840801394360. * CHL3_11 / LambdaNP2 * CHW / LambdaNP2
            -3303466457194. * CHL3_11 / LambdaNP2 * CHWB / LambdaNP2
            -15718197048. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +815713130. * CeH_33r / LambdaNP2 * CeH_33r / LambdaNP2
            +1787711513. * CeH_33r / LambdaNP2 * CuH_22r / LambdaNP2
            -674635336. * CeH_33r / LambdaNP2 * CuH_33r / LambdaNP2
            +846500334. * CeH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +2908087411. * CeH_33r / LambdaNP2 * CHD / LambdaNP2
            +1015185343169. * CeH_33r / LambdaNP2 * CHB / LambdaNP2
            +309601548568. * CeH_33r / LambdaNP2 * CHW / LambdaNP2
            -555604930506. * CeH_33r / LambdaNP2 * CHWB / LambdaNP2
            -2655703454. * CeH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +982269617. * CuH_22r / LambdaNP2 * CuH_22r / LambdaNP2
            -775279437. * CuH_22r / LambdaNP2 * CuH_33r / LambdaNP2
            +918435000. * CuH_22r / LambdaNP2 * CdH_33r / LambdaNP2
            +3341211565. * CuH_22r / LambdaNP2 * CHD / LambdaNP2
            +1166633704511. * CuH_22r / LambdaNP2 * CHB / LambdaNP2
            +355788826107. * CuH_22r / LambdaNP2 * CHW / LambdaNP2
            -638491722405. * CuH_22r / LambdaNP2 * CHWB / LambdaNP2
            -3050463265. * CuH_22r / LambdaNP2 * CLL_1221 / LambdaNP2
            +255914011. * CuH_33r / LambdaNP2 * CuH_33r / LambdaNP2
            -290846045. * CuH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            -2201679820. * CuH_33r / LambdaNP2 * CHD / LambdaNP2
            -770194323138. * CuH_33r / LambdaNP2 * CHB / LambdaNP2
            -234886522689. * CuH_33r / LambdaNP2 * CHW / LambdaNP2
            +421522795086. * CuH_33r / LambdaNP2 * CHWB / LambdaNP2
            +2005593436. * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +227138666. * CdH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +1254894771. * CdH_33r / LambdaNP2 * CHD / LambdaNP2
            +437662580431. * CdH_33r / LambdaNP2 * CHB / LambdaNP2
            +133474161702. * CdH_33r / LambdaNP2 * CHW / LambdaNP2
            -239530140207. * CdH_33r / LambdaNP2 * CHWB / LambdaNP2
            -1147259967. * CdH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +4735398687. * CHD / LambdaNP2 * CHD / LambdaNP2
            +3313068504426. * CHD / LambdaNP2 * CHB / LambdaNP2
            +1010388050204. * CHD / LambdaNP2 * CHW / LambdaNP2
            -1813222785915. * CHD / LambdaNP2 * CHWB / LambdaNP2
            -8627361765. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +579490834877631. * CHB / LambdaNP2 * CHB / LambdaNP2
            +353455181491632. * CHB / LambdaNP2 * CHW / LambdaNP2
            -634303809067317. * CHB / LambdaNP2 * CHWB / LambdaNP2
            -3017999431667. * CHB / LambdaNP2 * CLL_1221 / LambdaNP2
            +53896695945874. * CHW / LambdaNP2 * CHW / LambdaNP2
            -193443929102058. * CHW / LambdaNP2 * CHWB / LambdaNP2
            -920400697180. * CHW / LambdaNP2 * CLL_1221 / LambdaNP2
            +173575360463754. * CHWB / LambdaNP2 * CHWB / LambdaNP2
            +1651733228597. * CHWB / LambdaNP2 * CLL_1221 / LambdaNP2
            +3929549262. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );
    
}
      
double NPSMEFTd6::GammaHmumuRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHmumuint + eHmumupar;

    width += deltaGammaHmumuRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHmumuRatio2();
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHmumuRatio1() const
{
    return ( +121249. * CHbox / LambdaNP2
            -199794752. * CeH_22r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
        
}

double NPSMEFTd6::deltaGammaHmumuRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return (  +3675338322. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -7350676644. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            -12112470680846. * CHbox / LambdaNP2 * CeH_22r / LambdaNP2
            -1837669161. * CHbox / LambdaNP2 * CHD / LambdaNP2
            +3675338322. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +3675338322. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            +12112470680846. * CHL3_11 / LambdaNP2 * CeH_22r / LambdaNP2
            +1837669161. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            -3675338322. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +9979485773431880. * CeH_22r / LambdaNP2 * CeH_22r / LambdaNP2
            +3028117670211. * CeH_22r / LambdaNP2 * CHD / LambdaNP2
            -6056235340423. * CeH_22r / LambdaNP2 * CLL_1221 / LambdaNP2
            +229708645. * CHD / LambdaNP2 * CHD / LambdaNP2
            -918834581. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +918834581. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );
    
}

double NPSMEFTd6::GammaHtautauRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHtautauint + eHtautaupar;

    width += deltaGammaHtautauRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHtautauRatio2();            
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHtautauRatio1() const
{
    return ( +121249. * CHbox / LambdaNP2
            -11880769. * CeH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
        
}

double NPSMEFTd6::deltaGammaHtautauRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3675338322. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -7350676644. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            -720266502504. * CHbox / LambdaNP2 * CeH_33r / LambdaNP2
            -1837669161. * CHbox / LambdaNP2 * CHD / LambdaNP2
            +3675338322. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +3675338322. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            +720266502504. * CHL3_11 / LambdaNP2 * CeH_33r / LambdaNP2
            +1837669161. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            -3675338322. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +35288168676963. * CeH_33r / LambdaNP2 * CeH_33r / LambdaNP2
            +180066625626. * CeH_33r / LambdaNP2 * CHD / LambdaNP2
            -360133251252. * CeH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +229708645. * CHD / LambdaNP2 * CHD / LambdaNP2
            -918834581. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +918834581. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );
    
}

double NPSMEFTd6::GammaHccRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHccint + eHccpar;

    width += deltaGammaHccRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHccRatio2();            
        }
    
    return width;
    
}

double NPSMEFTd6::deltaGammaHccRatio1() const
{
    return ( +121249. * CHbox / LambdaNP2
            -16420490. * CuH_22r / LambdaNP2
            -1001.52 * CuH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
}

double NPSMEFTd6::deltaGammaHccRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3675338322. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -7350676644. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            -995485129059. * CHbox / LambdaNP2 * CuH_22r / LambdaNP2
            -60716517. * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
            -1837669161. * CHbox / LambdaNP2 * CHD / LambdaNP2
            +3675338322. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +3675338322. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            +995485129059. * CHL3_11 / LambdaNP2 * CuH_22r / LambdaNP2
            +60716517. * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
            +1837669161. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            -3675338322. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +67403475035188. * CuH_22r / LambdaNP2 * CuH_22r / LambdaNP2
            +8290982561. * CuH_22r / LambdaNP2 * CuH_33r / LambdaNP2
            +248871282265. * CuH_22r / LambdaNP2 * CHD / LambdaNP2
            -497742564529. * CuH_22r / LambdaNP2 * CLL_1221 / LambdaNP2
            +15179129. * CuH_33r / LambdaNP2 * CHD / LambdaNP2
            -30358259. * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +229708645. * CHD / LambdaNP2 * CHD / LambdaNP2
            -918834581. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +918834581. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );
    
}

double NPSMEFTd6::GammaHbbRatio() const
{
      // SM (1) + intrinsic + parametric theory relative errors (free pars)
    double width = 1.0 + eHbbint + eHbbpar;
    
    width += deltaGammaHbbRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            width += deltaGammaHbbRatio2();        
        }
    
    return width;
}

double NPSMEFTd6::deltaGammaHbbRatio1() const
{    
    return ( +121249. * CHbox / LambdaNP2
            -562.029 * CuH_33r / LambdaNP2
            -5026895. * CdH_33r / LambdaNP2
            -30312.3 * CHD / LambdaNP2
            -60624.6 * DeltaGF() / v() / v() );
}

double NPSMEFTd6::deltaGammaHbbRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    return ( +3675338322. * CHbox / LambdaNP2 * CHbox / LambdaNP2
            -7350676644. * CHbox / LambdaNP2 * CHL3_11 / LambdaNP2
            -34072751. * CHbox / LambdaNP2 * CuH_33r / LambdaNP2
            -304753324752. * CHbox / LambdaNP2 * CdH_33r / LambdaNP2
            -1837669161. * CHbox / LambdaNP2 * CHD / LambdaNP2
            +3675338322. * CHbox / LambdaNP2 * CLL_1221 / LambdaNP2
            +3675338322. * CHL3_11 / LambdaNP2 * CHL3_11 / LambdaNP2
            +34072751. * CHL3_11 / LambdaNP2 * CuH_33r / LambdaNP2
            +304753324752. * CHL3_11 / LambdaNP2 * CdH_33r / LambdaNP2
            +1837669161. * CHL3_11 / LambdaNP2 * CHD / LambdaNP2
            -3675338322. * CHL3_11 / LambdaNP2 * CLL_1221 / LambdaNP2
            +1419189033. * CuH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +8518188. * CuH_33r / LambdaNP2 * CHD / LambdaNP2
            -17036375. * CuH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +6317281455355. * CdH_33r / LambdaNP2 * CdH_33r / LambdaNP2
            +76188331188. * CdH_33r / LambdaNP2 * CHD / LambdaNP2
            -152376662376. * CdH_33r / LambdaNP2 * CLL_1221 / LambdaNP2
            +229708645. * CHD / LambdaNP2 * CHD / LambdaNP2
            -918834581. * CHD / LambdaNP2 * CLL_1221 / LambdaNP2
            +918834581. * CLL_1221 / LambdaNP2 * CLL_1221 / LambdaNP2 );        
    
}

double NPSMEFTd6::Br_H_exo() const
{    
    return BrHexo;
}

double NPSMEFTd6::Br_H_inv() const
{    
    return BrHinv;
}


///////////////////////SPECIAL OBSERVABLES/////////////////////////

double NPSMEFTd6::muttHZbbboost(const double sqrt_s) const
{
    /* Ratios of BR with the SM*/
    double BrHbbrat = BrHbbRatio();
    double BrZbbSM = (trueSM.GammaZ(quarks[BOTTOM]))/trueSM.Gamma_Z();
    double BrZbbrat = BR_Zf(quarks[BOTTOM])/BrZbbSM;

    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-mtpole / v());    
    double dkt = dKappa_t.real();
    
    double dgV = deltaGV_f(quarks[TOP]);
    double dgA = deltaGA_f(quarks[TOP]);
    double gLSM = quarks[TOP].getIsospin() 
    - (quarks[TOP].getCharge())*sW2_tree;
    double gRSM = - (quarks[TOP].getCharge())*sW2_tree;
    
    double dgL = 0.5*(dgV + dgA)/gLSM;
    double dgR = 0.5*(dgV - dgA)/gRSM;

    /* VERY CRUDE APPROX. */
    double dsigmarat;
    
    dsigmarat = 1.0 + 
            2.0 * dkt -
            2.0 * (gLSM*gLSM*dgL + gRSM*gRSM*dgR)/(gLSM*gLSM + gRSM*gRSM);
    
    return dsigmarat * (BrHbbrat / BrZbbrat);
    
}

double NPSMEFTd6::muttHbb(const double sqrt_s) const
{   
    return muttH(sqrt_s) * BrHbbRatio();
    
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::deltag1ZNP() const
{
      double NPdirect, NPindirect;
      
      /*    From own calculations. Agrees with with LHCHXWG-INT-2015-001 for common interactions */
      NPdirect = sW_tree / sqrt( 4.0 * M_PI * ale );
      NPdirect = - NPdirect * (Mz * Mz / v () / v() ) * CDHW * v2_over_LambdaNP2;
      
      NPindirect = - 1.0 / (cW2_tree-sW2_tree);
      
      NPindirect = NPindirect * (sW_tree * CHWB / cW_tree 
              + 0.25 * CHD ) * v2_over_LambdaNP2
              + 0.5 * NPindirect * DeltaGF() ;
      
      return NPdirect + NPindirect;
}
      
double NPSMEFTd6::deltaKgammaNP() const
{
      double NPdirect;

      /*    Translate from LHCHXWG-INT-2015-001: Checked with own calculations  */
      NPdirect = sqrt( 4.0 * M_PI * ale ) / 4.0 / sW2_tree;
      
      NPdirect = NPdirect * ( (4.0 * sW_tree * cW_tree / sqrt( 4.0 * M_PI * ale ) ) * CHWB 
              - sW_tree * CDHW 
              - cW_tree * CDHB ) * v2_over_LambdaNP2;
      
      return NPdirect;
}
      
double NPSMEFTd6::lambdaZNP() const
{
      double NPdirect;

      /*    Translate from LHCHXWG-INT-2015-001: Checked with own calculations  */
      NPdirect = - (3.0 / 2.0) * (sqrt( 4.0 * M_PI * ale ) / sW_tree) * CW * v2_over_LambdaNP2;

      return NPdirect;
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::dxseeWWdcos(const double sqrt_s, const double cos) const
{
    double s = sqrt_s * sqrt_s;
    double cos2 = cos * cos;
    double sin2 = 1.0 - cos2;
    double sin = sqrt(sin2);
    
    double topb = 0.3894*1000000000.0;
    
//  NC and CC couplings
    double gLe, gRe;
    gslpp::complex Uenu;
    
    gLe = -0.5 + sW2_tree + deltaGL_f(leptons[ELECTRON]);
    gRe = sW2_tree + deltaGR_f(leptons[ELECTRON]);
    
    Uenu = deltaGL_Wff(leptons[NEUTRINO_1], leptons[ELECTRON]);
    Uenu = 1.0 + Uenu;

//  Wigner functions
    double d1pp[2],d1mm[2],d1p0[2],d1m0[2],d10p[2],d10m[2],d100[2];
    
    d1pp[0]=sqrt((1.0 - cos2)/2.0);
    d1pp[1]=-sqrt((1.0 - cos2)/2.0);
    
    d1mm[0]=d1pp[0];
    d1mm[1]=d1pp[1];
    
    d1p0[0]=(1.0 - cos)/2.0;
    d1p0[1]=(1.0 + cos)/2.0;
    
    d1m0[0]=d1p0[1];
    d1m0[1]=d1p0[0];
    
    d10p[0]=d1p0[1];
    d10p[1]=d1p0[0];
    
    d10m[0]=d1p0[0];
    d10m[1]=d1p0[1];
    
    d100[0]=d1pp[0];
    d100[1]=d1pp[1];
    
    gslpp::matrix<double> d1LH(3, 3, 0.0);
    
    gslpp::matrix<double> d1RH(3, 3, 0.0);
    
    d1LH.assign(0,0, d1pp[0]);
    d1LH.assign(0,1, d1p0[0]);
    d1LH.assign(0,2, 0.0);
    
    d1LH.assign(1,0, d10p[0]);
    d1LH.assign(1,1, d100[0]);
    d1LH.assign(1,2, d10m[0]);
    
    d1LH.assign(2,0, 0.0);
    d1LH.assign(2,1, d1m0[0]);
    d1LH.assign(2,2, d1mm[0]);
    
    d1RH.assign(0,0, d1pp[1]);
    d1RH.assign(0,1, d1p0[1]);
    d1RH.assign(0,2, 0.0);
    
    d1RH.assign(1,0, d10p[1]);
    d1RH.assign(1,1, d100[1]);
    d1RH.assign(1,2, d10m[1]);
    
    d1RH.assign(2,0, 0.0);
    d1RH.assign(2,1, d1m0[1]);
    d1RH.assign(2,2, d1mm[1]);
    
//  TGC parameterization
    double g1Z,g1ga,kZ,kga,lambZ,lambga,g4Z,g4ga,g5Z,g5ga,ktZ,ktga,lambtZ,lambtga;

//  TGC present in the SM     
    g1Z=1.0 + deltag1ZNP();
    g1ga=1.0;
    kZ=1.0 + deltag1ZNP() - (sW2_tree/cW2_tree) * deltaKgammaNP();
    kga=1.0 + deltaKgammaNP();
//  TGC not present in the SM
    lambZ=lambdaZNP(); //Check normalization
    lambga=lambZ;
    g4Z=0.0;
    g4ga=0.0;
    g5Z=0.0;
    g5ga=0.0;
    ktZ=0.0;
    ktga=0.0;
    lambtZ=0.0;
    lambtga=0.0;
    
    double f3Z, f3ga;
    
    f3Z = g1Z + kZ + lambZ;    
    f3ga = g1ga + kga + lambga;
    
 // Kinematic factors
    double beta, gamma, gamma2;
    
    beta = sqrt(1.0 - 4.0 * Mw_tree() * Mw_tree() / s);
    gamma = sqrt_s/(2.0*Mw_tree());
    gamma2= gamma*gamma;
    
//  J=1 Subamplitudes: Z
    gslpp::complex AZpp, AZmm, AZp0, AZm0, AZ0p, AZ0m, AZ00;
    
    AZpp = gslpp::complex(g1Z + 2.0* gamma2* lambZ, (ktZ + lambtZ - 2.0*lambtZ)/beta , false);
    AZmm = gslpp::complex(g1Z + 2.0* gamma2* lambZ, -(ktZ + lambtZ - 2.0*lambtZ)/beta , false);
    AZp0 = gslpp::complex(f3Z + beta * g5Z , -g4Z + (ktZ-lambtZ)/beta , false);
    AZp0 = gamma * AZp0;
    AZm0 = gslpp::complex(f3Z - beta * g5Z , -g4Z - (ktZ-lambtZ)/beta , false);
    AZm0 = gamma * AZm0;
    AZ0p = gslpp::complex(f3Z - beta * g5Z , g4Z + (ktZ-lambtZ)/beta , false);
    AZ0p = gamma * AZ0p;
    AZ0m = gslpp::complex(f3Z + beta * g5Z , g4Z - (ktZ-lambtZ)/beta , false);
    AZ0m = gamma * AZ0m;
    AZ00 = gslpp::complex( g1Z + 2.0*gamma2*kZ, 0.0 , false);
    
//  Collect in matrices and separate LH and RH
    gslpp::matrix<gslpp::complex> AmpZLH(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> AmpZRH(3, 3, 0.0);
    
    AmpZLH.assign(0,0, AZpp *  d1LH(0,0) );
    AmpZLH.assign(0,1, AZp0 *  d1LH(0,1));
    AmpZLH.assign(0,2, 0.0);
    
    AmpZLH.assign(1,0, AZ0p *  d1LH(1,0));
    AmpZLH.assign(1,1, AZ00 *  d1LH(1,1));
    AmpZLH.assign(1,2, AZ0m *  d1LH(1,2));
    
    AmpZLH.assign(2,0, 0.0);
    AmpZLH.assign(2,1, AZm0 *  d1LH(2,1));
    AmpZLH.assign(2,2, AZmm *  d1LH(2,2));
    
    AmpZLH = AmpZLH * beta * s/(s-Mz*Mz);
    
//  Add the correct Zff coupling
    AmpZLH = AmpZLH * gLe / sW2_tree;
    
    AmpZRH.assign(0,0, AZpp *  d1RH(0,0) );
    AmpZRH.assign(0,1, AZp0 *  d1RH(0,1));
    AmpZRH.assign(0,2, 0.0);
    
    AmpZRH.assign(1,0, AZ0p *  d1RH(1,0));
    AmpZRH.assign(1,1, AZ00 *  d1RH(1,1));
    AmpZRH.assign(1,2, AZ0m *  d1RH(1,2));
    
    AmpZRH.assign(2,0, 0.0);
    AmpZRH.assign(2,1, AZm0 *  d1RH(2,1));
    AmpZRH.assign(2,2, AZmm *  d1RH(2,2));
    
    AmpZRH = AmpZRH * beta * s/(s-Mz*Mz);

//  Add the correct Zff coupling    
    AmpZRH = AmpZRH * gRe / sW2_tree;

//  J=1 Subamplitudes: gamma
    gslpp::complex Agapp, Agamm, Agap0, Agam0, Aga0p, Aga0m, Aga00;
    
    Agapp = gslpp::complex(g1ga + 2.0* gamma2* lambga, (ktga + lambtga - 2.0*lambtga)/beta , false);
    Agamm = gslpp::complex(g1ga + 2.0* gamma2* lambga, -(ktga + lambtga - 2.0*lambtga)/beta , false);
    Agap0 = gslpp::complex(f3ga + beta * g5ga , -g4ga + (ktga-lambtga)/beta , false);
    Agap0 = gamma * Agap0;
    Agam0 = gslpp::complex(f3ga - beta * g5ga , -g4ga - (ktga-lambtga)/beta , false);
    Agam0 = gamma * Agam0;
    Aga0p = gslpp::complex(f3ga - beta * g5ga , g4ga + (ktga-lambtga)/beta , false);
    Aga0p = gamma * Aga0p;
    Aga0m = gslpp::complex(f3ga + beta * g5ga , g4ga - (ktga-lambtga)/beta , false);
    Aga0m = gamma * Aga0m;
    Aga00 = gslpp::complex( g1ga + 2.0*gamma2*kga, 0.0 , false);

//  Collect in matrices. Here LH = RH, except for the Wigner functions  
    gslpp::matrix<gslpp::complex> AmpgaLH(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> AmpgaRH(3, 3, 0.0);
    
    AmpgaLH.assign(0,0, Agapp * d1LH(0,0));
    AmpgaLH.assign(0,1, Agap0 * d1LH(0,1));
    AmpgaLH.assign(0,2, 0.0);
    
    AmpgaLH.assign(1,0, Aga0p * d1LH(1,0));
    AmpgaLH.assign(1,1, Aga00 * d1LH(1,1));
    AmpgaLH.assign(1,2, Aga0m * d1LH(1,2));
    
    AmpgaLH.assign(2,0, 0.0);
    AmpgaLH.assign(2,1, Agam0 * d1LH(2,1));
    AmpgaLH.assign(2,2, Agamm * d1LH(2,2));
    
    AmpgaRH.assign(0,0, Agapp * d1RH(0,0));
    AmpgaRH.assign(0,1, Agap0 * d1RH(0,1));
    AmpgaRH.assign(0,2, 0.0);
    
    AmpgaRH.assign(1,0, Aga0p * d1RH(1,0));
    AmpgaRH.assign(1,1, Aga00 * d1RH(1,1));
    AmpgaRH.assign(1,2, Aga0m * d1RH(1,2));
    
    AmpgaRH.assign(2,0, 0.0);
    AmpgaRH.assign(2,1, Agam0 * d1RH(2,1));
    AmpgaRH.assign(2,2, Agamm * d1RH(2,2));
    
    AmpgaLH = -beta * AmpgaLH;
    AmpgaRH = -beta * AmpgaRH;
    
//  J=1 Subamplitudes: neutrino
    gslpp::complex Bpp, Bmm, Bp0, Bm0, B0p, B0m, B00;
    gslpp::complex Cpp, Cmm, Cp0, Cm0, C0p, C0m, C00;
    
    Bpp = gslpp::complex(1.0 , 0.0 , false);
    Bmm = Bpp;
    Bp0 = gslpp::complex( 2.0 * gamma, 0.0 , false);
    Bm0 = Bp0;
    B0p = Bp0;
    B0m = Bp0;
    B00 = gslpp::complex( 2.0 * gamma2, 0.0 , false);
    
    Cpp = gslpp::complex(1.0/gamma2 , 0.0 , false);
    Cmm = Cpp;
    Cp0 = gslpp::complex( 2.0 * (1.0 + beta)/gamma, 0.0 , false);
    Cm0 = gslpp::complex( 2.0 * (1.0 - beta)/gamma, 0.0 , false);
    C0p = Cm0;
    C0m = Cp0;
    C00 = gslpp::complex( 2.0 / gamma2, 0.0 , false);
    
//  Collect in matrices. Here LH = RH    
    gslpp::matrix<gslpp::complex> Bnu(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> Cnu(3, 3, 0.0);
    
    Bnu.assign(0,0, Bpp * d1LH(0,0));
    Bnu.assign(0,1, Bp0 * d1LH(0,1));
    Bnu.assign(0,2, 0.0);
    
    Bnu.assign(1,0, B0p * d1LH(1,0));
    Bnu.assign(1,1, B00 * d1LH(1,1));
    Bnu.assign(1,2, B0m * d1LH(1,2));
    
    Bnu.assign(2,0, 0.0);
    Bnu.assign(2,1, Bm0 * d1LH(2,1));
    Bnu.assign(2,2, Bmm * d1LH(2,2));
    
    Cnu.assign(0,0, Cpp * d1LH(0,0));
    Cnu.assign(0,1, Cp0 * d1LH(0,1));
    Cnu.assign(0,2, 0.0);
    
    Cnu.assign(1,0, C0p * d1LH(1,0));
    Cnu.assign(1,1, C00 * d1LH(1,1));
    Cnu.assign(1,2, C0m * d1LH(1,2));
    
    Cnu.assign(2,0, 0.0);
    Cnu.assign(2,1, Cm0 * d1LH(2,1));
    Cnu.assign(2,2, Cmm * d1LH(2,2));
 
//  The matrix with the total J=1 neutrino amplitude (only LH neutrinos)  
    gslpp::matrix<gslpp::complex> Ampnu1(3, 3, 0.0);
    
    Ampnu1 = Bnu - Cnu/(1.0 + beta*beta - 2.0 * beta * cos);
    
    Ampnu1 = Uenu * Uenu.conjugate() * Ampnu1 / (2.0 * beta * sW2_tree);
    
    gslpp::matrix<gslpp::complex> Ampnu2(3, 3, 0.0);

    Ampnu2.assign(0,2, (1.0 - cos)/2.0 );
    Ampnu2.assign(1,1, 0.0);
    Ampnu2.assign(2,0, -(1.0 + cos)/2.0);
    
    Ampnu2 = (8.0 * M_PI * ale / sW2_tree)* Uenu * Uenu.conjugate() * Ampnu2 * sin / (1.0 + beta*beta - 2.0*beta*cos);
    
//  Total amplitudes 
    gslpp::matrix<gslpp::complex> MRH(3, 3, 0.0);
    gslpp::matrix<gslpp::complex> MLH(3, 3, 0.0);

    MRH = sqrt(2.0) * 4.0 * M_PI * ale * (AmpZRH + AmpgaRH);
    MLH = - sqrt(2.0) * 4.0 * M_PI * ale * (AmpZLH + AmpgaLH + Ampnu1) + Ampnu2;
    
//  Total amplitude squared and differential cross section (in pb)
    gslpp::matrix<double> M2(3, 3, 0.0);
    double dxsdcos;
    
    dxsdcos = 0.0;
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            M2.assign(i,j, (MRH(i,j)* (MRH(i,j).conjugate()) 
                    + MLH(i,j)* (MLH(i,j).conjugate())).real() );
            
            dxsdcos = dxsdcos + M2(i,j);
        }
    }
    
//  Differential cross section in pb
    dxsdcos = (topb * beta / 32.0 / M_PI / s) * dxsdcos;

    return dxsdcos;
}

double NPSMEFTd6::dxseeWWdcosBin(const double sqrt_s, const double cos1, const double cos2) const
{
    gsl_integration_cquad_workspace * w_WW;/**< Gsl integral variable */
    w_WW = gsl_integration_cquad_workspace_alloc(100);
    
    double xsWWbin;/**< Gsl integral variable */
//    double errWW;/**< Gsl integral variable */
    
//    gsl_function FR;/**< Gsl integral variable */
    
//    FR = convertToGslFunction(boost::bind(&NPSMEFTd6::dxseeWWdcos,&(*this), sqrt_s, _1));
    
//    FR.function(&NPSMEFTd6::dxsWWdcos);
    
//    gsl_integration_cquad(&FR, cos1, cos2, 1.e-5, 1.e-4, w_WW, &xsWWbin, &errWW, NULL);
    
//  Simple integration for testing
    double cosx;
    
    xsWWbin = 0.0;
    
    for (int i=1; i<100; i++){
        cosx = cos1 +  i*(cos2-cos1)/100;
        xsWWbin = xsWWbin + dxseeWWdcos(sqrt_s, cosx);
    }
    
    xsWWbin = xsWWbin + 0.5 * (dxseeWWdcos(sqrt_s, cos1) + dxseeWWdcos(sqrt_s, cos2));
    
    xsWWbin = xsWWbin * (cos2-cos1)/100;
    
    return xsWWbin;
}

double NPSMEFTd6::xseeWW(const double sqrt_s) const
{    
    return dxseeWWdcosBin(sqrt_s, -1.0, 1.0);
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::kappaGeff() const
{
      return sqrt(GammaHggRatio());
}

double NPSMEFTd6::kappaZeff() const
{
      return sqrt(GammaHZZRatio());
}

double NPSMEFTd6::kappaWeff() const
{
      return sqrt(GammaHWWRatio());
}

double NPSMEFTd6::kappaAeff() const
{
      return sqrt(GammaHgagaRatio());
}

double NPSMEFTd6::kappaZAeff() const
{
      return sqrt(GammaHZgaRatio());
}

///////////////////////////////////////////////////////////////////////////////

double NPSMEFTd6::CLL_mu() const
{
    return (CLL_1122 + CLL_2211 + CLL_1221 + CLL_2112);
}

double NPSMEFTd6::CLL_tau() const
{
    return (CLL_1133 + CLL_3311 + CLL_1331 + CLL_3113);
}

double NPSMEFTd6::CLL_up() const
{
    return (CLQ1_1111-CLQ3_1111);
}

double NPSMEFTd6::CLL_down() const
{
    return (CLQ1_1111+CLQ3_1111);
}

double NPSMEFTd6::CLL_charm() const
{
    return (CLQ1_1122+CLQ1_2211-CLQ3_1122-CLQ3_2211);
}

double NPSMEFTd6::CLL_strange() const
{
    return (CLQ1_1122+CLQ1_2211+CLQ3_1122+CLQ3_2211);
}

double NPSMEFTd6::CLL_bottom() const
{
    return (CLQ1_1133+CLQ1_3311+CLQ3_1133+CLQ3_3311);
}

double NPSMEFTd6::CLR_mu() const
{
    return (CLe_1122+CLe_2211);
}

double NPSMEFTd6::CLR_tau() const
{
    return (CLe_1133+CLe_3311);
}

double NPSMEFTd6::CLR_up() const
{
    return (CLu_1111);
}

double NPSMEFTd6::CLR_down() const
{
    return (CLd_1111);
}

double NPSMEFTd6::CLR_charm() const
{
    return (CLu_1122+CLu_2211);
}

double NPSMEFTd6::CLR_strange() const
{
    return (CLd_1122+CLd_2211);
}

double NPSMEFTd6::CLR_bottom() const
{
    return (CLd_1133+CLd_3311);
}

double NPSMEFTd6::CRL_mu() const
{
    return (CLe_1122+CLe_2211);
}

double NPSMEFTd6::CRL_tau() const
{
    return (CLe_1133+CLe_3311);
}

double NPSMEFTd6::CRL_up() const
{
    return (CQe_1111);
}

double NPSMEFTd6::CRL_down() const
{
    return (CQe_1111);
}

double NPSMEFTd6::CRL_charm() const
{
    return (CQe_1122+CQe_2211);
}

double NPSMEFTd6::CRL_strange() const
{
    return (CQe_1122+CQe_2211);
}

double NPSMEFTd6::CRL_bottom() const
{
    return (CQe_1133+CQe_3311);
}

double NPSMEFTd6::CRR_mu() const
{
    return (Cee_1122+Cee_2211);
}

double NPSMEFTd6::CRR_tau() const
{
    return (Cee_1133+Cee_3311);
}


double NPSMEFTd6::CRR_up() const
{
    return (Ceu_1111);
}

double NPSMEFTd6::CRR_down() const
{
    return (Ced_1111);
}

double NPSMEFTd6::CRR_charm() const
{
    return (Ceu_1122+Ceu_2211);
}

double NPSMEFTd6::CRR_strange() const
{
    return (Ced_1122+Ced_2211);
}

double NPSMEFTd6::CRR_bottom() const
{
    return (Ced_1133+Ced_3311);
}
