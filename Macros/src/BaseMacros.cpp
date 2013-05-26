/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TROOT.h>
#include "BaseMacros.h"

void BaseMacros::DefineNewColours() 
{
    TColor *color[200];
    color[0] = new TColor(1300, 1.0, 1.0, 1.0, "");
    color[1] = new TColor(1301, 0.0, 0.0, 0.0, "");
    color[2] = new TColor(1302, 1.0, 0.0, 0.0, "");
    color[3] = new TColor(1303, 0.0, 1.0, 0.0, "");
    color[4] = new TColor(1304, 0.0, 0.0, 1.0, "");
    color[11] = new TColor(1311, 0.6, 0.4, 1.0, "");
    color[12] = new TColor(1312, 1.0, 0.5, 0.5, "");
    color[13] = new TColor(1313, 0.5, 1.0, 0.4, "");
    color[14] = new TColor(1314, 0.8, 0.6, 1.0, "");
    color[15] = new TColor(1315, 1.0, 0.90, 0.82, "");
    color[16] = new TColor(1316, 1.0, 0.75, 0.67, "");
    color[17] = new TColor(1317, 0.5, 1.0, 0.5, "");
    color[18] = new TColor(1318, 0.95, 1.0, 1.0, "");
    color[19] = new TColor(1319, 0.7, 1.0, 0.8, "");
    color[20] = new TColor(1320, 0.97, 0.87, 0.97, "");
    color[21] = new TColor(1321, 0.85, 0.44, 0.84, "");
    color[22] = new TColor(1322, 0.73, 0.33, 0.83, "");
    color[23] = new TColor(1323, 0.60, 0.20, 0.80, "");
    color[24] = new TColor(1324, 0.58, 0.00, 0.83, "");
    color[25] = new TColor(1325, 0.54, 0.17, 0.89, "");
    color[26] = new TColor(1326, 0.63, 0.13, 0.94, "");
    color[27] = new TColor(1327, 0.58, 0.44, 0.86, "");
    color[28] = new TColor(1328, 0.95, 0.85, 0.95, "");
    color[29] = new TColor(1329, 0.93, 0.51, 0.93, "");
    color[30] = new TColor(1330, 1.00, 0.98, 0.98, "");
    color[31] = new TColor(1331, 0.97, 0.97, 1.00, "");
    color[32] = new TColor(1332, 1.00, 1.00, 1.00, "");
    color[33] = new TColor(1333, 0.52, 0.44, 1.00, "");
    color[34] = new TColor(1334, 0.12, 0.56, 1.00, "");
    color[35] = new TColor(1335, 0.00, 0.75, 1.00, "");
    color[36] = new TColor(1336, 0.88, 1.00, 1.00, "");
    color[37] = new TColor(1337, 0.37, 0.62, 0.63, "");
    color[38] = new TColor(1338, 0.40, 0.80, 0.67, "");
    color[39] = new TColor(1339, 0.50, 1.00, 0.83, "");
    color[40] = new TColor(1340, 0.00, 0.39, 0.00, "");
    color[41] = new TColor(1341, 0.33, 0.42, 0.18, "");
    color[42] = new TColor(1342, 0.56, 0.74, 0.56, "");
    color[43] = new TColor(1343, 0.18, 0.55, 0.34, "");
    color[44] = new TColor(1344, 0.24, 0.70, 0.44, "");
    color[45] = new TColor(1345, 0.13, 0.70, 0.67, "");
    color[46] = new TColor(1346, 0.60, 0.98, 0.60, "");
    color[47] = new TColor(1347, 0.00, 1.00, 0.50, "");
    color[48] = new TColor(1348, 0.49, 0.99, 0.00, "");
    color[49] = new TColor(1349, 0.00, 1.00, 0.00, "");
    color[50] = new TColor(1350, 0.50, 1.00, 0.00, "");
    color[51] = new TColor(1351, 0.00, 0.98, 0.60, "");
    color[52] = new TColor(1352, 0.68, 1.00, 0.18, "");
    color[53] = new TColor(1353, 0.20, 0.80, 0.20, "");
    color[54] = new TColor(1354, 0.60, 0.80, 0.20, "");
    color[55] = new TColor(1355, 0.13, 0.55, 0.13, "");
    color[56] = new TColor(1356, 0.42, 0.56, 0.14, "");
    color[57] = new TColor(1357, 0.74, 0.72, 0.42, "");
    color[58] = new TColor(1358, 0.94, 0.90, 0.55, "");
    color[59] = new TColor(1359, 0.93, 0.91, 0.67, "");
    color[60] = new TColor(1360, 0.99, 0.99, 0.75, "");
    color[61] = new TColor(1361, 1.00, 0.70, 0.00, "");
    color[62] = new TColor(1362, 1.00, 1.00, 0.00, "");
    color[63] = new TColor(1363, 1.00, 0.84, 0.00, "");
    color[64] = new TColor(1364, 1.00, 0.95, 0.60, "");
    color[65] = new TColor(1365, 0.85, 0.65, 0.13, "");
    color[66] = new TColor(1366, 0.72, 0.53, 0.04, "");
    color[67] = new TColor(1367, 0.74, 0.56, 0.56, "");
    color[68] = new TColor(1368, 0.80, 0.36, 0.36, "");
    color[69] = new TColor(1369, 0.55, 0.27, 0.07, "");
    color[70] = new TColor(1370, 0.63, 0.32, 0.18, "");
    color[71] = new TColor(1371, 0.80, 0.52, 0.25, "");
    color[72] = new TColor(1372, 0.99, 0.85, 0.70, "");
    color[73] = new TColor(1373, 0.96, 0.96, 0.86, "");
    color[74] = new TColor(1374, 0.99, 0.90, 0.75, "");
    color[75] = new TColor(1375, 0.96, 0.64, 0.38, "");
    color[76] = new TColor(1376, 0.82, 0.71, 0.55, "");
    color[77] = new TColor(1377, 0.82, 0.41, 0.12, "");
    color[78] = new TColor(1378, 0.70, 0.13, 0.13, "");
    color[79] = new TColor(1379, 0.65, 0.16, 0.16, "");
    color[80] = new TColor(1380, 0.91, 0.59, 0.48, "");
    color[81] = new TColor(1381, 0.98, 0.50, 0.45, "");
    color[82] = new TColor(1382, 1.00, 0.85, 0.80, "");
    color[83] = new TColor(1383, 1.00, 0.65, 0.00, "");
    color[84] = new TColor(1384, 1.00, 0.55, 0.00, "");
    color[85] = new TColor(1385, 1.00, 0.50, 0.31, "");
    color[86] = new TColor(1386, 0.94, 0.50, 0.50, "");
    color[87] = new TColor(1387, 1.00, 0.39, 0.28, "");
    color[88] = new TColor(1388, 1.00, 0.27, 0.00, "");
    color[89] = new TColor(1389, 1.00, 0.00, 0.00, "");

    color[90] = new TColor(1390, 0.90, 0.60, 0.60, ""); //red
    color[91] = new TColor(1391, 0.70, 0.25, 0.25, "");
    color[92] = new TColor(1392, 0.87, 0.87, 0.91, ""); //blue
    color[93] = new TColor(1393, 0.59, 0.58, 0.91, "");
    color[94] = new TColor(1394, 0.65, 0.55, 0.85, ""); //violet
    color[95] = new TColor(1395, 0.49, 0.26, 0.64, "");
    color[96] = new TColor(1396, 0.95, 0.95, 0.45, ""); // yellow
    color[97] = new TColor(1397, 0.95, 0.95, 0.05, "");
    color[98] = new TColor(1398, 0.75, 0.92, 0.68, ""); //green
    color[99] = new TColor(1399, 0.36, 0.57, 0.30, "");
    color[100] = new TColor(1400, 0.97, 0.50, 0.09, ""); // orange
    color[101] = new TColor(1401, 0.76, 0.34, 0.09, "");
    color[102] = new TColor(1402, 0.97, 0.52, 0.75, ""); // pink 
    color[103] = new TColor(1403, 0.76, 0.32, 0.51, "");
    color[104] = new TColor(1404, 0.49, 0.60, 0.82, ""); // dark blue
    color[105] = new TColor(1405, 0.43, 0.48, 0.52, "");
    color[106] = new TColor(1406, 0.70, 0.70, 0.70, ""); // black
    color[107] = new TColor(1407, 0.40, 0.40, 0.40, "");
    color[108] = new TColor(1408, 0.10, 0.60, 0.10, ""); // dark green
    color[109] = new TColor(1409, 0.10, 0.40, 0.10, "");
    color[110] = new TColor(1410, 0.00, 0.00, 0.00, "");

    /* 
     *   cf. white   1,1,1
     *       black   0,0,0
     *       red     1,0,0
     *       green   0,1,0
     *       blue    0,0,1
     *       yellow  1,1,0
     *       magenta 1,0,1
     *       cyan    0,1,1
     */
    color[152] = new TColor(2152, 100./255., 149./255., 237./255.); // corn flower blue #6495ED
    //
    color[153] = new TColor(2153, 255./255., 140./255., 0.); // dark orange #FF8C00
    //
    color[154] = new TColor(2154, 0., 1., 0.); // lime #00FF00
    color[155] = new TColor(2155, 173./255., 1., 47./255.); // green yellow #ADFF2F
    color[156] = new TColor(2156, 0., 100./255., 0.); // dark green #006400
    //
    color[157] = new TColor(2157, 192./255., 192./255., 192./255.); // silver
    //
    color[158] = new TColor(2158, 72./255., 61./255., 139./255.); // dark slate blue #483D8B
    color[159] = new TColor(2159, 106./255., 90./255., 205./255.); // slate blue #6A5ACD
    //
    color[160] = new TColor(2160, 1.00, 0.00, 1.00); // magenta
    color[161] = new TColor(2161, 1., 20./255., 147./255.); // deep pink #FF1493
    color[162] = new TColor(2162, 238./255., 130./255., 238./255.); // violet #EE82EE
    //
    color[163] = new TColor(2163, 0., 206./255., 209./255.); // dark turquoise #00CED1
    color[164] = new TColor(2164, 64./255., 224./255., 208./255.); // turquoise #40E0D0
    color[165] = new TColor(2165, 175./255., 238./255., 238./255.); // pale turquoise #AFEEEE
    //
    color[166] = new TColor(2166, 135./255., 206./255., 250./255.); // light sky blue #87CEFA
    color[167] = new TColor(2167, 224./255., 1., 1.); // light cyan #E0FFFF
}


TString BaseMacros::ConvertTitle(TString orig) const 
{
    TString str;

    if (strncmp(orig, "sin2b", 5) == 0) str = "sin2#beta";
    else if (strncmp(orig, "sin2bpg", 7) == 0) str = "sin(2#beta+#gamma)";
    else if (strncmp(orig, "sin2a", 5) == 0) str = "sin2#alpha";
    else if (strncmp(orig, "alpha", 5) == 0) str = "#alpha[^{o}]";
    else if (strncmp(orig, "beta", 4) == 0) str = "#beta[^{o}]";
    else if (strncmp(orig, "gamma", 5) == 0) str = "#gamma[^{o}]";
    else if (strncmp(orig, "2bpg", 4) == 0) str = "2#beta+#gamma[^{o}]";
    else if (strncmp(orig, "rho", 3) == 0) str = "#bar{#rho}";
    else if (strncmp(orig, "eta", 3) == 0) str = "#bar{#eta}";
    else if (strncmp(orig, "rhof", 4) == 0) str = "#bar{#rho}";
    else if (strncmp(orig, "etaf", 4) == 0) str = "#bar{#eta}";
    else if (strncmp(orig, "Lambda", 6) == 0) str = "#lambda";
    else if (strncmp(orig, "Relamt", 6) == 0) str = "Re#lambda_{t}[10^{-3}]";
    else if (strncmp(orig, "Imlamt", 6) == 0) str = "Im#lambda_{t}[10^{-5}]";
    else if (strncmp(orig, "Rb", 2) == 0) str = "R_{b}";
    else if (strncmp(orig, "Rt", 2) == 0) str = "R_{t}";
    else if (strncmp(orig, "Vub", 3) == 0) str = "V_{ub}";
    else if (strncmp(orig, "Vcb", 3) == 0) str = "V_{cb}";
    else if (strncmp(orig, "Vtd", 3) == 0) str = "V_{td}#times10^{3}";
    else if (strncmp(orig, "DeltaMs", 7) == 0) str = "#Delta m_{s}[ps^{-1}]";
    else if (strncmp(orig, "FbSqrtb", 7) == 0) str = "F_{B}#sqrt{B}";
    else if (strncmp(orig, "Bk", 2) == 0) str = "B_{k}";
    else if (strncmp(orig, "Xi", 2) == 0) str = "#xi";
    else if (strncmp(orig, "Ftt", 3) == 0) str = "F_{tt}";
    else if (strncmp(orig, "epsk", 4) == 0) str = "#epsilon_{K}";
    else if (strncmp(orig, "btaunu", 6) == 0) str = "BR(B#rightarrow#tau#nu)";

    //-- Model parameters --
    else if (strncmp(orig, "AlsMz", 5) == 0) str = "#alpha_{#lower[-0.2]{s}}(M^{#kern[0.2]{2}}_{#kern[0.2]{#lower[-0.2]{Z}}})";
    else if (strncmp(orig, "dAle5Mz", 7) == 0) str = "#Delta#alpha_{#lower[-0.2]{had}}^{(5)}(M^{#kern[0.2]{2}}_{#kern[0.2]{#lower[-0.2]{Z}}})";
    else if (strncmp(orig, "Mz", 2) == 0) str = "M_{#kern[0.2]{#lower[-0.2]{Z}}}#kern[0.1]{[GeV]}";
    //else if (strncmp(orig, "mtop", 4) == 0) str = "M_{#kern[0.7]{#lower[-0.2]{t}}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "mtop", 4) == 0) str = "m_{#kern[0.5]{#lower[-0.2]{t}}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "mHl", 3) == 0) str = "m_{#kern[0.2]{#lower[-0.2]{h}}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "obliqueS", 8) == 0) str = "S" ;
    else if (strncmp(orig, "obliqueT", 8) == 0) str = "T" ;
    else if (strncmp(orig, "obliqueU", 8) == 0) str = "U" ;
    else if (strncmp(orig, "epsilon_1", 9) == 0) str = "#varepsilon_{#lower[-0.2]{1}}" ;
    else if (strncmp(orig, "epsilon_2", 9) == 0) str = "#varepsilon_{#lower[-0.2]{2}}" ;    
    else if (strncmp(orig, "epsilon_3", 9) == 0) str = "#varepsilon_{#lower[-0.2]{3}}" ;
    else if (strncmp(orig, "epsilon_b", 9) == 0) str = "#varepsilon_{#lower[-0.2]{b}}" ;
    else if (strncmp(orig, "epsilon1", 8) == 0) str = "#varepsilon_{#lower[-0.2]{1}}" ;
    else if (strncmp(orig, "epsilon2", 8) == 0) str = "#varepsilon_{#lower[-0.2]{2}}" ;    
    else if (strncmp(orig, "epsilon3", 8) == 0) str = "#varepsilon_{#lower[-0.2]{3}}" ;
    else if (strncmp(orig, "epsilonb", 8) == 0) str = "#varepsilon_{#lower[-0.2]{b}}" ;
    else if (strncmp(orig, "#varepsilon_{1}", 15) == 0) str = "#varepsilon_{#lower[-0.2]{1}}" ;
    else if (strncmp(orig, "#varepsilon_{2}", 15) == 0) str = "#varepsilon_{#lower[-0.2]{2}}" ;
    else if (strncmp(orig, "#varepsilon_{3}", 15) == 0) str = "#varepsilon_{#lower[-0.2]{3}}" ;
    else if (strncmp(orig, "#varepsilon_{b}", 15) == 0) str = "#varepsilon_{#lower[-0.2]{b}}" ;
    
    //-- EW precision observables --
    else if (strncmp(orig, "GammaZ", 6) == 0) str = "#Gamma_{#lower[-0.2]{Z}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "#Gamma_{Z}", 10) == 0) str = "#Gamma_{#lower[-0.2]{Z}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "sigmaHadron", 11) == 0) str = "#sigma^{0}_{#lower[-0.2]{h}}#kern[0.1]{[nb]}";
    else if (strncmp(orig, "#sigma_{had}", 12) == 0) str = "#sigma^{0}_{#lower[-0.2]{h}}#kern[0.1]{[nb]}";
    else if (strncmp(orig, "Rlepton", 7) == 0) str = "R^{0}_{#kern[0.4]{#lower[-0.2]{l}}}";
    else if (strncmp(orig, "R_{l}", 5) == 0) str = "R^{0}_{#kern[0.4]{#lower[-0.2]{l}}}";
    else if (strncmp(orig, "AFBlepton", 9) == 0) str = "A_{#lower[-0.2]{FB}}^{0,l}";
    else if (strncmp(orig, "A_{FB}^{l}", 10) == 0) str = "A_{#lower[-0.2]{FB}}^{0,l}";
    else if (strncmp(orig, "PtauPol", 7) == 0) str = "P_{#lower[-0.2]{#tau}}^{Pol}";
    else if (strncmp(orig, "P_{#tau}^{Pol}", 14) == 0) str = "P_{#lower[-0.2]{#tau}}^{Pol}";
    else if (strncmp(orig, "sin2thetaEff", 12) == 0) str = "sin^{2}#theta_{#lower[-0.3]{eff}}^{lept}";
    else if (strncmp(orig, "sin^{2}#theta_{eff}", 19) == 0) str = "sin^{2}#theta_{#lower[-0.3]{eff}}^{lept}";
    else if (strncmp(orig, "Alepton", 7) == 0) str = "A_{#kern[0.2]{#lower[-0.2]{l}}}";
    else if (strncmp(orig, "A_{l}", 5) == 0) str = "A_{#kern[0.2]{#lower[-0.2]{l}}}";
    else if (strncmp(orig, "Rbottom", 7) == 0) str = "R^{0}_{#lower[-0.2]{b}}";
    else if (strncmp(orig, "R_{b}", 5) == 0) str = "R^{0}_{#lower[-0.2]{b}}";
    else if (strncmp(orig, "Rcharm", 6) == 0) str = "R^{0}_{#lower[-0.2]{c}}";
    else if (strncmp(orig, "R_{c}", 5) == 0) str = "R^{0}_{#lower[-0.2]{c}}";
    else if (strncmp(orig, "AFBbottom", 9) == 0) str = "A_{#lower[-0.2]{FB}}^{0,b}";
    else if (strncmp(orig, "A_{FB}^{b}", 10) == 0) str = "A_{#lower[-0.2]{FB}}^{0,b}";
    else if (strncmp(orig, "AFBcharm", 8) == 0) str = "A_{#lower[-0.2]{FB}}^{0,c}";
    else if (strncmp(orig, "A_{FB}^{c}", 10) == 0) str = "A_{#lower[-0.2]{FB}}^{0,c}";
    else if (strncmp(orig, "Abottom", 7) == 0) str = "A_{#lower[-0.2]{b}}";
    else if (strncmp(orig, "A_{b}", 5) == 0) str = "A_{#lower[-0.2]{b}}";
    else if (strncmp(orig, "Acharm", 6) == 0) str = "A_{#lower[-0.2]{c}}";
    else if (strncmp(orig, "A_{c}", 5) == 0) str = "A_{#lower[-0.2]{c}}";
    else if (strncmp(orig, "Mw", 2) == 0) str = "M_{#kern[0.2]{#lower[-0.2]{W}}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "M_{W}", 5) == 0) str = "M_{#kern[0.2]{#lower[-0.2]{W}}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "GammaW", 6) == 0) str = "#Gamma_{#lower[-0.2]{W}}#kern[0.1]{[GeV]}";
    else if (strncmp(orig, "#Gamma_{W}", 10) == 0) str = "#Gamma_{#lower[-0.2]{W}}#kern[0.1]{[GeV]}";
    //
    else if (strncmp(orig, "AllForST", 8) == 0) str = "All (incl. #Gamma_{#lower[-0.2]{W}}, #sigma^{0}_{#lower[-0.2]{h}}, R^{0}_{#lower[-0.5]{f}})";
    else if (strncmp(orig, "AllForEPS", 9) == 0) str = "All (incl. #sigma^{0}_{#lower[-0.2]{h}}, R^{0}_{#lower[-0.5]{f}})";
    else if (strncmp(orig, "AllAsymmetries", 14) == 0) str = "sin^{2}#theta_{#lower[-0.3]{eff}}^{lept}, P_{#lower[-0.3]{#tau}}^{Pol}, A_{#kern[0.2]{#lower[-0.3]{f}}}, A_{#lower[-0.3]{FB}}^{0,f}";    
    else if (strncmp(orig, "EPS2EPSBSM", 9) == 0) str = "#varepsilon_{#lower[-0.2]{2}}=#varepsilon_{#lower[-0.2]{2}}^{SM}, #varepsilon_{#lower[-0.2]{b}}=#varepsilon_{#lower[-0.2]{b}}^{SM}" ;    
    //
    else if (strncmp(orig, "#deltag_{V}^{b}", 15) == 0) str = "#delta#kern[0.1]{g_{#kern[0.1]{#lower[-0.2]{V}}}^{#kern[0.2]{b}}}";
    else if (strncmp(orig, "#deltag_{A}^{b}", 15) == 0) str = "#delta#kern[0.1]{g_{#kern[0.1]{#lower[-0.2]{A}}}^{#kern[0.2]{b}}}";
    else if (strncmp(orig, "#deltag_{R}^{b}", 15) == 0) str = "#delta#kern[0.1]{g_{#kern[0.1]{#lower[-0.2]{R}}}^{#kern[0.2]{b}}}";
    else if (strncmp(orig, "#deltag_{L}^{b}", 15) == 0) str = "#delta#kern[0.1]{g_{#kern[0.1]{#lower[-0.2]{L}}}^{#kern[0.2]{b}}}";
    else if (strncmp(orig, "#delta#rho_{Z}^{b}", 18) == 0) str = "#delta#kern[0.1]{#rho_{#kern[0.1]{#lower[-0.2]{Z}}}^{#kern[0.2]{b}}}";
    else if (strncmp(orig, "#delta#kappa_{Z}^{b}", 20) == 0) str = "#delta#kern[0.1]{#kappa_{#kern[0.1]{#lower[-0.1]{Z}}}^{#kern[0.2]{b}}}";
    
    else str = orig;

    str.ReplaceAll("WHITESPACE", " ");
    
    return str;
}


