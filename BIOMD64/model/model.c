#ifdef SIZE_DEFINITIONS
#define N_METABS 26
#define N_ODE_METABS 0
#define N_INDEP_METABS 13
#define N_COMPARTMENTS 2
#define N_GLOBAL_PARAMS 15
#define N_KIN_PARAMS 70
#define N_REACTIONS 17

#define N_ARRAY_SIZE_P  96	// number of parameters
#define N_ARRAY_SIZE_X  13	// number of initials
#define N_ARRAY_SIZE_Y  4	// number of assigned elements
#define N_ARRAY_SIZE_XC 13	// number of x concentration
#define N_ARRAY_SIZE_PC 9	// number of p concentration
#define N_ARRAY_SIZE_YC 4	// number of y concentration
#define N_ARRAY_SIZE_DX 13	// number of ODEs 
#define N_ARRAY_SIZE_CT 1	// number of conserved totals

#endif // SIZE_DEFINITIONS

#ifdef TIME
#define T  <set here a user name for the time variable> 
#endif // TIME

#ifdef NAME_ARRAYS
const char* p_names[] = {"Extracellular Glucose", "Glycogen", "Trehalose", "CO2", "Succinate", "Ethanol", "Glycerol", "sum of AXP conc", "F2,6P", "extracellular", "cytosol", "gR", "KmPFKF6P", "KmPFKATP", "Lzero", "CiPFKATP", "KiPFKATP", "CPFKAMP", "KPFKAMP", "CPFKF26BP", "KPFKF26BP", "CPFKF16BP", "KPFKF16BP", "CPFKATP", "AK eq constant", "TPI eq constant", "KeqGLK", "KmGLKADP", "KmGLKATP", "KmGLKG6P", "KmGLKGLCi", "VmGLK", "KeqPGI_2", "KmPGIF6P_2", "KmPGIG6P_2", "VmPGI_2", "v", "v", "VmPFK", "KeqALD", "KmALDDHAP", "KmALDF16P", "KmALDGAP", "KmALDGAPi", "VmALD", "KmGAPDHBPG", "KmGAPDHGAP", "KmGAPDHNAD", "KmGAPDHNADH", "VmGAPDHf", "VmGAPDHr", "KeqPGK", "KmPGKADP", "KmPGKATP", "KmPGKBPG", "KmPGKP3G", "VmPGK", "KeqPGM", "KmPGMP2G", "KmPGMP3G", "VmPGM", "KeqENO", "KmENOP2G", "KmENOPEP", "VmENO", "KeqPYK", "KmPYKADP", "KmPYKATP", "KmPYKPEP", "KmPYKPYR", "VmPYK", "KmPDCPYR", "VmPDC", "nPDC", "KSUCC", "KeqGLT", "KmGLTGLCi", "KmGLTGLCo", "VmGLT", "KeqADH", "KiADHACE", "KiADHETOH", "KiADHNAD", "KiADHNADH", "KmADHACE", "KmADHETOH", "KmADHNAD", "KmADHNADH", "VmADH", "KeqG3PDH", "KmG3PDHDHAP", "KmG3PDHGLY", "KmG3PDHNAD", "KmG3PDHNADH", "VmG3PDH", "KATPASE",  "" };
const char* x_names[] = {"High energy phosphates", "Glucose 6 Phosphate", "Triose-phosphate", "NAD", "Acetaldehyde", "2-phosphoglycerate", "1,3-bisphosphoglycerate", "Fructose 6 Phosphate", "Glucose in Cytosol", "Phosphoenolpyruvate", "Pyruvate", "Fructose-1,6 bisphosphate", "3-phosphoglycerate",  "" };
const char* y_names[] = {"NADH", "ATP concentration", "ADP concentration", "AMP concentration",  "" };
const char* xc_names[] = {"High energy phosphates", "Glucose 6 Phosphate", "Triose-phosphate", "NAD", "Acetaldehyde", "2-phosphoglycerate", "1,3-bisphosphoglycerate", "Fructose 6 Phosphate", "Glucose in Cytosol", "Phosphoenolpyruvate", "Pyruvate", "Fructose-1,6 bisphosphate", "3-phosphoglycerate",  "" };
const char* pc_names[] = {"Extracellular Glucose", "Glycogen", "Trehalose", "CO2", "Succinate", "Ethanol", "Glycerol", "sum of AXP conc", "F2,6P",  "" };
const char* yc_names[] = {"NADH", "ATP concentration", "ADP concentration", "AMP concentration",  "" };
const char* dx_names[] = {"ODE High energy phosphates", "ODE Glucose 6 Phosphate", "ODE Triose-phosphate", "ODE NAD", "ODE Acetaldehyde", "ODE 2-phosphoglycerate", "ODE 1,3-bisphosphoglycerate", "ODE Fructose 6 Phosphate", "ODE Glucose in Cytosol", "ODE Phosphoenolpyruvate", "ODE Pyruvate", "ODE Fructose-1,6 bisphosphate", "ODE 3-phosphoglycerate",  "" };
const char* ct_names[] = {"CT NADH",  "" };
#endif // NAME_ARRAYS

#ifdef INITIAL
x[0] = 6.31;	//metabolite 'High energy phosphates': reactions
x[1] = 2.45;	//metabolite 'Glucose 6 Phosphate': reactions
x[2] = 0.96;	//metabolite 'Triose-phosphate': reactions
x[3] = 1.2;	//metabolite 'NAD': reactions
x[4] = 0.17;	//metabolite 'Acetaldehyde': reactions
x[5] = 0.12;	//metabolite '2-phosphoglycerate': reactions
x[6] = 0;	//metabolite '1,3-bisphosphoglycerate': reactions
x[7] = 0.62;	//metabolite 'Fructose 6 Phosphate': reactions
x[8] = 0.087;	//metabolite 'Glucose in Cytosol': reactions
x[9] = 0.07;	//metabolite 'Phosphoenolpyruvate': reactions
x[10] = 1.85;	//metabolite 'Pyruvate': reactions
x[11] = 5.51;	//metabolite 'Fructose-1,6 bisphosphate': reactions
x[12] = 0.9;	//metabolite '3-phosphoglycerate': reactions
#endif /* INITIAL */

#ifdef FIXED
ct[0] = 1.5899999999999999;	//ct[0] conserved total for 'NADH'
p[0] = 50;	//metabolite 'Extracellular Glucose': fixed
p[1] = 0;	//metabolite 'Glycogen': fixed
p[2] = 0;	//metabolite 'Trehalose': fixed
p[3] = 1;	//metabolite 'CO2': fixed
p[4] = 0;	//metabolite 'Succinate': fixed
p[5] = 50;	//metabolite 'Ethanol': fixed
p[6] = 0.15;	//metabolite 'Glycerol': fixed
p[7] = 4.1;	//metabolite 'sum of AXP conc': fixed
p[8] = 0.02;	//metabolite 'F2,6P': fixed
p[9] = 1;	//compartment 'extracellular':fixed
p[10] = 1;	//compartment 'cytosol':fixed
p[11] = 5.12;	//global quantity 'gR':fixed
p[12] = 0.1;	//global quantity 'KmPFKF6P':fixed
p[13] = 0.71;	//global quantity 'KmPFKATP':fixed
p[14] = 0.66;	//global quantity 'Lzero':fixed
p[15] = 100;	//global quantity 'CiPFKATP':fixed
p[16] = 0.65;	//global quantity 'KiPFKATP':fixed
p[17] = 0.0845;	//global quantity 'CPFKAMP':fixed
p[18] = 0.0995;	//global quantity 'KPFKAMP':fixed
p[19] = 0.0174;	//global quantity 'CPFKF26BP':fixed
p[20] = 0.000682;	//global quantity 'KPFKF26BP':fixed
p[21] = 0.397;	//global quantity 'CPFKF16BP':fixed
p[22] = 0.111;	//global quantity 'KPFKF16BP':fixed
p[23] = 3;	//global quantity 'CPFKATP':fixed
p[24] = 0.45;	//global quantity 'AK eq constant':fixed
p[25] = 0.045;	//global quantity 'TPI eq constant':fixed
p[26] = 3800;	//reaction 'Hexokinase':  kinetic parameter 'KeqGLK'
p[27] = 0.23;	//reaction 'Hexokinase':  kinetic parameter 'KmGLKADP'
p[28] = 0.15;	//reaction 'Hexokinase':  kinetic parameter 'KmGLKATP'
p[29] = 30;	//reaction 'Hexokinase':  kinetic parameter 'KmGLKG6P'
p[30] = 0.08;	//reaction 'Hexokinase':  kinetic parameter 'KmGLKGLCi'
p[31] = 226.452;	//reaction 'Hexokinase':  kinetic parameter 'VmGLK'
p[32] = 0.314;	//reaction 'Glucose-6-phosphate isomerase':  kinetic parameter 'KeqPGI_2'
p[33] = 0.3;	//reaction 'Glucose-6-phosphate isomerase':  kinetic parameter 'KmPGIF6P_2'
p[34] = 1.4;	//reaction 'Glucose-6-phosphate isomerase':  kinetic parameter 'KmPGIG6P_2'
p[35] = 339.677;	//reaction 'Glucose-6-phosphate isomerase':  kinetic parameter 'VmPGI_2'
p[36] = 6;	//reaction 'Glycogen synthesis':  kinetic parameter 'v'
p[37] = 2.4;	//reaction 'Trehalose 6-phosphate synthase':  kinetic parameter 'v'
p[38] = 182.903;	//reaction 'Phosphofructokinase':  kinetic parameter 'VmPFK'
p[39] = 0.069;	//reaction 'Aldolase':  kinetic parameter 'KeqALD'
p[40] = 2.4;	//reaction 'Aldolase':  kinetic parameter 'KmALDDHAP'
p[41] = 0.3;	//reaction 'Aldolase':  kinetic parameter 'KmALDF16P'
p[42] = 2;	//reaction 'Aldolase':  kinetic parameter 'KmALDGAP'
p[43] = 10;	//reaction 'Aldolase':  kinetic parameter 'KmALDGAPi'
p[44] = 322.258;	//reaction 'Aldolase':  kinetic parameter 'VmALD'
p[45] = 0.0098;	//reaction 'Glyceraldehyde 3-phosphate dehydrogenase':  kinetic parameter 'KmGAPDHBPG'
p[46] = 0.21;	//reaction 'Glyceraldehyde 3-phosphate dehydrogenase':  kinetic parameter 'KmGAPDHGAP'
p[47] = 0.09;	//reaction 'Glyceraldehyde 3-phosphate dehydrogenase':  kinetic parameter 'KmGAPDHNAD'
p[48] = 0.06;	//reaction 'Glyceraldehyde 3-phosphate dehydrogenase':  kinetic parameter 'KmGAPDHNADH'
p[49] = 1184.52;	//reaction 'Glyceraldehyde 3-phosphate dehydrogenase':  kinetic parameter 'VmGAPDHf'
p[50] = 6549.8;	//reaction 'Glyceraldehyde 3-phosphate dehydrogenase':  kinetic parameter 'VmGAPDHr'
p[51] = 3200;	//reaction 'Phosphoglycerate kinase':  kinetic parameter 'KeqPGK'
p[52] = 0.2;	//reaction 'Phosphoglycerate kinase':  kinetic parameter 'KmPGKADP'
p[53] = 0.3;	//reaction 'Phosphoglycerate kinase':  kinetic parameter 'KmPGKATP'
p[54] = 0.003;	//reaction 'Phosphoglycerate kinase':  kinetic parameter 'KmPGKBPG'
p[55] = 0.53;	//reaction 'Phosphoglycerate kinase':  kinetic parameter 'KmPGKP3G'
p[56] = 1306.45;	//reaction 'Phosphoglycerate kinase':  kinetic parameter 'VmPGK'
p[57] = 0.19;	//reaction 'Phosphoglycerate mutase':  kinetic parameter 'KeqPGM'
p[58] = 0.08;	//reaction 'Phosphoglycerate mutase':  kinetic parameter 'KmPGMP2G'
p[59] = 1.2;	//reaction 'Phosphoglycerate mutase':  kinetic parameter 'KmPGMP3G'
p[60] = 2525.81;	//reaction 'Phosphoglycerate mutase':  kinetic parameter 'VmPGM'
p[61] = 6.7;	//reaction 'Enolase':  kinetic parameter 'KeqENO'
p[62] = 0.04;	//reaction 'Enolase':  kinetic parameter 'KmENOP2G'
p[63] = 0.5;	//reaction 'Enolase':  kinetic parameter 'KmENOPEP'
p[64] = 365.806;	//reaction 'Enolase':  kinetic parameter 'VmENO'
p[65] = 6500;	//reaction 'Pyruvate kinase':  kinetic parameter 'KeqPYK'
p[66] = 0.53;	//reaction 'Pyruvate kinase':  kinetic parameter 'KmPYKADP'
p[67] = 1.5;	//reaction 'Pyruvate kinase':  kinetic parameter 'KmPYKATP'
p[68] = 0.14;	//reaction 'Pyruvate kinase':  kinetic parameter 'KmPYKPEP'
p[69] = 21;	//reaction 'Pyruvate kinase':  kinetic parameter 'KmPYKPYR'
p[70] = 1088.71;	//reaction 'Pyruvate kinase':  kinetic parameter 'VmPYK'
p[71] = 4.33;	//reaction 'Pyruvate decarboxylase':  kinetic parameter 'KmPDCPYR'
p[72] = 174.194;	//reaction 'Pyruvate decarboxylase':  kinetic parameter 'VmPDC'
p[73] = 1.9;	//reaction 'Pyruvate decarboxylase':  kinetic parameter 'nPDC'
p[74] = 21.4;	//reaction 'Succinate synthesis':  kinetic parameter 'KSUCC'
p[75] = 1;	//reaction 'Glucose transport':  kinetic parameter 'KeqGLT'
p[76] = 1.1918;	//reaction 'Glucose transport':  kinetic parameter 'KmGLTGLCi'
p[77] = 1.1918;	//reaction 'Glucose transport':  kinetic parameter 'KmGLTGLCo'
p[78] = 97.264;	//reaction 'Glucose transport':  kinetic parameter 'VmGLT'
p[79] = 6.9e-05;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KeqADH'
p[80] = 1.1;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KiADHACE'
p[81] = 90;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KiADHETOH'
p[82] = 0.92;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KiADHNAD'
p[83] = 0.031;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KiADHNADH'
p[84] = 1.11;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KmADHACE'
p[85] = 17;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KmADHETOH'
p[86] = 0.17;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KmADHNAD'
p[87] = 0.11;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'KmADHNADH'
p[88] = 810;	//reaction 'Alcohol dehydrogenase':  kinetic parameter 'VmADH'
p[89] = 4300;	//reaction 'Glycerol 3-phosphate dehydrogenase':  kinetic parameter 'KeqG3PDH'
p[90] = 0.4;	//reaction 'Glycerol 3-phosphate dehydrogenase':  kinetic parameter 'KmG3PDHDHAP'
p[91] = 1;	//reaction 'Glycerol 3-phosphate dehydrogenase':  kinetic parameter 'KmG3PDHGLY'
p[92] = 0.93;	//reaction 'Glycerol 3-phosphate dehydrogenase':  kinetic parameter 'KmG3PDHNAD'
p[93] = 0.023;	//reaction 'Glycerol 3-phosphate dehydrogenase':  kinetic parameter 'KmG3PDHNADH'
p[94] = 70.15;	//reaction 'Glycerol 3-phosphate dehydrogenase':  kinetic parameter 'VmG3PDH'
p[95] = 33.7;	//reaction 'ATPase activity':  kinetic parameter 'KATPASE'
#endif /* FIXED */

#ifdef ASSIGNMENT
y[0] = ct[0]-x[3];	//metabolite 'NADH': reactions
y[1] = (x_c[0]-y_c[2])/2.00000000000000000 * p[10];	//model entity 'ATP concentration':assignment
y[2] = (p_c[7]-pow((pow(x_c[0],2.00000000000000000)*(1.00000000000000000-4.00000000000000000*p[24])+2.00000000000000000*p_c[7]*x_c[0]*(4.00000000000000000*p[24]-1.00000000000000000)+pow(p_c[7],2.00000000000000000)),0.50000000000000000))/(1.00000000000000000-4.00000000000000000*p[24]) * p[10];	//model entity 'ADP concentration':assignment
y[3] = p_c[7]-y_c[1]-y_c[2] * p[10];	//model entity 'AMP concentration':assignment
x_c[0] = x[0]/p[10];	//concentration of metabolite 'High energy phosphates': reactions
x_c[1] = x[1]/p[10];	//concentration of metabolite 'Glucose 6 Phosphate': reactions
x_c[2] = x[2]/p[10];	//concentration of metabolite 'Triose-phosphate': reactions
x_c[3] = x[3]/p[10];	//concentration of metabolite 'NAD': reactions
x_c[4] = x[4]/p[10];	//concentration of metabolite 'Acetaldehyde': reactions
x_c[5] = x[5]/p[10];	//concentration of metabolite '2-phosphoglycerate': reactions
x_c[6] = x[6]/p[10];	//concentration of metabolite '1,3-bisphosphoglycerate': reactions
x_c[7] = x[7]/p[10];	//concentration of metabolite 'Fructose 6 Phosphate': reactions
x_c[8] = x[8]/p[10];	//concentration of metabolite 'Glucose in Cytosol': reactions
x_c[9] = x[9]/p[10];	//concentration of metabolite 'Phosphoenolpyruvate': reactions
x_c[10] = x[10]/p[10];	//concentration of metabolite 'Pyruvate': reactions
x_c[11] = x[11]/p[10];	//concentration of metabolite 'Fructose-1,6 bisphosphate': reactions
x_c[12] = x[12]/p[10];	//concentration of metabolite '3-phosphoglycerate': reactions
y_c[0] = y[0]/p[10];	//concentration of metabolite 'NADH': reactions
y_c[1] = y[1]/p[10];	//concentration of metabolite 'ATP concentration': assignment
y_c[2] = y[2]/p[10];	//concentration of metabolite 'ADP concentration': assignment
y_c[3] = y[3]/p[10];	//concentration of metabolite 'AMP concentration': assignment
p_c[0] = p[0]/p[9];	//concentration of metabolite 'Extracellular Glucose': fixed
p_c[1] = p[1]/p[10];	//concentration of metabolite 'Glycogen': fixed
p_c[2] = p[2]/p[10];	//concentration of metabolite 'Trehalose': fixed
p_c[3] = p[3]/p[10];	//concentration of metabolite 'CO2': fixed
p_c[4] = p[4]/p[10];	//concentration of metabolite 'Succinate': fixed
p_c[5] = p[5]/p[10];	//concentration of metabolite 'Ethanol': fixed
p_c[6] = p[6]/p[10];	//concentration of metabolite 'Glycerol': fixed
p_c[7] = p[7]/p[10];	//concentration of metabolite 'sum of AXP conc': fixed
p_c[8] = p[8]/p[10];	//concentration of metabolite 'F2,6P': fixed
#endif /* ASSIGNMENT */

#ifdef FUNCTIONS_HEADERS
double FunctionForHexokinase(double modif_0, double modif_1, double prod_0, double sub_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5); 
double FunctionForGlucose_6_phosphateIsomerase(double prod_0, double sub_0, double param_0, double param_1, double param_2, double param_3); 
double ConstantFlux_irreversible_(double param_0); 
double R_PFK(double varb_0, double varb_1, double varb_2, double varb_3, double varb_4); 
double L_PFK(double varb_0, double varb_1, double varb_2, double varb_3, double varb_4, double varb_5, double varb_6, double varb_7, double varb_8, double varb_9, double varb_10, double varb_11, double varb_12); 
double T_PFK(double varb_0, double varb_1, double varb_2); 
double FunctionForPhosphofructokinase(double modif_0, double modif_1, double param_0, double param_1, double param_2, double param_3, double param_4, double prod_0, double modif_2, double sub_0, double param_5, double param_6, double param_7, double param_8, double param_9, double param_10, double param_11, double param_12, double param_13); 
double FunctionForAldolase(double sub_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5, double prod_0, double param_6); 
double FunctionForGlyceraldehyde3_phosphateDehydrogenase(double prod_0, double param_0, double param_1, double param_2, double param_3, double param_4, double sub_0, double prod_1, double sub_1, double param_5, double param_6); 
double FunctionForPhosphoglycerateKinase(double modif_0, double modif_1, double sub_0, double param_0, double param_1, double param_2, double param_3, double param_4, double prod_0, double param_5); 
double FunctionForPhosphoglycerateMutase(double param_0, double param_1, double param_2, double prod_0, double sub_0, double param_3); 
double FunctionForEnolase(double param_0, double param_1, double param_2, double sub_0, double prod_0, double param_3); 
double FunctionForPyruvateKinase(double modif_0, double modif_1, double param_0, double param_1, double param_2, double param_3, double param_4, double sub_0, double prod_0, double param_5); 
double FunctionForPyruvateDecarboxylase(double param_0, double sub_0, double param_1, double param_2); 
double FunctionForSuccinateSynthesis(double sub_0, double param_0); 
double FunctionForGlucoseTransport(double prod_0, double sub_0, double param_0, double param_1, double param_2, double param_3); 
double FunctionForAlcoholDehydrogenase(double sub_0, double prod_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5, double param_6, double param_7, double param_8, double prod_1, double sub_1, double param_9, double volume_0); 
double FunctionForGlycerol3_phosphateDehydrogenase(double prod_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5, double prod_1, double sub_0, double sub_1, double param_6); 
double FunctionForATPaseActivity(double modif_0, double param_0); 
#endif /* FUNCTIONS_HEADERS */

#ifdef FUNCTIONS
double FunctionForHexokinase(double modif_0, double modif_1, double prod_0, double sub_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5) 	//Function for Hexokinase
{return  param_5/(param_4*param_2)*(sub_0*modif_1-prod_0*modif_0/param_0)/((1.00000000000000000+sub_0/param_4+prod_0/param_3)*(1.00000000000000000+modif_1/param_2+modif_0/param_1));} 
double FunctionForGlucose_6_phosphateIsomerase(double prod_0, double sub_0, double param_0, double param_1, double param_2, double param_3) 	//Function for Glucose-6-phosphate isomerase
{return  param_3/param_2*(sub_0-prod_0/param_0)/(1.00000000000000000+sub_0/param_2+prod_0/param_1);} 
double ConstantFlux_irreversible_(double param_0) 	//Constant flux (irreversible)
{return  param_0;} 
double R_PFK(double varb_0, double varb_1, double varb_2, double varb_3, double varb_4) 	//R_PFK
{return  1.00000000000000000+varb_4/varb_0+varb_3/varb_1+varb_2*(varb_4/varb_0)*(varb_3/varb_1);} 
double L_PFK(double varb_0, double varb_1, double varb_2, double varb_3, double varb_4, double varb_5, double varb_6, double varb_7, double varb_8, double varb_9, double varb_10, double varb_11, double varb_12) 	//L_PFK
{return  varb_0*pow(((1.00000000000000000+varb_1*(varb_9/varb_2))/(1.00000000000000000+varb_9/varb_2)),2.00000000000000000)*pow(((1.00000000000000000+varb_3*(varb_10/varb_4))/(1.00000000000000000+varb_10/varb_4)),2.00000000000000000)*pow(((1.00000000000000000+varb_5*varb_12/varb_6+varb_7*varb_11/varb_8)/(1.00000000000000000+varb_12/varb_6+varb_11/varb_8)),2.00000000000000000);} 
double T_PFK(double varb_0, double varb_1, double varb_2) 	//T_PFK
{return  1.00000000000000000+varb_0*(varb_2/varb_1);} 
double FunctionForPhosphofructokinase(double modif_0, double modif_1, double param_0, double param_1, double param_2, double param_3, double param_4, double prod_0, double modif_2, double sub_0, double param_5, double param_6, double param_7, double param_8, double param_9, double param_10, double param_11, double param_12, double param_13) 	//Function for Phosphofructokinase
{return  param_12*param_13*(sub_0/param_10)*(modif_1/param_9)*R_PFK(param_10,param_9,param_13,modif_1,sub_0)/(pow(R_PFK(param_10,param_9,param_13,modif_1,sub_0),2.00000000000000000)+L_PFK(param_11,param_4,param_8,param_0,param_5,param_3,param_7,param_2,param_6,modif_1,modif_0,prod_0,modif_2)*pow(T_PFK(param_1,param_9,modif_1),2.00000000000000000));} 
double FunctionForAldolase(double sub_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5, double prod_0, double param_6) 	//Function for Aldolase
{return  param_6/param_3*(sub_0-param_1/(1.00000000000000000+param_1)*prod_0*(1.00000000000000000/(1.00000000000000000+param_1))*prod_0/param_0)/(1.00000000000000000+sub_0/param_3+param_1/(1.00000000000000000+param_1)*prod_0/param_4+1.00000000000000000/(1.00000000000000000+param_1)*prod_0/param_2+param_1/(1.00000000000000000+param_1)*prod_0*(1.00000000000000000/(1.00000000000000000+param_1))*prod_0/(param_4*param_2)+sub_0*(param_1/(1.00000000000000000+param_1))*prod_0/(param_5*param_3));} 
double FunctionForGlyceraldehyde3_phosphateDehydrogenase(double prod_0, double param_0, double param_1, double param_2, double param_3, double param_4, double sub_0, double prod_1, double sub_1, double param_5, double param_6) 	//Function for Glyceraldehyde 3-phosphate dehydrogenase
{return  (param_5*(param_0/(1.00000000000000000+param_0))*sub_1*sub_0/(param_2*param_3)-param_6*prod_0*prod_1/(param_1*param_4))/((1.00000000000000000+param_0/(1.00000000000000000+param_0)*sub_1/param_2+prod_0/param_1)*(1.00000000000000000+sub_0/param_3+prod_1/param_4));} 
double FunctionForPhosphoglycerateKinase(double modif_0, double modif_1, double sub_0, double param_0, double param_1, double param_2, double param_3, double param_4, double prod_0, double param_5) 	//Function for Phosphoglycerate kinase
{return  param_5/(param_4*param_2)*(param_0*sub_0*modif_0-prod_0*modif_1)/((1.00000000000000000+sub_0/param_3+prod_0/param_4)*(1.00000000000000000+modif_1/param_2+modif_0/param_1));} 
double FunctionForPhosphoglycerateMutase(double param_0, double param_1, double param_2, double prod_0, double sub_0, double param_3) 	//Function for Phosphoglycerate mutase
{return  param_3/param_2*(sub_0-prod_0/param_0)/(1.00000000000000000+sub_0/param_2+prod_0/param_1);} 
double FunctionForEnolase(double param_0, double param_1, double param_2, double sub_0, double prod_0, double param_3) 	//Function for Enolase
{return  param_3/param_1*(sub_0-prod_0/param_0)/(1.00000000000000000+sub_0/param_1+prod_0/param_2);} 
double FunctionForPyruvateKinase(double modif_0, double modif_1, double param_0, double param_1, double param_2, double param_3, double param_4, double sub_0, double prod_0, double param_5) 	//Function for Pyruvate kinase
{return  param_5/(param_3*param_1)*(sub_0*modif_0-prod_0*modif_1/param_0)/((1.00000000000000000+sub_0/param_3+prod_0/param_4)*(1.00000000000000000+modif_1/param_2+modif_0/param_1));} 
double FunctionForPyruvateDecarboxylase(double param_0, double sub_0, double param_1, double param_2) 	//Function for Pyruvate decarboxylase
{return  param_1*(pow(sub_0,param_2)/pow(param_0,param_2))/(1.00000000000000000+pow(sub_0,param_2)/pow(param_0,param_2));} 
double FunctionForSuccinateSynthesis(double sub_0, double param_0) 	//Function for Succinate synthesis
{return  param_0*sub_0;} 
double FunctionForGlucoseTransport(double prod_0, double sub_0, double param_0, double param_1, double param_2, double param_3) 	//Function for Glucose transport
{return  param_3/param_2*(sub_0-prod_0/param_0)/(1.00000000000000000+sub_0/param_2+prod_0/param_1+0.91000000000000003*sub_0*prod_0/(param_2*param_1));} 
double FunctionForAlcoholDehydrogenase(double sub_0, double prod_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5, double param_6, double param_7, double param_8, double prod_1, double sub_1, double param_9, double volume_0) 	//Function for Alcohol dehydrogenase
{return  (-volume_0)*(param_9/(param_3*param_6)*(prod_1*prod_0-sub_1*sub_0/param_0)/(1.00000000000000000+prod_1/param_3+param_7*prod_0/(param_3*param_6)+param_8*sub_0/(param_4*param_5)+sub_1/param_4+prod_1*prod_0/(param_3*param_6)+param_8*prod_1*sub_0/(param_3*param_4*param_5)+param_7*prod_0*sub_1/(param_3*param_6*param_4)+sub_1*sub_0/(param_4*param_5)+prod_1*prod_0*sub_0/(param_3*param_6*param_1)+prod_0*sub_1*sub_0/(param_2*param_4*param_5)))/volume_0;} 
double FunctionForGlycerol3_phosphateDehydrogenase(double prod_0, double param_0, double param_1, double param_2, double param_3, double param_4, double param_5, double prod_1, double sub_0, double sub_1, double param_6) 	//Function for Glycerol 3-phosphate dehydrogenase
{return  param_6/(param_2*param_5)*(1.00000000000000000/(1.00000000000000000+param_1)*sub_1*sub_0-prod_0*prod_1/param_0)/((1.00000000000000000+1.00000000000000000/(1.00000000000000000+param_1)*sub_1/param_2+prod_0/param_3)*(1.00000000000000000+sub_0/param_5+prod_1/param_4));} 
double FunctionForATPaseActivity(double modif_0, double param_0) 	//Function for ATPase activity
{return  param_0*modif_0;} 
#endif /* FUNCTIONS */

#ifdef ODEs
dx[0] = -FunctionForHexokinase(y_c[2], y_c[1], x_c[1], x_c[8], p[26], p[27], p[28], p[29], p[30], p[31])*p[10]-ConstantFlux_irreversible_(p[36])*p[10]-ConstantFlux_irreversible_(p[37])*p[10]-FunctionForPhosphofructokinase(y_c[3], y_c[1], p[17], p[23], p[21], p[19], p[15], x_c[11], p_c[8], x_c[7], p[18], p[22], p[20], p[16], p[13], p[12], p[14], p[38], p[11])*p[10]+FunctionForPhosphoglycerateKinase(y_c[2], y_c[1], x_c[6], p[51], p[52], p[53], p[54], p[55], x_c[12], p[56])*p[10]+FunctionForPyruvateKinase(y_c[2], y_c[1], p[65], p[66], p[67], p[68], p[69], x_c[9], x_c[10], p[70])*p[10]-4*FunctionForSuccinateSynthesis(x_c[4], p[74])*p[10]-FunctionForATPaseActivity(y_c[1], p[95])*p[10];
dx[1] = FunctionForHexokinase(y_c[2], y_c[1], x_c[1], x_c[8], p[26], p[27], p[28], p[29], p[30], p[31])*p[10]-FunctionForGlucose_6_phosphateIsomerase(x_c[7], x_c[1], p[32], p[33], p[34], p[35])*p[10]-ConstantFlux_irreversible_(p[36])*p[10]-2*ConstantFlux_irreversible_(p[37])*p[10];
dx[2] = 2*FunctionForAldolase(x_c[11], p[39], p[25], p[40], p[41], p[42], p[43], x_c[2], p[44])*p[10]-FunctionForGlyceraldehyde3_phosphateDehydrogenase(x_c[6], p[25], p[45], p[46], p[47], p[48], x_c[3], y_c[0], x_c[2], p[49], p[50])*p[10]-FunctionForGlycerol3_phosphateDehydrogenase(p_c[6], p[89], p[25], p[90], p[91], p[92], p[93], x_c[3], y_c[0], x_c[2], p[94])*p[10];
dx[3] = -FunctionForGlyceraldehyde3_phosphateDehydrogenase(x_c[6], p[25], p[45], p[46], p[47], p[48], x_c[3], y_c[0], x_c[2], p[49], p[50])*p[10]-3*FunctionForSuccinateSynthesis(x_c[4], p[74])*p[10]+FunctionForAlcoholDehydrogenase(x_c[4], p_c[5], p[79], p[80], p[81], p[82], p[83], p[84], p[85], p[86], p[87], x_c[3], y_c[0], p[88], p[10])*p[10]+FunctionForGlycerol3_phosphateDehydrogenase(p_c[6], p[89], p[25], p[90], p[91], p[92], p[93], x_c[3], y_c[0], x_c[2], p[94])*p[10];
dx[4] = FunctionForPyruvateDecarboxylase(p[71], x_c[10], p[72], p[73])*p[10]-2*FunctionForSuccinateSynthesis(x_c[4], p[74])*p[10]-FunctionForAlcoholDehydrogenase(x_c[4], p_c[5], p[79], p[80], p[81], p[82], p[83], p[84], p[85], p[86], p[87], x_c[3], y_c[0], p[88], p[10])*p[10];
dx[5] = FunctionForPhosphoglycerateMutase(p[57], p[58], p[59], x_c[5], x_c[12], p[60])*p[10]-FunctionForEnolase(p[61], p[62], p[63], x_c[5], x_c[9], p[64])*p[10];
dx[6] = FunctionForGlyceraldehyde3_phosphateDehydrogenase(x_c[6], p[25], p[45], p[46], p[47], p[48], x_c[3], y_c[0], x_c[2], p[49], p[50])*p[10]-FunctionForPhosphoglycerateKinase(y_c[2], y_c[1], x_c[6], p[51], p[52], p[53], p[54], p[55], x_c[12], p[56])*p[10];
dx[7] = FunctionForGlucose_6_phosphateIsomerase(x_c[7], x_c[1], p[32], p[33], p[34], p[35])*p[10]-FunctionForPhosphofructokinase(y_c[3], y_c[1], p[17], p[23], p[21], p[19], p[15], x_c[11], p_c[8], x_c[7], p[18], p[22], p[20], p[16], p[13], p[12], p[14], p[38], p[11])*p[10];
dx[8] = -FunctionForHexokinase(y_c[2], y_c[1], x_c[1], x_c[8], p[26], p[27], p[28], p[29], p[30], p[31])*p[10]+FunctionForGlucoseTransport(x_c[8], p_c[0], p[75], p[76], p[77], p[78]);
dx[9] = FunctionForEnolase(p[61], p[62], p[63], x_c[5], x_c[9], p[64])*p[10]-FunctionForPyruvateKinase(y_c[2], y_c[1], p[65], p[66], p[67], p[68], p[69], x_c[9], x_c[10], p[70])*p[10];
dx[10] = FunctionForPyruvateKinase(y_c[2], y_c[1], p[65], p[66], p[67], p[68], p[69], x_c[9], x_c[10], p[70])*p[10]-FunctionForPyruvateDecarboxylase(p[71], x_c[10], p[72], p[73])*p[10];
dx[11] = FunctionForPhosphofructokinase(y_c[3], y_c[1], p[17], p[23], p[21], p[19], p[15], x_c[11], p_c[8], x_c[7], p[18], p[22], p[20], p[16], p[13], p[12], p[14], p[38], p[11])*p[10]-FunctionForAldolase(x_c[11], p[39], p[25], p[40], p[41], p[42], p[43], x_c[2], p[44])*p[10];
dx[12] = FunctionForPhosphoglycerateKinase(y_c[2], y_c[1], x_c[6], p[51], p[52], p[53], p[54], p[55], x_c[12], p[56])*p[10]-FunctionForPhosphoglycerateMutase(p[57], p[58], p[59], x_c[5], x_c[12], p[60])*p[10];
#endif /* ODEs */
