function xdot = Teusink2000_Glycolysis(time,x)
%  synopsis:
%     xdot = Teusink2000_Glycolysis (time, x)
%     x0 = Teusink2000_Glycolysis
%
%  Teusink2000_Glycolysis can be used with odeN functions as follows:
%
%  x0 = Teusink2000_Glycolysis;
%  [t,x] = ode23s(@Teusink2000_Glycolysis, [0 100], Teusink2000_Glycolysis);
%  plot (t,x);
%
%  where 100 is the end time
%
%  When Teusink2000_Glycolysis is used without any arguments it returns a vector of
%  the initial concentrations of the 19 floating species.
%  Otherwise Teusink2000_Glycolysis should be called with two arguments:
%  time and x.  time is the current time. x is the vector of the
%  concentrations of the 19 floating species.
%  When these parameters are supplied Teusink2000_Glycolysis returns a vector of
%  the rates of change of the concentrations of the 19 floating species.
%
%  the following table shows the mapping between the vector
%  index of a floating species and the species name.
%
%  NOTE for compartmental models
%  matlab translator generates code that when simulated in matlab,
%  produces results which have the units of species amounts. Users
%  should divide the results for each species with the volume of the
%  compartment it resides in, in order to obtain concentrations.
%
%  Indx      Name
%  x(1)        GLCi
%  x(2)        G6P
%  x(3)        F6P
%  x(4)        F16P
%  x(5)        TRIO
%  x(6)        BPG
%  x(7)        P3G
%  x(8)        P2G
%  x(9)        PEP
%  x(10)        PYR
%  x(11)        ACE
%  x(12)        P
%  x(13)        NAD
%  x(14)        NADH
%  x(15)        ATP
%  x(16)        ADP
%  x(17)        AMP
%  x(18)        SUM_P
%  x(19)        F26BP

xdot = zeros(19, 1);

% List of Compartments
vol__extracellular = 1;		%extracellular
vol__cytosol = 1;		%cytosol

% Global Parameters
g_p1 = 5.12;		% gR
g_p2 = 0.1;		    % KmPFKF6P
g_p3 = 0.71;		% KmPFKATP
g_p4 = 0.66;		% Lzero
g_p5 = 100;		    % CiPFKATP
g_p6 = 0.65;		% KiPFKATP
g_p7 = 0.0845;		% CPFKAMP
g_p8 = 0.0995;		% KPFKAMP
g_p9 = 0.0174;		% CPFKF26BP
g_p10 = 0.000682;	% KPFKF26BP
g_p11 = 0.397;		% CPFKF16BP
g_p12 = 0.111;		% KPFKF16BP
g_p13 = 3;		    % CPFKATP
g_p14 = 0.45;		% KeqAK
g_p15 = 0.045;		% KeqTPI

% Boundary Conditions
g_p16 = 0;		% Glyc = Glycogen[Concentration]
g_p17 = 0;		% Trh = Trehalose[Concentration]
g_p18 = 1;		% CO2 = CO2[Concentration]
g_p19 = 0;		% SUCC = Succinate[Concentration]
g_p20 = 50;		% GLCo = Extracellular Glucose[Concentration]
g_p21 = 50;		% ETOH = Ethanol[Concentration]
g_p22 = 0.15;	% GLY = Glycerol[Concentration]

% Local Parameters
par_r1_p1 = 226.452;	% [vGLK, VmGLK]
par_r1_p2 = 0.08;		% [vGLK, KmGLKGLCi]
par_r1_p3 = 0.15;		% [vGLK, KmGLKATP]
par_r1_p4 = 3800;		% [vGLK, KeqGLK]
par_r1_p5 = 30;		    % [vGLK, KmGLKG6P]
par_r1_p6 = 0.23;		% [vGLK, KmGLKADP]
par_r2_p1 = 339.677;	% [vPGI, VmPGI_2]
par_r2_p2 = 1.4;		% [vPGI, KmPGIG6P_2]
par_r2_p3 = 0.314;		% [vPGI, KeqPGI_2]
par_r2_p4 = 0.3;		% [vPGI, KmPGIF6P_2]
par_r3_p1 = 6;		    % [vGLYCO, KGLYCOGEN_3]
par_r4_p1 = 2.4;		% [vTreha, KTREHALOSE]
par_r5_p1 = 182.903;	% [vPFK, VmPFK]
par_r6_p1 = 322.258;	% [vALD, VmALD]
par_r6_p2 = 0.3;		% [vALD, KmALDF16P]
par_r6_p3 = 0.069;		% [vALD, KeqALD]
par_r6_p4 = 2;		    % [vALD, KmALDGAP]
par_r6_p5 = 2.4;		% [vALD, KmALDDHAP]
par_r6_p6 = 10;		    % [vALD, KmALDGAPi]
par_r7_p1 = 1184.52;	% [vGAPDH, VmGAPDHf]
par_r7_p2 = 0.21;		% [vGAPDH, KmGAPDHGAP]
par_r7_p3 = 0.09;		% [vGAPDH, KmGAPDHNAD]
par_r7_p4 = 6549.8;		% [vGAPDH, VmGAPDHr]
par_r7_p5 = 0.0098;		% [vGAPDH, KmGAPDHBPG]
par_r7_p6 = 0.06;		% [vGAPDH, KmGAPDHNADH]
par_r8_p1 = 1306.45;	% [vPGK, VmPGK]
par_r8_p2 = 0.53;		% [vPGK, KmPGKP3G]
par_r8_p3 = 0.3;		% [vPGK, KmPGKATP]
par_r8_p4 = 3200;		% [vPGK, KeqPGK]
par_r8_p5 = 0.003;		% [vPGK, KmPGKBPG]
par_r8_p6 = 0.2;		% [vPGK, KmPGKADP]
par_r9_p1 = 2525.81;	% [vPGM, VmPGM]
par_r9_p2 = 1.2;		% [vPGM, KmPGMP3G]
par_r9_p3 = 0.19;		% [vPGM, KeqPGM]
par_r9_p4 = 0.08;		% [vPGM, KmPGMP2G]
par_r10_p1 = 365.806;	% [vENO, VmENO]
par_r10_p2 = 0.04;		% [vENO, KmENOP2G]
par_r10_p3 = 6.7;		% [vENO, KeqENO]
par_r10_p4 = 0.5;		% [vENO, KmENOPEP]
par_r11_p1 = 1088.71;	% [vPYK, VmPYK]
par_r11_p2 = 0.14;		% [vPYK, KmPYKPEP]
par_r11_p3 = 0.53;		% [vPYK, KmPYKADP]
par_r11_p4 = 6500;		% [vPYK, KeqPYK]
par_r11_p5 = 21;		% [vPYK, KmPYKPYR]
par_r11_p6 = 1.5;		% [vPYK, KmPYKATP]
par_r12_p1 = 174.194;	% [vPDC, VmPDC]
par_r12_p2 = 1.9;		% [vPDC, nPDC]
par_r12_p3 = 4.33;		% [vPDC, KmPDCPYR]
par_r13_p1 = 21.4;		% [vSUC, KSUCC]
par_r14_p1 = 97.264;	% [vGLT, VmGLT]
par_r14_p2 = 1.1918;	% [vGLT, KmGLTGLCo]
par_r14_p3 = 1;		    % [vGLT, KeqGLT]
par_r14_p4 = 1.1918;	% [vGLT, KmGLTGLCi]
par_r15_p1 = 810;		% [vADH, VmADH]
par_r15_p2 = 0.92;		% [vADH, KiADHNAD]
par_r15_p3 = 17;		% [vADH, KmADHETOH]
par_r15_p4 = 6.9e-005;	% [vADH, KeqADH]
par_r15_p5 = 0.17;		% [vADH, KmADHNAD]
par_r15_p6 = 0.11;		% [vADH, KmADHNADH]
par_r15_p7 = 0.031;		% [vADH, KiADHNADH]
par_r15_p8 = 1.11;		% [vADH, KmADHACE]
par_r15_p9 = 1.1;		% [vADH, KiADHACE]
par_r15_p10 = 90;		% [vADH, KiADHETOH]
par_r16_p1 = 70.15;		% [vG3PDH, VmG3PDH]
par_r16_p2 = 0.4;		% [vG3PDH, KmG3PDHDHAP]
par_r16_p3 = 0.023;		% [vG3PDH, KmG3PDHNADH]
par_r16_p4 = 4300;		% [vG3PDH, KeqG3PDH]
par_r16_p5 = 1;		    % [vG3PDH, KmG3PDHGLY]
par_r16_p6 = 0.93;		% [vG3PDH, KmG3PDHNAD]
par_r17_p1 = 33.7;		% [vATP, KATPASE]

if (nargin == 0)

    % set initial conditions
   xdot(1) = 0.087*vol__cytosol;		% GLCi = Glucose in Cytosol [Concentration]
   xdot(2) = 2.45*vol__cytosol;		    % G6P = Glucose 6 Phosphate [Concentration]
   xdot(3) = 0.62*vol__cytosol;		    % F6P = Fructose 6 Phosphate [Concentration]
   xdot(4) = 5.51*vol__cytosol;		    % F16P = Fructose-1,6 bisphosphate [Concentration]
   xdot(5) = 0.96*vol__cytosol;		    % TRIO = Triose-phosphate [Concentration]
   xdot(6) = 0*vol__cytosol;		    % BPG = 1,3-bisphosphoglycerate [Concentration]
   xdot(7) = 0.9*vol__cytosol;		    % P3G = 3-phosphoglycerate [Concentration]
   xdot(8) = 0.12*vol__cytosol;		    % P2G = 2-phosphoglycerate [Concentration]
   xdot(9) = 0.07*vol__cytosol;		    % PEP = Phosphoenolpyruvate [Concentration]
   xdot(10) = 1.85*vol__cytosol;		% PYR = Pyruvate [Concentration]
   xdot(11) = 0.17*vol__cytosol;		% ACE = Acetaldehyde [Concentration]
   xdot(12) = 6.31*vol__cytosol;		% P = High energy phosphates [Concentration]
   xdot(13) = 1.2*vol__cytosol;		    % NAD = NAD [Concentration]
   xdot(14) = 0.39*vol__cytosol;		% NADH = NADH [Concentration]
   xdot(15) = 0;		                % ATP = ATP concentration [Amount]
   xdot(16) = 0;		                % ADP = ADP concentration [Amount]
   xdot(17) = 0;		                % AMP = AMP concentration [Amount]
   xdot(18) = 4.1*vol__cytosol;		    % SUM_P = sum of AXP conc [Concentration]
   xdot(19) = 0.02*vol__cytosol;		% F26BP = F2,6P [Concentration]

else

    % listOfRules
   x(16) = ((x(18))-pow(pow((x(12)),2)*(1-4*g_p14)+2*(x(18))*(x(12))*(4*g_p14-1)+pow((x(18)),2),0.5))/(1-4*g_p14);
   x(15) = ((x(12))-(x(16)))/2;
   x(17) = (x(18))-(x(15))-(x(16));

    % calculate rates of change
   R0 = vol__cytosol*par_r1_p1/(par_r1_p2*par_r1_p3)*((x(1))*(x(15))-(x(2))*(x(16))/par_r1_p4)/((1+(x(1))/par_r1_p2+(x(2))/par_r1_p5)*(1+(x(15))/par_r1_p3+(x(16))/par_r1_p6));
   R1 = vol__cytosol*par_r2_p1/par_r2_p2*((x(2))-(x(3))/par_r2_p3)/(1+(x(2))/par_r2_p2+(x(3))/par_r2_p4);
   R2 = vol__cytosol*par_r3_p1;
   R3 = vol__cytosol*par_r4_p1;
   R4 = vol__cytosol*par_r5_p1*g_p1*((x(3))/g_p2)*((x(15))/g_p3)*R_PFK(g_p2,g_p3,g_p1,(x(15)),(x(3)))/(pow(R_PFK(g_p2,g_p3,g_p1,(x(15)),(x(3))),2)+L_PFK(g_p4,g_p5,g_p6,g_p7,g_p8,g_p9,g_p10,g_p11,g_p12,(x(15)),(x(17)),(x(4)),(x(19)))*pow(T_PFK(g_p13,g_p3,(x(15))),2));
   R5 = vol__cytosol*par_r6_p1/par_r6_p2*((x(4))-g_p15/(1+g_p15)*(x(5))*(1/(1+g_p15))*(x(5))/par_r6_p3)/(1+(x(4))/par_r6_p2+g_p15/(1+g_p15)*(x(5))/par_r6_p4+1/(1+g_p15)*(x(5))/par_r6_p5+g_p15/(1+g_p15)*(x(5))*(1/(1+g_p15))*(x(5))/(par_r6_p4*par_r6_p5)+(x(4))*(g_p15/(1+g_p15))*(x(5))/(par_r6_p6*par_r6_p2));
   R6 = vol__cytosol*(par_r7_p1*(g_p15/(1+g_p15))*(x(5))*(x(13))/(par_r7_p2*par_r7_p3)-par_r7_p4*(x(6))*(x(14))/(par_r7_p5*par_r7_p6))/((1+g_p15/(1+g_p15)*(x(5))/par_r7_p2+(x(6))/par_r7_p5)*(1+(x(13))/par_r7_p3+(x(14))/par_r7_p6));
   R7 = vol__cytosol*par_r8_p1/(par_r8_p2*par_r8_p3)*(par_r8_p4*(x(6))*(x(16))-(x(7))*(x(15)))/((1+(x(6))/par_r8_p5+(x(7))/par_r8_p2)*(1+(x(15))/par_r8_p3+(x(16))/par_r8_p6));
   R8 = vol__cytosol*par_r9_p1/par_r9_p2*((x(7))-(x(8))/par_r9_p3)/(1+(x(7))/par_r9_p2+(x(8))/par_r9_p4);
   R9 = vol__cytosol*par_r10_p1/par_r10_p2*((x(8))-(x(9))/par_r10_p3)/(1+(x(8))/par_r10_p2+(x(9))/par_r10_p4);
   R10 = vol__cytosol*par_r11_p1/(par_r11_p2*par_r11_p3)*((x(9))*(x(16))-(x(10))*(x(15))/par_r11_p4)/((1+(x(9))/par_r11_p2+(x(10))/par_r11_p5)*(1+(x(15))/par_r11_p6+(x(16))/par_r11_p3));
   R11 = vol__cytosol*par_r12_p1*(pow((x(10)),par_r12_p2)/pow(par_r12_p3,par_r12_p2))/(1+pow((x(10)),par_r12_p2)/pow(par_r12_p3,par_r12_p2));
   R12 = vol__cytosol*par_r13_p1*(x(11));
   R13 = par_r14_p1/par_r14_p2*((g_p20)-(x(1))/par_r14_p3)/(1+(g_p20)/par_r14_p2+(x(1))/par_r14_p4+0.91*(g_p20)*(x(1))/(par_r14_p2*par_r14_p4));
   R14 = -vol__cytosol*(par_r15_p1/(par_r15_p2*par_r15_p3)*((x(13))*(g_p21)-(x(14))*(x(11))/par_r15_p4)/(1+(x(13))/par_r15_p2+par_r15_p5*(g_p21)/(par_r15_p2*par_r15_p3)+par_r15_p6*(x(11))/(par_r15_p7*par_r15_p8)+(x(14))/par_r15_p7+(x(13))*(g_p21)/(par_r15_p2*par_r15_p3)+par_r15_p6*(x(13))*(x(11))/(par_r15_p2*par_r15_p7*par_r15_p8)+par_r15_p5*(g_p21)*(x(14))/(par_r15_p2*par_r15_p3*par_r15_p7)+(x(14))*(x(11))/(par_r15_p7*par_r15_p8)+(x(13))*(g_p21)*(x(11))/(par_r15_p2*par_r15_p3*par_r15_p9)+(g_p21)*(x(14))*(x(11))/(par_r15_p10*par_r15_p7*par_r15_p8)));
   R15 = vol__cytosol*par_r16_p1/(par_r16_p2*par_r16_p3)*(1/(1+g_p15)*(x(5))*(x(14))-(g_p22)*(x(13))/par_r16_p4)/((1+1/(1+g_p15)*(x(5))/par_r16_p2+(g_p22)/par_r16_p5)*(1+(x(14))/par_r16_p3+(x(13))/par_r16_p6));
   R16 = vol__cytosol*par_r17_p1*(x(15));

   xdot = [
      - R0 + R13
      + R0 - R1 - R2 - 2*R3
      + R1 - R4
      + R4 - R5
      + 2*R5 - R6 - R15
      + R6 - R7
      + R7 - R8
      + R8 - R9
      + R9 - R10
      + R10 - R11
      + R11 - 2*R12 - R14
      - R0 - R2 - R3 - R4 + R7 + R10 - 4*R12 - R16
      - R6 - 3*R12 + R14 + R15
      + R6 + 3*R12 - R14 - R15
        0
        0
        0
        0
        0
   ];
end;


function z = pow (x, y)
    z = x^y;


function z = sqr (x)
    z = x*x;


function z = piecewise(varargin)
		numArgs = nargin;
		result = 0;
		foundResult = 0;
		for k=1:2: numArgs-1
			if varargin{k+1} == 1
				result = varargin{k};
				foundResult = 1;
				break;
			end
		end
		if foundResult == 0
			result = varargin{numArgs};
		end
		z = result;


function z = gt(a,b)
   if a > b
   	  z = 1;
   else
      z = 0;
   end


function z = lt(a,b)
   if a < b
   	  z = 1;
   else
      z = 0;
   end


function z = geq(a,b)
   if a >= b
   	  z = 1;
   else
      z = 0;
   end


function z = leq(a,b)
   if a <= b
   	  z = 1;
   else
      z = 0;
   end


function z = neq(a,b)
   if a ~= b
   	  z = 1;
   else
      z = 0;
   end


function z = and(varargin)
		result = 1;
		for k=1:nargin
		   if varargin{k} ~= 1
		      result = 0;
		      break;
		   end
		end
		z = result;


function z = or(varargin)
		result = 0;
		for k=1:nargin
		   if varargin{k} ~= 0
		      result = 1;
		      break;
		   end
		end
		z = result;


function z = xor(varargin)
		foundZero = 0;
		foundOne = 0;
		for k = 1:nargin
			if varargin{k} == 0
			   foundZero = 1;
			else
			   foundOne = 1;
			end
		end
		if foundZero && foundOne
			z = 1;
		else
		  z = 0;
		end



function z = not(a)
   if a == 1
   	  z = 0;
   else
      z = 1;
   end


function z = root(a,b)
	z = a^(1/b);

