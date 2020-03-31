import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import solve_ivp
from collections import namedtuple

# set matplotlib backend
matplotlib.use('TkAgg')

# configure pandas to show full dataframe
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

# create named tuple type for parameters so we can refer to them by name and not a number in an array
param_named_tuple = namedtuple(
    'params',
    [
        # volumes
        'extracellular', 'cytosol',
        # kinetic parameters
        'KeqAK', 'KeqGLK', 'KmGLKADP', 'KmGLKATP', 'KmGLKG6P', 'KmGLKGLCi',
        'VmGLK', 'KeqPGI_2', 'KmPGIF6P_2', 'KmPGIG6P_2', 'VmPGI_2', 'vGLYCO_v', 'vTreha_v', 'CPFKAMP', 'CPFKATP',
        'CPFKF16BP', 'CPFKF26BP', 'CiPFKATP', 'KPFKAMP', 'KPFKF16BP', 'KPFKF26BP', 'KiPFKATP', 'KmPFKATP', 'KmPFKF6P',
        'Lzero', 'gR', 'VmPFK', 'KeqTPI', 'KeqALD', 'KmALDDHAP', 'KmALDF16P', 'KmALDGAP', 'KmALDGAPi', 'VmALD',
        'KmGAPDHBPG', 'KmGAPDHGAP', 'KmGAPDHNAD', 'KmGAPDHNADH', 'VmGAPDHf', 'VmGAPDHr', 'KeqPGK', 'KmPGKADP',
        'KmPGKATP', 'KmPGKBPG', 'KmPGKP3G', 'VmPGK', 'KeqPGM', 'KmPGMP2G', 'KmPGMP3G', 'VmPGM', 'KeqENO', 'KmENOP2G',
        'KmENOPEP', 'VmENO', 'KeqPYK', 'KmPYKADP', 'KmPYKATP', 'KmPYKPEP', 'KmPYKPYR', 'VmPYK', 'KmPDCPYR', 'VmPDC',
        'nPDC', 'KSUCC', 'KeqGLT', 'KmGLTGLCi', 'KmGLTGLCo', 'VmGLT', 'KeqADH', 'KiADHACE', 'KiADHETOH', 'KiADHNAD',
        'KiADHNADH', 'KmADHACE', 'KmADHETOH', 'KmADHNAD', 'KmADHNADH', 'VmADH', 'KeqG3PDH', 'KmG3PDHDHAP', 'KmG3PDHGLY',
        'KmG3PDHNAD', 'KmG3PDHNADH', 'VmG3PDH', 'KATPASE',
        # fixed parameters
        'Glyc', 'Trh', 'CO2', 'SUCC', 'ETOH', 'GLY', 'SUM_P', 'F26BP', 'GLCo',
    ])



def teusink2000(t, y0, p):
    """
    Callable function that can be integrated with scipy.integrate functions.

    Args:
        t (vector, list, np.array): integration time points
        y0 (vector, list, np.array): initial concentration parameters
        p (named tuple): all other model parameters

    Returns:

    """
    # Collect parameters again in a namedTuple
    p = param_named_tuple(*p)
    GLCi, G6P, F6P, F16P, TRIO, BPG, P3G, P2G, PEP, PYR, ACE, P, NAD, NADH = y0

    # assignment rules
    adp = (p.SUM_P - (P ** 2 * (1 - 4 * p.KeqAK)
                    + 2 * p.SUM_P * P * (4 * p.KeqAK - 1)
                    + p.SUM_P ** 2) ** 0.5) / (1 - 4 * p.KeqAK)
    atp = (P - adp) / 2
    amp = p.SUM_P - atp - adp

    # reaction functions
    def glycogen_synthesis():
        return p.vGLYCO_v

    def trehalose6p_synthesis():
        return p.vTreha_v

    def alcohol_dehydrogenase():
        numerator = (p.VmADH / (p.KiADHNAD * p.KmADHETOH)) * (NAD * p.ETOH - (NADH * ACE / p.KeqADH))
        denom = 1 + (NAD / p.KiADHNAD) \
                + p.KmADHNAD * p.ETOH / (p.KiADHNAD * p.KmADHETOH) \
                + p.KmADHNADH * ACE / (p.KiADHNADH * p.KmADHACE) \
                + NADH / p.KiADHNADH \
                + NAD * p.ETOH / (p.KiADHNAD * p.KmADHETOH) \
                + p.KmADHNADH * NAD * ACE / (p.KiADHNAD * p.KiADHNADH * p.KmADHACE) \
                + p.KmADHNAD * p.ETOH * NADH / (p.KiADHNAD * p.KmADHETOH * p.KiADHNADH) \
                + NADH * ACE / (p.KiADHNADH * p.KmADHACE) \
                + NAD * p.ETOH * ACE / (p.KiADHNAD * p.KmADHETOH * p.KiADHACE) \
                + p.ETOH * NADH * ACE / (p.KiADHETOH * p.KiADHNADH * p.KmADHACE)
        return -numerator / denom

    def glycerol_3_phosphate_dehydrogenase():
        numerator = (p.VmG3PDH / (p.KmG3PDHDHAP * p.KmG3PDHNADH)) * ((1 / (1 + p.KeqTPI)) * TRIO * NADH
                                                                     - p.GLY * NAD / p.KeqG3PDH)

        denom = (1 + (1 / (1 + p.KeqTPI)) * TRIO / p.KmG3PDHDHAP + p.GLY / p.KmG3PDHGLY) * (
                1 + NADH / p.KmG3PDHNADH + NAD / p.KmG3PDHNAD)
        return numerator / denom

    def r_pfk():
        return 1 + F6P / p.KmPFKF6P + atp / p.KmPFKATP + p.gR * (F6P / p.KmPFKF6P) * (atp / p.KmPFKATP)

    def t_pfk():
        return 1 + p.CPFKATP * (atp / p.KmPFKATP)

    def l_pfk():
        return p.Lzero * ((1 + p.CiPFKATP * (atp / p.KiPFKATP)) / (1 + atp / p.KiPFKATP)) ** 2 \
               * ((1 + p.CPFKAMP * (amp / p.KPFKAMP)) / (1 + amp / p.KPFKAMP)) ** 2 \
               * ((1 + p.CPFKF26BP * p.F26BP / p.KPFKF26BP + p.CPFKF16BP * F16P / p.KPFKF16BP) / (
                1 + p.F26BP / p.KPFKF26BP + F16P / p.KPFKF16BP)) ** 2

    def hexokinase():
        return (p.VmGLK / (p.KmGLKGLCi * p.KmGLKATP)) * (GLCi * atp - G6P * adp / p.KeqGLK) / (
                (1 + GLCi / p.KmGLKGLCi + G6P / p.KmGLKG6P) * (1 + atp / p.KmGLKATP + adp / p.KmGLKADP))

    def glucose_6_phosphate_isomerase():
        return (p.VmPGI_2 / p.KmPGIG6P_2) * (G6P - F6P / p.KeqPGI_2) / (
                1 + G6P / p.KmPGIG6P_2 + F6P / p.KmPGIF6P_2)

    def aldolase():
        return (p.VmALD / p.KmALDF16P) * (
                F16P - (p.KeqTPI / (1 + p.KeqTPI)) * TRIO * (1 / (1 + p.KeqTPI)) * TRIO / p.KeqALD) / (
                       1 + F16P / p.KmALDF16P + (p.KeqTPI / (1 + p.KeqTPI)) * TRIO / p.KmALDGAP + (
                       1 / (1 + p.KeqTPI)) * TRIO / p.KmALDDHAP + (p.KeqTPI / (1 + p.KeqTPI)) * TRIO * (
                               1 / (1 + p.KeqTPI)) * TRIO / (p.KmALDGAP * p.KmALDDHAP) + F16P * (
                               p.KeqTPI / (1 + p.KeqTPI)) * TRIO / (p.KmALDGAPi * p.KmALDF16P))

    def phosphoglycerate_kinase():
        return (p.VmPGK / (p.KmPGKP3G * p.KmPGKATP)) * (p.KeqPGK * BPG * adp - P3G * atp) / (
                (1 + BPG / p.KmPGKBPG + P3G / p.KmPGKP3G) * (1 + atp / p.KmPGKATP + adp / p.KmPGKADP))

    def glyceraldehyde_3_phosphate_dehydrogenase():
        numerator = p.VmGAPDHf * (p.KeqTPI / (1 + p.KeqTPI)) * TRIO * NAD / (p.KmGAPDHGAP * p.KmGAPDHNAD) \
                    - p.VmGAPDHr * BPG * NADH / (p.KmGAPDHBPG * p.KmGAPDHNADH)
        denom = (1 + (p.KeqTPI / (1 + p.KeqTPI)) * TRIO / p.KmGAPDHGAP + BPG / p.KmGAPDHBPG) \
                * (1 + NAD / p.KmGAPDHNAD + NADH / p.KmGAPDHNADH)
        return numerator / denom

    def enolase():
        return (p.VmENO / p.KmENOP2G) * (P2G - PEP / p.KeqENO) / (1 + P2G / p.KmENOP2G + PEP / p.KmENOPEP)

    def pyruvate_decarboxylase():
        x = PYR ** p.nPDC / (p.KmPDCPYR ** p.nPDC)
        return p.VmPDC * x / (1 + x)

    def succinate_synthesis():
        return p.KSUCC * ACE

    def glucose_transport():
        return (p.VmGLT / p.KmGLTGLCo) * (p.GLCo - GLCi / p.KeqGLT) / (
                1 + p.GLCo / p.KmGLTGLCo + GLCi / p.KmGLTGLCi + 0.91 * p.GLCo * GLCi / (p.KmGLTGLCo * p.KmGLTGLCi))

    def phosphoglycerate_mutase():
        return (p.VmPGM / p.KmPGMP3G) * (P3G - P2G / p.KeqPGM) / (1 + P3G / p.KmPGMP3G + P2G / p.KmPGMP2G)

    def pyruvate_kinase():
        return (p.VmPYK / (p.KmPYKPEP * p.KmPYKADP)) * (PEP * adp - PYR * atp / p.KeqPYK) / (
                (1 + PEP / p.KmPYKPEP + PYR / p.KmPYKPYR) * (1 + atp / p.KmPYKATP + adp / p.KmPYKADP))

    def atpase_activity():
        return p.KATPASE * atp

    def phosphofructokinase():
        return p.VmPFK * p.gR * (F6P / p.KmPFKF6P) * (atp / p.KmPFKATP) * r_pfk() / (
                r_pfk() ** 2 + l_pfk() * t_pfk() ** 2)

    # list of ODEs
    return [
        # GLCi, glucose in p.cytosol
        -p.cytosol * hexokinase()
        + glucose_transport(),

        # G6P, glucose 6 phosphase
        + p.cytosol * hexokinase()
        - p.cytosol * glucose_6_phosphate_isomerase()
        - p.cytosol * glycogen_synthesis()
        - 2 * p.cytosol * trehalose6p_synthesis(),

        # F6P, fructose 6 phosphate
        + p.cytosol * glucose_6_phosphate_isomerase()
        - p.cytosol * phosphofructokinase(),

        # F16BP, fructose-1-6-bisphosphate
        + p.cytosol * phosphofructokinase()
        - p.cytosol * aldolase(),

        # TRIO Triose-phosphate
        - p.cytosol * glycerol_3_phosphate_dehydrogenase()
        + 2 * p.cytosol * aldolase()
        - p.cytosol * glyceraldehyde_3_phosphate_dehydrogenase(),

        # BPG, 1,3-bisphosphoglycerate
        + p.cytosol * glyceraldehyde_3_phosphate_dehydrogenase()
        - p.cytosol * phosphoglycerate_kinase(),

        # P3G, 3-phosphoglycerate
        + p.cytosol * phosphoglycerate_kinase()
        - p.cytosol * phosphoglycerate_mutase(),

        # P2G, 2-phosphoglycerate
        + p.cytosol * phosphoglycerate_mutase()
        - p.cytosol * enolase(),

        # PEP, Phosphoenolpyruvate
        - p.cytosol * pyruvate_kinase()
        + p.cytosol * enolase(),

        # PYR, Pyruvate
        + p.cytosol * pyruvate_kinase()
        - p.cytosol * pyruvate_decarboxylase(),

        # ACE, Acetaldehyde
        + p.cytosol * pyruvate_decarboxylase()
        - 2 * p.cytosol * succinate_synthesis()
        - p.cytosol * alcohol_dehydrogenase(),

        # P high energy phosphates
        - p.cytosol * hexokinase()
        + p.cytosol * pyruvate_kinase()
        - 4 * p.cytosol * succinate_synthesis()
        - p.cytosol * atpase_activity()
        - p.cytosol * glycogen_synthesis()
        - p.cytosol * trehalose6p_synthesis()
        - p.cytosol * phosphofructokinase()
        + p.cytosol * phosphoglycerate_kinase(),

        # NAD,
        - 3 * p.cytosol * succinate_synthesis()
        + p.cytosol * alcohol_dehydrogenase()
        + p.cytosol * glycerol_3_phosphate_dehydrogenase()
        - p.cytosol * glyceraldehyde_3_phosphate_dehydrogenase(),

        # NADH,
        + 3 * p.cytosol * succinate_synthesis()
        - p.cytosol * alcohol_dehydrogenase()
        - p.cytosol * glycerol_3_phosphate_dehydrogenase()
        + p.cytosol * glyceraldehyde_3_phosphate_dehydrogenase()
    ]

# initial concentration parameters
GLCi = 0.087
G6P = 2.45
F6P = 0.62
F16P = 5.51
TRIO = 0.96
BPG = 0
P3G = 0.9
P2G = 0.12
PEP = 0.07
PYR = 1.85
ACE = 0.17
P = 6.31
NAD = 1.2
NADH = 0.39

# collect initial concentration parameters into list
y0 = [GLCi, G6P, F6P, F16P, TRIO, BPG, P3G, P2G, PEP, PYR, ACE, P, NAD, NADH]

# create instance of named tuple for kinetic parameters
parameters = param_named_tuple(
    # Compartment
    extracellular=1,
    cytosol=1,

    # parameters used in assignments
    KeqAK=0.45,

    # Hexokinase reaction parameters
    KeqGLK=3800,
    KmGLKADP=0.23,
    KmGLKATP=0.15,
    KmGLKG6P=30,
    KmGLKGLCi=0.08,
    VmGLK=226.452,

    # Glucose-6-phosphate isomerase parameters
    KeqPGI_2=0.314,
    KmPGIF6P_2=0.3,
    KmPGIG6P_2=1.4,
    VmPGI_2=339.677,

    # Glycogen synthesis parameters
    vGLYCO_v=6,

    # Trehalose 6-phosphate synthase parameters
    vTreha_v=2.4,

    # Phosphofructokinase parameters
    CPFKAMP=0.0845,
    CPFKATP=3,
    CPFKF16BP=0.397,
    CPFKF26BP=0.0174,
    CiPFKATP=100,
    KPFKAMP=0.0995,
    KPFKF16BP=0.111,
    KPFKF26BP=0.000682,
    KiPFKATP=0.65,
    KmPFKATP=0.71,
    KmPFKF6P=0.1,
    Lzero=0.66,
    gR=5.12,
    VmPFK=182.903,

    # Aldolase parameters
    KeqTPI=0.045,  # note parameter also used in Glycer 3p dehy
    KeqALD=0.069,
    KmALDDHAP=2.4,
    KmALDF16P=0.3,
    KmALDGAP=2,
    KmALDGAPi=10,
    VmALD=322.258,

    # Glyceraldehyde 3-phosphate dehydrogenase parameters
    KmGAPDHBPG=0.0098,
    KmGAPDHGAP=0.21,
    KmGAPDHNAD=0.09,
    KmGAPDHNADH=0.06,
    VmGAPDHf=1184.52,
    VmGAPDHr=6549.8,

    # Phosphoglycerate kinase parameters
    KeqPGK=3200,
    KmPGKADP=0.2,
    KmPGKATP=0.3,
    KmPGKBPG=0.003,
    KmPGKP3G=0.53,
    VmPGK=1306.45,

    # Phosphoglycerate mutase parameters
    KeqPGM=0.19,
    KmPGMP2G=0.08,
    KmPGMP3G=1.2,
    VmPGM=2525.81,

    # Enolase parameters
    KeqENO=6.7,
    KmENOP2G=0.04,
    KmENOPEP=0.5,
    VmENO=365.806,

    # Pyruvate kinase parameters
    KeqPYK=6500,
    KmPYKADP=0.53,
    KmPYKATP=1.5,
    KmPYKPEP=0.14,
    KmPYKPYR=21,
    VmPYK=1088.71,

    # Pyruvate decarboxylase parameters
    KmPDCPYR=4.33,
    VmPDC=174.194,
    nPDC=1.9,

    # Succinate synthesis parameters
    KSUCC=21.4,

    # Glucose transport parameters
    KeqGLT=1,
    KmGLTGLCi=1.1918,
    KmGLTGLCo=1.1918,
    VmGLT=97.264,

    # Alcohol dehydrogenase parameters
    KeqADH=6.9e-05,
    KiADHACE=1.1,
    KiADHETOH=90,
    KiADHNAD=0.92,
    KiADHNADH=0.031,
    KmADHACE=1.11,
    KmADHETOH=17,
    KmADHNAD=0.17,
    KmADHNADH=0.11,
    VmADH=810,

    # Glycerol 3-phosphate dehydrogenase parameters
    KeqG3PDH=4300,
    KmG3PDHDHAP=0.4,
    KmG3PDHGLY=1,
    KmG3PDHNAD=0.93,
    KmG3PDHNADH=0.023,
    VmG3PDH=70.15,

    # ATPase activity parameters
    KATPASE=33.7,

    # These metabolites are sinks. Not sure why CO2 has to be 1 though.
    # Should be 0? does it make a difference?
    Glyc=0,
    Trh=0,
    CO2=1,
    SUCC=0,
    ETOH=50,
    GLY=0.15,
    SUM_P=4.1,
    F26BP=0.02,
    GLCo=50
)
# list of column names
cols = ['GLCi','G6P', 'F6P', 'F16P', 'TRIO', 'BPG', 'P3G', 'P2G', 'PEP', 'PYR', 'ACE', 'P', 'NAD', 'NADH']

# integrate the model
data = solve_ivp(teusink2000, t_span=(0, 10), y0=y0, args=(parameters,), method="LSODA")
df = pd.DataFrame(data.y, columns=data.t, index=cols).transpose()

# plot data
for i in df.columns:
    plt.figure()
    plt.plot(df.index, df[i])
    plt.xlabel('Time')
    plt.ylabel(i)

# after 10 time steps the model is at steady state
print(df.iloc[-1])

'''
Steady state calculated by scipy integration
--------------------------------------------
GLCi    0.098759
G6P     1.033246
F6P     0.112813
F16P    0.601908
TRIO    0.777524
BPG     0.000330
P3G     0.356484
P2G     0.044844
PEP     0.073617
PYR     8.523153
ACE     0.170114
P       6.308882
NAD     1.545560
NADH    0.044440


Steady state calculated by tellurium
--------------------------------------
            0
GLCi  0.098759
G6P   1.033246
F6P   0.112813
F16P  0.601908
TRIO  0.777524
BPG   0.000330
P3G   0.356484
P2G   0.044844
PEP   0.073617
PYR   8.523152
ACE   0.170114
P     6.308882
NAD   1.545560
NADH  0.044440



Steady state calculated by copasi steady state task
---------------------------------------------------
GLCi  0.09875869199169003       
G6P   1.033245613681812         
F6P   0.1128128145855018        
F16P  0.6019076395836982        
TRIO  0.7775235367088307        
BPG   0.0003295738869195608     
P3G   0.3564840378592988        
P2G   0.04484371111903231       
PEP   0.07361684247063681       
PYR   8.52315246355188          
ACE   0.17011445161350183       
P     6.308881637770621         
NAD   1.5455597670249495        
NADH  0.044440232975050405      



'''

