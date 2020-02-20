import os
# import site
# site.addsitedir(r'D:\tellurium')
import tellurium as te
# from tellurium.utils.misc import ODEExtractor
import tesedml as libsedml
import matplotlib.pyplot as plt

# # find the sbml and sedml directories
# sbml_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sbml')
#
# if not os.path.isdir(sbml_dir):
#     raise ValueError('cannot find sbml directory {}'.format(sbml_dir))
#
# # find the sbml and sedml files
# teusink2000 = os.path.join(sbml_dir, 'model.xml')

# if not os.path.isfile(teusink2000):
#     raise FileNotFoundError(teusink2000)
#
# r = te.loadSBMLModel(teusink2000)

# print(r.getAntimony())


# y0 = [
#     GLCi, G6P, F6P, F16P, TRIO,
#     BPG, P3G, P2G, PEP, PYR, ACE, P,
#     NAD, NADH, Glyc, Trh, CO2, SUCC,
#     GLCo, ETOH, GLY, SUM_P, F26BP]

'''
vGLK
vPGI
vGLYCO
vTreha
vPFK
vALD
vGAPDH
vPGK
vPGM
vENO
vPYK
vPDC
vSUC
vGLT
vADH
vG3PDH
vATP
'''

# vGLK
# vPGI
# vGLYCO
# vTreha
# vPFK
# vALD
# vGAPDH
# vPGK
# vPGM
# vENO
# vPYK
# vPDC
# vSUC
# vGLT
# vADH
# vG3PDH
# vATP

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
NAD = 1.2
NADH = 0.39
GLCo = 50
ETOH = 50
GLY = 0.15
F26BP = 0.02

y0 = [GLCi, G6P, F6P, F16P, TRIO, BPG, P3G, P2G, PEP, PYR, ACE, NAD, NADH, GLCo, ETOH, GLY, F26BP, F26BP ]


def teusink2000(y0, t):
    # Unpack Species
    GLCi, G6P, F6P, F16P, TRIO, BPG, P3G, P2G, PEP, PYR, ACE, NAD, NADH, GLCo, ETOH, GLY, F26BP, F26BP = y0

    P = 6.31
    SUM_P = 4.1
    Glyc = 0
    Trh = 0
    CO2 = 1
    SUCC = 0
    # Compartment
    extracellular = 1
    cytosol = 1

    # Variable
    KeqAK = 0.45
    gR = 5.12
    KmPFKF6P = 0.1
    KmPFKATP = 0.71
    Lzero = 0.66
    CiPFKATP = 100
    KiPFKATP = 0.65
    CPFKAMP = 0.0845
    KPFKAMP = 0.0995
    CPFKF26BP = 0.0174
    KPFKF26BP = 0.000682
    CPFKF16BP = 0.397
    KPFKF16BP = 0.111
    CPFKATP = 3
    KeqTPI = 0.045

    KeqGLK = 3800
    KmGLKADP = 0.23
    KmGLKATP = 0.15
    KmGLKG6P = 30
    KmGLKGLCi = 0.08
    VmGLK = 226.452
    KeqPGI_2 = 0.314
    KmPGIF6P_2 = 0.3
    KmPGIG6P_2 = 1.4
    VmPGI_2 = 339.677
    vGLYCO_v = 6
    vTreha_v = 2.4
    VmPFK = 182.903
    KeqALD = 0.069
    KmALDDHAP = 2.4
    KmALDF16P = 0.3
    KmALDGAP = 2
    KmALDGAPi = 10
    VmALD = 322.258
    KmGAPDHBPG = 0.0098
    KmGAPDHGAP = 0.21
    KmGAPDHNAD = 0.09
    KmGAPDHNADH = 0.06
    VmGAPDHf = 1184.52
    VmGAPDHr = 6549.8
    KeqPGK = 3200
    KmPGKADP = 0.2
    KmPGKATP = 0.3
    KmPGKBPG = 0.003
    KmPGKP3G = 0.53
    VmPGK = 1306.45
    KeqPGM = 0.19
    KmPGMP2G = 0.08
    KmPGMP3G = 1.2
    VmPGM = 2525.81
    KeqENO = 6.7
    KmENOP2G = 0.04
    KmENOPEP = 0.5
    VmENO = 365.806
    KeqPYK = 6500
    KmPYKADP = 0.53
    KmPYKATP = 1.5
    KmPYKPEP = 0.14
    KmPYKPYR = 21
    VmPYK = 1088.71
    KmPDCPYR = 4.33
    VmPDC = 174.194
    nPDC = 1.9
    KSUCC = 21.4
    KeqGLT = 1
    KmGLTGLCi = 1.1918
    KmGLTGLCo = 1.1918
    VmGLT = 97.264
    KeqADH = 6.9e-05
    KiADHACE = 1.1
    KiADHETOH = 90
    KiADHNAD = 0.92
    KiADHNADH = 0.031
    KmADHACE = 1.11
    KmADHETOH = 17
    KmADHNAD = 0.17
    KmADHNADH = 0.11
    VmADH = 810
    KeqG3PDH = 4300
    KmG3PDHDHAP = 0.4
    KmG3PDHGLY = 1
    KmG3PDHNAD = 0.93
    KmG3PDHNADH = 0.023
    VmG3PDH = 70.15
    KATPASE = 33.7

    # assignment rules
    ADP = (SUM_P - (P ** (2 * (1 - 4 * KeqAK)) + 2 * SUM_P * P * (4 * KeqAK - 1) + SUM_P ** 2) ** 0.5) / (1 - 4 * KeqAK)
    ATP = (P - ADP) / 2
    AMP = SUM_P - ATP - ADP

    # functions

    def Alcohol_dehydrogenase():
        return -cytosol * ((VmADH / (KiADHNAD * KmADHETOH)) * (NAD * ETOH - NADH * ACE / KeqADH) / (
                1 + NAD / KiADHNAD + KmADHNAD * ETOH / (KiADHNAD * KmADHETOH) + KmADHNADH * ACE / (
                KiADHNADH * KmADHACE) + NADH / KiADHNADH + NAD * ETOH / (
                        KiADHNAD * KmADHETOH) + KmADHNADH * NAD * ACE / (
                        KiADHNAD * KiADHNADH * KmADHACE) + KmADHNAD * ETOH * NADH / (
                        KiADHNAD * KmADHETOH * KiADHNADH) + NADH * ACE / (
                        KiADHNADH * KmADHACE) + NAD * ETOH * ACE / (
                        KiADHNAD * KmADHETOH * KiADHACE) + ETOH * NADH * ACE / (
                        KiADHETOH * KiADHNADH * KmADHACE))) / cytosol

    def Glycerol_3_phosphate_dehydrogenase():
        return (VmG3PDH / (KmG3PDHDHAP * KmG3PDHNADH)) * ((1 / (1 + KeqTPI)) * TRIO * NADH - GLY * NAD / KeqG3PDH) / (
                (1 + (1 / (1 + KeqTPI)) * TRIO / KmG3PDHDHAP + GLY / KmG3PDHGLY) * (
                1 + NADH / KmG3PDHNADH + NAD / KmG3PDHNAD))

    def R_PFK():
        return 1 + F6P / KmPFKF6P + ATP / KmPFKATP + gR * (F6P / KmPFKF6P) * (ATP / KmPFKATP)

    def T_PFK():
        return 1 + CPFKATP * (ATP / KmPFKATP)

    def L_PFK():
        return Lzero * ((1 + CiPFKATP * (ATP / KiPFKATP)) / (1 + ATP / KiPFKATP)) ** 2 * (
                (1 + CPFKAMP * (AMP / KPFKAMP)) / (1 + AMP / KPFKAMP)) ** 2 * (
                       (1 + CPFKF26BP * F26BP / KPFKF26BP + CPFKF26BP * F16P / KPFKF16BP) / (
                       1 + F26BP / KPFKF26BP + F16P / KPFKF16BP)) ** 2

    def Hexokinase():
        return (VmGLK / (KmGLKGLCi * KmGLKATP)) * (GLCi * ATP - G6P * ADP / KeqGLK) / (
                (1 + GLCi / KmGLKGLCi + G6P / KmGLKG6P) * (1 + ATP / KmGLKATP + ADP / KmGLKADP))

    def Glucose_6_phosphate_isomerase():
        return (VmPGI_2 / KmPGIG6P_2) * (G6P - F6P / KeqPGI_2) / (1 + G6P / KmPGIG6P_2 + F6P / KmPGIF6P_2)

    def Aldolase():
        return (VmALD / KmALDF16P) * (F16P - (KeqTPI / (1 + KeqTPI)) * TRIO * (1 / (1 + KeqTPI)) * TRIO / KeqALD) / (
                1 + F16P / KmALDF16P + (KeqTPI / (1 + KeqTPI)) * TRIO / KmALDGAP + (
                1 / (1 + KeqTPI)) * TRIO / KmALDDHAP + (KeqTPI / (1 + KeqTPI)) * TRIO * (
                        1 / (1 + KeqTPI)) * TRIO / (KmALDGAP * KmALDDHAP) + F16P * (
                        KeqTPI / (1 + KeqTPI)) * TRIO / (KmALDGAPi * KmALDF16P))

    def Phosphoglycerate_kinase():
        return (VmPGK / (KmPGKP3G * KmPGKATP)) * (KeqPGK * BPG * ADP - P3G * ATP) / (
                (1 + BPG / KmPGKBPG + P3G / KmPGKP3G) * (1 + ATP / KmPGKATP + ADP / KmPGKADP))

    def Glyceraldehyde_3_phosphate_dehydrogenase():
        return (VmGAPDHf * (KeqTPI / (1 + KeqTPI)) * TRIO * NAD / (KmGAPDHGAP * KmGAPDHNAD) - VmGAPDHr * BPG * NADH / (
                KmGAPDHBPG * KmGAPDHNADH)) / ((1 + (KeqTPI / (1 + KeqTPI)) * TRIO / KmGAPDHGAP + BPG / KmGAPDHBPG) * (
                1 + NAD / KmGAPDHNAD + NADH / KmGAPDHNADH))

    def Enolase():
        return (VmENO / KmENOP2G) * (P2G - PEP / KeqENO) / (1 + P2G / KmENOP2G + PEP / KmENOPEP)

    def Pyruvate_decarboxylase():
        return VmPDC * (PYR ** nPDC / KmPDCPYR ** nPDC) / (1 + PYR ** nPDC / KmPDCPYR ** nPDC)

    def Succinate_synthesis():
        return KSUCC * ACE

    def Glucose_transport():
        return (VmGLT / KmGLTGLCo) * (GLCo - GLCi / KeqGLT) / (
                1 + GLCo / KmGLTGLCo + GLCi / KmGLTGLCi + 0.91 * GLCo * GLCi / (KmGLTGLCo * KmGLTGLCi))

    def Phosphoglycerate_mutase():
        return (VmPGM / KmPGMP3G) * (P3G - P2G / KeqPGM) / (1 + P3G / KmPGMP3G + P2G / KmPGMP2G)

    def Pyruvate_kinase():
        return (VmPYK / (KmPYKPEP * KmPYKADP)) * (PEP * ADP - PYR * ATP / KeqPYK) / (
                (1 + PEP / KmPYKPEP + PYR / KmPYKPYR) * (1 + ATP / KmPYKATP + ADP / KmPYKADP))

    def ATPase_activity():
        return KATPASE * ATP

    # Phosphofructokinase(AMP, ATP, CPFKAMP, CPFKATP, CPFKF16BP, CPFKF26BP, CiPFKATP, F16P, F26BP, F6P, KPFKAMP, KPFKF16BP, KPFKF26BP, KiPFKATP, KmPFKATP, KmPFKF6P, Lzero, VmPFK, gR)
    def Phosphofructokinase():
        return VmPFK * gR * (F6P / KmPFKF6P) * (ATP / KmPFKATP) * R_PFK() / (R_PFK() ** 2 + L_PFK() * T_PFK() ** 2)

    dydt = [
        # vGLK
        cytosol * Hexokinase(),
        # vPGI
        cytosol * Glucose_6_phosphate_isomerase(),
        # vGLYCO
        cytosol * vGLYCO_v,
        # vTreha
        cytosol * vTreha_v,
        # vPFK
        cytosol * Phosphofructokinase(),
        # vALD
        cytosol * Aldolase(),
        # vGAPDH
        cytosol * Glyceraldehyde_3_phosphate_dehydrogenase(),
        # vPGK
        cytosol * Phosphoglycerate_kinase(),
        # vPGM
        cytosol * Phosphoglycerate_mutase(),
        # vENO
        cytosol * Enolase(),
        # vPYK
        cytosol * Pyruvate_kinase(),
        # vPDC
        cytosol * Pyruvate_decarboxylase(),
        # vSUC
        cytosol * Succinate_synthesis(),
        # vGLT
        Glucose_transport(),
        # vADH
        cytosol * Alcohol_dehydrogenase(),
        # vG3PDH
        cytosol * Glycerol_3_phosphate_dehydrogenase(),
        # vATP
        cytosol * ATPase_activity(),
    ]

    return dydt


import numpy as np
from scipy.integrate import odeint

t = np.linspace(0, 10, 101)
sol = odeint(teusink2000, y0, t)


import pandas as pd
import matplotlib
matplotlib.use('TkAgg')

data = pd.DataFrame(sol)
print(data)

for i in data:
    plt.figure()
    plt.plot(data.index, data[i])

plt.show()





