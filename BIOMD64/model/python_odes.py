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
GLCi
G6P
F6P
F16BP
TRIO
BPG
P3G
P2G
PEP
PYR
ACE
P
NAD
NADH
'''

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

y0 = [GLCi, G6P, F6P, F16P, TRIO, BPG, P3G, P2G, PEP, PYR, ACE, P, NAD, NADH]


def teusink2000(y0, t):
    # Unpack Species
    GLCi, G6P, F6P, F16P, TRIO, BPG, P3G, P2G, PEP, PYR, ACE, P, NAD, NADH = y0

    Glyc = 0  # fixed
    Trh = 0  # fixed
    CO2 = 1  # fixed
    SUCC = 0  # fixed

    ETOH = 50  # fixed
    GLCo = 50  # fixed
    GLY = 0.15  # fixed
    SUM_P = 4.1  # fixed
    F26BP = 0.02  # fixed

    # Compartment
    extracellular = 1
    cytosol = 1

    # parameters
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

    def adp():
        return (SUM_P
                - (P ** (2 * (1 - 4 * KeqAK))
                   + 2 * SUM_P * P * (4 * KeqAK - 1)
                   + SUM_P ** 2) ** 0.5) / (1 - 4 * KeqAK)

    def atp():
        return (P - adp()) / 2

    def amp():
        return SUM_P - atp() - adp()

    # functions

    def alcohol_dehydrogenase():
        return -cytosol * ((VmADH / (KiADHNAD * KmADHETOH)) * (NAD * ETOH - NADH * ACE / KeqADH) / (
                1 + NAD / KiADHNAD + KmADHNAD * ETOH / (KiADHNAD * KmADHETOH) + KmADHNADH * ACE / (
                KiADHNADH * KmADHACE) + NADH / KiADHNADH + NAD * ETOH / (
                        KiADHNAD * KmADHETOH) + KmADHNADH * NAD * ACE / (
                        KiADHNAD * KiADHNADH * KmADHACE) + KmADHNAD * ETOH * NADH / (
                        KiADHNAD * KmADHETOH * KiADHNADH) + NADH * ACE / (
                        KiADHNADH * KmADHACE) + NAD * ETOH * ACE / (
                        KiADHNAD * KmADHETOH * KiADHACE) + ETOH * NADH * ACE / (
                        KiADHETOH * KiADHNADH * KmADHACE))) / cytosol

    def glycerol_3_phosphate_dehydrogenase():
        return (VmG3PDH / (KmG3PDHDHAP * KmG3PDHNADH)) * ((1 / (1 + KeqTPI)) * TRIO * NADH - GLY * NAD / KeqG3PDH) / (
                (1 + (1 / (1 + KeqTPI)) * TRIO / KmG3PDHDHAP + GLY / KmG3PDHGLY) * (
                1 + NADH / KmG3PDHNADH + NAD / KmG3PDHNAD))

    def r_pfk():
        return 1 + F6P / KmPFKF6P + atp() / KmPFKATP + gR * (F6P / KmPFKF6P) * (atp() / KmPFKATP)

    def t_pfk():
        return 1 + CPFKATP * (atp() / KmPFKATP)

    def l_pfk():
        return Lzero * ((1 + CiPFKATP * (atp() / KiPFKATP)) / (1 + atp() / KiPFKATP)) ** 2 * (
                (1 + CPFKAMP * (amp() / KPFKAMP)) / (1 + amp() / KPFKAMP)) ** 2 * (
                       (1 + CPFKF26BP * F26BP / KPFKF26BP + CPFKF26BP * F16P / KPFKF16BP) / (
                       1 + F26BP / KPFKF26BP + F16P / KPFKF16BP)) ** 2

    def hexokinase():
        return (VmGLK / (KmGLKGLCi * KmGLKATP)) * (GLCi * atp() - G6P * adp() / KeqGLK) / (
                (1 + GLCi / KmGLKGLCi + G6P / KmGLKG6P) * (1 + atp() / KmGLKATP + adp() / KmGLKADP))

    def glucose_6_phosphate_isomerase():
        return (VmPGI_2 / KmPGIG6P_2) * (G6P - F6P / KeqPGI_2) / (1 + G6P / KmPGIG6P_2 + F6P / KmPGIF6P_2)

    def aldolase():
        return (VmALD / KmALDF16P) * (F16P - (KeqTPI / (1 + KeqTPI)) * TRIO * (1 / (1 + KeqTPI)) * TRIO / KeqALD) / (
                1 + F16P / KmALDF16P + (KeqTPI / (1 + KeqTPI)) * TRIO / KmALDGAP + (
                1 / (1 + KeqTPI)) * TRIO / KmALDDHAP + (KeqTPI / (1 + KeqTPI)) * TRIO * (
                        1 / (1 + KeqTPI)) * TRIO / (KmALDGAP * KmALDDHAP) + F16P * (
                        KeqTPI / (1 + KeqTPI)) * TRIO / (KmALDGAPi * KmALDF16P))

    def phosphoglycerate_kinase():
        return (VmPGK / (KmPGKP3G * KmPGKATP)) * (KeqPGK * BPG * adp() - P3G * atp()) / (
                (1 + BPG / KmPGKBPG + P3G / KmPGKP3G) * (1 + atp() / KmPGKATP + adp() / KmPGKADP))

    def glyceraldehyde_3_phosphate_dehydrogenase():
        return (VmGAPDHf * (KeqTPI / (1 + KeqTPI)) * TRIO * NAD / (KmGAPDHGAP * KmGAPDHNAD) - VmGAPDHr * BPG * NADH / (
                KmGAPDHBPG * KmGAPDHNADH)) / ((1 + (KeqTPI / (1 + KeqTPI)) * TRIO / KmGAPDHGAP + BPG / KmGAPDHBPG) * (
                1 + NAD / KmGAPDHNAD + NADH / KmGAPDHNADH))

    def enolase():
        return (VmENO / KmENOP2G) * (P2G - PEP / KeqENO) / (1 + P2G / KmENOP2G + PEP / KmENOPEP)

    def pyruvate_decarboxylase():
        return VmPDC * (PYR ** nPDC / KmPDCPYR ** nPDC) / (1 + PYR ** nPDC / KmPDCPYR ** nPDC)

    def succinate_synthesis():
        return KSUCC * ACE

    def glucose_transport():
        return (VmGLT / KmGLTGLCo) * (GLCo - GLCi / KeqGLT) / (
                1 + GLCo / KmGLTGLCo + GLCi / KmGLTGLCi + 0.91 * GLCo * GLCi / (KmGLTGLCo * KmGLTGLCi))

    def phosphoglycerate_mutase():
        return (VmPGM / KmPGMP3G) * (P3G - P2G / KeqPGM) / (1 + P3G / KmPGMP3G + P2G / KmPGMP2G)

    def pyruvate_kinase():
        return (VmPYK / (KmPYKPEP * KmPYKADP)) * (PEP * adp() - PYR * atp() / KeqPYK) / (
                (1 + PEP / KmPYKPEP + PYR / KmPYKPYR) * (1 + atp() / KmPYKATP + adp() / KmPYKADP))

    def ATPase_activity():
        return KATPASE * atp()

    def phosphofructokinase():
        return VmPFK * gR * (F6P / KmPFKF6P) * (atp() / KmPFKATP) * r_pfk() / (r_pfk() ** 2 + l_pfk() * t_pfk() ** 2)

    def glycogen_synthesis():
        return vGLYCO_v

    def trehalose6p_synthesis():
        return vTreha_v

    return [
        # GLCi, glucose in cytosol
        -cytosol * hexokinase()
        + glucose_transport(),

        # G6P, glucose 6 phosphase
        + cytosol * hexokinase()
        - cytosol * glucose_6_phosphate_isomerase()
        - cytosol * glycogen_synthesis()
        - 2 * cytosol * trehalose6p_synthesis(),

        # F6P, fructose 6 phosphate
        + cytosol * glucose_6_phosphate_isomerase()
        - cytosol * phosphofructokinase(),

        # F16BP, fructose-1-6-bisphosphate
        + cytosol * phosphofructokinase()
        - cytosol * aldolase(),

        # TRIO Triose-phosphate
        - cytosol * glycerol_3_phosphate_dehydrogenase()
        + 2 * cytosol * aldolase()
        - cytosol * glyceraldehyde_3_phosphate_dehydrogenase(),

        # BPG, 1,3-bisphosphoglycerate
        + cytosol * glyceraldehyde_3_phosphate_dehydrogenase()
        - cytosol * phosphoglycerate_kinase(),

        # P3G, 3-phosphoglycerate
        + cytosol * phosphoglycerate_kinase()
        - cytosol * phosphoglycerate_mutase(),

        # P2G, 2-phosphoglycerate
        + cytosol * phosphoglycerate_mutase()
        - cytosol * enolase(),

        # PEP, Phosphoenolpyruvate
        - cytosol * pyruvate_kinase()
        + cytosol * enolase(),

        # PYR, Pyruvate
        + cytosol * pyruvate_kinase()
        - cytosol * pyruvate_decarboxylase(),

        # ACE, Acetaldehyde
        - cytosol * pyruvate_decarboxylase()
        + 2 * cytosol * succinate_synthesis()
        - cytosol * alcohol_dehydrogenase(),

        # P high energy phosphates
        - cytosol * hexokinase()
        + cytosol * pyruvate_kinase()
        - 4 * cytosol * succinate_synthesis()
        - cytosol * ATPase_activity()
        - cytosol * glycogen_synthesis()
        - cytosol * trehalose6p_synthesis()
        - cytosol * phosphofructokinase()
        + cytosol * phosphoglycerate_kinase(),

        # NAD,
        - 3 * cytosol * succinate_synthesis()
        + cytosol * alcohol_dehydrogenase()
        + cytosol * glycerol_3_phosphate_dehydrogenase()
        - cytosol * glyceraldehyde_3_phosphate_dehydrogenase(),

        # NADH,
        + 3 * cytosol * succinate_synthesis()
        - cytosol * alcohol_dehydrogenase()
        - cytosol * glycerol_3_phosphate_dehydrogenase()
        + cytosol * glyceraldehyde_3_phosphate_dehydrogenase()
    ]


import numpy as np
from scipy.integrate import odeint

t = np.linspace(0, 10, 11)
sol = odeint(teusink2000, y0, t)

import pandas as pd
import matplotlib

matplotlib.use('TkAgg')

data = pd.DataFrame(sol)
data.columns = ['GLCi', 'G6P', 'F6P', 'F16P', 'TRIO', 'BPG', 'P3G', 'P2G', 'PEP', 'PYR', 'ACE', 'P', 'NAD', 'NADH']
print(data)

for i in data:
    plt.figure()
    plt.plot(data.index, data[i])

plt.show()
