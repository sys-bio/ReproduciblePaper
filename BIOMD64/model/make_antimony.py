import os
# import site
# site.addsitedir(r'D:\tellurium')
import tellurium as te
# from tellurium.utils.misc import ODEExtractor
import tesedml as libsedml
import matplotlib.pyplot as plt

# find the sbml and sedml directories
sbml = os.path.join(os.path.dirname(__file__), 'model.xml')

r = te.loadSBMLModel(sbml)

x = r.getRatesOfChange()
print(x)

#                                       GLCi,        G6P ,       F6P ,           F16P,           TRIO,       BPG ,          P3G    ,         P2G    ,   PEP ,           PYR ,           ACE ,        P,              NAD ,          NADH

te_rates_of_changes_t0             = [7.86928894, 46.29403461, -108.83652571, -135.55092418, 370.81078843, 559.34800552, -608.50001547, -68.15681063, 79.81121704, 133.30623792, -1486.13822923, -593.81613491,  1413.3163239, -1413.3163239]
scipy_rates_of_changes_with_odeint = [7.86928894, 46.29403461, -108.83652571, -135.55092418, 370.81078843, 559.34800552, -608.50001547, -68.15681063, 79.81121704, 133.30623792, -1486.13822923, -593.81613491,  1413.3163239, -1413.3163239]

scipy_rates_of_changes_with_ode    = [7.86928894, 46.29403461, -108.83652571, -135.55092418, 370.81078843, 559.34800552, -608.50001547, -68.15681063, 79.81121704, 133.30623792,  1529.34728874, -593.81613491, -1602.1691940,  1602.16919407]



# ant = r.getAntimony()
#
#
#  
# r = te.loada(ant)
#
# with open('teusink2000.xml', 'wb') as f:
#     f.write(r.getSBML().encode('utf8'))
#
# tellurium
# [62.81013265336725, -8.646809100044877, -128.2652241725271, -116.12222571736251, 370.81078843055485, 223.6360824107639, -272.78809236505265, -68.15681062504174, -44.395453513908336, 257.51290847338646, -1529.347288740996, -45.74448995255109, 1413.3163238997802, -1413.3163238997802]
# [62.81013265336725, -8.646809100044877, -128.2652241725271, -116.12222571736251, 370.81078843055485, 223.6360824107639, -272.78809236505265, -68.15681062504174, -44.395453513908336, 257.51290847338646, -1529.347288740996, -45.74448995255109, 1413.3163238997802, -1413.3163238997802]
#
#
# # # te
# #  [[        0.087,    2.45,     0.62,     5.51,     0.96,           0,      0.9,      0.12,      0.07,    1.85,     0.17,    6.31,     1.2,      0.39],
# #   [    0.0949711, 1.10836, 0.125191,  0.62835, 0.795584, 0.000380475, 0.384457, 0.0487723, 0.0832136, 9.59589, 0.193613, 6.44476,  1.5503, 0.0396956],
# #   [    0.0986913, 1.03437,    0.113, 0.602272, 0.777776, 0.000330413, 0.356996, 0.0449151,  0.073787, 8.54612, 0.170734, 6.31121,  1.5457, 0.0442966],
# #   [    0.0987575, 1.03327, 0.112816, 0.601914, 0.777528, 0.000329589, 0.356493,  0.044845, 0.0736198, 8.52355, 0.170125, 6.30892, 1.54556, 0.0444377],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736169, 8.52316, 0.170115, 6.30888, 1.54556, 0.0444402],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736169, 8.52316, 0.170115, 6.30888, 1.54556, 0.0444402],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736169, 8.52316, 0.170115, 6.30888, 1.54556, 0.0444402],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736168, 8.52315, 0.170114, 6.30888, 1.54556, 0.0444402],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736168, 8.52315, 0.170114, 6.30888, 1.54556, 0.0444402],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736168, 8.52315, 0.170114, 6.30888, 1.54556, 0.0444402],
# #   [    0.0987587, 1.03325, 0.112813, 0.601908, 0.777524, 0.000329574, 0.356484, 0.0448437, 0.0736168, 8.52315, 0.170114, 6.30888, 1.54556, 0.0444402]]
#
# '''
# 	GLCi	G6P	F6P	F16P	TRIO	BPG	P3G	P2G	PEP	PYR	ACE	P	NAD	NADH
# 0.087	2.45	0.62	5.51	0.96	0	0.9	0.12	0.07	1.85	0.17	6.31	1.2	0.39
# 2.218186445	0.117881532	0.014716619	4.622480825	2.576573598	7.90E-06	0.037374387	0.004297382	0.005337802	1.861034095	-0.001547563	2.58290736	0.184380562	1.405619438
# 1.52926424e-315	5.91E-222	1.09E-202	7.13E+160	4.90E+252	4.98E+151	5.02E+223	1.17E+171	5e-324	1.060997896e-313	9.91120574086165e-310	1.13E+217	1.14E+217	-1.54E-154
# 1.41E+294	-2.01018158106456e-310	-3.18E-294	6.93E+173	2.86E+289	6.41E+290	9.71E+221	2.62E+289	4.04E+202	-3.26E-241	3.12E+289	4.04E+202	6.71E+202	5.38E+173
# 1.12E-303	3.34E+289	3.10E+289	3.38E+183	3.09E+289	3.39E-111	3.33E-111	3.68E+180	8.52E+223	9.30E+242	1.06E+248	5.67E-302	2.05E-289	3.63E-217
# 6.88571056076e-313	1.76E-292	4.27E-306	8.91E-307	4.23E-307	5.99973186104874e-310	1.28E-288	4.99E-307	-4.46E+127	2.65E+180	3.977582017853e-312	3.12E-308	1.39E+188	5.19723188e-314
# 1.21E+98	1.03E+243	8.04E+165	5.33E+233	9.92E+247	1.23E+171	4.6960284415216e-310	4.63E-303	1.43E+241	4.5873305009129e-310	5.3995336897e-313	3.91E-305	1.33E+241	5.43E+247
# 9.45E+218	3.35E-157	2.21649091603e-313	1.27E-306	1.33E+241	4.75029977766064e-310	9.39E+93	2.35E-306	1.33E+241	0.00E+00	3.28E-157	8.28E+170	4.33E-307	1.33E+241
# -2.38E-22	4.7628858652e-313	0.00E+00	7.29E-304	4.70E-294	1.00E-303	8.37E-299	1.21E+98	1.03E+243	8.04E+165	5.33E+233	-2.0642926032517e-310	4.889078303396e-311	5.19E-303
# 2.22E-301	2.50E-306	3.73E-301	3.26E-292	2.13E-287	7.02E-292	1.42E-248	3.12E-282	2.15E-282	0.00E+00	5.57E-307	4.60E-287	2.49E-306	4.37E-210
# 3.94E-268	3.97E-263	2.60E-258	1.71E-253	1.71E-253	4.60E-287	2.49E-306	1.74E+137	-1.81E+60	1.64E-287	2.13E-287	4.86E-239	2.02E-233	2.34E-258
#
# '''
#
# # x = r.simulate(0, 10, 11)
# # print(x)
#
# print(r.ATP)
# print(r.ADP)
# print(r.AMP)
#
# #te
# ATP ADP AMP
# 2.5091904762788237 1.2916190474423517 0.2991904762788242
# 2.5091904762788237 1.2916190474423517 0.2991904762788242
