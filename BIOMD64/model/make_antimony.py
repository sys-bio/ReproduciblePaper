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

te_rates_of_changes_t0 = [7.86928894, 46.29403461, -108.83652571, -135.55092418, 370.81078843, 559.34800552,
                          -608.50001547, -68.15681063, 79.81121704, 133.30623792, -1486.13822923, -593.81613491,
                          1413.3163239, -1413.3163239]
scipy_rates_of_changes_with_odeint = [7.86928893, 46.29403461, -108.83652571, -135.55092418, 370.81078843, 559.34800552,
                                      -608.50001547, -68.15681063, 79.81121704, 133.30623792, 1529.34728874,
                                      -593.81613491, -1602.1691940, 1602.16919407]


# ant = r.getAntimony()
#
#
#  
# r = te.loada(ant)
#
# with open('teusink2000.xml', 'wb') as f:
#     f.write(r.getSBML().encode('utf8'))
#
