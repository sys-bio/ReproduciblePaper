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
ant = r.getAntimony()

r = te.loada(ant)

with open('teusink2000.xml', 'wb') as f:
    f.write(r.getSBML().encode('utf8'))


