import os
import tellurium as te
import tesedml as libsedml
import matplotlib.pyplot as plt

sbml_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sbml')
sedml_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sedml')
if not os.path.isdir(sbml_dir):
    raise ValueError('cannot find sbml directory')
if not os.path.isdir(sedml_dir):
    raise ValueError('cannot find sedml directory')
teusink2000 = os.path.join(sbml_dir, 'model.xml')
sedml_file = os.path.join(sedml_dir, 'simulation.xml')
print(teusink2000)
r = te.loadSBMLModel(teusink2000)

sedml_doc = libsedml.readSedML(sedml_file)
n_errors = sedml_doc.getErrorLog().getNumFailsWithSeverity(libsedml.LIBSEDML_SEV_ERROR)
print('Read SED-ML file, number of errors: {}'.format(n_errors))
if n_errors > 0:
    print(sedml_doc.getErrorLog().toString())

with open(sedml_file, 'r') as f:
    sedml = f.read()

te.executeSEDML(sedml, workingDir=os.path.dirname(__file__))

