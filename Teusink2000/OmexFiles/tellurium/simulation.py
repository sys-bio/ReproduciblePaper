import tellurium as te
import os

# get current directory
current_directory = os.path.dirname(__file__)

# get antimony string
antimony_file = os.path.join(current_directory, 'teusink2000.ant')

with open(antimony_file) as f:
    antimony = f.read()

# get roadrunner model using tellurium front end
rr = te.loada(antimony)

# set conserved moiety analysis to True
rr.conservedMoietyAnalysis = True

# compute steady state
rr.steadyState()

# collect steady state values and names
steadystate_values = rr.getSteadyStateValues()
steadystate_names = rr.getSteadyStateSelectionStrings()

# collect into convenient data structure (python dict)
ss = dict(zip(steadystate_names, steadystate_values))

# display steady state
for name, value in ss.items():
    print(name, value)











