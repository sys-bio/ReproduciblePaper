import os
# import site
# site.addsitedir(r'D:\tellurium')
import tellurium as te
# from tellurium.utils.misc import ODEExtractor
import tesedml as libsedml
import matplotlib.pyplot as plt

# find the sbml and sedml directories
# sbml = os.path.join(os.path.dirname(__file__), 'teusink2000.xml')

r = te.loada(
    """
    model new_model
        R1: A => B ; k1*A;
        R2: B => A ; k2*B;
        k1 = 0.1;
        k2 = 0.1; 
        A = 100;
        B = 0
    end
    """
)

r.conservedMoityAnalysis = True

ss = r.steadyState()















