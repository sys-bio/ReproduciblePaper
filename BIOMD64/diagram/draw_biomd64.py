import os
import site
import matplotlib.pyplot as plt

site.addsitedir(r'D:\libsbml-draw\src\python')
from libsbml_draw import SBMLlayout, Style

model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'model')
fname = os.path.join(model_dir, 'teusink2000.xml')

style = Style()

s = SBMLlayout(fname)

print(s.drawNetwork(scaling_factor=0.4))



plt.show()