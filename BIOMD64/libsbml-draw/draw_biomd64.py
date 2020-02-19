import os
from libsbml_draw import SBMLlayout

s = SBMLlayout('BIOMD0000000064.xml')

print(s.drawNetwork())