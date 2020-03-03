import os
import site
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

site.addsitedir(r'D:\libsbml-draw\src\python')
from libsbml_draw import SBMLlayout, Style

model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'model')
fname = os.path.join(model_dir, 'teusink2000.xml')
fname = r'C:\Users\cwelsh\Downloads\graphfavb\BIO64_layout.xml'
fname = r'D:\ReproduciblePaper\BIOMD64\diagram\Teusink_layout.xml'
fname = r'D:\ReproduciblePaper\BIOMD64\model\JWSOnline_Teusink.xml'
fnameimg = r'D:\ReproduciblePaper\BIOMD64\model\image.png'

if os.path.isfile(fnameimg):
    os.remove(fnameimg)

style = Style()
style.font.size = 35
style.edge.edgecolor = '#808080'
style.edge.fillcolor = '#4d4d4d'
style.edge.width = 8


style.node.edgecolor = '#ff8c1a'
style.node.fillcolor = '#ffe6cc'
style.node.edgewidth = 8

style.arrow.scale = 75

s = SBMLlayout(fname, style=style)


for i in range(1):
    fnameimg = fr'D:\ReproduciblePaper\BIOMD64\model\image{i}.pdf'
    s.regenerateLayout()
    # s.setReactionColor('all', 'black')
    # s.setReactionCurveWidth('all', 25)
    s.apply_style()
    s.drawNetwork(show=False, save_file_name=fnameimg)
    s.writeSBML()






