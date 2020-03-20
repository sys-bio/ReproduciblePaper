# semsim only available on a very specific version of cygwin
import semsim
from semsim import *
fname = r'D:\ReproduciblePaper\BIOMD64\model\model.xml'
with open(fname) as f:
    sbml = f.read()


importer = SBMLImporter(sbml)
model = importer.getSBMLModel()


c = model.getComponentForId('cytosol')
c.setAnnotation(c.getAnnotation().makeComposite(semsim.PhysicalProperty(semsim.OPB.get(523))))
c.getCompositeAnnotation().addTerm(bqb.isPartOf, CL.get(169))

rdf = model.getRDF(fname, 'turtle')
rdf_file = r'/BIOMD64/model/annotation_turtle.txt'
with open(rdf_file, 'w') as f:
    f.write(rdf) 










