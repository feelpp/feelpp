import feelpp as feelpp 
import pytest
from feelpp.toolboxes.fluid import *


collision_cases = [
    ("circle",'fluid/moving_body/gravity/collisions/circle.cfg',["Circle"],["Cir"]),
    ("ellipse",'fluid/moving_body/gravity/collisions/ellipse.cfg',["Ellipse"],["Ell"])
]

def remesh(f,required_facets,required_elts):
    new_mesh, cpt = feelpp.remesh(mesh=f.mesh(), metric="gradedls({},{})".format(0.01, 0.2), 
    required_elts=required_elts, required_facets=required_facets, params='{"remesh":{ "verbose":-1}}')    
    f.applyRemesh(f.mesh(),new_mesh)


#d
# @pytest.mark.skip(reason="This test is being skipped for now.")
@pytest.mark.parametrize("case,casefile,required_facets,required_elts", collision_cases)
def test_collision(case,casefile,required_facets,required_elts):
    
    feelpp.Environment.setConfigFile(casefile)
    f = fluid(dim=2, orderVelocity=2, orderPressure=1)
    f.init()
    
    remesh(f,required_facets,required_elts)

    radius = 0.125
    forceRange = 0.03

    f.addContactForceModel()
    f.addContactForceResModel()
    f.startTimeStep()
    while not f.timeStepBase().isFinished():
        
        if (f.timeStepBase().iteration() % 5 == 0):
            remesh(f,required_facets,required_elts)
            f.addContactForceModel()
            f.addContactForceResModel()
        
        f.solve()
        f.exportResults()

        if (case == "circle"):
            pos = f.postProcessMeasures().values().get('Quantities_body_Circle.mass_center_1')
            assert(pos < forceRange + radius)
            assert(pos > radius)
                
        if (case == "ellipse"):
            pos = f.postProcessMeasures().values().get('Quantities_body_Ellipse.mass_center_1')
            assert(pos > 0.)
                         
        f.updateTimeStep()
    
    if (case == "ellipse"):
        assert(f.postProcessMeasures().values().get('Quantities_body_Ellipse.rigid_rotation_angles') < 0)
    

#test_collision("circle",'fluid/moving_body/gravity/collisions/circle.cfg',["Circle"],["Cir"])
#test_collision("ellipse",'fluid/moving_body/gravity/collisions/ellipse.cfg',["Ellipse"],["Ell"])


