import feelpp.core as fppc
import feelpp.core.meshmover as mm
from feelpp.core.measure import measure

def run(m, geo):
    
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy = geo
    mesh = fppc.load(m, mesh_name, 0.1)

    Xh = fppc.functionSpace(mesh=mesh,space="Pchv")
    u = Xh.element()
    M=measure(range=fppc.elements(mesh))

    if mesh.dimension()==2:
        assert(abs(M-2)<1e-10)
    else:
        assert(abs(M-1)<1e-10)
    
    u.on(range=fppc.elements(mesh), expr=fppc.expr("{x,y}:x:y" if dim ==2 else "{x,y,z}:x:y:z",row=dim)) 

    newmesh = mm.meshMove(mesh,u)
    M=measure(range=fppc.elements(newmesh))   

    # Equal to 8 in 2d and 3d
    assert(abs(M-8)<1e-10)

def test_meshmove(init_feelpp):
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/meshmove/")
    geo = {
        '2': fppc.create_rectangle(),# [0,1]x[0,2]
        '3': fppc.create_box(), # [0,1]x[0,2]x[0,0.5]
    }
    run(fppc.mesh(dim=2), geo['2'])
    run(fppc.mesh(dim=3, realdim=3), geo['3'])