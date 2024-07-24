#!/usr/bin/env python

__docformat__ = 'restructuredtext'

#from fipy.tools.parser import parse

#numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 2500)
#numberOfSteps = parse('--numberOfSteps', action = 'store', type = 'int', default = 10)

#### this py code should reproduce the 1D one-mode results in paper PRE 70/ 051605 (2004) , i.e. :
#### psi = psi_offest + A sin qx 
#### r > -0.25 : periodic state for offset psi as psi_offset < sqrt(-r/3)
####		   : constant state for offset psi as psi_offset > sqrt(-r/3)
#### r < -0.25 : coexistence of liquid and solid as shown in plot Fig.2 
#### r = -0.03: sqrt (-r/3) = 0.1; r = -0.12: sqrt(-r/3) = 0.2

#### also try to subtract due to Eq. 31 the free energy corresponding to psi_offset:
#### omega_0 psi_offset^2/2 + psi_offset^4/4

from math import pi as PI

import fipy.tools.numerix as numerix
steps = 1000
nx = 100
dx = 0.02*PI
L = dx * nx


r =  -0.5
psi_offset = 0.35

from fipy.meshes.grid1D import Grid1D
mesh = Grid1D(dx, nx)

from fipy.variables.cellVariable import CellVariable
from fipy.tools.numerix import random

var = CellVariable(name = "A",
                   mesh = mesh,
                   value = psi_offset*2.*random.random(nx))
					#value = sin(nx*dx*Pi/L)
					
				   
faceVar = var.getArithmeticFaceValue()

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.transientTerm import TransientTerm

Eq28_diffTerm2 = ImplicitDiffusionTerm(coeff = (r+1.))
Eq28_diffTerm4 = ImplicitDiffusionTerm(coeff = (2.,1.))
Eq28_diffTerm6 = ImplicitDiffusionTerm(coeff = (1.,1.,1.))
Eq28_nonlinDiffTerm = ImplicitDiffusionTerm(coeff = 3.*faceVar*faceVar)
Eq28_constChemPotential = 0.5*r*psi_offset*psi_offset + 0.25*psi_offset*psi_offset*psi_offset*psi_offset

eqch = TransientTerm() - Eq28_diffTerm2 - Eq28_diffTerm4 - Eq28_diffTerm6 - Eq28_nonlinDiffTerm + Eq28_constChemPotential

from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
##solver = LinearLUSolver(tolerance = 1e-15,steps = 1000) ## unsolved problem to implement LinearLUSolver
solver = LinearPCGSolver(tolerance = 1e-10,steps = 1000)

from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition

#BCs = (
	#FixedValue(faces=mesh.getFacesRight(), value=1),
    #FixedValue(faces=mesh.getFacesLeft(), value=.5),
    #NthOrderBoundaryCondition(faces=mesh.getFacesLeft(), value=0, order=2),
    #NthOrderBoundaryCondition(faces=mesh.getFacesRight(), value=0, order=2))

	
BCs = (FixedValue(mesh.getFacesRight(), psi_offset),
       FixedValue(mesh.getFacesLeft(), psi_offset),
       NthOrderBoundaryCondition(mesh.getFacesLeft(), 0, 2),
       NthOrderBoundaryCondition(mesh.getFacesRight(), 0, 2),
	   NthOrderBoundaryCondition(mesh.getFacesLeft(), 0, 4),
       NthOrderBoundaryCondition(mesh.getFacesRight(), 0, 4)  )


if __name__ == '__main__':

    import fipy.viewers
	
    viewer = fipy.viewers.make(vars = var, limits = {'datamin': 0., 'datamax': 1.0})
    viewer.plot()
	
    dexp=-5

for step in range(steps):
    dt = numerix.exp(dexp)
    dt = min(100, dt)
    dexp += 0.01
    var.updateOld()
    eqch.solve(var, boundaryConditions = BCs, solver = solver, dt = dt)

    if __name__ == '__main__':
        viewer.plot()
        print 'step',step,'dt',dt
	
def _run():
    pass
            
