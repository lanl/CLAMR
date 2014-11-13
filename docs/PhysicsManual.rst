==============
Physics Manual
==============

CLAMR is a cell-based adaptive mesh refinement (AMR) hydrodynamic
mini-app that simulates the long range propagation of waves. CLAMR
is being developed at LANL as a test-bed for algorithm development
for next-generation architectures. The computational model uses the
shallow water equations to simulate fluid flow and harnesses the three
conservation laws of mass, x momentum, and y momentum, that are inherent
in the shallow water equations {[}16{]} as shown in Equations \ref{eq:mass}-\ref{eq:y-momentum}. 

.. math::

   \begin{eqnarray*}
   \frac{\partial h}{\partial t}+\frac{\partial(hu)}{\partial x}+\frac{\partial(hv)}{\partial y} & = & 0\hspace{1 cm}\text{(Conservation of mass)}\\
   \frac{\partial(hu)}{\partial t}+\frac{\partial}{\partial x}\left(hu^{2}+\frac{1}{2}gh^{2}\right)+\frac{\partial}{\partial y}(huv) & = & 0\hspace{1 cm}\text{(Conservation of $x$-momentum)}\\
   \frac{\partial(hv)}{\partial t}+\frac{\partial}{\partial x}(hvu)+\frac{\partial}{\partial y}\left(hv^{2}+\frac{1}{2}gh^{2}\right) & = & 0\hspace{1 cm}\text{(Conservation of $y$-momentum)}
   \end{eqnarray*}

where h is the total fluid column height and the vector (u, v) represents
the fluid\textquoteright{}s horizontal velocities over time. The constant
of gravitational acceleration :math:`g` is approximately :math:`9.8` meters
per :math:`s^{2}`. Simplifying assumptions include that the flow in the
vertical direction is negligible and that the fluid bottom is flat.
At each time step the state variables for height and momentum, in
both the x and y directions, are updated. Due to the incompressibility
of water, the density of the water can be treated as constant, and
thus the height of the water column is essentially the total mass
of the column. Figure 1 shows a conceptual diagram of the shallow
water simulation. 

In this study the standard test problem of a circular dam break is
used. For this problem a cylindrical pulse is created at the center
of the mesh and the shallow-water equations are used to calculate
the wave moving outward. The circular dam break is considered to be
a challenging problem due to the sharp rise at the dam break location,
the need to maintain cylindrical symmetry, and the near zero, or dry,
condition that results at the center of the simulation caused by the
flow moving symmetrically outward. As the shock reverberates off of
the boundaries, the waves generated begin to dampen. This dampening
eventually leads to a steady-state where the water settles to a uniform
average pool height. The conservation of mass becomes crucially important
in validating the consistency of the algorithm as it progresses forward
in simulation time, and across processors in calculations done in
parallel. Because of the importance of mass conservation, both theoretically
and operationally, it is checked frequently during the run.

CLAMR is representative of a large class of high performance computing
applications based on conservation laws with a finite difference or
finite volume formulation. CLAMR uses the shallow water equations
because they are one of the simpler sets of equations to work with.
This simpler equation set reduces the effort needed for exploring
new computational and programming methods. Slightly more complex are
the Euler equations, which are used for compressible fluid flow modeling.
The Euler equations add an energy equation to the mass and momentum
equations increasing the number of equations from three to four and
adding some additional terms. With total energy conserved, another
conservation check is available and necessary for checking for a valid
simulation. 

