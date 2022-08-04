# Internal Fluid Flow

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This is a very short introduction to an incompressible and inviscid fluid steady internal flow and to the Internal FluidFlow Toolbox for Scilab.

Our focus here is a small set of equations that described the phenomenon and are required to solve problems on internal fluid flow. Fluid mechanics is a pretty extensive topic in fluid dynamics and there are a lot of important and interesting observations related to the topic that are not taken into account in this text, because they have do no direct impact the computation performed by the functions in this toolbox.

This text is divided in two main sections: The Theory and The Internal Fluid Flow Toolbox.

## The Theory

### The Bernoulli Equation

The Bernoulli equation is an expression of the mechanical energy balance for a very particular situation:

- internal steady flow of an
- incompressible inviscid fluid, where
- friction effects and tube fittings can be neglected.

For such a case, the mechanical energy is conserved, and for any two points 1 and 2 we have

$$
{\rho v_2^2 \over 2} + \rho g z_2 + p_2 =
{\rho v_1^2 \over 2} + \rho g z_1 + p_1
$$

or

$$
{v_2^2 \over 2g}+z_2+{p_2 \over \rho g}=
{v_1^2 \over 2g}+z_1+{p_1 \over \rho g}
$$

where

- *&rho;* is the fluid's density,
- *v* is the speed flow,
- *g* is the gravitational acceleration,
- *z* is the elevation, and
- *p* is the static pressure.

### Head Loss

The flow of viscous fluids is accompanied of energy dispersion, which can be measured as a pressure drop or, equivalently, as a head loss *h*, by the Darcy-Weisbach equation,

$$
h=f{v^2 \over 2g} {L \over D}
$$

where *f* is the Darcy friction factor, *L* is the pipe's length and *D* is the pipe's hydraulic diameter,

$$
D={4A \over P}
$$

where *A* is the cross-sectional area of the flow and *P* is the wet perimeter of the cross-section. *f* is described as a function of the Reynolds number,

$$
Re={\rho vg \over \mu}
$$

and the pipe's relative roughness,

$$
\varepsilon={k \over D}
$$

where

- *&mu;* is the fluid's dynamic viscosity and
- *k* is the pipe's[ internal surface] roughness.

The Reynolds number *Re*, the Darcy friction factor *f*, and the relative roughness *&epsilon;* completely describe the internal flow of incompressible viscous fluids, for both laminar and turbulent regimes. Usually, *f* is given as a function of *Re* and *&epsilon;*.

The simplest problems on internal fluid flow consist on computing one of them given the two other. More complex situations arise when only one or none of those variables is known. Instead, dimensional variables involved are given. However not always, in most cases iterative computation is required.

### Laminar Flow and Turbulent Flow

For laminar flow, the Darcy friction factor is given by the Poiseuille condition,

$$
f={64 \over Re}
$$

For turbulent flow, the Darcy friction factor is given implicitly by the Colebrook-White equation,

$$
{1 \over \sqrt{f}}=2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon} + {2.51 \over {Re \sqrt{f}}}}
$$

## The Internal Fluid Flow for Scilab Toolbox

This package provides a set of functions designed to solve problems of internal fluid flow. All functions are based on the Poiseuille condition for laminar flow, the Colebrook-White equation for turbulent flow, and the Darcy-Weisbach equation for head loss. The simplest problems on internal flow consist in computing either the Reynolds number or the Darcy friction factor given the other and the relative roughness. For those cases, this package provides functions epsRe2fD and epsfD2Re, respectively. More elaborated problems consist in computing both the Reynolds number and the Darcy friction factor given the head loss, the tube length, the fluid's density and dynamic viscosity, the gravitational acceleration, the relative roughness and either the dynamic diameter or the linear velocity or the volumetric flow. For those cases, this package provides functions hDeps2fRe, hveps2fRe and hQeps2fRe, respectively. A slightly more elaborate situation arises when roughness is given instead of relative roughness along with the linear velocity or the volumetric flow. For those cases, this package provides functions hvthk2fRe and hQthk2fRe, respectively. All function in this package offer the option of plotting the solution on a schematic Moody diagram.

Internal Fluid Flow Toolbox provides the following functions:

- epsRe2fD
- epsfD2Re
- hDeps2fDRe
- hveps2fDRe
- hQeps2fDRe
- hvthk2fDRe
- hQthk2fDRe

### epsRe2fD

epsRe2fD computes the Darcy friction factor *f* given the relative roughness *&epsilon;* and the Reynolds number *Re*. If given *Re* < 2500, then flow is assumed to be laminar and *f* is computed using of the Poiseuille condition. Otherwise, flow is assumed to be turbulent and *f* is computed using the Colebrook-White equation.

**Syntax:**

``[f]=epsRe2fD(Re,[eps[,s]])``

*e.g.* Compute the Darcy friction factor *f* given the Reynolds number *Re*=2.5e4 and the relative roughness *&epsilon;*=0.0044:

``--> f=epsRe2fD(2.5e4,0.0044,%f)``

or

``--> Re=2.5e4,eps=0.0044,f=epsRe2fD(Re,eps)``

*e.g.* Compute the Darcy friction factor *f* given the Reynolds number *Re*=2.5e4 and the relative roughness *&epsilon;*=0.0044 and plot a representation of the solution on a schematic Moody diagram:

``--> f=epsRe2fD(2.5e4,0.0044,%t)``

*e.g.* Compute the Darcy friction factor *f* given the Reynolds number *Re*=2.5e4:

``--> f=epsRe2fD(2.5e4)``

### epsfD2Re

espfD2Re computes the Reynolds number *Re* given the relative roughness *&epsilon;* and the Darcy friction factor *f*. Depending on the inputs, solution may be laminar or turbulent flow, or either for smooth pipes with higher friction, or none for lower friction and rough pipes. If the Poiseuille condition produces Re < 2500, laminar solution is accepted. If given *f* is possible for turbulent flow,

$$
{1 \over \sqrt{f}}=2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon}}
$$

(which is Colebrook-White equation for for elevated *Re*) the turbulent solution is accepted. If both solutions are accepted, espfD2Re returns both answers. If neither laminar or turbulent solutions are accepted, espfD2Re returns an empty matrix. If given *&epsilon;* > 0.05, execution is aborted.

**Syntax:**

``[Re]=epsfD2Re(f[,eps[,s]])``

*e.g.* Compute the Reynolds *Re* number given the Darcy friction factor *f*=0.033 and the relative roughness *&epsilon;*=0.0044:

``--> Re=epsfD2Re(0.033,0.0044,%f)``

or

``--> f=0.033,eps=0.0044,Re=epsfD2Re(f,eps)``

*e.g.* Compute the Reynolds *Re* number given the Darcy friction factor *f*=0.033 and the relative roughness *&epsilon;*=0.0044 and plot a representation of the solution on a schematic Moody diagram:

``--> Re=epsfD2Re(0.033,0.0044,%t)``

*e.g.* Compute the Reynolds number factor *f* given the Darcy friction *f*=0.033:

``--> Re=epsfD2Re(0.033)``

### hDeps2fDRe

hDeps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L*, relative roughness *&epsilon;* and hydraulic diameter *D*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing speed flow *v* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
Re^2 f={2gh\rho^2D^3 \over {\mu^2 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively, however an analytic solution is possible in this case.

**Syntax:**

``--> [Re,f]=hDeps2fDRe(h,g,mu,rho,D,L,eps[,s])``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h*=40 (cm),

the gravitational acceleration *h*=981 (cm/s/s),

the fluid's the dynamic viscosity *&mu;*=0.0089 (g/cm/s) and density *&rho;*=0.98 (g/cu.cm), and

the pipe's hydraulic diameter *D*=10 (cm), length *L*=2500 (cm) and relative roughness *&epsilon;*=0.0025:

``--> [Re,fD]=hDeps2fDRe(40,981,0.0089,0.98,10,2500,0.0025,%f)``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

``--> [Re,f]=hDeps2fDRe(40,981,0.0089,0.98,10,2500,0.0025,%t)``

### hveps2fDRe

hveps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness *&epsilon;*, the speed flow *v*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``--> [Re,f]=hveps2fDRe(h,g,mu,rho,v,L,eps[,s])``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h*=40 (cm),

the gravitational acceleration *g*=981 (cm/s/s),

the fluid's the dynamic viscosity *&mu;*=0.0089 (g/cm/s) and density *&rho;*=0.98 (g/cu.cm),

the speed flow *v*=110 (cm/s), and

the pipe's length *L*=2500 (cm) and relative roughness *&epsilon;*=0.0025:

``--> [Re,f]=hveps2fDRe(40,981,0.0089,0.98,110,2500,0.0025,%f)``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

``--> [Re,f]=hveps2fDRe(40,981,0.0089,0.98,110,2500,0.0025,%t)``

### hQeps2fDRe

hQeps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness *&epsilon;*, the volumetric flow rate Q, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``--> [Re,f]=hQeps2fDRe(h,g,mu,rho,Q,L,eps[,s])``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h*=40 (cm),

the gravitational acceleration *g*=981 (cm/s/s),

the fluid's the dynamic viscosity *&mu;*=0.0089 (g/cm/s) and density *&rho;*=0.98 (g/cu.cm),

the volumetric flow rate *Q*=8666 (cu. cm/s), and

the pipe's length *L*=2500 (cm) and relative roughness *&epsilon;*=0.0025:

``--> [Re,f]=hQeps2fDRe(40,981,0.0089,0.98,8666,2500,0.0025,%f)``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

``--> [Re,f]=hQeps2fDRe(40,981,0.0089,0.98,8666,2500,0.0025,%t)``

### hvthk2fDRe

hvthk2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the speed flow *v*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``--> [Re,f]=hvthk2fDRe(h,g,mu,rho,v,L,thk[,s])``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h*=40 (cm),

the gravitational acceleration *g*=981 (cm/s/s),

the fluid's the dynamic viscosity *&mu;*=0.0089 (g/cm/s) and density *&rho;*=0.98 (g/cu.cm),

the speed flow *v*=110 (cm/s), and

the pipe's length *L*=2500 (cm) and roughness *k*=0.025 (cm):

``--> [Re,f]=hvthk2fDRe(40,981,0.0089,0.98,110,2500,0.025,%f)``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

``--> [Re,f]=hvthk2fDRe(40,981,0.0089,0.98,110,2500,0.025,%t)``

### hQthk2fDRe

hQthk2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the volumetric flow rate Q, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``--> [Re,f]=hQthk2fDRe(h,g,mu,rho,Q,L,thk[,s])``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h*=40 (cm),

the gravitational acceleration *g*=981 (cm/s/s),

the fluid's the dynamic viscosity *&mu;*=0.0089 (g/cm/s) and density *&rho;*=0.98 (g/cu.cm),

the volumetric flow rate *Q*=8666 (cu. cm/s), and

the pipe's length *L*=2500 (cm) and roughness *k*=0.025:

``--> [Re,f]=hQthk2fDRe(40,981,0.0089,0.98,8666,2500,0.025,%f)``

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

``--> [Re,f]=hQthk2fDRe(40,981,0.0089,0.98,8666,2500,0.025,%t)``

Copyright &copy; 2022 Alexandre Umpierre

email: aumpierre@gmail.com
