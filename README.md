# Internal Fluid Flow

[![DOI](https://zenodo.org/badge/509430202.svg)](https://zenodo.org/badge/latestdoi/509430202)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/aumpierre-unb/Internal-Fluid-Flow-for-Scilab)

## Install and load Internal Fluid Flow Toolbox for Scilab

1\. Download ```v0.2.1.zip``` file from [Releases](https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-Scilab/releases).

2\. Set Scilab current directory to the download folder.

3\. To install and load the toolbox, type:

```dotnetcli
atomsInstall('v0.2.1.zip'),..
atomsLoad('inflowlib')
```

## Citation of Internal Fluid Flow Toolbox for Scilab

You can cite all versions (both released and pre-released), by using
[DOI 105281/zenodo.6927994](https://doi.org/10.5281/zenodo.6927994).

This DOI represents all versions, and will always resolve to the latest one.

To cite the last released version, please check
https://zenodo.org/account/settings/github/repository/aumpierre-unb/Internal-Fluid-Flow-for-Scilab.

---

The following is a very short introduction to the steady internal flow of an incompressible and inviscid fluid and to the Internal Fluid Flow Toolbox for Scilab.

Internal flow is a pretty extensive topic in fluid mechanics and there are a lot of important and interesting observations related to it that are not taken into account in this text, because they have no direct impact the computation performed by the functions in this toolbox. Our focus here is a small set of equations that described the phenomenon and are required to solve problems on internal fluid flow.

This text is divided in two main sections: The Theory and The Internal Fluid Flow Toolbox for Scilab.

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
- *v* is the flow speed,
- *g* is the gravitational acceleration,
- *z* is the elevation, and
- *p* is the static pressure.

### Head Loss

The flow of viscous fluids is accompanied of energy dispersion, which can be measured as pressure drop or, equivalently, as head loss *h*, by the Darcy-Weisbach equation,

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

For laminar flow, *Re* < 2,300 (typically), the Darcy friction factor is given by the Poiseuille condition,

$$
f={64 \over Re}
$$

For turbulent flow, *Re* > 2,300 (typically), the Darcy friction factor is given implicitly by the Colebrook-White equation,

$$
{1 \over \sqrt{f}}=2\ \mathrm{log} {1 \over\displaystyle {\varepsilon \over 3.7} + {2.51 \over {Re \sqrt{f}}}}
$$

## The Internal Fluid Flow Toolbox for Scilab

Internal Fluid Flow provides the following functions:

- epsRe2fD
- epsfD2Re
- hDeps2fDRe
- hveps2fDRe
- hQeps2fDRe
- hvthk2fDRe
- hQthk2fDRe

### epsRe2fD

epsRe2fD computes the Darcy friction factor *f* given the relative roughness *&epsilon;* and the Reynolds number *Re*. If given *Re* < 2,500, then flow is assumed to be laminar and *f* is computed using of the Poiseuille condition. Otherwise, flow is assumed to be turbulent and *f* is computed using the Colebrook-White equation.

**Syntax:**

```dotnetcli
[f]=epsRe2fD(Re,[eps[,s]])
```

*e.g.* Compute the Darcy friction factor *f* given the Reynolds number *Re*=25,000 and the relative roughness *&epsilon;* = 0.0044:

```dotnetcli
f=epsRe2fD(2.5e4,4.4e-3,%f)
```

or

```dotnetcli
Re=2.5e4,eps=4.4e-3,f=epsRe2fD(Re,eps)
```

*e.g.* Compute the Darcy friction factor *f* given the Reynolds number *Re*=25,000 and the relative roughness *&epsilon;* = 0.0044 and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
f=epsRe2fD(2.5e4,4.4e-3,%t)
```

*e.g.* Compute the Darcy friction factor *f* given the Reynolds number *Re*=25,000:

```dotnetcli
f=epsRe2fD(2.5e4)
```

### epsfD2Re

espfD2Re computes the Reynolds number *Re* given the relative roughness *&epsilon;* and the Darcy friction factor *f*. Depending on the inputs, solution may be laminar or turbulent flow, or either for smooth pipes with higher friction, or none for lower friction and rough pipes. If the Poiseuille condition produces *Re* < 2,500, laminar solution is accepted. If given *f* is possible for turbulent flow,

$$
{1 \over \sqrt{f}}<2\ \mathrm{log} {1 \over\displaystyle {\varepsilon \over 3.7}}
$$

(which is Colebrook-White equation for for elevated *Re*) the turbulent solution is accepted. If both solutions are accepted, espfD2Re returns both answers. If neither laminar or turbulent solutions are accepted, espfD2Re returns an empty matrix. If given *&epsilon;* > 0.05, execution is aborted.

**Syntax:**

```dotnetcli
[Re]=epsfD2Re(f[,eps[,s]])
```

*e.g.* Compute the Reynolds *Re* number given the Darcy friction factor *f* = 0.033 and the relative roughness *&epsilon;* = 0.0044:

```dotnetcli
Re=epsfD2Re(3.3e-2,4.4e-3,%f)
```

or

```dotnetcli
f=3.3e-2,eps=4.4e-3,Re=epsfD2Re(f,eps)
```

*e.g.* Compute the Reynolds *Re* number given the Darcy friction factor *f* = 0.033 and the relative roughness *&epsilon;* = 0.0044 and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
Re=epsfD2Re(3.3e-2,4.4e-3,%t)
```

*e.g.* Compute the Reynolds number factor *f* given the Darcy friction *f* = 0.033:

```dotnetcli
Re=epsfD2Re(3.3e-2)
```

### hDeps2fDRe

hDeps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L*, relative roughness *&epsilon;* and hydraulic diameter *D*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing flow speed *v* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
Re^2 f={2gh\rho^2D^3 \over {\mu^2 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively, however an analytic solution is possible in this case.

**Syntax:**

```dotnetcli
[Re,f]=hDeps2fDRe(h,g,mu,rho,D,L,eps[,s])
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h* = 0.40 m,

the gravitational acceleration *h* = 9.81 m/s/s,

the fluid's the dynamic viscosity *&mu;* = 0.89 cP and density *&rho;* = 0.98 kg/L, and

the pipe's hydraulic diameter *D*=10 (cm), length *L* = 25 m and relative roughness *&epsilon;* = 0.0025:

```dotnetcli
[Re,fD]=hDeps2fDRe(40,981,8.9e-3,0.98,10,2.5e3,2.5e-3,%f)
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
[Re,f]=hDeps2fDRe(40,981,8.9e-3,0.98,10,2.5e3,2.5e-3,%t)
```

### hveps2fDRe

hveps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness *&epsilon;*, the flow speed *v*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hveps2fDRe(h,g,mu,rho,v,L,eps[,s])
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h* = 0.40 m,

the gravitational acceleration *g* = 9.81 m/s/s,

the fluid's the dynamic viscosity *&mu;* = 0.89 cP and density *&rho;* = 0.98 kg/L,

the flow speed *v* = 1.1 m/s, and

the pipe's length *L* = 25 m and relative roughness *&epsilon;* = 0.0025:

```dotnetcli
[Re,f]=hveps2fDRe(40,981,8.9e-3,0.98,1.1e2,2.5e3,2.5e-3,%f)
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
[Re,f]=hveps2fDRe(40,981,8.9e-3,0.98,1.1e2,2.5e3,2.5e-3,%t)
```

### hQeps2fDRe

hQeps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness *&epsilon;*, the volumetric flow rate *Q*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hQeps2fDRe(h,g,mu,rho,Q,L,eps[,s])
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h* = 0.40 m,

the gravitational acceleration *g* = 9.81 m/s/s,

the fluid's the dynamic viscosity *&mu;* = 0.89 cP and density *&rho;* = 0.98 kg/L,

the volumetric flow rate *Q* = 8.6 L/s, and

the pipe's length *L* = 25 m and relative roughness *&epsilon;* = 0.0025:

```dotnetcli
[Re,f]=hQeps2fDRe(40,981,8.9e-3,0.98,8.6e3,2.5e3,2.5e-3,%f)
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
[Re,f]=hQeps2fDRe(40,981,8.9e-3,0.98,8.6e3,2.5e3,2.5e-3,%t)
```

### hvthk2fDRe

hvthk2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the flow speed *v*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hvthk2fDRe(h,g,mu,rho,v,L,thk[,s])
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h* = 0.40 m,

the gravitational acceleration *g* = 9.81 m/s/s,

the fluid's the dynamic viscosity *&mu;* = 0.89 cP and density *&rho;* = 0.98 kg/L,

the flow speed *v* = 1.1 m/s, and

the pipe's length *L* = 25 m and roughness *k* = 0.25 mm:

```dotnetcli
[Re,f]=hvthk2fDRe(40,981,8.9e-3,0.98,1.1e2,2.5e3,2.5e-2,%f)
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
[Re,f]=hvthk2fDRe(40,981,8.9e-3,0.98,1.1e2,2.5e3,2.5e-2,%t)
```

### hQthk2fDRe

hQthk2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the volumetric flow rate *Q*, the gravitational acceleration *g*, and the fluid's density *&rho;* and dynamic viscosity *&mu;*. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hQthk2fDRe(h,g,mu,rho,Q,L,thk[,s])
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given

the head loss *h* = 0.40 m,

the gravitational acceleration *g* = 9.81 m/s/s,

the fluid's the dynamic viscosity *&mu;* = 0.89 cP and density *&rho;* = 0.98 kg/L,

the volumetric flow rate *Q* = 8.6 L/s, and

the pipe's length *L* = 25 m and roughness *k* = 0.25 mm:

```dotnetcli
[Re,f]=hQthk2fDRe(40,981,8.9e-3,0.98,8.6e3,2.5e3,2.5e-2,%f)
```

*e.g.* Compute the Reynolds number *Re* and the Darcy friction factor *f* given the same inputs and plot a representation of the solution on a schematic Moody diagram:

```dotnetcli
[Re,f]=hQthk2fDRe(40,981,8.9e-3,0.98,8.6e3,2.5e3,2.5e-2,%t)
```

Copyright &copy; 2022 2023 Alexandre Umpierre

email: aumpierre@gmail.com
