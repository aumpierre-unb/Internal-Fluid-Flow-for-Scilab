// Copyright (C) 2022 2023 Alexandre Umpierre
// This file is part of internal-fluid-flow Toolbox.
// internal-fluid-flow Toolbox is free software:
// you can redistribute it and/or modify it under the terms
// of the GNU General Public License (GPL) version 3
// as published by the Free Software Foundation.
// internal-fluid-flow Toolbox is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the
// GNU General Public License along with this program.
// It is also available at www.gnu.org/licenses/.

function [Re,fD]=hDeps2fDRe(h,g,mu,rho,D,L,eps,varargin)
    // hDeps2fDRe computes
    // the Reynolds number and
    // the Darcy friction factor based on
    // the pipe's hydraulic diameter and
    // the pipe's relative roughness
    //
    // Syntax
    // [Re,fD]=hDeps2fDRe(h,g,mu,rho,D,L,eps[,fig])
    //
    // Parameters
    // h: head loss
    // g: gravitational acceleration
    // mu: fluid's dynamic viscosity
    // rho: fluid's density
    // D: pipe's hydraulic diameter
    // L: pipe's length
    // eps: pipe's relative roughness
    // fig: optional, boolean for display plot (default is fig=%f)
    // Re: Reynolds number
    // fD: Darcy friction factor
    //
    // Description
    // hDeps2fDRe computes
    // the Reynolds number and
    // the Darcy friction factor for a internal fluid flow given 
    // the head loss h, 
    // the gravitational acceleration g, 
    // the fluid's dynamic viscosity mu and density rho, and 
    // the pipe's hydraulic diameter D, length L and relative roughness eps.
    // Inputs are to be given in a consistent system of units.
    //
    // Examples
    // // Compute the Reynolds number Re and
    // // the Darcy friction factor fD given
    // // the head loss h=40 cm,
    // // the gravitational acceleration g=981 cm/s/s,
    // // the fluid's the dynamic viscosity mu=8.9e-3 g/cm/s and
    // // density rho=0.98 g/cc, and
    // // the pipe's hydraulic diameter D=7 cm,
    // // length L=2.5e3 cm and
    // // relative roughness eps=2.5e-3:
    // //
    // // This call computes Re e fD:
    // [Re,fD]=hDeps2fDRe(40,981,8.9e-3,0.98,7,2.5e3,2.5e-3,%f)
    // // Alternatively:
    // h=40;.. //head loss (cm)
    // g=981;.. //gravitational acceleration (cm/s/s)
    // mu=8.9e-3;.. //fluid's dynamic viscosity (g/cm/s)
    // rho=0.98;.. //fluid's density (g/cc)
    // D=7;.. //pipe's hydraulic diameter (cm)
    // L=2.5e3;.. //pipe's length (cm)
    // eps=2.5e-3;.. //pipe's relative roughness
    // [Re,fD]=hDeps2fDRe(h,g,mu,rho,D,L,eps)
    // // This call computes Re e fD
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // [Re,fD]=hDeps2fDRe(40,981,8.9e-3,0.98,7,2.5e3,2.5e-3,%t)
    // // Compute the Reynolds number Re and
    // // the Darcy friction factor fD given
    // // the head loss h=40 cm,
    // // the gravitational acceleration g=981 cm/s/s,
    // // the fluid's the dynamic viscosity mu=8.9e-3 g/cm/s and
    // // density rho=0.98 g/cc, and
    // // the pipe's hydraulic diameter D=0.7 cm,
    // // length L=2.5e3 cm and
    // // relative roughness eps=2.5e-3:
    // //
    // // This call computes Re e fD
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // [Re,fD]=hDeps2fDRe(40,981,8.9e-3,0.98,0.7,2.5e3,2.5e-3,%t)
    //
    // See also
    //  epsfD2Re
    //  epsRe2fD
    //  hveps2fDRe
    //  hvthk2fDRe
    //  hQeps2fDRe
    //  hQthk2fDRe
    //
    // Authors
    //  Alexandre Umpierre

    K=2*g*h*rho^2*D^3/mu^2/L
    function y=foo(fD)
        y=1/fD^.5+2*log10(eps/3.7+2.51/(K/fD)^.5/fD^.5)
    endfunction
    fD=root(foo,1d-2,1e-1,1d-4)
    Re=(K/fD)^.5
    if Re>2.3e3
        islam=%f
    else
        Re=K/64
        fD=64/Re
        islam=%t
    end
    if argn(2)==8 && varargin(1)
        if winsid()==[]
            scf(0)
        else
            scf(max(winsid())+1)
        end
        x=[5e-2 2.5e-2 1e-2 3e-3 1e-3 3e-4 1e-4]
        for i=1:length(x)
            turbulent(x(i),"k")
        end
        rough("-.b")
        if eps~=0
            smoothline("-.b")
        end
        if islam
            laminar("r")
        else
            laminar("k")
            turbulent(eps,"r")
        end
        loglog(Re,fD,"rd")
        loglog([(K/1e-2)^.5 (K/1e-1)^.5],[1d-2 1d-1],"--r")
        xgrid(33,1,7)
        xlabel("$Re={ {\rho v D} \over\displaystyle \mu }$",..
               "fontsize",4)
        ylabel("$f={{2 g h \rho^2 D^3} \over\displaystyle {\mu^2 L}} Re^{-2}$",..
               "fontsize",4)
        gca().tight_limits=["on","on"]
        gca().data_bounds=[1d2 1d8 6e-3 1d-1]
        gca().grid=[1,1]
        gca().grid_style=[9,9]
        gcf().figure_size=[600,600]
    end
endfunction
