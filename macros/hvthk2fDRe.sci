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

function [Re,fD]=hvthk2fDRe(h,g,mu,rho,v,L,thk,varargin)
    // hvthk2fDRe computes
    // the Reynolds number and
    // the Darcy friction factor based on
    // the flow speed and
    // the pipe's roughness
    //
    // Syntax
    // [Re,fD]=hvthk2fDRe(h,g,mu,rho,v,L,thk[,fig])
    //
    // Parameters
    // h: head loss
    // g: gravitational acceleration
    // mu: fluid's dynamic viscosity
    // rho: fluid's density
    // v: flow speed
    // L: pipe's length
    // thk: pipe's roughness
    // fig: optional, boolean for display plot (default is fig=%f)
    // Re: Reynolds number
    // fD: Darcy friction factor
    //
    // Description
    // hQthk2fDRe computes
    // the Reynolds number and
    // the Darcy friction factor for a internal fluid flow given
    // the head loss h,
    // the gravitational acceleration g,
    // the fluid's dynamic viscosity mu and density rho, and
    // the flow speed v, and
    // the pipe's length L and roughness thk.
    // Inputs are to be given in a consistent system of units.
    // hvthk2fDRe is a main function of
    // the internal-fluid-flow toolbox for Scilab.
    //
    // Examples
    // // Compute the Reynolds number Re and
    // // the Darcy friction factor fD given
    // // the head loss h=40 cm,
    // // the gravitational acceleration g=981 cm/s/s,
    // // the fluid's the dynamic viscosity mu=8.9e-3 g/cm/s and
    // // density rho=0.98 g/cc,
    // // the flow speed v=1.1e2 cm/s and
    // // the pipe's length L=2.5e3 cm and
    // // roughness thk=2.5e-2 cm:
    // //
    // // This call computes Re e fD:
    // [Re,fD]=hvthk2fDRe(40,981,8.9e-3,0.98,1.1e2,2.5e3,2.5e-2,%f)
    // // Alternatively:
    // h=40;.. //head loss
    // g=981;.. //gravitational acceleration (cm/s/s)
    // mu=8.9e-3;.. //fluid's dynamic viscosity (g/cm/s)
    // rho=0.98;.. //fluid's density (g/cc)
    // v=1.1e2;.. //flow speed (cm/s)
    // L=2.5e3;.. //pipe's length (cm)
    // thk=2.5e-2;.. //pipe's roughness (cm)
    // [Re,fD]=hvthk2fDRe(h,g,mu,rho,v,L,thk)
    // // This call computes Re e fD
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // [Re,fD]=hvthk2fDRe(40,981,8.9e-3,0.98,1.1e2,2.5e3,2.5e-2,%t)
    // // Compute the Reynolds number Re and
    // // the Darcy friction factor fD given
    // // the head loss h=40 cm,
    // // the gravitational acceleration g=981 cm/s/s,
    // // the fluid's the dynamic viscosity mu=8.9e-3 g/cm/s and
    // // density rho=0.98 g/cc,
    // // the flow speed v=20 cm/s and
    // // the pipe's length L=2.5e3 cm and
    // // roughness thk=2.5e-3 cm:
    // //
    // // This call computes Re e fD
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // [Re,fD]=hvthk2fDRe(40,981,8.9e-3,0.98,20,2.5e3,2.5e-3,%t)
    // // Compute the Reynolds number Re and
    // // the Darcy friction factor fD given
    // // the head loss h=30 cm,
    // // the gravitational acceleration g=981 cm/s/s,
    // // the fluid's the dynamic viscosity mu=8.9e-3 g/cm/s and
    // // density rho=0.98 g/cc,
    // // the flow speed v=26 cm/s and
    // // the pipe's length L=2.5e3 cm and
    // // roughness thk=2.5e-2 cm:
    // //
    // // This call computes Re e fD
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // [Re,fD]=hvthk2fDRe(40,981,8.9e-3,0.98,26,2.5e3,2.5e-2,%t)
    //
    // See also
    //  epsfD2Re
    //  epsRe2fD
    //  hDeps2fDRe
    //  hveps2fDRe
    //  hQeps2fDRe
    //  hQthk2fDRe
    //
    // Authors
    //  Alexandre Umpierre

    Re=[]
    fD=[]
    M=2*g*mu*h/v^3/rho/L
    function y=foo(fD)
        y=1/fD^.5+2*log10(..
        thk/(fD/M*mu/rho/v).. // this is eps
        /3.7+2.51/(fD/M)/fD^.5)
    endfunction
    f_=root(foo,1d-2,1e-1,1d-4)
    Re_=f_/M
    isturb=%f
    if Re_>2.3e3
        Re=[Re;Re_]
        fD=[fD;f_]
        D=Re_*mu/rho/v
        eps=thk/D
        isturb=%t
    end
    Re_=(64/M)^.5
    if Re_<2.3e3
        Re=[Re;Re_]
        fD=[fD;64/Re_]
    end
    if ~isturb
        D = Re_ * mu / rho / v
        eps = thk / D
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
        if thk~=0
            smoothline("-.b")
        end
        if min(Re)<2.3e3
            laminar("r")
        else
            laminar("k")
        end
        if max(Re)>2.3e3
            turbulent(eps,"r")
        else
            turbulent(eps,"k")
        end
        loglog(Re,fD,"rd")
        loglog([1e-2/M 1e-1/M],[1d-2 1d-1],"--r")
        xgrid(33,1,7)
        xlabel("$Re={ {\rho v D} \over\displaystyle \mu }$",..
               "fontsize",4)
        ylabel("$f={{2 g h \mu} \over\displaystyle {\rho v^3 L}} Re$",..
               "fontsize",4)
        gca().tight_limits=["on","on"];
        gca().data_bounds=[1d2 1d8 6e-3 1d-1]
        gca().grid=[1,1]
        gca().grid_style=[9,9]
        gcf().figure_size=[600,600]
    end
endfunction
