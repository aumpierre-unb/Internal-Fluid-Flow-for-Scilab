/*
Copyright (C) 2022 Alexandre Umpierre

This file is part of Internal Fluid Flow Toolbox.
Internal Fluid Flow Toolbox is free software:
you can redistribute it and/or modify it under the terms
of the GNU General Public License (GPL) version 3
as published by the Free Software Foundation.

Internal Fluid Flow Toolbox is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the
GNU General Public License along with this program
(license.txt).
It is also available at www.gnu.org/licenses/.
*/

function [Re,fD]=hvthk2fDRe(h,g,mu,rho,v,L,thk,varargin)
    // Solution for internal fluid flow based on flow speed and pipe's roughness
    //
    // Syntax
    // [Re,fD]=hvthk2fDRe(h,g,mu,rho,v,L,thk[,s])
    //
    // Parameters
    // h: head loss
    // g: gravitational acceleration
    // mu: fluid's dynamic viscosity
    // rho: fluid's density
    // v: flow speed
    // L: pipe's length
    // thk: pipe's roughness
    // s: optional, boolean for display plot (default is s=%f)
    // Re: Reynolds number
    // fD: Darcy friction factor
    //
    // Description
    // hQthk2fDRe computes the Reynolds number and 
    // the Darcy friction factor for a internal fluid flow, given 
    // the head loss h, 
    // the gravitational acceleration g, 
    // the fluid's dynamic viscosity mu and density rho, and 
    // the flow speed v, and
    // the pipe's length L and roughness thk.
    // Inputs are to be given in a consistent system of units.
    //
    // Examples
    // ..// e.g. Compute the Reynolds number Re and
    // ..// the Darcy friction factor fD given
    // ..// the head loss h=40 cm,
    // ..// the gravitational acceleration g=981 cm/s/s,
    // ..// the fluid's the dynamic viscosity mu=0.0089 g/cm/s and
    // ..// density rho=0.98 g/cc,
    // ..// the flow speed v=100 cm/s and
    // ..// the pipe's length L=2500 cm and
    // ..// roughness thk=0.025 cm:
    // ..//
    // ..// This call computes Re e fD:
    // [Re,fD]=hvthk2fDRe(40,981,0.0089,0.98,100,2500,0.025,%f)
    // ..// Alternatively:
    // h=40;.. //head loss
    // g=981;.. //gravitational acceleration (cm/s/s)
    // mu=0.0089;.. //fluid's dynamic viscosity (g/cm/s)
    // rho=0.98;.. //fluid's density (g/cc)
    // v=100;.. //flow speed (cm/s)
    // L=2500;.. //pipe's length (cm)
    // thk=0.025;.. //pipe's roughness (cm)
    // [Re,fD]=hvthk2fDRe(h,g,mu,rho,v,L,thk)
    // ..// This call computes Re e fD
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // [Re,fD]=hvthk2fDRe(40,981,0.0089,0.98,100,2500,0.025,%t)
    // ..// e.g. Compute the Reynolds number Re and
    // ..// the Darcy friction factor fD given
    // ..// the head loss h=40 cm,
    // ..// the gravitational acceleration g=981 cm/s/s,
    // ..// the fluid's the dynamic viscosity mu=0.0089 g/cm/s and
    // ..// density rho=0.98 g/cc,
    // ..// the flow speed v=20 cm/s and
    // ..// the pipe's length L=2500 cm and
    // ..// roughness thk=0.0025 cm:
    // ..//
    // ..// This call computes Re e fD
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // [Re,fD]=hvthk2fDRe(40,981,0.0089,0.98,20,2500,0.0025,%t)
    // ..// e.g. Compute the Reynolds number Re and
    // ..// the Darcy friction factor fD given
    // ..// the head loss h=30 cm,
    // ..// the gravitational acceleration g=981 cm/s/s,
    // ..// the fluid's the dynamic viscosity mu=0.0089 g/cm/s and
    // ..// density rho=0.98 g/cc,
    // ..// the flow speed v=26 cm/s and
    // ..// the pipe's length L=2500 cm and
    // ..// roughness thk=0.0025 cm:
    // ..//
    // ..// This call computes Re e fD
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // [Re,fD]=hvthk2fDRe(40,981,0.0089,0.98,26,2500,0.0025,%t)
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

    M=2*g*mu*h/v^3/rho/L
    r=%t
    Re=1e4
    fD=M*Re
    D=Re*mu/v/rho
    eps=thk/D
    fD=epsRe2fD(Re,eps)
    while abs(fD-Re*M)/fD>5e-3
        if fD-Re*M>0 Re=Re*1.02
        else
            Re=Re*0.98
            if Re<2.3e3
                Re=sqrt(64/M)
                fD=epsRe2fD(Re,eps)
                D=Re*mu/v/rho
                eps=thk/D
                r=%f
                break
            end
        end
        fD=M*Re
        D=Re*mu/v/rho
        eps=thk/D
        fD=epsRe2fD(Re,eps)
    end
    if r && sqrt(64/M)<2.3e3
        Re=[sqrt(64/M);Re]
        fD=[64/sqrt(64/M);fD]
    end
    if size(Re,1)==2
        warning("Laminar and turbulent solutions found.")
    end
    if (argn(2)==8 && varargin(1))
        if winsid()==[] scf(0)
        elseif scf(max(winsid())+1) end
        if min(Re)<2.3e3 laminar("r")
        else laminar("k") end
        if max(Re)>2.3e3 turb(eps,"r")
        else turb(eps,"k") end
        if eps*3<5e-2 turb(eps*3,"k")
        else turb(eps/2,"k") end
        if eps*10<5e-2 turb(eps*10,"k")
        else turb(eps/7,"k") end
        turb(eps/3,"k")
        turb(eps/10,"k")
        rough("b")
        smooth("b")
        loglog(Re,fD,"rd")
        loglog([1e-2/M 1e-1/M],[1d-2 1d-1],"--r")
        xgrid(33,1,7)
        xlabel("$Re={ {\rho v D} \over\displaystyle \mu }$",..
               "fontsize",4)
        ylabel("$f={{2 g h \mu} \over\displaystyle {\rho v^3 L}} Re$",..
               "fontsize",4)
        gca().data_bounds=[1d2 1d8 1d-2 1d-1]
        gca().grid=[1,1]
        gca().grid_style=[9,9]
        gcf().figure_size=[600,600]
    end
endfunction

function laminar(t)
    Re=[5e2 4e3]
    f=64 ./ Re
    loglog(Re,f,t)
endfunction

function turb(eps,t)
    N=51
    for i=1:N
        w=log10(2d3)+(i-1)*(log10(1d8)-log10(2d3))/(N-1)
        Re(i)=10^w
        function y=foo(fD)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re(i)/sqrt(fD))
        endfunction
        f(i)=root(foo,6e-4,1e-1,1e-4)
    end
    loglog(Re,f,t)
endfunction

function smooth(t)
    N=51
    for i=1:N
        w=log10(2d3)+(i-1)*(log10(1d7)-log10(2d3))/(N-1)
        Re(i)=10^w
        function y=foo(fD)
            y=1/sqrt(fD)+2*log10(2.51/Re(i)/sqrt(fD))
        endfunction
        f(i)=root(foo,6e-3,1e-1,1e-4)
    end
    loglog(Re,f,t)
end

function rough(t)
    eps=[]
    f=[]
    Re=[]
    N=31
    for i=1:N
        w=log10(4e-5)+(i-1)*(log10(5e-2)-log10(4e-5))/(N-1)
        eps=[eps;10^w]
        f=[f;1.01*(2*log10(3.7/eps($)))^-2]
        z=epsfD2Re(f($),eps($))
        Re=[Re;z($)]
    end
    loglog(Re,f,t)
end

function x2=root(f,x1,x2,tol)
    while abs(f(x2))>tol
        x=(x1+x2)/2
        if f(x)*f(x1)>0 x1=x
        else x2=x end
    end
end
