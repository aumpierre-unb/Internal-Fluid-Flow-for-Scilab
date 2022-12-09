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

function [Re,fD]=hDeps2fDRe(h,g,mu,rho,D,L,eps,varargin)
    // Solution for internal fluid flow based on pipe's hydraulic diameter and relative roughness
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
    // hDeps2fDRe computes the Reynolds number and 
    // the Darcy friction factor for a internal fluid flow, given 
    // the head loss h, 
    // the gravitational acceleration g, 
    // the fluid's dynamic viscosity mu and density rho, and 
    // the pipe's hydraulic diameter D, length L and relative roughness eps.
    // Inputs are to be given in a consistent system of units.
    //
    // Examples
    // ..// e.g. Compute the Reynolds number Re and
    // ..// the Darcy friction factor fD given
    // ..// the head loss h=40 cm,
    // ..// the gravitational acceleration g=981 cm/s/s,
    // ..// the fluid's the dynamic viscosity mu=0.0089 g/cm/s and
    // ..// density rho=0.98 g/cc, and
    // ..// the pipe's hydraulic diameter D=7 cm,
    // ..// length L=2500 cm and
    // ..// relative roughness eps=0.0025:
    // ..//
    // ..// This call computes Re e fD:
    // [Re,fD]=hDeps2fDRe(40,981,0.0089,0.98,7,2500,0.0025,%f)
    // ..// Alternatively:
    // h=40;.. //head loss (cm)
    // g=981;.. //gravitational acceleration (cm/s/s)
    // mu=0.0089;.. //fluid's dynamic viscosity (g/cm/s)
    // rho=0.98;.. //fluid's density (g/cc)
    // D=7;.. //pipe's hydraulic diameter (cm)
    // L=2500;.. //pipe's length (cm)
    // eps=0.0025;.. //pipe's relative roughness
    // [Re,fD]=hDeps2fDRe(h,g,mu,rho,D,L,eps)
    // ..// This call computes Re e fD
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // [Re,fD]=hDeps2fDRe(40,981,0.0089,0.98,7,2500,0.0025,%t)
    // ..// e.g. Compute the Reynolds number Re and
    // ..// the Darcy friction factor fD given
    // ..// the head loss h=40 cm,
    // ..// the gravitational acceleration g=981 cm/s/s,
    // ..// the fluid's the dynamic viscosity mu=0.0089 g/cm/s and
    // ..// density rho=0.98 g/cc, and
    // ..// the pipe's hydraulic diameter D=0.7 cm,
    // ..// length L=2500 cm and
    // ..// relative roughness eps=0.0025:
    // ..//
    // ..// This call computes Re e fD
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // [Re,fD]=hDeps2fDRe(40,981,0.0089,0.98,0.7,2500,0.0025,%t)
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
    islam=%t
    Re=K/64
    fD=64/Re
    if Re>2.3e3
        islam=%f
        Re=1e4
        fD=epsRe2fD(Re,eps)
        while abs(fD-K/Re^2)/fD>5e-3
            if fD-K/Re^2<0 Re=Re*1.02
            else
                Re=Re*0.98
                if Re<2.3e3
                    islam=%t
                    Re=K/64
                    fD=64/Re
                    warning("Solution found in extended laminar range.")
                    break
                end
            end
            fD=epsRe2fD(Re,eps)
        end
    end
    if (argn(2)==8 && varargin(1))
        if winsid()==[] scf(0)
        else scf(max(winsid())+1) end
        if islam
            laminar("r")
            turb(eps,"k")
        else
            laminar("k")
            turb(eps,"r")
        end
        if eps<1e-4, turb(1e-5,'k')
        else turb(eps/3,'k') end
        if eps<1e-4, turb(1e-4,'k')
        else turb(eps/10,'k') end
        if eps<1e-4, turb(1e-3,'k')
        elseif eps*3>5e-2, turb(5e-2,'k')
        else turb(eps*3,'k') end
        if eps<1e-4, turb(5e-3,'k')
        elseif eps*10>5e-2, turb(eps/6,'k')
        else turb(eps*10,'k') end
        rough("b")
        smooth("b")
        loglog(Re,fD,"rd")
        loglog([(K/1e-2)^(1/2) (K/1e-1)^(1/2)],[1d-2 1d-1],"--r")
        xgrid(33,1,7)
        xlabel("$Re={ {\rho v D} \over\displaystyle \mu }$",..
               "fontsize",4)
        ylabel("$f={{2 g h \rho^2 D^3} \over\displaystyle {\mu^2 L}} Re^{-2}$",..
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
