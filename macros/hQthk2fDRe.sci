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

function [Re,fD]=hQthk2fDRe(h,g,mu,rho,Q,L,thk,varargin)
    // Solution for internal fluid flow based on volumetric flow rate and pipe's roughness
    //
    // Syntax
    // [Re,fD]=hQthk2fDRe(h,g,mu,rho,Q,L,thk[,s])
    //
    // Parameters
    // h: head loss
    // g: gravitational acceleration
    // mu: fluid's dynamic viscosity
    // rho: fluid's density
    // Q: volumetric flow rate
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
    // the volumetric flow rate Q, and
    // the pipe's length L and roughness thk.
    // Inputs are to be given in a consistent system of units.
    //
    // Examples
    // ..// e.g. Compute the Reynolds number Re and
    // ..// the Darcy friction factor fD given
    // ..// the head loss h=40 (cm),
    // ..// the gravitational acceleration g=981 (cm/s/s),
    // ..// the fluid's the dynamic viscosity mu=0.0089 (g/cm/s) and
    // ..// density rho=0.98 (g/cu.cm),
    // ..// the volumetric flow rate Q=8666 (cu.cm/s) and
    // ..// the pipe's length L=2500 (cm) and
    // ..// roughness thk=0.025 (cm):
    // ..//
    // ..// This call computes Re e fD:
    // [Re,fD]=hQthk2fDRe(40,981,0.0089,0.98,8666,2500,0.025,%f)
    // ..// Alternatively:
    // h=40;.. //head loss (cm)
    // g=981;.. //gravitational acceleration (cm/s/s)
    // mu=0.0089;.. //fluid's dynamic viscosity (g/cm/s)
    // rho=0.98;.. //fluid's density (g/cu.cm)
    // Q=8666;.. //volumetric flow rate (cu.cm/s)
    // L=2500;.. //pipe's length (cm)
    // thk=0.025;.. //pipe's roughness (cm)
    // [Re,fD]=hQthk2fDRe(h,g,mu,rho,Q,L,thk)
    // ..// This call computes Re e fD
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // [Re,fD]=hQthk2fDRe(40,981,0.0089,0.98,8666,2500,0.025,%t)
    //
    // See also
    //  epsfD2Re
    //  epsRe2fD
    //  hDeps2fDRe
    //  hveps2fDRe
    //  hvthk2fDRe
    //  hQeps2fDRe
    //
    // Authors
    //  Alexandre Umpierre

    P=2*g*h*Q^3/(%pi/4)^3/(mu/rho)^5/L
    Re=(P/64)^(1/4)
    fD=64/Re
    if Re>25e2
        Re=1e4
        fD=P/Re^5
        D=rho*Q/Re/mu/(%pi/4)
        eps=thk/D
        fD=epsRe2fD(Re,eps)
        while abs(fD-P/Re^5)/fD>5e-3
            if fD-P/Re^5<0 Re=Re*1.02
            else Re=Re*0.98 end
            fD=P/Re^5
            D=rho*Q/Re/mu/(%pi/4)
            eps=thk/D
            fD=epsRe2fD(Re,eps)
        end
    end
    if (argn(2)==8 && varargin(1))
        if winsid()==[] scf(0)
        elseif scf(max(winsid())+1) end
        laminar()
        turb(eps)
        turb(eps*5)
        turb(eps*10)
        turb(eps/5)
        turb(eps/10)
        rough()
        xgrid(33,1,7)
        plot(Re,fD,"rd")
        plot([Re/10 Re*10],[P/(Re/10)^5 P/(Re*10)^5],"--r")
        xlabel("$Re=\frac{\rho vD}{\mu }$",..
               "fontsize",4)
        ylabel("$f=\frac{h}{\frac{v^{2}}{2g}\frac{L}{D}}$",..
               "fontsize",4)
        gca().data_bounds=[1d3 1d8 1d-2 1d-1]
        gca().grid=[1,1]
        gca().grid_style=[9,9]
        gcf().figure_size=[600,600]
    end
endfunction

function laminar()
    Re=[5e2 4e3]
    f=64 ./ Re
    plot2d("ll",Re,f)
endfunction

function turb(eps)
    N=50
    for i=1:N
        w=log10(2d3)+i*(log10(1d8)-log10(2d3))/N
        Re(i)=10^w
        function y=foo(fD)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re(i)/sqrt(fD))
        endfunction
        f(i)=root(foo,6e-3,1e-1,1e-4)
    end
    plot2d("ll",Re,f)
endfunction

function rough()
    eps=[]
    f=[]
    Re=[]
    N=30
    for i=1:N
        w=log10(4e-5)+i*(log10(5e-2)-log10(4e-5))/N
        eps=[eps;10^w]
        f=[f;1.01*(2*log10(3.7/eps($)))^-2]
        z=epsfD2Re(f($),eps($))
        Re=[Re;z($)]
    end
    plot2d("ll",Re,f,2)
end

function x2=root(f,x1,x2,tol)
    while abs(f(x2))>tol
        x=(x1+x2)/2
        if f(x)*f(x1)>0 x1=x
        else x2=x end
    end
end
