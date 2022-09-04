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

function [Re]=epsfD2Re(fD,varargin)
    // Computes the Reynolds number given the Darcy friction factor and the relative roughness
    //
    // Syntax
    // [Re]=epsfD2Re(fD[,eps[,fig]])
    //
    // Parameters
    // fD: Darcy friction factor
    // eps: optional, relative roughness (default is eps=2e-3)
    // fig: optional, boolean for display plot (default is fig=%f)
    // Re: Reynolds number
    //
    // Description
    // epsfD2Re computes the Reynolds number, given
    // the Darcy friction factor and
    // the relative roughness for for laminar regime and,
    // when possible, also for turbulent regime.
    // By default, tube is assumed to be smooth, eps=0. 
    // If eps>5e-2, eps is reset to 5e-2. 
    // Computation is based on the Colebrooke-White equation 
    // for turbulent flow and the Poiseuille condition 
    // for laminar flow. 
    // Inputs are to be given in a consistent system of units.
    //
    // Examples
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.033 and
    // ..// the relative roughness eps=0.00044.
    // ..// In this case,
    // ..// both laminar and turbulent solutions are acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.033 and eps=0.0044:
    // Re=epsfD2Re(0.033,0.00044,%f)
    // ..// Alternatively:
    // fD=0.033;..
    // eps=0.0044;..
    // Re=epsfD2Re(fD,eps)
    // ..// This call computes Re
    // ..// for fD=0.033 and eps=0.0044
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // Re=epsfD2Re(0.033,0.0044,%t)
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.055 and
    // ..// the relative roughness eps=0.00044.
    // ..// In this case, due to higher friction,
    // ..// only the laminar solution is acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.055 and eps=0.00044
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // Re=epsfD2Re(0.055,0.00044,%t)
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.00044 and
    // ..// the relative roughness eps=0.0022.
    // ..// In this case, due to lower friction,
    // ..// only the turbulent solution is acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.022 and eps=0.00044
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // Re=epsfD2Re(0.022,0.00044,%t)
    // ..// This call computes Re
    // ..// for fD=0.022 and
    // ..// the default smooth condition, eps=0:
    // Re=epsfD2Re(0.022)
    //
    // See also
    //  epsRe2fD
    //  hDeps2fDRe
    //  hveps2fDRe
    //  hvthk2fDRe
    //  hQeps2fDRe
    //  hQthk2fDRe
    //
    // Authors
    //  Alexandre Umpierre

    if argn(2)==1
        eps=0
        warning("Relative roughness assined to eps=0")
    else
        eps=varargin(1)
        if eps>5e-2
            eps=5e-2
            warning("Relative roughness reassined to eps=5e-2.")
        end
    end
    Re=[]
    f=[]
    if 64/fD<2.3e3
        Re=[Re;64/fD]
        f=[f;fD]
    end
    if fD>(2*log10(3.7/eps))^-2
        function y=foo(Re)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re/sqrt(fD))
        endfunction
        re=root(foo,1e3,1e8,1e-4)
        if re>2.3e3
            Re=[Re;re]
            f=[f;fD]
        end
    end
    if isempty(f)
        printf("No solution found.\n")
        abort
    end
    if ~isempty(f) && argn(2)==3 && varargin(2)
        if winsid()==[] scf(0)
        elseif scf(max(winsid())+1) end
        if min(Re)<2.3e3 laminar("r")
        else laminar("k") end
        if max(Re)>2.3e3 turb(eps,"r")
        else turb(eps,"k") end
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
        if ~eps==0 smooth("b") end
        loglog(Re,f,"rd")
        loglog([1d2 1d8],[fD fD],"--r")
        xgrid(33,1,7)
        xlabel("$Re$","fontsize",4)
        ylabel("$f$","fontsize",4)
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
