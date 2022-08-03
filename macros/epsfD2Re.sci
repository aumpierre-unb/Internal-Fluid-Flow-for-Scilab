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
    // [Re]=epsfD2Re(fD[,eps[,s]])
    //
    // Parameters
    // fD: Darcy friction factor
    // eps: optional, relative roughness (default is eps=2e-3)
    // s: optional, boolean for display plot (default is s=%f)
    // Re: Reynolds number
    //
    // Description
    // epsfD2Re computes the Reynolds number, given
    // the Darcy friction factor and
    // the relative roughness for for laminar regime and,
    // when possible, also for turbulent regime.
    // By default, relative roughness is eps=2e-3. 
    // If eps>5e-2, execution is aborted. 
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
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.055 and
    // ..// the relative roughness eps=0.00044.
    // ..// In this case, due to higher friction,
    // ..// only the laminar solution is acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.055 and eps=0.00044:
    // Re=epsfD2Re(0.055,0.00044,%f)
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.00044 and
    // ..// the relative roughness eps=0.0022.
    // ..// In this case, due to lower friction,
    // ..// only the turbulent solution is acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.022 and eps=0.00044:
    // Re=epsfD2Re(0.022,0.00044,%f)
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.033 and
    // ..// the relative roughness eps=0.044.
    // ..// In this case, due to higher roughness,
    // ..// only the laminar solution is acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.033 and eps=0.044:
    // Re=epsfD2Re(0.033,0.044,%f)
    // ..// e.g. Compute the Reynolds number Re given
    // ..// the Darcy friction factor fD=0.033 and
    // ..// the relative roughness eps=0.0044.
    // ..// In this case,
    // ..// both laminar and turbulent solutions are acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.033 and eps=0.0044:
    // Re=epsfD2Re(0.033,0.0044,%f)
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
    // ..// the Darcy friction factor fD=0.033.
    // ..// When not assigned, relative roughness is
    // ..// the default eps=2e-3.
    // ..// In this case, 
    // ..// both laminar and turbulent solutions are acceptable:
    // ..//
    // ..// This call computes Re
    // ..// for fD=0.033 and the default eps=2e-3:
    // Re=epsfD2Re(0.033)
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
        eps=2e-3
    else
        eps=varargin(1)
        if eps>5e-2 abort end
    end
    Re=[]
    f=[]
    if 64/fD<2.5e3
        Re=[Re;64/fD]
        f=[f;fD]
    end
    if fD>(2*log10(3.7/eps))^-2
        function y=foo(Re)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re/sqrt(fD))
        endfunction
        re=root(foo,1e3,1e8,1e-4)
        if re>2.5e3
            Re=[Re;re]
            f=[f;fD]
        end
    end
    if ~isempty(f) && argn(2)==3 && varargin(2)
        if winsid()==[] scf(0)
        elseif scf(max(winsid())+1) end
        laminar()
        turb(eps)
        turb(eps*3)
        turb(eps*10)
        turb(eps/3)
        turb(eps/10)
        rough()
        smooth()
        xgrid(33,1,7)
        loglog(Re,f,"rd")
//        xlabel("$Re={ {\rho v D} \over\displaystyle \mu }$",..
//               "fontsize",4)
        xlabel("$Re$","fontsize",4)
//        ylabel("$f={h \over\displaystyle {v^2 \over\displaystyle 2g}{L \over\displaystyle D}}$",..
//               "fontsize",4)
        ylabel("$f$","fontsize",4)
        gca().data_bounds=[1d2 1d8 1d-2 1d-1]
        gca().grid=[1,1]
        gca().grid_style=[9,9]
        gcf().figure_size=[600,600]
    end
endfunction

function laminar()
    Re=[5e2 4e3]
    f=64 ./ Re
    loglog(Re,f,"k")
endfunction

function turb(eps)
    N=50
    for i=1:N
        w=log10(2d3)+i*(log10(1d8)-log10(2d3))/N
        Re(i)=10^w
        function y=foo(fD)//Colebrooke-White equation
            y=1/sqrt(fD)+2*log10(eps/3.7+..
             +2.51/Re(i)/sqrt(fD))
        endfunction
        function y=foobar(fD)//Tkachenko-Mileikovskyi equation
            y=fD+..
             -1/(0.8284*log10(eps/4.913+10.31/Re(i)))^2
        endfunction
        f(i)=root(foo,6e-3,1e-1,1e-4)
    end
    loglog(Re,f,"k")
endfunction

function smooth()
    N=50
    for i=1:N
        w=log10(2d3)+i*(log10(1d7)-log10(2d3))/N
        Re(i)=10^w
        function y=foo(fD)
            y=1/sqrt(fD)+2*log10(2.51/Re(i)/sqrt(fD))
        endfunction
        f(i)=root(foo,6e-3,1e-1,1e-4)
    end
    loglog(Re,f,"--b")
end

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
    loglog(Re,f,"--b")
end

function x2=root(f,x1,x2,tol)
    while abs(f(x2))>tol
        x=(x1+x2)/2
        if f(x)*f(x1)>0 x1=x
        else x2=x end
    end
end
