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

function [fD]=epsRe2fD(Re,varargin)
    // Computes the Darcy friction factor given the Reynolds number and the relative roughness
    //
    // Syntax
    // [fD]=epsRe2fD(Re[,eps[,s]])
    //
    // Parameters
    // Re: Reynolds number
    // eps: optional, relative roughness (default is eps=2e-3)
    // s: optional, boolean for display plot (default is s=%f)
    // fD: Darcy friction factor
    //
    // Description
    // epsRe2fD computes the Reynolds number, given
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
    // ..// e.g. Compute the Darcy friction factor fD given
    // ..// the Reynolds number Re=2.5e4 and
    // ..// the relative roughness eps=0.0044:
    // ..//
    // ..// This call computes fD
    // ..// for Re=2.5e4 and eps=0.0044:
    // fD=epsRe2fD(2.5e4,0.0044,%f)
    // ..// Alternatively:
    // Re=2.5e4;..
    // eps=0.0044;..
    // fD=epsRe2fD(Re,eps)
    // ..// This call computes fD
    // ..// for Re=2.5e4 and eps=0.0044
    // ..// and plots a representation of the solution
    // ..// on a schematic Moody diagram:
    // fD=epsRe2fD(2.5e4,0.0044,%t)
    // ..// e.g. Compute the Darcy friction factor fD given
    // ..// the Reynolds number Re=2.5e4:
    // ..//
    // ..// This call computes fD
    // ..// for the default eps=2e-3:
    // fD=epsRe2fD(2.5e4)
    //
    // See also
    //  epsfD2Re
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
    if Re<3e3
        fD=64/Re
    else
        function y=foo(fD)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re/sqrt(fD))
        endfunction
        fD=root(foo,1d-2,1e-1,1d-4)
    end
    if (argn(2)==3 && varargin(2))
        if winsid()==[] scf(0)
        elseif scf(max(winsid())+1) end
        laminar()
        turb(eps)
        turb(eps*5)
        turb(eps*10)
        turb(eps/5)
        turb(eps/10)
        xgrid(33,1,7)
        plot(Re,fD,"rd")
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
