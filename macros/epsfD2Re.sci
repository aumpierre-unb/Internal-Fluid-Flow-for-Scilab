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
    // Re=epsfD2Re(fD[,eps[,s]])
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
    // // Compute the Reynolds number Re given 
    // // the relative roughness eps=0.0044 and 
    // // the Darcy friction factor fD=0.033:
    //
    // Re=epsfD2Re(0.033,0.0044,%f)
    // // This call computes Re only. Alternatively:
    // fD=0.033,..
    // eps=0.0044,..
    // Re=epsfD2Re(fD,eps)
    //
    // Re=epsfD2Re(0.033,0.0044,%t)
    // // This call computes Re 
    // // for eps=0.0044 and fD=0.033 and 
    // // plots a representation of the solution 
    // // on a schematic Moody diagram.
    //
    // // Compute the Reynolds number Re given 
    // // the Darcy friction factor fD=0.033 only:
    //
    // Re=epsfD2Re(0.033)
    // // This call computes Re 
    // // for the default eps=2e-3 and 
    // // plots a representation of the solution 
    // // on a schematic Moody diagram.
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
    if 64/fD<3e3
        Re=[Re;64/fD]
        f=[f;fD]
    end
    if fD>(2*log10(3.7/eps))^-2
        function y=foo(Re)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re/sqrt(fD))
        endfunction
        Re=[Re;root(foo,1e3,1e8,1e-4)]
        f=[f;fD]
    end
    if argn(2)==3 && varargin(2)
        if winsid()==[] scf(0)
        elseif scf(max(winsid())+1) end
        laminar()
        turb(eps)
        turb(eps*3)
        turb(eps*10)
        turb(eps/3)
        turb(eps/10)
        rough()
        xgrid(33,1,7)
        plot(Re,f,"rd")
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
        f=[f;1.02*(2*log10(3.7/eps($)))^-2]
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
