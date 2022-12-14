// Copyright (C) 2022 Alexandre Umpierre
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

function [fD]=epsRe2fD(Re,varargin)
    // epsRe2fD computes the Darcy friction factor given
    // the Reynolds number and
    // the pipe"s relative roughness
    //
    // Syntax
    // [fD]=epsRe2fD(Re[,eps[,fig]])
    //
    // Parameters
    // Re: Reynolds number
    // eps: optional, relative roughness (default is eps=2e-3)
    // fig: optional, boolean for display plot (default is fig=%f)
    // fD: Darcy friction factor
    //
    // Description
    // epsRe2fD computes the Darcy friction factor given
    // the Reynolds number and
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
    // // Compute the Darcy friction factor fD given
    // // the Reynolds number Re=2.5e4 and
    // // the relative roughness eps=4.4e-3:
    // //
    // // This call computes fD
    // // for Re=2.5e4 and eps=4.4e-3:
    // fD=epsRe2fD(2.5e4,4.4e-3,%f)
    // // Alternatively:
    // Re=2.5e4;..
    // eps=4.4e-3;..
    // fD=epsRe2fD(Re,eps)
    // // This call computes fD
    // // for Re=2.5e4 and eps=4.4e-3
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // fD=epsRe2fD(2.5e4,4.4e-3,%t)
    // // Compute the Darcy friction factor fD given
    // // the Reynolds number Re=2.5e4:
    // //
    // // This call computes fD
    // // for Re=2.5e4 and
    // // the default smooth condition, eps=0:
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
        eps=0
        warning("Relative roughness assined to eps=0")
    else
        eps=varargin(1)
        if eps>5e-2
            eps=5e-2
            warning("Relative roughness reassined to eps=5e-2.")
        end
    end
    if Re<2.3e3
        fD=64/Re
    else
        function y=foo(fD)
            y=1/sqrt(fD)+2*log10(eps/3.7..
             +2.51/Re/sqrt(fD))
        endfunction
        fD=root(foo,1d-2,1e-1,1d-4)
    end
    if argn(2)==3 && varargin(2)
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
        if Re<2.3e3
            laminar("r")
            turbulent(eps,"k")
        else
            laminar("k")
            turbulent(eps,"r")
        end
        loglog(Re,fD,"rd")
        loglog([Re Re],[1e-2 1e-1],"--r")
        xgrid(33,1,7)
        xlabel("$Re$","fontsize",4)
        ylabel("$f$","fontsize",4)
        gca().data_bounds=[1d2 1d8 1d-2 1d-1]
        gca().grid=[1,1]
        gca().grid_style=[9,9]
        gcf().figure_size=[600,600]
    end
endfunction
