//  Copyright (C) 2022 2023 Alexandre Umpierre
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

function [Re]=epsfD2Re(fD,varargin)
    // epsfD2Re computes the Reynolds number given
    // the Darcy friction factor and
    // the pipe's relative roughness
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
    // epsfD2Re computes the Reynolds number given
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
    // // Compute the Reynolds number Re given
    // // the Darcy friction factor fD=3.3e-2 and
    // // the relative roughness eps=4.4e-4.
    // // In this case,
    // // both laminar and turbulent solutions are acceptable:
    // //
    // // This call computes Re
    // // for fD=3.3e-2 and eps=4.4e-3:
    // Re=epsfD2Re(3.3e-2,4.4e-4,%f)
    // // Alternatively:
    // fD=3.3e-2;..
    // eps=4.4e-3;..
    // Re=epsfD2Re(fD,eps)
    // // This call computes Re
    // // for fD=3.3e-2 and eps=4.4e-3
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // Re=epsfD2Re(3.3e-2,4.4e-3,%t)
    // // Compute the Reynolds number Re given
    // // the Darcy friction factor fD=5.5e-2 and
    // // the relative roughness eps=4.4e-4.
    // // In this case, due to higher friction,
    // // only the laminar solution is acceptable:
    // //
    // // This call computes Re
    // // for fD=5.5e-2 and eps=4.4e-4
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // Re=epsfD2Re(5.5e-2,4.4e-4,%t)
    // // Compute the Reynolds number Re given
    // // the Darcy friction factor fD=4.4e-4 and
    // // the relative roughness eps=2.2e-3.
    // // In this case, due to lower friction,
    // // only the turbulent solution is acceptable:
    // //
    // // This call computes Re
    // // for fD=2.2e-2 and eps=4.4e-4
    // // and plots a representation of the solution
    // // on a schematic Moody diagram:
    // Re=epsfD2Re(2.2e-2,4.4e-4,%t)
    // // This call computes Re
    // // for fD=2.2e-2 and
    // // the default smooth condition, eps=0:
    // Re=epsfD2Re(2.2e-2)
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
    Re_=64/fD
    if Re_<2.3e3
        Re=[Re;Re_]
        f=[f;fD]
    end
    if fD>(2*log10(3.7/eps))^-2
        Re_=2.51/(10^(1/fD^.5/-2)-eps/3.7)/fD^.5
        if Re_>2.3e3
            Re=[Re;Re_]
            f=[f;fD]
        end
    end
    if isempty(f)
        printf("No solution found.\n")
        abort
    end
    if ~isempty(f) && argn(2)==3 && varargin(2)
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
