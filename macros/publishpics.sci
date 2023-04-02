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

function publishpics(K,eps)
    function y=foo(fD)
        y=1/fD^.5+2*log10(eps/3.7+2.51/(K/fD)^.5/fD^.5)
    endfunction
    fD=root(foo,1d-2,1e-1,1d-4)
    Re=(K/fD)^.5
    if Re>2.3e3
        islam=%f
    else
        Re=K/64
        fD=64/Re
        islam=%t
    end
    if winsid()==[]
        scf(0)
    else
        scf(max(winsid())+1)
    end
    x=[5e-2 2.5e-2 1e-2 3e-3 1e-3 3e-4 1e-4]
    for i=1:length(x)
        turbulent(x(i),"k")
    end
    //loglog(Re,fD,"rd")
    rough("-.b")
    smoothline("-.b")
    laminar_2("k")
    //loglog([(K/1e-2)^.5 (K/1e-1)^.5],[1d-2 1d-1],"--r")
    xgrid(33,1,7)
    xlabel("$Re$","fontsize",4)
    ylabel("$f$","fontsize",4)
    gca().data_bounds=[1d2 1d8 1d-2 1d-1]
    gca().grid=[1,1]
    gca().grid_style=[9,9]
    gcf().figure_size=[600,600]
endfunction

function laminar_2(t)
    Re=[5e2 2.3e3]
    f=64 ./ Re
    loglog(Re,f,t)
endfunction
