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

function smoothline(t)
    f=[]
    Re=[]
    N=30
    for n=1:N
        w=log10(2d3)+(n-1)*(log10(1d7)-log10(2d3))/(N-1)
        Re=[Re;10^w]
        function y=foo(fD)
            y=1/fD^.5+2*log10(2.51/10^w/fD^.5)
        endfunction
        f=[f;root(foo,6e-3,1e-1,1e-4)]
    end
    loglog(Re,f,t)
end
