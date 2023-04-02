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

function rough(t)
    f=[]
    Re=[]
    N=20
    for n=1:N
        w=log10(4e-5)+(n-1)*(log10(5e-2)-log10(4e-5))/(N-1)
        eps=10^w
        f=[f;1.01*(2*log10(3.7/eps))^-2]
        z=epsfD2Re(f($),eps)
        Re=[Re;z($)]
    end
    loglog(Re,f,t)
end
