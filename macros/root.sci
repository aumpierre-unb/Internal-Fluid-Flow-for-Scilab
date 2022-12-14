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

function x2=root(f,x1,x2,tol)
    while abs(f(x2))>tol
        x=(x1+x2)/2
        if f(x)*f(x1)>0 x1=x
        else x2=x end
    end
end
