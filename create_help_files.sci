// Copyright (C) 2022 2023 Alexandre Umpierre
// This file is part of psyhrometrics Toolbox.
// psyhrometrics Toolbox is free software:
// you can redistribute it and/or modify it under the terms
// of the GNU General Public License (GPL) version 3
// as published by the Free Software Foundation.
// internal-fluid-flow Toolbox is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the
// GNU General Public License along with this program.
// It is also available at www.gnu.org/licenses/.

// create_help_files is part of
// the internal-fluid-flow toolbox for Scilab.

path=get_absolute_file_path()
chdir(path+'macros')
mkdir(path+'help')
help_from_sci(path+'macros',path+'help')

chdir(path)
chdir(path+'help')
xmltojar(path+'help','Internal Fluid Flow','en_US')

add_help_chapter("Internal Fluid Flow",path)
