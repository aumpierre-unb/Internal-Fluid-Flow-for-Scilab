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

path=get_absolute_file_path();
tbx_builder_macros(path);
tbx_builder_src(path);
tbx_builder_gateway(path);
tbx_build_localization(path);
tbx_builder_help(path);
tbx_build_loader(path);
tbx_build_cleaner(path);

