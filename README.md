# Internal-Fluid-Flow-for-Scilab
A toolbox for internal fluid flow for Scilab</br>
This package provides a set of functions 
designed to solve problems of internal fluid flow 
on Scilab under GNU GPLv3. All functions are based 
on the Poiseuille condition for laminar flow, 
the Colebrooke-White equation for turbulent flow, 
and the Darcy-Weisbach equation for head loss.</br>
Basic functions allow the user to compute either 
the Reynolds number or the Darcy friction factor, 
given one of them is given along with 
the relative roughness.</br>
More advanced functions allow 
to compute both the Reynolds number and 
the Darcy friction factor given the head loss and 
other variables that describe the fluid flow.</br>
Enjoy!

TODO:

inflowlib v0.0.5 (2022-07-22)
=============================
* Spell of "Weisbach" is now correct
* espfD2Re now computes Reynolds number for laminar flow and/or turbulent flow
when possible
* The second input of epsfD2Re and epsRe2fD now is optional with default value
eps=2e-3
* Plots now include a complete turbulence line for reference
 
inflowlib v0.0.3 (2022-07-09)
=============================
* First release
* Module under construction       
