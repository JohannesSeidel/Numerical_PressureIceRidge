@echo off

md .\debug\bin
md .\debug\obj
md .\debug\mod 
md .\output\save
md .\output\vtk
md .\output\csv 
md .\output\ridge

gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_ContactData.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_AABB.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_Quaternion.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_Ridge.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_Iceblock.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_Element.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\type_Wall.f08

gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Constants.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Functions.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Global.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Output.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Neighbour.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_OverlapComputation.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Particle.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_Force.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_ImportExport.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\mod_PredCorr.f08
gfortran -fdump-core -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c .\sourcecde\Main.f08

gfortran *.o -Og -fdump-core -fcheck=all -fbounds-check -fbacktrace -o NumRidge

copy NumRidge.exe ./debug/bin
move *.mod ./debug/mod
move *.o ./debug/obj


