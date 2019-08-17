#!/bin/bash

clear

mkdir -p ./optimize/bin ./optimize/obj ./optimize/mod ./output/save ./output/vtk ./output/csv ./output/ridge

gfortran -O3 -c ./sourcecode/type_ContactData.f08
gfortran -O3 -c ./sourcecode/type_AABB.f08
gfortran -O3 -c ./sourcecode/type_Quaternion.f08
gfortran -O3 -c ./sourcecode/type_Ridge.f08
gfortran -O3 -c ./sourcecode/type_Iceblock.f08
gfortran -O3 -c ./sourcecode/type_Element.f08
gfortran -O3 -c ./sourcecode/type_Wall.f08

gfortran -O3 -c ./sourcecode/mod_Constants.f08
gfortran -O3 -c ./sourcecode/mod_Functions.f08
gfortran -O3 -c ./sourcecode/mod_Global.f08
gfortran -O3 -c ./sourcecode/mod_Output.f08
gfortran -O3 -c ./sourcecode/mod_Neighbour.f08
gfortran -O3 -c ./sourcecode/mod_OverlapComputation.f08
gfortran -O3 -c ./sourcecode/mod_Force.f08
gfortran -O3 -c ./sourcecode/mod_Particle.f08
gfortran -O3 -c ./sourcecode/mod_ImportExport.f08
gfortran -O3 -c ./sourcecode/mod_PredCorr.f08
gfortran -O3 -c ./sourcecode/Main.f08

gfortran *.o -O3 -o NumRidge

cp NumRidge ./optimize/bin
mv *.mod ./optimize/mod
mv *.o ./optimize/obj
