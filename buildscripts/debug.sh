#!/bash/bin

clear

mkdir -p ./debug/bin ./debug/obj ./debug/mod ./output/save ./output/vtk ./output/csv ./output/ridge

gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_ContactData.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_AABB.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_Quaternion.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_Ridge.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_Iceblock.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_Element.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/type_Wall.f08

gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Constants.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Functions.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Global.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Output.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Neighbour.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_OverlapComputation.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Particle.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_Force.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_ImportExport.f08
gfortran -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/mod_PredCorr.f08
gfortran -fdump-core -fcheck=all -fbounds-check -fbacktrace -Og -Wall -Wextra -pedantic -c ./sourcecode/Main.f08

gfortran *.o -Og -fdump-core -fcheck=all -fbounds-check -fbacktrace -o NumRidge

cp NumRidge ./debug/bin
mv *.mod ./debug/mod
mv *.o ./debug/obj
