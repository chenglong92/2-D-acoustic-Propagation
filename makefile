objects=CAA_1_main.F90 DRP7_sub.f90 F_new.F90 initialize_sub.f90 \
	LDDRK_sub.f90 Partial_Deriv_sub.f90
two_dimensional_acoustics_transmition.exe: $(objects)
	gfortran $(objects) -o two_dimensional_acoustics_transmition.exe
