This program calculate the many-body properties of a unitary Bose gas with renormalized interaction. The renormalization function is universal function that can be found in my paper. The major part of the program solves the renormalized GP equation (both time-independent and time-dependent) using the split operator method in a finite difference basis. It also incorporates several subroutines to calculate the many-body wave function using Thomas-Fermi approximation(thomasfermi.f90) and hyper spherical method(hypersub.f90), and to calculate the breathing mode excitation using Bogoliubov approximation (bogoliubovsub3.f90). 

Please read the Makefile to see how to compile the code. The entire code contains the following source files:


solgpnew[4,7,8].f90, solgptimedep2.f90: These are the main code of the program. When solving time-independent GP, compile solgpnew[4,7,8].f90. When solving time-dependent GP, compile solgptimedep2.f90.
solgpnew4: complete time independent calculation including Thomas-Fermi and hyperspherical method, typically used in time independent calculation by scanning both the scattering length and the number of particles. 
solgpnew7: in addition to solgpnew4, added bogoliubov approximaton in calculating elementary excitations.
solgpnew8: in addition to solgpnew4, added the calculation of two body contact.
solgptimedep2.f90: realtime propagaion of the TDGP, quenching the scattering legnth from aprog*2/3 to aprog, the matrix size is ajustable in parameter matdim. 
bogoliubovsub3.f90: This subroutine solves the coupled Bogoliubov equations to get the breathing mode excitation using linear mean-field solution basis, with the ground state wave function as the interaction core.
thomasfermi.f90: This subroutine calculates the mean-field wave function using Thomas-Fermi approximation with the inclusion of renormalized interaction. 
hypersub.f90: this subroutine solves the hyper-radial Schrodinger equation by finding the minimum of the hyper-radial potential curve. 





Input files: When solving time-independent GP, use parameter7.inp. When solving time-dependent GP, use paratimedep.inp.
parameter7.inp: This is the input file specifying the parameters when solving the time-indepnedent GP.
	a0initial,a0final,napower,numapoints: parameters that determines how to scan the scattering length in units of harmonic oscillator length;
	numparticle: not used;
	dr: finite difference grids spacing in units of harmonic oscillator length; the number of grid points is always 1000;
	nmin, nmax, nnpower, nnpoints: parameters that determines how to scan the number of particles;
	
paratimedep.inp: 
	a0initial, a0final, napower, numapoints: parameters that determines how to scan the scattering length in units of harmonic oscillator length;
	matdim: number of finite difference basis; 
	dr: finite difference grids spacing in units of harmonic oscillator length;
	nmin, nmax, nnpoints, nnpower:  parameters that determines how to scan the number of particles;



