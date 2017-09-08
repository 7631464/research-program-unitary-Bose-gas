FC = ifort

CMPLFLAGS =  -O3 #-xW -tpp7

LIB_CARTER = -llapack -lblas #-ffree-line-length-none

LIBRARIES = $(LIB_CARTER)

SOURCE = solgpnew7.f90 bogoliubovsub3.f90 thomasfermi.f90 hypersub.f90

OBJECT=$(SOURCE) 

Rmatrix_calc: $(OBJECT)
	$(FC) $(CMPLFLAGS) -o beccalc $(OBJECT) $(LIBRARIES)

.f90.o:              
	$(FC) $(CMPLFLAGS) -c $< 
.f.o: 
	$(FC) $(CMPLFLAGS)  -c $<

clean:
	rm $(DEL)
