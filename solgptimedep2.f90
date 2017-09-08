program main
implicit double precision (a-h,o-z)
!use split operator method, H=TA+TB+V
!the wavefunction is calculated in y(r)=r*psi(r)
!this program is to test the energy conservation of GP time propagation without renormalization, we ramep the external potential from 1/2 r^2 to
!1/2 (r/2)^2.
double precision,allocatable::expV(:)   !imaginary time evolution of V
double precision,allocatable::rint(:)    !grids of r
double precision,allocatable::ytrial(:)  !initial wavefunction to do the imaginary time diffusion,allocatable
double precision,allocatable::y(:)   !diffusion,allocatable wavefunction solution, unnormalized
double precision,allocatable::y1(:)  !diffuison wave function solution, normalized
double precision,allocatable::Tdiag(:)   !diagnal part of kinetic matrix
double precision,allocatable::Tnondiag(:)   !nondiagal part of kinetic matrix
double precision,allocatable::Vho(:)    !harmonic potential array
double precision,allocatable::Vint1(:)   !interaction potential
double precision,allocatable::Vint2(:)   !interaction potential
double precision,allocatable::kf(:)
double precision,allocatable::rho(:)
double precision,allocatable::exptad(:),exptasd(:)
double precision,allocatable::exptbd(:),exptbsd(:)
double precision::hmtave       !average value of H, the orbital energy
double precision::hmt2ave      !average value of H^2
double precision::dnorm        !wavefunction normalization constant
double precision::cnorm(50000)      !wavefunction normalization at each time step
double precision,allocatable::hmtdiag(:)   !diagonal part of hamiltonian matrix
double precision::e0     !total energy in each time step
double precision::orbit  !orbital energy in each time step
double precision::error     !error in each imaginary time diffusion
double precision::errortrack(500)
complex(8)::dtc     !propagation time step dtc=i*t
complex(8),allocatable::yc(:)    !wave function in realtime propagaton, unnomalized
complex(8),allocatable::y1c(:)  !wave function in realtime propagaton, nomalized
complex(8),allocatable::expTadc(:),exptasdc(:)
complex(8),allocatable::expTbdc(:),exptbsdc(:)
complex(8),allocatable::expVc(:)
complex(8)::ddnorm
complex(8)::cj      !imaginary unit
complex(8)::yk(2000)    !wavefunction in k space
complex(8)::ykrho(2000)
double precision::num   !total number of particles
double precision::kint(2000)   !momentum grids
double precision::krho(2000)    !transverse momentum grids
double precision::F(2000)   !  density |yk|^2
double precision::fkrho(2000)   !probability distribution of transverse momentum
integer::ipiv(20)
double precision::cs(2000),dz(2000),we(2000),cssav(2000),wesav(2000)!cs,cssav are points of Gauss-Legendre quadrature,we,wesav are the weights
!double precision::avesize(1000) !average size of the wavefunction
!double precision::deltae0 !the change of total energy after each time propagation step
double precision:: raverecord(20)
double precision:: nmin, nmax, nnstep

interface
function zeta(x,matdim)    !renormalization function
double precision::x(matdim),zeta(matdim)
end function
function dzeta(x,matdim)       !derivatives of renormalization function
double precision::x(matdim),dzeta(matdim)
end function
function trimatmult(a,b,c,matdim)
double precision::a(matdim),b(matdim-1),c(matdim),trimatmult(matdim)
end function
function trimatmultc(a,b,c,matdim)
complex(8)::a(matdim),b(matdim-1),c(matdim),trimatmultc(matdim)
end function
function dotprod(a,b,matdim)
double precision::a(matdim),b(matdim),dotprod
end function
function dotprodc(a,b,matdim)
complex(8)::a(matdim),b(matdim),dotprodc
end function
end interface

open(unit=10,file='timedep.dat')
open(unit=20,file='wavefunction.txt')
open(unit=30,file='wfreal.txt')
open(unit=40,file='wfimag.txt')
!open(unit=50,file='frho.txt')
!open(unit=60,file='wkimag2000.txt')
!open(unit=70,file='avesize2000.txt')
open(13,file='runtime.dat')
open(18,file='breathfreq.dat')

!write(10,*)'thi program solve the time depended GPE, intial a=0.01 then a is ramped to&
!&a=100, total number of particle is num=2000. finite difference spacing dr=0.05. I calculate the momentum distribution and &
!&integrate it over kz to get the transverse momentum distribution, as measured experimentally.&
!& the toal propagation time is 10, with time step dt=0.001.'

!All of the values are in oscillator units
!number of particles N=50000
!initial scattering length a=0.00218
!final scattering length a=21.8
!experimental Thomas-Fermi density <n>=200osu=5.5*10^12 cm^-3
!experimental Thomas-Fermi radius r0=4.39osu=15 um
!experimental Fermi momentum kF=11osu=6.9um^-1
!the ground state enery at a=0.01 is 7750.1388734280845
!the orbital energy at a=0.01 is  5.2224803759827783
!the ground state enery at a=100 is 30553.924918197881
!the orbital energy at a=100 is 18.075979790899172
!the ground state enery at a=0.00218 is 357073.58289817191
!the orbital energy at a=0.00218 is 9.8718940886520592
!the ground state enery at a=21.8 is 2218811.5831129011
!the orbital energy at a=21.8 is 52.610414342602304
!the total propagating time in experiment is 0.04398=700 us
!the propagating time in the momentum distribution measurement is 0.0106809=170 us

cj=(0d0,1d0)
pi=3.14159265358979d0



read(*,*)
read(*,*) ainitial, afinal, numpoints,  matdim, dr
read(*,*)

read(*,*)
read(*,*)nmin, nmax, numNpoints, nnpower
read(*,*)

read(*,*)
read(*,*) napower

allocate(expV(matdim),rint(matdim),y(matdim),y1(matdim), ytrial(matdim),Tdiag(matdim) ,Tnondiag(matdim-1) , &
vho(matdim), vint1(matdim),vint2(matdim),  &
kf(matdim),rho(matdim), exptad(matdim),exptasd(matdim-1),exptbd(matdim),exptbsd(matdim-1),hmtdiag(matdim) , &
yc(matdim), y1c(matdim), expTadc(matdim),exptasdc(matdim-1), expTbdc(matdim),exptbsdc(matdim-1), &
expVc(matdim) )


!ainitial=1d-6
!afinal=1d0
!napower=3
astep=(afinal**(1d0/dfloat(napower))-ainitial**(1d0/dfloat(napower)))/dfloat(numpoints-1)
nnstep=(nmax**(1d0/dfloat(nnpower))-nmin**(1d0/dfloat(nnpower)))/dfloat(numNpoints-1)


call cpu_time(time1)
do inn=1,1              !scan nn
!if(inn==1) then
!num=1d4
!else
!num=5d4
!endif
num=(nmin**(1d0/dfloat(nnpower))+nnstep*dfloat(inn-1))**nnpower

do ia=1,numpoints             !scan a    (quench from aprog/2 to aprog)
aprog=(ainitial**(1d0/dfloat(napower))+astep*dfloat(ia-1))**napower

egtrack=0d0
Eg=0d0
istep=1
!a0=0.0418d0
a0=aprog*2d0/3d0
rhoave=0d0
kfave=0d0
!num=1d3
dt=0.03d0
eigpre=0d0
iexit=0
!dr=0.03d0
istep=1
e0=0d0
orbit=0d0
error=0d0
avesize=0d0
call ZELEGA(2000,CSsav,DZ)
call WELEGA(2000,CSsav,DZ,WEsav)

!call cpu_time(time1)

dlambda=1d0+(2d0/pi)**(1d0/10d0)*(num*a0)**(1d0/5d0)

!phi\[Lambda] = 1/\[Lambda]^(3/2)/\[Pi]^(3/4) *Exp[-(r/\[Lambda])^2/2]
!rint(1000) is the grid point array, ytrial(1000) is the trial wave function evaluated at each grid point.
do i=1,matdim
rint(i)=dr*dfloat(i)
ytrial(i)=1d0/(dlambda)**(3d0/2d0)/pi**(3d0/4d0)*exp(-(rint(i)/dlambda)**2/2d0)*rint(i) 
!ytrial(i)=1d0/(4.6d0)**(3d0/2d0)/pi**(3d0/4d0)*exp(-(rint(i)/4.6d0)**2/2d0)*rint(i)
end do
do i =1,2000
kint(i)=2d0*pi/400d0*dfloat(i-1)
end do 

do i=1,2000
krho(i)=2d0*pi/400d0*dfloat(i-1)
end do 
!write(*,*) 'matdim=',matdim
!write(*,*) 'dr=',dr
!write(*,*) 'rint=',rint
!write(*,*) zeta(rint,matdim)
!stop
!------------------------------------------------------------------------------------------------
!construct expTa and expTb
exptad=cosh(dt/(4d0*dr**2))
exptasd=0d0
exptasd(1:madim-1:2)=sinh(dt/(4d0*dr**2))


exptbd(1)=1d0
exptbd(matdim)=1d0
exptbd(2:matdim-1)=cosh(dt/(2d0*dr**2))
exptbsd=0d0
exptbsd(2:matdim-2:2)=sinh(dt/(2d0*dr**2))

!------------------------------------------------------------------------------------------------
!construct Tdiag,Tdiagmat,Tnondiag,Vho
Tdiag=1d0/dr**2
Tnondiag=-0.5d0/dr**2
Vho=rint**2/2d0

!------------------------------------------------------------------------------------------------
!Start imaginary time diffusion to converge the BEC ground state
y1=ytrial
step=1
error=0d0
errortrack=0d0
do while(.true.)

rho=num*y1**2/rint**2
kf=(6*pi**2*num*(y1/rint)**2)**(1d0/3d0)

Vint1=num*4*pi/3d0*(dzeta(kf*a0,matdim)*a0*((y1/rint)**2)-zeta(kf*a0,matdim)*((y1/rint)**4)**(1d0/3d0)/(6*pi**2*num)**(1d0/3d0))
Vint2=num*4*pi*zeta(kf*a0,matdim)/(6*pi**2*num)**(1d0/3d0)*((y1/rint)**4)**(1d0/3d0)


!here is right after the construction of everything for hamiltonian


!write(*,*)e0(istep)
hmtdiag=Tdiag+Vho+Vint2

hmt2ave=4*pi*dotprod(y1,trimatmult(hmtdiag,tnondiag,trimatmult(hmtdiag,tnondiag,y1,matdim),matdim),matdim)*dr
hmtave=4*pi*dotprod(y1,trimatmult(hmtdiag,tnondiag,y1,matdim),matdim)*dr
orbit=hmtave
error=sqrt(hmt2ave-hmtave**2)
errortrack(1:499)=errortrack(2:500)
errortrack(500)=error
if (error<1d-3 .or. istep>30000 .or. dt<1d-10) exit
if (istep>500) then
if (abs(error-errortrack(200))/error<1d-3 ) then

dt=dt/2d0

exptad=cosh(dt/(4d0*dr**2))
exptasd=0d0
exptasd(1:matdim-1:2)=sinh(dt/(4d0*dr**2))


exptbd(1)=1d0
exptbd(matdim)=1d0
exptbd(2:matdim-1)=cosh(dt/(2d0*dr**2))
exptbsd=0d0
exptbsd(2:matdim-2:2)=sinh(dt/(2d0*dr**2))

end if
end if
expv=exp(-(vint1+vint2+Tdiag+Vho)*dt/2d0)
!here is right before time diffusion

y=expv*trimatmult(exptad,exptasd,trimatmult(exptbd,exptbsd,trimatmult(exptad,exptasd,expv*y1,matdim),matdim),matdim)

dnorm=0d0
do i=1,matdim
dnorm=dnorm+y(i)**2
end do
dnorm=dnorm*4*pi*dr
y1=y/sqrt(dnorm)

!here is right after time diffusion

istep=istep+1

end do
e0=4*pi*dr*dotprod(y1,trimatmult((Tdiag+Vho)*num+Vint2*num/2d0,Tnondiag*num,y1,matdim),matdim)
write(*,*)'the ground state energy is'
write(*,*) e0
write(*,*)'the orbit energy is'
write(*,*) orbit
write(*,*)'the error is'
write(*,*) error
!write(10,*)'the wave function is'
!write(20,*) y1(1:1000)
write(*,*) 'the total step is'
write(*,*) istep
write(*,*) 'the final imaginary time step is'
write(*,*) dt
write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!------------------------------------------------------------------------------------------------------------------------
!Start real time propagating
e0=0d0
orbit=0d0
error=0d0
!a0=0.1d0
a0=aprog
istep=1
y1c=dcmplx(y1)
dtc=(0d0,0.00005d0)
dt=0.00005d0
expvc=(0d0,0d0)
!ddnorm=4d0*pi*matmul(transpose(conjg(y1cs)),y1cs)*dr
!write(*,*)ddnorm(1,1)
!y1c=y1cs(:,1)

!construct expTac and expTbc
exptadc=cos(dt/(4d0*dr**2))
exptasdc=(0d0,0d0)
exptasdc(1:matdim-1:2)=cj*sin(dt/(4d0*dr**2))

exptbdc(1)=(1d0,0d0)
exptbdc(matdim)=(1d0,0d0)
exptbdc(2:matdim-1)=cos(dt/(2d0*dr**2))
exptbsdc=0d0
exptbsdc(2:matdim-2:2)=cj*sin(dt/(2d0*dr**2))

iir=1
ravepre=0d0
ravepre2=0d0
!start real time propagation
do while(istep<=200001)

rho=num*abs(y1c)**2/rint**2
kf=(6*pi**2*num*abs(y1c/rint)**2)**(1d0/3d0)
!vint1=0d0
Vint1=num*4*pi/3d0*(dzeta(kf*a0,matdim)*a0*(abs(y1c/rint)**2)-zeta(kf*a0,matdim)*(abs(y1c/rint)**4)**(1d0/3d0)/(6*pi**2*num)**(1d0/3d0))
Vint2=num*4*pi*zeta(kf*a0,matdim)/(6*pi**2*num)**(1d0/3d0)*(abs(y1c/rint)**4)**(1d0/3d0)



expvc=exp(-(vint1+vint2+Tdiag+Vho)*dtc/2d0)

!here is right after the construnction of hamiltonian


hmtdiag=Tdiag+Vho+Vint1+Vint2

!----------------calculate average orbital energy and average total energy and average r-----------
orbit=real(4*pi*dotprodc(y1c,trimatmultc(dcmplx(hmtdiag),dcmplx(tnondiag),y1c,matdim),matdim)*dr)

e0=abs(4*pi*dr*dotprodc(y1c,trimatmultc(dcmplx((Tdiag+Vho)*num+Vint2*num/2d0),dcmplx(Tnondiag*num),y1c,matdim),matdim))
!test=aimag(4*pi*dr*dotprodc(y1c,trimatmultc(dcmplx((Tdiag+Vho)*num+Vint2*num/2d0),dcmplx(Tnondiag*num),y1c)))
rave=0d0
do i=1,matdim
rave=rave+4d0*pi*dr*abs(y1c(i))**2*rint(i)
end do 

if (rave<ravepre .and. ravepre>ravepre2 .and. istep>2) then
raverecord(iir)=dt*dfloat(istep-2)
iir=iir+1
end if

if (rave>ravepre .and. ravepre<ravepre2 .and. istep>2) then
raverecord(iir)=dt*dfloat(istep-2)
iir=iir+1
end if

ravepre2=ravepre
ravepre=rave


go to 100

yk=0d0
f=0d0
do m=1,1000
yk(1)=yk(1)+4*pi*y1c(m)*rint(m)*dr
end do
do i=2,2000
do m=1,1000
yk(i)=yk(i)+4*pi*y1c(m)*sin(kint(i)*rint(m))/kint(i)*dr
end do
end do
f=abs(yk)**2

fkrho=0d0
do i=1,2000
dkzmax=sqrt(100d0*pi**2-krho(i)**2)
cs=dkzmax/2d0*(cssav+1d0)
we=dkzmax/2d0*wesav
iki=2
do m=1,2000
dummyk=sqrt(krho(i)**2+cs(m)**2)
do ik=iki,2000
if (dummyk<kint(ik)) then
fkrho(i)=fkrho(i)+((f(ik)-f(ik-1))/(2d0*pi/400d0)*(dummyk-kint(ik-1))+f(ik-1))*we(m)
iki=ik
exit
end if
end do

end do
end do

fkrho=2d0*fkrho

100 continue

if (mod(istep,100)==1) then
!write(10,'(10e22.12)') dt*dfloat(istep),e0, orbit,rave
!write(20,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!write(20,*)abs(y1c(1:1000))
!write(30,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!write(30,*)real(y1c(1:1000))
!write(40,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!write(40,*)aimag(y1c(1:1000))
!write(50,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!write(50,*)fkrho(1:2000)
!write(60,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!write(60,*)aimag(fkrho(1:2000))
end if

!here is right before the time propagation

y1c=expvc*trimatmultc(exptadc,exptasdc,trimatmultc(exptbdc,exptbsdc,trimatmultc(exptadc,exptasdc,expvc*y1c,matdim),matdim),matdim)

!ddnorm=4*pi*dotprodc(y1c,y1c)*dr
!if (aimag(ddnorm)>1d-8) write(*,*) 'there is an error in ddnorm'
!cnorm(istep)=real(ddnorm)
!if (abs(1d0-cnorm(istep))>1d-10) write(*,*)'there is error in normalization'
!here is right after the time propagation

istep=istep+1

end do


write(18,'(10e20.10)') num, aprog, raverecord(1:4)

end do   ! end of scanning a

end do    ! end of scanning nn



call cpu_time(time2)
write(13,*)'runtime=',time2-time1
!write(10,*)'the error is'
!write(10,*) errhmt(1,1)
!write(10,*)'the normalization at each time step is'
!write(10,*) (cnorm(i),i=1,istep-1,10)
!write(10,*)'the average size of the wave functino at istep/50 is'
!write(10,*) (avesize(i),i=1,500)
end program

!****************************************************************************************
function zeta(x,matdim)
double precision::x(matdim),zeta(matdim)
!zeta=0.395255838d0-1.138d0*atan(0.362d0-0.994d0*x)
zeta=0.395255838d0-5.272769443287473d-10-1.138d0*atan(0.362d0-0.994d0*x/1.0001131697794252d0)
end function

function dzeta(x,matdim)
double precision x(matdim),dzeta(matdim)
!dzeta=1.131172d0/(1d0+(0.362d0-0.994d0*x)**2)
dzeta=1.131044d0/(1d0+(0.362d0-0.9938875219683654d0*x)**2)
end function

function zetan(x)
double precision::x,zetan
!zetan=0.395255838d0-1.138d0*atan(0.362d0-0.994d0*x)
zetan=0.395255838d0-5.272769443287473d-10-1.138d0*atan(0.362d0-0.994d0*x/1.0001131697794252d0)
end function

function dzetan(x)
double precision x,dzetan
!dzetan=1.131172d0/(1d0+(0.362d0-0.994d0*x)**2)
dzetan=1.131044d0/(1d0+(0.362d0-0.9938875219683654d0*x)**2)
end function
!****************************************************************************************
double precision function fint(f,krho,r)
double precision::f(500),krho(500),r
do i=2,500
if (r<krho(i)) then
fint=(f(i)-f(i-1))/(2d0*pi/100d0)*(r-krho(i-1))+f(i-1)
exit
end if
end do
end function
!****************************************************************************************
function trimatmult(a,b,c,matdim)
integer:: matdim
! the matrix is unitary and symmetric
!a(n) is the diagonal elements
!b(n-1) is the 1 and -1 diagonal elements
double precision::a(matdim),b(matdim-1),c(matdim),ac(matdim),bc(matdim-1),bc2(matdim-1),trimatmult(matdim)
ac=a*c
bc=b*c(2:matdim)
bc2=b*c(1:matdim-1)
trimatmult(1)=ac(1)+bc(1)
trimatmult(matdim)=ac(matdim)+bc2(matdim-1)
do i=2,matdim-1
trimatmult(i)=ac(i)+bc(i)+bc2(i-1)
end do
end function
!****************************************************************************************
function trimatmultc(a,b,c,matdim)
! the matrix is unitary and symmetric
!a(n) is the diagonal elements
!b(n-1) is the 1 and -1 diagonal elements
integer:: matdim
complex(8)::a(matdim),b(matdim-1),c(matdim),ac(matdim),bc(matdim-1),bc2(matdim-1),trimatmultc(matdim)
ac=a*c
bc=b*c(2:matdim)
bc2=b*c(1:matdim-1)
trimatmultc(1)=ac(1)+bc(1)
trimatmultc(matdim)=ac(matdim)+bc2(matdim-1)
do i=2,matdim-1
trimatmultc(i)=ac(i)+bc(i)+bc2(i-1)
end do
end function
!****************************************************************************************
function diagmatmult(a,c,matdim)
integer:: matdim
!the matrix is diagonal
!a(n) is the diagonal elements
!c(n) is the multiplied vector
complex(8)::a(matdim),c(matdim),diagmatmult(matdim)
do i=1,matdim
diagmatmult(i)=a(i)*c(i)
end do
end function
!****************************************************************************************
function dotprod(a,b,matdim)
integer::matdim
double precision::a(matdim),b(matdim),dotprod
dotprod=0d0
do i=1,matdim
dotprod=dotprod+a(i)*b(i)
end do
end function
!****************************************************************************************
function dotprodc(a,b,matdim)
integer:: matdim
complex(8)::a(matdim),b(matdim),dotprodc
dotprodc=(0d0,0d0)
do i=1,matdim
dotprodc=dotprodc+conjg(a(i))*b(i)
end do
end function

!----------------------------------------------------------------------------------------
SUBROUTINE VALEPO(N,X,Y,DY,D2Y)
!**************************************************************
!*   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N
!*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
!*   N  = DEGREE OF THE POLYNOMIAL
!*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
!*   Y  = VALUE OF THE POLYNOMIAL IN X
!*   DY = VALUE OF THE FIRST DERIVATIVE IN X
!*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
!**************************************************************
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

Y   = 1.D0
DY  = 0.D0
D2Y = 0.D0
IF (N == 0) RETURN

Y   = X
DY  = 1.D0
D2Y = 0.D0
IF(N == 1) RETURN

YP   = 1.D0
DYP  = 0.D0
D2YP = 0.D0

DO I=2,N
C1 = DFLOAT(I)
C2 = 2.D0*C1-1.D0
C4 = C1-1.D0
YM = Y
Y  = (C2*X*Y-C4*YP)/C1
YP = YM
DYM  = DY
DY   = (C2*X*DY-C4*DYP+C2*YP)/C1
DYP  = DYM
D2YM = D2Y
D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1
D2YP = D2YM
END DO

RETURN
END
!----------------------------------------------------------------------------------------
SUBROUTINE ZELEGA(N,CS,DZ)
!***************************************************************
!*   COMPUTES THE ZEROES OF THE LEGENDRE POLYNOMIAL OF DEGREE N
!*   N  = THE NUMBER OF ZEROES
!*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
!*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
!***************************************************************
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION CS(N), DZ(N)
IF (N ==0) RETURN

CS(1) = 0.D0
DZ(1) = 1.D0

IF(N == 1) RETURN

N2 = N/2
IN = 2*N-4*N2-1
PH = 1.57079632679489661923D0
C  = PH/(2.D0*DFLOAT(N)+1.D0)

DO I=1,N2
DI  = DFLOAT(I)
CSX = DCOS(C*(4.D0*DI-1.D0))

DO IT=1,8
CALL VALEPO(N,CSX,Y,DY,D2Y)
CSX = CSX-Y/DY
END DO

CS(I) = -CSX
CS(N-I+1) = CSX
DZ(I) = DY*DFLOAT(IN)
DZ(N-I+1) = DY
END DO

IF(IN == -1) RETURN
CSX = 0.D0
CS(N2+1) = CSX
CALL VALEPO(N,CSX,Y,DY,D2Y)
DZ(N2+1) = DY

END
!-------------------------------------------------------------------------------------
SUBROUTINE WELEGA(N,CS,DZ,WE)
!*****************************************************************
!*   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
!*   N  = ORDER OF THE FORMULA
!*   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
!*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
!*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
!*****************************************************************
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION CS(N), DZ(N), WE(N)
IF (N == 0) RETURN

N2 = N/2

DO I=1,N2
X  = CS(I)
DY = DZ(I)
WEX = 2.D0/((1.D0-X*X)*DY*DY)
WE(I) = WEX
WE(N-I+1) = WEX
END DO

IF(N == 2*N2) RETURN
DY = DZ(N2+1)
WE(N2+1) = 2.D0/(DY*DY)

RETURN
END
