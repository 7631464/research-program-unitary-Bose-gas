program main
implicit double precision (a-h,o-z)
!use split operator method, H=TA+TB+V
!the wavefunction is calculated in y(r)=r*psi(r)
double precision::expV(1000)   !imaginary time evolution of V
double precision::rint(1000)    !grids of r
double precision::ytrial(1000)  !initial wavefunction to do the imaginary time diffusion
double precision::y(1000)   !diffusion wavefunction solution, unnormalized
double precision::y1(1000)  !diffuison wave function solution, normalized
double precision::Tdiag(1000)   !diagnal part of kinetic matrix
double precision::Tnondiag(999)   !nondiagal part of kinetic matrix
double precision::Vho(1000)    !harmonic potential array
double precision::Vint1(1000)   !interaction potential
double precision::Vint2(1000)   !interaction potential
double precision::kf(1000)
double precision::rho(1000)
double precision::exptad(1000),exptasd(999)
double precision::exptbd(1000),exptbsd(999)
double precision::hmtave       !average value of H, the orbital energy
double precision::hmt2ave      !average value of H^2
double precision::dnorm        !wavefunction normalization constant
double precision::cnorm(50000)      !wavefunction normalization at each time step
double precision::hmtdiag(1000)   !diagonal part of hamiltonian matrix
double precision::e0     !total energy in each time step
double precision::orbit  !orbital energy in each time step
double precision::error     !error in each imaginary time diffusion
double precision::errortrack(500)
complex(8)::dtc     !propagation time step dtc=i*t
complex(8)::yc(1000)    !wave function in realtime propagaton, unnomalized
complex(8)::y1c(1000)  !wave function in realtime propagaton, nomalized
complex(8)::expTadc(1000),exptasdc(999)
complex(8)::expTbdc(1000),exptbsdc(999)
complex(8)::expVc(1000)
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
double precision::avesize(1000) !average size of the wavefunction
double precision:: nmin, nmax, nstep
double precision:: dy1dr(1000)


interface
function zeta(x)    !renormalization function
double precision::x(1000),zeta(1000)
end function
function dzeta(x)       !derivatives of renormalization function
double precision::x(1000),dzeta(1000)
end function
function trimatmult(a,b,c)
double precision::a(1000),b(999),c(1000),trimatmult(1000)
end function
function trimatmultc(a,b,c)
complex(8)::a(1000),b(999),c(1000),trimatmultc(1000)
end function
function dotprod(a,b)
double precision::a(1000),b(1000),dotprod
end function
function dotprodc(a,b)
complex(8)::a(1000),b(1000),dotprodc
end function
end interface

!open(unit=10,file='solgprealtime.txt')
open(unit=20,file='wavefunction.dat')
!open(unit=30,file='wfreal.txt')
!open(unit=40,file='wfimag.txt')
!open(unit=50,file='frho.txt')
!open(unit=60,file='wkimag2000.txt')
!open(unit=70,file='avesize2000.txt')
open(12,file='energyvsa.dat')
open(13,file='runtime.dat')



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
egtrack=0d0
Eg=0d0
istep=1
!a0=0.00218d0
rhoave=0d0
kfave=0d0
!num=50000d0
dt=0.016d0
eigpre=0d0
iexit=0
!dr=0.04d0
istep=1
e0=0d0
orbit=0d0
error=0d0
avesize=0d0



read(*,*)
read(*,*) a0initial, a0final, numpoints,  numparticle, dr
read(*,*)

read(*,*)
read(*,*)nmin, nmax, numNpoints, nnpower
read(*,*)

read(*,*)
read(*,*) napower

num=dfloat(numparticle)


do i=1,1000
rint(i)=dr*dfloat(i)
ytrial(i)=1d0/(4.6d0)**(3d0/2d0)/pi**(3d0/4d0)*exp(-(rint(i)/4.6d0)**2/2d0)*rint(i)
end do

!a0step=(a0final-a0initial)/dfloat(numpoints-1)

!a0step=(a0final**(1d0/3d0)-a0initial**(1d0/3d0))/dfloat(numpoints-1)
a0step=(a0final**(1d0/dfloat(napower))-a0initial**(1d0/dfloat(napower)) )/dfloat(numpoints-1)

!num=dfloat(numparticle)

call ZELEGA(2000,CSsav,DZ)
call WELEGA(2000,CSsav,DZ,WEsav)

!phi\[Lambda] = 1/\[Lambda]^(3/2)/\[Pi]^(3/4) *Exp[-(r/\[Lambda])^2/2]
!rint(1000) is the grid point array, ytrial(1000) is the trial wave function evaluated at each grid point.
do i=1,1000
rint(i)=dr*dfloat(i)
ytrial(i)=1d0/(4.6d0)**(3d0/2d0)/pi**(3d0/4d0)*exp(-(rint(i)/4.6d0)**2/2d0)*rint(i)
end do
do i =1,2000
kint(i)=2d0*pi/400d0*dfloat(i-1)
end do 

do i=1,2000
krho(i)=2d0*pi/400d0*dfloat(i-1)
end do 

!------------------------------------------------------------------------------------------------
!construct expTa and expTb
exptad=cosh(dt/(4d0*dr**2))
exptasd=0d0
exptasd(1:999:2)=sinh(dt/(4d0*dr**2))


exptbd(1)=1d0
exptbd(1000)=1d0
exptbd(2:999)=cosh(dt/(2d0*dr**2))
exptbsd=0d0
exptbsd(2:998:2)=sinh(dt/(2d0*dr**2))

!------------------------------------------------------------------------------------------------
!construct Tdiag,Tdiagmat,Tnondiag,Vho
Tdiag=1d0/dr**2
Tnondiag=-0.5d0/dr**2
Vho=rint**2/2d0

!------------------------------------------------------------------------------------------------
!Start imaginary time diffusion to converge the BEC ground state
call cpu_time(time1)

!nmin=1d2
!nmax=1d6
!nnpower=4
!numNpoints=11
nstep=(nmax**(1d0/dfloat(nnpower))-nmin**(1d0/dfloat(nnpower)))/dfloat(numNpoints-1)

do inn=1, numNpoints
num=(nmin**(1d0/dfloat(nnpower))+nstep*dfloat(inn-1))**nnpower
write(*,*) 'num=',num

do ia=1, numpoints    ! scanning a (scattering length )
!a0=a0initial+a0step*dfloat(ia-1)    ! in ho units
!a0=(a0initial**(1d0/3d0)+a0step*dfloat(ia-1))**3    ! in ho units
a00=(a0initial**(1d0/dfloat(napower))+a0step*dfloat(ia-1))**napower    ! in ho units
ainterval=a00/2d3
do iter=1,3

if (iter==1) then
a0=a00-ainterval 
elseif (iter==2) then
a0=a00+ainterval
else
a0=a00
end if

dlambda=1d0+(2d0/pi)**(1d0/10d0)*(num*a0)**(1d0/5d0)

dt=0.03d0


write(*,*)'dlambda=',dlambda

do i=1,1000
rint(i)=dr*dfloat(i)
ytrial(i)=1d0/(dlambda)**(3d0/2d0)/pi**(3d0/4d0)*exp(-(rint(i)/dlambda)**2/2d0)*rint(i)   !y is the full 3D rescaled wave function(y=psi*r) including the angular part Y00(although it's a constant 1/sqrt(4pi)) 
end do

y1=ytrial
!step=1
error=0d0
errortrack=0d0
istep=1
e0=0d0
orbit=0d0
error=0d0
avesize=0d0


do while(.true.)     ! start the imaginary time defussion iteration  

rho=num*y1**2/rint**2
kf=(6*pi**2*num*(y1/rint)**2)**(1d0/3d0)

Vint1=num*4*pi/3d0*(dzeta(kf*a0)*a0*((y1/rint)**2)-zeta(kf*a0)*((y1/rint)**4)**(1d0/3d0)/(6*pi**2*num)**(1d0/3d0))
Vint2=num*4*pi*zeta(kf*a0)/(6*pi**2*num)**(1d0/3d0)*((y1/rint)**4)**(1d0/3d0)


!here is right after the construction of everything for hamiltonian


!write(*,*)e0(istep)
hmtdiag=Tdiag+Vho+Vint1+Vint2

hmt2ave=4*pi*dotprod(y1,trimatmult(hmtdiag,tnondiag,trimatmult(hmtdiag,tnondiag,y1)))*dr
hmtave=4*pi*dotprod(y1,trimatmult(hmtdiag,tnondiag,y1))*dr
orbit=hmtave
error=sqrt(hmt2ave-hmtave**2)
errortrack(1:499)=errortrack(2:500)
errortrack(500)=error
if (abs(error/orbit)<1d-5 .or. istep>100000 .or. dt<1d-10) exit
if (istep>500 .and. abs(error-errortrack(200))/error<1d-3 ) then

dt=dt/2d0

exptad=cosh(dt/(4d0*dr**2))
exptasd=0d0
exptasd(1:999:2)=sinh(dt/(4d0*dr**2))


exptbd(1)=1d0
exptbd(1000)=1d0
exptbd(2:999)=cosh(dt/(2d0*dr**2))
exptbsd=0d0
exptbsd(2:998:2)=sinh(dt/(2d0*dr**2))

end if
expv=exp(-(vint1+vint2+Tdiag+Vho)*dt/2d0)
!here is right before time diffusion

y=expv*trimatmult(exptad,exptasd,trimatmult(exptbd,exptbsd,trimatmult(exptad,exptasd,expv*y1)))

dnorm=0d0
do i=1,1000
dnorm=dnorm+y(i)**2
end do
dnorm=dnorm*4*pi*dr
y1=y/sqrt(dnorm)

!here is right after time diffusion

istep=istep+1

end do


!---calculate ground state energy e0
e0=4*pi*dr*dotprod(y1,trimatmult((Tdiag+Vho)*num+Vint2*num/2d0,Tnondiag*num,y1))

if (iter==1) then 
e0pre=e0
elseif (iter==2) then
e0post=e0
endif

end do       ! end the interation of iter

!----calculate the two body contact---
denda=(e0post-e0pre)/2d0/ainterval
c2=8d0*pi*a0**2*denda

!---calculate average density with T-F approximation
call thomasfermi(a0,num,Rtf, e0tf, orbittf)

densityave=num/(4d0*pi/3d0*Rtf**3)
rrave=Rtf/num**(1d0/3d0)

!----calculate <n>, <n^2/3> and <n^1/3> with the weight of phi^2-----
density1=0d0
density23=0d0
density13=0d0

do i=1,1000
density1=density1+4d0*pi*num*y1(i)**4/rint(i)**2*dr
density23=density23+4d0*pi*num**(2d0/3d0)*y1(i)**2*(y1(i)/rint(i))**(4d0/3d0)*dr
density13=density13+4d0*pi*num**(1d0/3d0)*y1(i)**2*(y1(i)/rint(i))**(2d0/3d0)*dr
end do 


!---calcualte with Hyperspherical method------
call hypersphere(a0, num, e0hh, rminhh, omegahh, depndrhh)



!---calculate condensate fraction
!sumall=0d0
!do i=1,1000
!sumall=sumall+y1(i)*phi0(rint(i),1d0)*sqrt(4d0*pi)*rint(i)*dr
!end do 
!condensefrac=sumall**2

!y1norm=0d0
!phinorm=0d0
!do i=1,1000
!y1norm=y1norm+y1(i)**2*dr*4d0*pi
!phinorm=phinorm+phi0(rint(i),1d0)**2*rint(i)**2*dr
!end do 
!write(*,*) 'y1norm=',y1norm, 'phinorm=', phinorm,'dr=',dr

!------calculate the LHY correction-----
go to 112
eklhy=0d0
epotlhy=0d0
eintlhy=0d0
ealllhy=0d0
dy1dr(1)=(y1(2)-y1(1))/dr
dy1dr(1000)=(y1(1000)-y1(999))/dr
do i=2,999
dy1dr(i)=(y1(i+1)-y1(i-1))/2d0/dr
end do 
do i=1,1000
eklhy=eklhy+4d0*pi*num*dy1dr(i)**2/2d0*dr
epotlhy=epotlhy+4d0*pi*num*y1(i)**2*rint(i)**2/2d0*dr 
eintlhy=eintlhy+num*(num-1d0)/2d0*(4d0*pi)**2*a0*( 1d0+128d0/(15d0*sqrt(pi))*sqrt(num*(y1(i)/rint(i))**2 * a0**3)  )* (y1(i)**4/rint(i)**2) *dr
enddo
ealllhy=eklhy+epotlhy+eintlhy
112 continue
!write(10,*)'the ground state energy is'
!write(10,*) e0
!write(10,*)'the orbit energy is'
!write(10,*) orbit
!write(10,*)'the error is'
!write(10,*) error
!write(10,*)'the wave function is'
!write(20,'(1000e20.10)') y1(1:1000)
!write(10,*) 'the total step is'
!write(10,*) istep
!write(10,*) 'the final imaginary time step is'
!write(10,*) dt
!write(10,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(12,'(100e22.12)') num, a0, e0,orbit, density1, density23, density13, densityave, rrave, Rtf, e0tf, orbittf, e0hh, rminhh, omegahh, c2,error,dt, dfloat(istep), depndrhh
end do                                  !----------------- end of scaning a (scattering length)

end do                                  !----------------- end of scanning inn (number of particles)

call cpu_time(time2)

write(13,*) 'runtime=',time2-time1
stop      ! this is the end of the solving the ground state program 

!------------------------------------------------------------------------------------------------------------------------
!Start real time propagating
e0=0d0
orbit=0d0
error=0d0
a0=0.00218d0
istep=1
y1c=dcmplx(y1)
dtc=(0d0,0.0001d0)
dt=0.0001d0
expvc=(0d0,0d0)
Vho=Vho/4d0
!ddnorm=4d0*pi*matmul(transpose(conjg(y1cs)),y1cs)*dr
!write(*,*)ddnorm(1,1)
!y1c=y1cs(:,1)

!construct expTac and expTbc
exptadc=cos(dt/(4d0*dr**2))
exptasdc=(0d0,0d0)
exptasdc(1:999:2)=cj*sin(dt/(4d0*dr**2))

exptbdc(1)=(1d0,0d0)
exptbdc(1000)=(1d0,0d0)
exptbdc(2:999)=cos(dt/(2d0*dr**2))
exptbsdc=0d0
exptbsdc(2:998:2)=cj*sin(dt/(2d0*dr**2))

!start real time propagation
do while(istep<=10000)

rho=num*abs(y1c)**2/rint**2
kf=(6*pi**2*num*abs(y1c/rint)**2)**(1d0/3d0)

Vint1=num*4*pi/3d0*(dzeta(kf*a0)*a0*(abs(y1c/rint)**2)-zeta(kf*a0)*(abs(y1c/rint)**4)**(1d0/3d0)/(6*pi**2*num)**(1d0/3d0))
Vint2=num*4*pi*zeta(kf*a0)/(6*pi**2*num)**(1d0/3d0)*(abs(y1c/rint)**4)**(1d0/3d0)

expvc=exp(-(vint1+vint2+Tdiag+Vho)*dtc/2d0)

!here is right after the construnction of hamiltonian

if (mod(istep,100)==1) then

hmtdiag=Tdiag+Vho+Vint1+Vint2
orbit=real(4*pi*dotprodc(y1c,trimatmultc(dcmplx(hmtdiag),dcmplx(tnondiag),y1c))*dr)

e0=real(4*pi*dr*dotprodc(y1c,trimatmultc(dcmplx((Tdiag+Vho)*num+Vint2*num/2d0),dcmplx(Tnondiag*num),y1c)))

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
icon=int(istep/100)+1
do i=1,1000
avesize(icon)=avesize(icon)+4*pi*abs(y1c(i))**2*rint(i)*dr
end do

write(10,*) e0, orbit
write(20,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(20,*)abs(y1c(1:1000))
write(30,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(30,*)real(y1c(1:1000))
write(40,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(40,*)aimag(y1c(1:1000))
write(50,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(50,*)fkrho(1:2000)
!write(60,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!write(60,*)aimag(fkrho(1:2000))
end if

!here is right before the time propagation

y1c=expvc*trimatmultc(exptadc,exptasdc,trimatmultc(exptbdc,exptbsdc,trimatmultc(exptadc,exptasdc,expvc*y1c)))

!ddnorm=4*pi*dotprodc(y1c,y1c)*dr
!if (aimag(ddnorm)>1d-8) write(*,*) 'there is an error in ddnorm'
!cnorm(istep)=real(ddnorm)
!if (abs(1d0-cnorm(istep))>1d-10) write(*,*)'there is error in normalization'
!here is right after the time propagation

istep=istep+1

end do

!write(10,*)'the error is'
!write(10,*) errhmt(1,1)
!write(10,*)'the normalization at each time step is'
!write(10,*) (cnorm(i),i=1,istep-1,10)
!write(10,*)'the average size of the wave functino at istep/50 is'
!write(10,*) (avesize(i),i=1,500)
end program

!****************************************************************************************
function zeta(x)
double precision::x(1000),zeta(1000)
!zeta=0.395255838d0-1.138d0*atan(0.362d0-0.994d0*x)
zeta=0.395255838d0-5.272769443287473d-10-1.138d0*atan(0.362d0-0.994d0*x/1.0001131697794252d0)
end function

function dzeta(x)
double precision x(1000),dzeta(1000)
!dzeta=1.131172d0/(1d0+(0.362d0-0.994d0*x)**2)
dzeta=1.131044d0/(1d0+(0.362d0-0.9938875219683654d0*x)**2)
end function

function zetan(x)
double precision::x,zetan
zetan=0.395255838d0-5.272769443287473d-10-1.138d0*atan(0.362d0-0.994d0*x/1.0001131697794252d0)
!zetan=0.395255838d0-1.138d0*atan(0.362d0-0.994d0*x)
end function

function dzetan(x)
double precision x,dzetan
!dzetan=1.131172d0/(1d0+(0.362d0-0.994d0*x)**2)
dzetan=1.131044d0/(1d0+(0.362d0-0.9938875219683654d0*x)**2)
end function

double precision function phi0(r,lho)
implicit double precision (a-h,o-z)
double precision:: lho
pi=dacos(-1d0)
phi0=sqrt(4d0/(sqrt(pi)*lho**3))*exp(-r**2/(2d0*lho**2))
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
function trimatmult(a,b,c)
! the matrix is unitary and symmetric
!a(n) is the diagonal elements
!b(n-1) is the 1 and -1 diagonal elements
double precision::a(1000),b(999),c(1000),ac(1000),bc(999),bc2(999),trimatmult(1000)
ac=a*c
bc=b*c(2:1000)
bc2=b*c(1:999)
trimatmult(1)=ac(1)+bc(1)
trimatmult(1000)=ac(1000)+bc2(999)
do i=2,999
trimatmult(i)=ac(i)+bc(i)+bc2(i-1)
end do
end function
!****************************************************************************************
function trimatmultc(a,b,c)
! the matrix is unitary and symmetric
!a(n) is the diagonal elements
!b(n-1) is the 1 and -1 diagonal elements
complex(8)::a(1000),b(999),c(1000),ac(1000),bc(999),bc2(999),trimatmultc(1000)
ac=a*c
bc=b*c(2:1000)
bc2=b*c(1:999)
trimatmultc(1)=ac(1)+bc(1)
trimatmultc(1000)=ac(1000)+bc2(999)
do i=2,999
trimatmultc(i)=ac(i)+bc(i)+bc2(i-1)
end do
end function
!****************************************************************************************
function diagmatmult(a,c)
!the matrix is diagonal
!a(n) is the diagonal elements
!c(n) is the multiplied vector
complex(8)::a(1000),c(1000),diagmatmult(1000)
do i=1,1000
diagmatmult(i)=a(i)*c(i)
end do
end function
!****************************************************************************************
function dotprod(a,b)
double precision::a(1000),b(1000),dotprod
dotprod=0d0
do i=1,1000
dotprod=dotprod+a(i)*b(i)
end do
end function
!****************************************************************************************
function dotprodc(a,b)
complex(8)::a(1000),b(1000),dotprodc
dotprodc=(0d0,0d0)
do i=1,1000
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
