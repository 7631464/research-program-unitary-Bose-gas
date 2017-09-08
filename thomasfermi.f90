
subroutine thomasfermi(a,nn,Rtf, e0, orbit)
implicit double precision (a-h,o-z)
double precision:: a, nn, Rtf, e0, orbit,kf, kf1,kfi,kfim1
double precision:: psi2arr(2000),rarr(2000)

pi=acos(-1d0)

do iR=1,1000
R1=0.1d0*dfloat(iR)
aval1=anorm(R1,a,nn,psi2max)-1d0

if (aval1>0d0 ) exit

end do
!write(*,*)'good here'
!stop
R1=0.1d0*dfloat(iR-1)
R2=0.1d0*dfloat(iR)

aval1=anorm(R1,a,nn,psi2max)-1d0
aval2=anorm(R2,a,nn,psi2max)-1d0
!write(*,*) 'good here'
!stop

istep=0
aval3=1d0
do while (abs(aval3)>1d-5 .and. istep<1000 )

istep=istep+1
R3=(R1+R2)/2d0
aval3=anorm(R3,a,nn,psi2max)-1d0

if (aval3*aval1<0d0) then
R2=R3
aval2=aval3
else if (aval3*aval2<0d0) then
R1=R3
aval1=aval3
endif

end do

write(*,*) 'a=',a,'nn=',nn, 'R3=',R3,'aval3=',aval3
!stop


Rtf=R3

psistep=psi2max/2000d0
do i=1,2000
psi2arr(i)=psistep*dfloat(i)
rarr(i)=sqrt(Rtf**2-fpsi(psi2arr(i),a,nn))
end do

dummy=0d0
do i=1,2000-1
dummy=dummy+4d0*pi*psi2arr(i)*rarr(i)**2*abs(rarr(i)-rarr(i+1))
end do

orbit=R3**2/2d0

!-----integrate using tripozoid rule to calculate vext and vint
e0=0d0
vint=0d0
vext=0d0

kf1=(6d0*pi**2*nn*psi2arr(1))**(1d0/3d0)
vint=vint+nn*(nn-1d0)/2d0*(4d0*pi)**2*(psi2arr(1)**2*rarr(1)**2*( zetan(kf1*a)/kf1) )/2d0*abs(Rtf-rarr(1))
vext=vext+nn*4d0*pi*(psi2arr(1)*rarr(1)**4)/2d0*abs(Rtf-rarr(1))/2d0
do i=2,2000
kfi= (6d0*pi**2*nn*psi2arr(i))**(1d0/3d0)
kfim1= (6d0*pi**2*nn*psi2arr(i-1))**(1d0/3d0)
vext=vext+nn*4d0*pi*(psi2arr(i)*rarr(i)**4+psi2arr(i-1)*rarr(i-1)**4)/2d0*abs(rarr(i)-rarr(i-1))/2d0
vint=vint+nn*(nn-1d0)/2d0*(4d0*pi)**2*(psi2arr(i)**2*rarr(i)**2*( zetan(kfi*a)/kfi)+psi2arr(i-1)**2*rarr(i-1)**2*( zetan(kfim1*a)/kfim1) )/2d0*abs(rarr(i)-rarr(i-1))
end do 

!do i=1,2000-1
!kf= (6d0*pi**2*nn*psi2arr(i))**(1d0/3d0)
!vext=vext+nn*4d0*pi*psi2arr(i)*rarr(i)**4*abs(rarr(i)-rarr(i+1))/2d0
!vint=vint+nn*(nn-1d0)/2d0*(4d0*pi)**2*psi2arr(i)**2*rarr(i)**2*( zetan(kf*a)/kf )*abs(rarr(i)-rarr(i+1))
!end do

e0=vint+vext

end subroutine



!========================================================================

double precision function anorm(R, a, nn, psi2max)
implicit double precision (a-h,o-z)
double precision:: a, nn, psi2max
double precision:: psi2arr(2000), rarr(2000)

!solve the equation for maximum psi at a certain R
pi=acos(-1d0)


do i=1,1000
x1=0.1d0*dfloat(i)
fval1=fpsi(x1,a,nn)-R**2
if (fval1>0d0) exit
enddo

x1=0.1d0*dfloat(i-1)
x2=0.1d0*dfloat(i)
fval1=fpsi(x1,a,nn)-R**2
fval2=fpsi(x2,a,nn)-R**2

istep=0
fval3=1d0
do while (abs(fval3)>1d-6 .and. istep<1000)!----

istep=istep+1
x3=(x1+x2)/2d0
fval3=fpsi(x3,a,nn)-R**2

if (fval1*fval3<0d0) then
x2=x3
fval2=fval3
elseif (fval2*fval3<0d0) then
x1=x3
fval1=fval3
end if

end do!----

if (fval3<0) then
psi2max=x3
else
psi2max=x1
end if

psistep=psi2max/2000d0

do i=1,2000
psi2arr(i)=psistep*dfloat(i)
rarr(i)=sqrt(R**2-fpsi(psi2arr(i),a,nn))
end do

!---integrate psi^2*r^2dr using tripozoid rule to get anorm
anorm=0d0
anorm=anorm+4d0*pi*psi2arr(1)*rarr(1)**2/2d0*abs(R-rarr(1))
do i=2,2000
anorm=anorm+4d0*pi*(psi2arr(i)*rarr(i)**2+psi2arr(i-1)*rarr(i-1)**2)/2d0*abs(rarr(i)-rarr(i-1))
end do 

!anorm=0d0
!do i=1,2000-1
!anorm=anorm+4d0*pi*psi2arr(i)*rarr(i)**2*abs(rarr(i)-rarr(i+1))
!end do

end function


double precision function fpsi(x,a,nn)     ! R^2-r^2 as a functio of psi
implicit double precision (a-h,o-z)
double precision:: nn,a,x   ! x=psi^2(x>=0), a=scattering length, nn=number of particles

pi=acos(-1d0)

fpsi=8d0*pi*(nn-1d0)*( 1d0/3d0*dzetan( (6d0*pi**2*nn)**(1d0/3d0)*x**(1d0/3d0)*a ) *a*x +2d0/3d0*zetan( (6d0*pi**2*nn)**(1d0/3d0)*x**(1d0/3d0)*a ) / (6d0*pi**2*nn)**(1d0/3d0) * x**(2d0/3d0)  )

end function




