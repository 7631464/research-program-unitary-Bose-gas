subroutine bogoliubov(a,nn,y1,orbit,eigsmall)
implicit double precision (a-h,o-z)
double precision:: nn
double precision::AAmat(100,100)!,BBmat(2000,2000)
double precision:: phi(1000),y1(1000)   !y1 is the rescaled wave function, phi=y1/r
!double precision::Tmat(1000,2),
double precision:: Vhomat(1000),Vintmat1(1000),Vintmat2(1000),rint(1000), Vintmat3(1000)
double precision:: wr(100), wi(100)  ! the real and imaginary part of the eigen values
double precision:: lvector(100,100),rvector(100,100), work1(10000),work2(1000) 
double precision:: umat(50,50),vmat(50,50),eigval(50)
double precision:: hmatband(2,1000)
double precision:: kfvec(1000)
double precision:: evalue(1000),evector(1000,1000)
double precision:: fmat(50,50),gmat(50,50),hhmat(50,50),eemat(50,50), Imat(50,50)
!(Tmat+Vhomat+Vintmat1-(epsilon+\hbar\omega))u+Vintmat2*v=0
!(Tmat+Vhomat+Vintmat1-(epsilon-\hbar\omega))v+Vintmat2*u=0
pi=acos(-1d0)

dr=0.02d0
do i=1,1000
rint(i)=dr*dfloat(i)
end do
phi=y1/rint

!write(4731,'(1000e20.10)') (phi(i),i=1,1000)
!stop
kfvec=(6d0*pi**2*nn*phi**2)**(1d0/3d0)
!write(*,*) 'good here'

!Tmat=0d0
!do i=1,1000
!Tmat(i,1)=-1d0/2d0*(-2d0)/dr**2
!Tmat(i,2)=-1d0/2d0*(1d0)/dr**2
!end do
Tdiag=-1d0/2d0*(-2d0)/dr**2
Tsubdiag=-1d0/2d0*(1d0)/dr**2

Vhomat=0d0
Vintmat1=0d0
Vintmat2=0d0
do i=1,1000
Vhomat(i)=rint(i)**2/2d0
Vintmat1(i)=2d0*cc0func(phi(i),nn,a)*phi(i)**2+cc1func(phi(i),nn,a)*abs(phi(i))**3
Vintmat2(i)=-(cc0func(phi(i),nn,a)*phi(i)**2+cc1func(phi(i),nn,a)*abs(phi(i))**3)
Vintmat3(i)=nn*4d0*pi/3d0*(dzetan(kfvec(i)*a)*a+2d0*zetan(kfvec(i)*a)/kfvec(i))*phi(i)**2
end do

!write(*,*) 'good in Vintmat3'
!stop
!if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
hmatband(1,1:1000)=Tsubdiag
hmatband(2,1:1000)=Tdiag+Vhomat+Vintmat3
call  DSBEV( 'V', 'U', 1000, 1, hmatband, 2, evalue, evector, 1000, WORK1,INFO )
if (abs(INFO)>1d-8) then 
write(*,*)'INFO=/=0',info
stop
end if
!write(*,*) 'orbit=',orbit
!do i=1,50
!write(*,*) evalue(i)
!end do 
!evector is the rescaled wave function 
! the normalization constant of phi_i is 1/sqrt(4\pi dr), that is, r*phi_i=evector(:,i)/sqrt(4\pi dr), phi_i is normalized in 3D 
fmat=0d0
gmat=0d0
hhmat=0d0
do i=1,50
do j=1,50

do m=1,1000
fmat(i,j)=fmat(i,j)+evector(m,i)*Vintmat1(m)*evector(m,j)
gmat(i,j)=gmat(i,j)+evector(m,i)*Vintmat2(m)*evector(m,j)
hhmat(i,j)=hhmat(i,j)+evector(m,i)*Vintmat3(m)*evector(m,j)
end do 

end do 
end do 

eemat=0d0
imat=0d0
do i=1,50
eemat(i,i)=evalue(i)
imat(i,i)=1d0
end do 
!-------------------------------------------------------

aamat=0d0
aamat(1:50,1:50)=eemat+fmat-hhmat-evalue(1)*imat
aamat(1:50,51:100)=gmat
aamat(51:100,1:50)=-gmat
aamat(51:100,51:100)=-(eemat+fmat-hhmat-evalue(1)*imat)

!-------------------------------------------------------
call cpu_time(time1)
call DGEEV( 'N', 'N', 100, aamat, 100, WR, WI, lvector, 100, rvector,100, WORK2, 1000, INFO )
if (abs(INFO)>1d-8) then 
write(*,*)'INFO=/=0',info
stop
end if
call cpu_time(time2)
write(*,*) 'time in DGGEV=',time2-time1

!do i=1,100
!write(*,*) wr(i),wi(i)
!end do 
!stop
!do i=1,100
!if (abs(wr(i))<3d0) then 
!write(*,*) 'freq=',wr(i)
!end if
!end do 

!do i=1,2000
!write(478,'(10e22.12)') wr(i),wi(i)
!write(479,'(2000e20.10)') (rvector(i,j), j=1,2000)
!end do 

istep=1
umat=0d0
vmat=0d0
eigval=0d0

do i=1,100

if (wr(i)>1d-4 .and. abs(wi(i))<1d-8) then
!umat(1:50,istep)=rvector(1:50,i)
!vmat(1:50,istep)=rvector(51:100,i)
eigval(istep)=wr(i)
istep=istep+1
endif

end do 

!----find the smallest value of eigval, which corresponds to the breathing mode frequency------
eigsmall=eigval(1)
do i=2,istep-1
if (eigval(i)<eigsmall) eigsmall=eigval(i)
end do 

!write(4701,'(100e20.10)') (eigval(i),i=1,istep-1) 
!do i=1,50
!write(4711,'(100e20.10)') (umat(i,j),j=1,istep-1)
!write(4721,'(100e20.10)') (vmat(i,j),j=1,istep-1)
!end do 


!stop

end subroutine



double precision function cc0func(x,nn,a)   !x=|phi|
implicit double precision (a-h,o-z)
double precision:: nn
pi=acos(-1d0)
! aa, bb, cc are the coefficients of the renormalization function zeta(x)=aa-bb*ArcTan(cc-dd*x)
aa=0.395255838d0-5.272769443287473d-10
bb=1.138d0
cc=0.362d0
dd=0.994d0/1.0001131697794252d0

cc0func=4d0/9d0*(nn-1d0)*pi*((3d0*a*bb* dd)/(1d0+(CC-6d0**(1d0/3d0)*a*dd* (x**2* nn)**(1d0/3d0)*pi**(2d0/3d0))**2) &
+ (6d0**(2d0/3d0)*(AA - BB *atan(CC - 6d0**(1d0/3d0)* a* DD* (x**2* nn)**(1d0/3d0) *pi**(2d0/3d0)) ) ) &
/((x**2 * nn)**(1d0/3d0) *pi**(2d0/3d0)) )

end function



double precision function cc1func(x,nn,a)   !x=|phi|
implicit double precision (a-h,o-z)
double precision:: nn
pi=acos(-1d0)
! aa, bb, cc are the coefficients of the renormalization function zeta(x)=aa-bb*ArcTan(cc-dd*x)
aa=0.395255838d0-5.272769443287473d-10
bb=1.138d0
cc=0.362d0
dd=0.994d0/1.0001131697794252d0

cc1func=4d0/3d0*(nn-1)*pi*2d0**(1d0/3d0)/(3d0*abs(x)*(3d0*pi)**(2d0/3d0)) * &
(-(6d0**(1d0/3d0)*AA)/(x**2* nn)**(1d0/3d0) +  &
(a* BB*DD*(6d0**(2d0/3d0)* (1d0 + CC**2) -6d0* a *CC*DD*( x**2* nn)**(1d0/3d0)*pi**(2d0/3d0)) *pi**(2d0/3d0))/ &
(1d0 + (CC - 6d0**(1d0/3d0)* a* DD*(x**2 *nn)**(1d0/3d0)*pi**(2d0/3d0))**2)**2 + &
(6d0**(1d0/3d0) *BB* atan( CC - 6d0**(1d0/3d0) *a* DD *(x**2 *nn)**(1d0/3d0) *pi**(2d0/3d0)) )/(x**2* nn)**(1d0/3d0))

end function


