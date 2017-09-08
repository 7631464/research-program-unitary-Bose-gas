subroutine hypersphere(a, nn, e0, rmin, omega, depndr)
implicit double precision (a-h,o-z)

double precision:: nn
double precision:: a
double precision:: e0, rmin, omega, depndr
!double precision:: ueffarray(10), Rarray(10),ueff2parray(10)
double precision:: fd1pcoeff(5),fd2pcoeff(5), fdueff(5), Rposition(5)

pi=dacos(-1d0)

!call cpu_time(time1)

!----------------------------- numerical parameters-----------------------------
fd1pcoeff=(/ 1d0/12d0, -2d0/3d0, 0d0, 2d0/3d0, -1d0/12d0 /)

fd2pcoeff=(/ -1d0/12d0, 4d0/3d0, -5d0/2d0, 4d0/3d0, -1d0/12d0 /)

!-------------------------------------------------------------------------
!-----------------------start calculation----------------------------

!ib=0
!ueffarray=0d0

iR=1
Rstep=1d-2
ueffpnpre=1d10

do while(.true.)  !scan the hyperradius R
R=1d0+Rstep*dfloat(iR-1)



rmax=4d0*R   ! Rmax is the final point of the integration range
nrrpoints=1000
dr=rmax/nrrpoints

duint=0d0
do irr=1,nrrpoints
rr=dr*dfloat(irr)
duint=duint+dr*(phi0(rr,R)/sqrt(4d0*pi))**(4d0-2d0/3d0)*zetan((6d0*pi**2*nn)**(1d0/3d0)*(phi0(rr,R)/sqrt(4d0*pi))**(2d0/3d0)*a)*rr**2*4d0*pi
end do

ueffpn=9d0/(8d0*R**2)+R**2/2d0+0.5d0*4d0*pi*nn**(2d0/3d0)/(6d0*pi**2)**(1d0/3d0)*duint

!ueffarray(1:9)=ueffarray(2:10)
!ueffarray(10)=ueffpn


if (ueffpn>ueffpnpre .or. iR>1000) then
Rstart=R-0.02d0
Rend=R
exit
end if

!if (ib>10 .or. iR>1000) exit
ueffpnpre=ueffpn

iR=iR+1

end do    ! end of scaning hyperradius R

if (iR>1000) then
write(*,*) 'iR >1000, does not reach the minimum of the potential'
stop
end if
write(*,*) 'Rstart=',Rstart, 'Rend=', Rend
!-------------------------end of first round of scaning ------------------------------------------

!----------------------------------------------------------------------------------------
! do a second scaning on R
fdueff=0d0
Rposition=0d0
iR=1
Rstep=1d-4
ueff1ppre=1d10


! scan the hyperradius after refining the grids,
do while(.true.)  !scan the hyperradius R
R=Rstart+Rstep*dfloat(iR-1)

rmax=4d0*R   ! Rmax is the final point of the integration range
nrrpoints=1000
dr=rmax/nrrpoints

duint=0d0
do irr=1,nrrpoints
rr=dr*dfloat(irr)
duint=duint+dr*(phi0(rr,R)/sqrt(4d0*pi))**(4d0-2d0/3d0)*zetan((6d0*pi**2*nn)**(1d0/3d0)*(phi0(rr,R)/sqrt(4d0*pi))**(2d0/3d0)*a)*rr**2*4d0*pi
end do

ueffpn=9d0/(8d0*R**2)+R**2/2d0+0.5d0*4d0*pi*nn**(2d0/3d0)/(6d0*pi**2)**(1d0/3d0)*duint

fdueff(1:4)=fdueff(2:5)
fdueff(5)=ueffpn
Rposition(1:4)=Rposition(2:5)
Rposition(5)=R


if (iR>5) then

ueff1p=0d0
ueff2p=0d0

do i=1,5
ueff1p=ueff1p+fd1pcoeff(i)*fdueff(i)
ueff2p=ueff2p+fd2pcoeff(i)*fdueff(i)
end do
ueff1p=ueff1p/Rstep
ueff2p=ueff2p/Rstep**2

if (abs(ueff1p)>abs(ueff1ppre) .and. iR>6) then
exit
end if

ueff1ppre=ueff1p

ueff2ppre=ueff2p

end if

iR=iR+1

end do    ! end of scaning hyperradius R

rmin=Rstart+Rstep*dfloat(iR-3-1)
energypn=fdueff(2)
e0=energypn*nn
omega=sqrt(ueff2ppre)
depndr=ueff1ppre

!---------------------------------------------------------
!---- interpolate the hyperradial potential function -----

!call spline(Rarray,ueffarray,10,1d40,1d40,ueff2parray)


!call splint(Rarray,ueffarray,ueff2array,10,x,y)

!---------------------------------------------------------


!write(12,'(100e22.12)') nn, a, energypn, rmin , omega, ueff1ppre


!-------------------------------------------------------------------------
!call cpu_time(time2)
!write(14,*)'runtime=',time2=time1

end subroutine

