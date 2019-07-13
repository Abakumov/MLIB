program modelmaker
implicit none

! Model
! nz: Dimension Model z-direction
! nx: Dimension Model x-direction
! ny: Dimension Model y-direction (this is always the third dimension !)

! Wavelet
! nts:      Number of timesteps of the wavelet
! nscomp:   Number of wavelet components
! nsrchist: Number of sources DEFINED IN THE WAVELET
! nloc:     Number of sources DEFINED IN SLOCATIONS
! nsrc:     Number of sources You really want to have (?)



parameter nz=101, nx=102, ny=103
parameter nts=20, nscomp=6, nsrchist=1
real pi
real fund
real hwert
character*6 tag
!real, dimension(nz,nx,ny,3) :: model
real, dimension(nz) :: fakemodel
real, dimension(nts,nscomp,nsrchist) :: wavelet
integer n

integer i,j,k
integer luout
integer irec
integer io_error,read_error

integer n1,n2,n3,n4

integer es,zs,ds,vs,fs,ss,ju

real vp,vshear,rho,c11,c44
real dt, dh

real vp0, vp1, vp2, vp3, vp4, vp5, vp6
real vs0, vs1, vs2, vs3, vs4, vs5, vs6

integer value

write(*,*) 'in modelmaker'

irec=1
luout=12

tag='moduli'

!velocities

vp=4000.
vs=2350.


rho=3000.

dh=10.
dt=0.8*dh/4000. !stabcrit. (mit vmax)
!dt=0.00001. !e-05 erfuellt das Kriterium
fund=40.0 !Central Frequency of the Source Signal [Hz]
pi=3.14159265

write(0,*) 'field ready'

open(convert='big_endian', unit = luout, file = tag , &
form = 'unformatted', status='replace', action = 'write', &
access = 'direct', recl=nz)

   fakemodel=0.

   do i=1,ny
   do j=1,nx
   do k=1,nz
    fakemodel(k)=vp*vp*rho
   enddo

   write(luout,rec = irec) fakemodel(:)
   irec=irec+1

   enddo
   enddo

   write(*,*) 'c11 finished'

   fakemodel=0.

   do i=1,ny
   do j=1,nx

   do k=1,nz
    fakemodel(k)=vs*vs*rho
   enddo

   write(luout,rec = irec) fakemodel(:)
   irec=irec+1

   enddo
   enddo

   write(*,*) 'c44 finished'

   fakemodel=0.00001

   do i=1,ny
   do j=1,nx
   do k=1,nz
    fakemodel(k)=rho
  enddo

   write(luout,rec = irec) fakemodel(:)
   irec=irec+1

   enddo
   enddo

   write(*,*) 'rho finished'

close(luout)

!!!!! COEFF !!!!

open(convert='big_endian', unit = luout, file = 'coeff' , &
form = 'unformatted', status='replace', action = 'write', &
access = 'direct', recl=1)

write(luout,rec = 1)  1.0
write(luout,rec = 2) -1.0

close(luout)

!!!!! SLOCATION !!!!

open(convert='big_endian', unit = luout, file = 'slocation' , &
form = 'unformatted', status='replace', action = 'write', &
access = 'direct', recl=1)

write(luout,rec = 1)     51.  !z-loc
write(luout,rec = 2)     51.  !x-loc
write(luout,rec = 3)     51.  !y-loc

close(luout)

!!!! GEOLOC := location of receivers using option wantarbseis !!!!

open(convert='big_endian', unit = luout, file = 'geoloc' , &
form = 'unformatted', status='replace', action = 'write', &
access = 'direct', recl=1)

! along x-axis
	write(luout,rec = 1)   51.   !z-loc
	write(luout,rec = 2)   10.0   !x-loc
	write(luout,rec = 3)   51.0     !y-loc

	write(luout,rec = 4)   51.   !z-loc
	write(luout,rec = 5)   91.   !x-loc
	write(luout,rec = 6)   51.     !y-loc

! along y-axis
	write(luout,rec = 7)   51.   !z-loc
	write(luout,rec = 8)   51.   !x-loc
	write(luout,rec = 9)   11.     !y-loc

	write(luout,rec = 10)   51.   !z-loc1
	write(luout,rec = 11)   51.   !x-loc
	write(luout,rec = 12)   91.     !y-loc

! along z-axis
	write(luout,rec = 13)   11.  !z-loc
	write(luout,rec = 14)   51.   !x-loc
	write(luout,rec = 15)   51.     !y-loc

	write(luout,rec = 16)   91.   !z-loc
	write(luout,rec = 17)   51.   !x-loc
	write(luout,rec = 18)   51.     !y-loc

!45 Degrees
	write(luout,rec = 19)   51.  !z-loc
        write(luout,rec = 20)   79.   !x-loc
        write(luout,rec = 21)   79.     !y-loc


	write(luout,rec = 22)   74.  !z-loc
        write(luout,rec = 23)   74.   !x-loc
        write(luout,rec = 24)   74.     !y-loc



!Geophon Linie
!write(luout,rec = 19)   400.   !z-loc
!write(luout,rec = 20)   400.0  !x-loc
!write(luout,rec = 21)   400.0  !y-loc

!write(luout,rec = 22)   410.   !z-loc
!write(luout,rec = 23)   410.0  !x-loc
!write(luout,rec = 24)   410.0  !y-loc

!write(luout,rec = 25)   420.   !z-loc
!write(luout,rec = 26)   420.0  !x-loc
!write(luout,rec = 27)   420.0  !y-loc

!write(luout,rec = 28)   430.   !z-loc
!write(luout,rec = 29)   430.0  !x-loc
!write(luout,rec = 30)   430.0  !y-loc

!write(luout,rec = 31)   440.   !z-loc
!write(luout,rec = 32)   440.0  !x-loc
!write(luout,rec = 33)   440.0  !y-loc

!write(luout,rec = 34)   450.   !z-loc
!write(luout,rec = 35)   450.0  !x-loc
!write(luout,rec = 36)   450.0  !y-loc

close(luout)

!!!!! Wavelet !!!!!



wavelet=0.
j=1
do i=-nts/2+1,nts/2,1
!hwert=-((real(i)*dt*2.0*pi*fund)**2.0)
!write(*,*) hwert
!hwert=1.0
!ricker1

! to test: s_xx
wavelet(j,2,1)= -200000000000000000000.*(0.)*(2.0*pi*fund)*(dt*real(i))*&
exp(-((dt*real(i))*2.0*pi*fund)**2.0)

! to test: s_yy
wavelet(j,4,1)= -200000000000000000000.*(0.)*(2.0*pi*fund)*(dt*real(i))*&
exp(-((dt*real(i))*2.0*pi*fund)**2.0)

! to test: s_zz
wavelet(j,1,1)= -200000000000000000000.*(0.)*(2.0*pi*fund)*(dt*real(i))*&
exp(-((dt*real(i))*2.0*pi*fund)**2.0)

! to test: s_zy
wavelet(j,5,1)= -200000000000000000000.*(0.)*(2.0*pi*fund)*(dt*real(i))*&
exp(-((dt*real(i))*2.0*pi*fund)**2.0)

! to test: s_zx
wavelet(j,3,1)= -200000000000000000000.*(0.)*(2.0*pi*fund)*(dt*real(i))*&
exp(-((dt*real(i))*2.0*pi*fund)**2.0)

! to test: s_xy
wavelet(j,6,1)= -200000000000000000000.*(-1.)*(2.0*pi*fund)*(dt*real(i))*&
exp(-((dt*real(i))*2.0*pi*fund)**2.0)

!wavelet(j,1,3)=(1-2*((real(i)*dt*2.0*pi*fund)**2.0)) * &
!  exp(-((real(i)*dt*2.0*pi*fund)**2.0))


j=j+1
enddo

open(convert='big_endian', unit = luout, file = 'wavelet' , &
form = 'unformatted', status='replace', action = 'write', &
access = 'direct', recl=1)

irec=1

   do i=1,nsrchist
   do j=1,nscomp
   do k=1,nts

   write(luout,rec = irec) wavelet(k,j,i)
   irec=irec+1

   enddo
   enddo
   enddo

close(luout)

open(unit = 25, file = 'waveletheader', status='replace')

write(25,*) 'SEPlib Headerfile: quick & dirty: please always check !!!'
write(25,*) 'so far: only n1,n2,n3 will be part of output'
write(25,*) ' '
write(25,*) 'sets next: in="wavelet"'
write(25,*) ' '
ju=nts
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n1=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
ju=nscomp
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n2=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
ju=nsrchist
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n3=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
write(25,*) 'n4=1'

close(25)



!!!!! HEADER FOR MODULI !!!!

n1=nz
n2=nx
n3=ny
n4=3

open(unit = 25, file = tag//'header', status='replace')

write(25,*) 'SEPlib Headerfile: quick & dirty: please always check !!!'
write(25,*) 'so far: only n1,n2,n3,n4 will be part of output'
write(25,*) ' '
write(25,*) 'sets next: in="./',tag,'"'
write(25,*) ' '
ju=n1
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n1=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
ju=n2
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n2=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
if(n3==1) then
ju=n4
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n3=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
write(25,*) 'n4=1'
else
ju=n3
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n3=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
ju=n4
es=modulo(ju,10)/1
zs=modulo(ju,100)/10
ds=modulo(ju,1000)/100
vs=modulo(ju,10000)/1000
fs=modulo(ju,100000)/10000
ss=modulo(ju,1000000)/100000
write(25,*) 'n4=',char(48+ss),char(48+fs),char(48+vs),char(48+ds),char(48+zs),char(48+es)
endif

close(25)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VISCOELASTIC STUFF


open(unit=20,file='parfile',status ='replace',action ='write',iostat=io_error)

write(20,*) 'parfile for Heidimod'
write(20,*)

write(20,*) 'Parameters of Model'
write(20,*) 'nz',nz
write(20,*) 'nx',nx
write(20,*) 'ny',ny
write(20,*) 'dz',dh
write(20,*) 'dx',dh
write(20,*) 'dy',dh
if (ny==1) then
write(20,*) 'ndim',2
else
write(20,*) 'ndim',3
endif
write(20,*) 'ncomp',3
! Parameters are in moduli
write(20,*)
write(20,*) 'Parameters of Proc-distribution'
write(20,*) 'procz', 20
write(20,*) 'procx', 4 
write(20,*) 'lap', 2
write(20,*) 
write(20,*) 'Parameters of spatial derivatives'
write(20,*) 'ncoef', 2
write(20,*)
write(20,*) 'Parameters of Periodic boundaries'
write(20,*) 'xcb',0
write(20,*) 'ycb',0
write(20,*)
write(20,*) 'Parameters of Damping boundaries'
write(20,*) 'ndamp',1
write(20,*) 'damping',0.000005
write(20,*)
write(20,*) 'Parameters of viscoelastic functions'
write(20,*) 'nrelaxvisc',0
! Parameters in yvisc
write(20,*)
write(20,*) 'Parameters of simulation length and output'
!write(20,*) 'nt', 200 ! Number of timesteps
write(20,*) 'nt', 130 ! Number of timesteps
write(20,*) 'snap_i', 100 ! Snapshotintervall
write(20,*) 'cpiat',0
write(20,*) 'restart',0
write(20,*) 'wantsnap',1
write(20,*) 'wantseis',0
write(20,*) 'wantstress',0
write(20,*) 'geo_depth',3
write(20,*) 'wantarbseis',1
write(20,*) 'narbgeoloc',8
! Locations are in geoloc
write(20,*)
write(20,*) 'Parameters of wavelet-input file'
write(20,*) 'nts',nts
write(20,*) 'nscomp',nscomp
write(20,*) 'nsrchist',nsrchist
write(20,*) 'dt',dt
! Source function in wavelet
write(20,*)
write(20,*) 'Parameters of source geometry'
! src_type=1 := explosion; src_type=4 := displacement src_type=5 :=body force src_type=6 :=tensor 
write(20,*) 'src_type',6
write(20,*) 'linesource',0 !3
write(20,*) 'posx',100
write(20,*) 'posy',100
write(20,*) 'posz',50
write(20,*) 'radii',50
write(20,*) 'src_1',333
write(20,*) 'src_depth',200
write(20,*) 'lsz0',0.008
write(20,*) 'lsz1',0.008
write(20,*) 'lsx0',0.003
write(20,*) 'lsx1',0.007
write(20,*) 'srcextend',4 !0
write(20,*) 'srcdamp',0.15
!
! if useloc==1 then sourcelocations are defined in slocations
! nloc is number of different source locations in slocations
!
! different situations:
! nloc > 1 and nsrchist=1 -> same source-time-function
! at different locations: status: o.k.; srcdelay defines delay ...
!
! nloc > 1 and nloc=nsrchist -> time reverse-modus
! please set srcdelay to 0
!
write(20,*) 'useloc',1
! Localtion are in slocations
write(20,*) 'nloc',1
! have to be >0 (number of sources)
write(20,*) 'srcdelay',0
write(20,*) 
write(20,*) 'Parameters of high perfomance computing'
write(20,*) 'xlargescale',1
write(20,*) 'ylargescale',ny
write(20,*) 'seisint',20
write(20,*) 'wantzsnap', 1
write(20,*) 'wantxsnap', 1
write(20,*) 'wantysnap', 1
write(20,*) 'wantzseis', 0
write(20,*) 'wantxseis', 0
write(20,*) 'wantyseis', 0
write(20,*) 'wantzaseis', 1
write(20,*) 'wantxaseis', 1
write(20,*) 'wantyaseis', 1
write(20,*) 
write(20,*) 'Parameters of time reverse modeling'
write(20,*) 'trm',0
write(20,*) 'trmi',40
write(20,*) 'wantpe',0
write(20,*) 'wantse',0
write(20,*) 'wantge',0
write(20,*) 'wantme',0
write(20,*) 'currenttrc',1
write(20,*) 'trmswitch',0
write(20,*) 'trm2zero',0






close(20)


end
