

! RANDOM NUMBER GNERATOR
! SCROLL DOWN FOR MAIN PROGRAM


module mt19937_64

  use, intrinsic :: iso_fortran_env
  implicit none

  private
  public :: init_genrand64
  public :: init_by_array64
  public :: genrand64_real1
  public :: genrand64_real2
  public :: genrand64_real3

! NOTE: genrand64_int64 is kept private, as it generates different numbers
!       compared to the original C code. This is because the original C code
!       uses unsigned integers, while Fortran relies on signed integers.
!       This, however, has no impact on the generation of real numbers
!       (they are identical to those produced by the original C code).
! public :: genrand64_int64

  integer, parameter :: r8 = real64
  integer, parameter :: i8 = int64

  integer(i8), parameter :: nn       = 312_i8
  integer(i8), parameter :: mm       = 156_i8
  integer(i8), parameter :: seed_def = 5489_i8
  integer(i8), parameter :: matrix_a = -5403634167711393303_i8
  integer(i8), parameter :: um       = -2147483648_i8 ! most significant 33 bits
  integer(i8), parameter :: lm       =  2147483647_i8 ! least significant 31 bits

  real(r8),    parameter :: pi253_1  = 1._r8/(2._r8**53 - 1._r8)
  real(r8),    parameter :: pi253    = 1._r8/(2._r8**53)
  real(r8),    parameter :: pi252    = 1._r8/(2._r8**52)

  integer(i8) :: mt(nn)       ! array for the state vector
  integer     :: mti = nn+1   ! mti==nn+1 means mt(nn) is not initialized


contains


  !-----------------------------------------------------------------------------
  ! Initializes mt(nn) with a seed

  subroutine init_genrand64(seed)
    implicit none
    integer(i8), intent(in) :: seed
    integer :: i

    mt(1) = seed
    do i = 1, nn-1
      mt(i+1) = 6364136223846793005_i8 * ieor(mt(i), ishft(mt(i), -62)) + i
    end do

    mti = nn

  end subroutine init_genrand64


  !-----------------------------------------------------------------------------
  ! Initializes by an array with array-length
  !   init_key is the array for initializing keys

  subroutine init_by_array64(init_key)
    implicit none
    integer(i8), intent(in) :: init_key(:)
    integer(i8), parameter  :: c1 = 3935559000370003845_i8
    integer(i8), parameter  :: c2 = 2862933555777941757_i8
    integer(i8) :: i, j, k, kk, key_length

    call init_genrand64(19650218_i8)
    key_length = size(init_key)
    i = 1_i8; j = 0_i8
    k = max(nn, key_length)

    do kk = 1, k
      mt(i+1) = ieor(mt(i+1), c1 * ieor(mt(i), ishft(mt(i), -62))) &
                  + init_key(j+1) + j
      i = i+1; j = j+1
      if(i >= nn) then
        mt(1) = mt(nn)
        i = 1
      end if
      if(j >= key_length) j = 0
    end do

    do kk = 1, nn-1
      mt(i+1) = ieor(mt(i+1), c2 * ieor(mt(i), ishft(mt(i), -62))) - i
      i = i+1
      if(i >= nn) then
        mt(1) = mt(nn)
        i = 1
      end if
    end do

    mt(1) = ishft(1_i8, 63)  ! MSB is 1; assuring non-zero initial array

  end subroutine init_by_array64


  !-----------------------------------------------------------------------------
  ! Generates a random number on [-2^63, 2^63-1]-interval

  integer(r8) function genrand64_int64()
    implicit none
    integer(i8) :: mag01(0:1) = (/0_i8, matrix_a/)
    integer(i8) :: x
    integer     :: i

    if(mti >= nn) then ! generate nn words at one time

      ! if init_genrand64() has not been called, a default initial seed is used
      if(mti == nn+1) call init_genrand64(seed_def)

      do i = 1, nn-mm
        x = ior(iand(mt(i),um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do

      do i = nn-mm+1, nn-1
        x = ior(iand(mt(i), um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm-nn), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do

      x = ior(iand(mt(nn), um), iand(mt(1), lm))
      mt(nn) = ieor(ieor(mt(mm), ishft(x, -1)), mag01(iand(x, 1_i8)))

      mti = 0

    end if

    mti = mti + 1
    x = mt(mti)

    x = ieor(x, iand(ishft(x,-29), 6148914691236517205_i8))
    x = ieor(x, iand(ishft(x, 17), 8202884508482404352_i8))
    x = ieor(x, iand(ishft(x, 37),   -2270628950310912_i8))
    x = ieor(x, ishft(x, -43))

    genrand64_int64 = x

  end function genrand64_int64


  !-----------------------------------------------------------------------------
  ! Generates a random number on [0,1]-real-interval

  real(r8) function genrand64_real1()
    implicit none

    genrand64_real1 = real(ishft(genrand64_int64(), -11), kind=r8) * pi253_1

  end function genrand64_real1


  !-----------------------------------------------------------------------------
  ! Generates a random number on [0,1)-real-interval

  real(r8) function genrand64_real2()
    implicit none

    genrand64_real2 = real(ishft(genrand64_int64(), -11), kind=r8) * pi253

  end function genrand64_real2


  !-----------------------------------------------------------------------------
  ! Generates a random number on (0,1)-real-interval

  real(r8) function genrand64_real3()
    implicit none

    genrand64_real3 = real(ishft(genrand64_int64(), -12), kind=r8)
    genrand64_real3 = (genrand64_real3 + 0.5_r8) * pi252

  end function genrand64_real3


end module mt19937_64






!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


! MAIN PROGRAM
! COMPILE WITH INTEL FORTRAN

!The parameters in this program correspond to the weakest aspirin dose in the paper, fig 4.
!For other scenarios, adjust parmeters as specified in the paper.  


program main
use mt19937_64
USE IFPORT
implicit none

 
real(8)     :: n1, n2, n3, n4, n5, n6
real(8)     :: R12, R14, R23, R25, R36, R45, R56
real(8)     :: gamma3, gamma4, gamma5, gamma6
real(8)     :: k1, k2, d, Ncrypt   
real(8)     :: g1, gg1, g2, g3, g4, g5, gg5, g6, g7, g8, g9, gg9, g10, g11, g12, g13, g14  
real(8)     :: aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, mm, nn, ss, sss
real(8)     :: xtime, t  
real(8)     :: ran0, xran, xr1        
integer     :: count, pathcount, patharr(1000000), iii
integer(8)  :: iran, iruns
logical     :: aspirin
character*30 :: arg
  
 
 
open(10, file='r-path8c-asp-30.dat', form='formatted' )  ! file contains information to calculate age incidence curves
open(20, file='s-path8c-asp-30.dat', form='formatted' )  ! file contains information to check evolutionary pathways to advanced adenoma 
 
call random_seed() 
call random_number(xran)
iran = int(1E6*xran)
call init_genrand64(iran)   



R12=0.000107972 
R14=9.94617E-7  
R23=0.000127717  
R25=1.59139E-6 
R36=3.73976E-6  
R45= 0.000382222  
R56=0.000452117  

k1=1039.12 
k2=171.13 

gamma3=0.2
gamma4=0.07
gamma5=0.07
gamma6=1.012    !0.8  !0.
d=0.05

Ncrypt = 1E7


 

iruns = 0 
DO WHILE (.true.)
iruns = iruns + 1
write(*,*) iruns

if (iruns >= 500000) exit
 
t = 0
n1 = Ncrypt*exp(-(R12+R14)*t)
n2 = (  Ncrypt*R12*( exp(-(R23+R25)*t) - exp(-(R12+R14)*t) )   ) / ( R12+R14-R23-R25)

n3 = 0
n4 = 0
n5 = 0
n6 = 0



R12=0.000107972 
R14=9.94617E-7  
R23=0.000127717  
R25=1.59139E-6 
R36=3.73976E-6  
R45= 0.000382222  
R56=0.000452117  

k1=1039.12 
k2=171.13 

gamma3=0.2
gamma4=0.07
gamma5=0.07
gamma6=1.012    !0.8  !0.
d=0.05

Ncrypt = 1E7






aspirin = .false. 
xtime = 0.

pathcount = 0
patharr = 0

do while (.true.)  
 
!   ASPIRIN

   !goto 888   !UNCOMMENT THIS IF YOU WANT TO RUN WITHOUT ASPIRIN
   
   ! AGE AT WHICH ASPIRIN TREATMENT IS STARTED
   if ( xtime >=30 .and. .not. aspirin   ) then
      aspirin = .true.
      gamma3 = gamma3*0.9   
      gamma4 = gamma4*0.9      
      gamma6 = gamma6*0.9    
      
      d = d*1.5
      
      R12=0.0000848292  
      R14=9.31122E-7  
      R23=0.000109836    
      R25=1.36859E-6 
      R36=3.21619E-6   
      R45=0.000328711   
      R56=0.000388821 
   endif 
   
   ! AGE AT WHICH ASPIRIN TREATMENT IS ENDED
   if ( xtime >= 40 ) then
      R12=0.000107972 
      R14=9.94617E-7  
      R23=0.000127717  
      R25=1.59139E-6 
      R36=3.73976E-6  
      R45= 0.000382222  
      R56=0.000452117  

      k1=1039.12 
      k2=171.13 

      gamma3=0.2
      gamma4=0.07
      gamma5=0.07
      gamma6=1.012    !0.8  !0.
      d=0.05

      Ncrypt = 1E7
   endif   
!888 continue   !UNCOMMENT THIS IF YOU WANT TO RUN WITHOUT ASPIRIN   
    
    
   sss = n3+n4+n5

   g1  = abs( gamma3*n3*(1-sss/k1) )
   gg1 = gamma3*n3*(1-sss/k1)
   g2 = d*n3
   g3 = R23*n2
   g4 = R36*n3

   
   g5 = abs( gamma4*n4*(1-sss/k2) ) 
   gg5 = gamma4*n4*(1-sss/k2)
   g6 = d*n4
   g7 = R14*n1
   g8 = R45*n4

   
   g9  = abs( gamma5*n5*(1-sss/k2) ) 
   gg9 = gamma5*n5*(1-sss/k2)
   g10 = d*n5
   g11 = R25*n2
   g12 = R56*n5   
       
   g13 = gamma6*n6
   g14 = d*n6
   
    
   
   ss = g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 
   if ( ss == 0 ) exit
   
   
   
   aa = g1/ss
   bb = g2/ss
   cc = g3/ss
   dd = g4/ss
   ee = g5/ss
   ff = g6/ss
   gg = g7/ss
   hh = g8/ss
   ii = g9/ss
   jj = g10/ss
   kk = g11/ss
   ll = g12/ss
   mm = g13/ss
   nn = g14/ss
 
   
   
   
   
   xran = ran0()
   
   if ( xran < aa) then
      !division
      if ( gg1 > 0 ) then      
         n3 = n3 + 1
      else if ( gg1 < 0 ) then
         n3 = n3 - 1
      endif
          
   else if ( xran < aa+bb ) then
      !death
      n3 = n3 - 1   
      
   else if ( xran < aa+bb+cc ) then
      !production
      n3 = n3 + 1      
      
   else if ( xran < aa+bb+cc+dd ) then
      !conversion
      n3 = n3 - 1      
      n6 = n6 + 1
      pathcount = pathcount + 1
      patharr(pathcount) = 1
   
      
      
      
   else if ( xran < aa+bb+cc+dd+ee) then
      !division
      if ( gg5 > 0 ) then      
         n4 = n4 + 1
      else if ( gg5 < 0 ) then
         n4 = n4 - 1
      endif   
      
   else if ( xran < aa+bb+cc+dd+ee+ff) then
      !death
      n4 = n4 - 1
      
   else if ( xran < aa+bb+cc+dd+ee+ff+gg) then
      !production
      n4 = n4 + 1
      
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh) then
      !conversion
      n4 = n4 - 1
      n5 = n5 + 1   
         
      
      
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh+ii ) then
      !division
      if ( gg9 > 0 ) then
         n5 = n5 + 1  
      else if ( gg9 < 0 ) then
         n5 = n5 - 1
      endif
   
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh+ii+jj ) then
      !death
      n5 = n5 - 1
      
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+kk ) then
      !production
      n5 = n5 + 1   
      
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+kk+ll ) then
      !conversion
      n5 = n5 - 1   
      n6 = n6 + 1
      pathcount = pathcount + 1
      patharr(pathcount) = 0
      
      
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+kk+ll+mm ) then
      !growth  
      n6 = n6 + 1   
          
   else if ( xran < aa+bb+cc+dd+ee+ff+gg+hh+ii+jj+kk+ll+mm+nn ) then
      !death  
      n6 = n6 - 1          
          
   endif   
      
      
      
   if ( n6 >=1E2) then
      write(10,'(1x,6(f20.10,3x))') xtime, n1, n2, n3, n4, n5
      write(20,'(1x, 1000000(i1,3x))') ( patharr(iii), iii=1,pathcount )
      goto 777
   endif    
  
   if ( xtime >= 80 .and. n6<1E2 ) then 
      write(10,'(1x,6(f20.10,3x))') -1., -1., -1., -1., -1., -1.
      write(20,'(1x,i2)') -1
      goto 777
   endif  
    
    
   xran = ran0()
   xtime = xtime + ( (-log(xran))/ss )  
   
   
   
   

   t = xtime
   n1 = Ncrypt*exp(-(R12+R14)*t)
   n2 = (  Ncrypt*R12*( exp(-(R23+R25)*t) - exp(-(R12+R14)*t) )   ) / ( R12+R14-R23-R25)

  

    
enddo


777 continue

ENDDO  !iruns
 



 
 
 
end program








 






real(8) function ran0()
use mt19937_64
implicit none

integer, save :: idum=0  
real(8)       :: ran3, xran
integer(8)    :: iran


ran0 = genrand64_real1()



return
end function ran0

















