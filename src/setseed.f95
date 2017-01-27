! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! This subroutine set the seed for the Fortran
! random number generator


subroutine setseed(seed)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
integer, intent(in) :: seed

! Internal arguments:
integer :: seedsize
integer, dimension(:), allocatable :: seedval

! initialize random seed
call random_seed()

! find out dimension of seed array
! (this is system dependent)
call random_seed(size=seedsize)

! allocate seedval array
allocate(seedval(seedsize))

! assign value to seedval
seedval = seed

! set user defined seed value
call random_seed(put=seedval)

end subroutine setseed
