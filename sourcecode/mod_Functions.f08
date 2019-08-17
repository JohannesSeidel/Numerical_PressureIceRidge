module mod_Functions
  implicit none

contains
  function Cross (a ,b) result(axb)
    ! berechnet das Kreuzprodukt von a und b
    implicit none
    double precision, intent (in) :: a(3), b(3)
    double precision :: axb(3)
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end function Cross

  subroutine Swap_vertices (v_ix, v_iy)
    !
    ! This subroutine swaps the values of two vertices
    !
    implicit none
    double precision, intent (in out) :: v_ix(:), v_iy(:)
    double precision, allocatable     :: dummy(:)
    allocate (dummy(size (v_ix)))
    dummy = v_iy
    v_iy  = v_ix
    v_ix  = dummy
    deallocate (dummy)
  end subroutine Swap_vertices

  subroutine Bubble_Sort (a, sort_axis)
    !
    ! Bubble Sort:
    ! This subroutine sorts array columns according
    ! to an user defined row.
    !
    implicit none
    double precision, intent (in out) :: a(:, :)
    integer,          intent (in)     :: sort_axis
    integer                           :: i, j
    double precision, allocatable     :: temp(:)
    logical :: swapped
    allocate (temp(size (a, 1)))
    do j = size (a, 2)-1, 1, -1
      swapped = .false.
      do i = 1, j
        if (a(sort_axis, i) > a(sort_axis, i+1)) then
          temp      = a(:, i)
          a(:, i)   = a(:, i+1)
          a(:, i+1) = temp
          swapped = .true.
        end if
      end do
      if (.not. swapped) exit
    end do
  end subroutine Bubble_sort

  subroutine M33INV (A, AINV, OK_FLAG)
    !***********************************************************************************************************************************
    !  M33INV  -  Compute the inverse of a 3x3 matrix.
    !  from: http://web.hku.hk/~gdli/UsefulFiles/matrix/m33inv_f90.txt
    !
    !  A       = input 3x3 matrix to be inverted
    !  AINV    = output 3x3 inverse of matrix A
    !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
    !***********************************************************************************************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
    DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
    LOGICAL, INTENT(OUT) :: OK_FLAG

    DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-14
    DOUBLE PRECISION :: DET
    DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


    DET = A(1,1)*A(2,2)*A(3,3)  &
      - A(1,1)*A(2,3)*A(3,2)  &
      - A(1,2)*A(2,1)*A(3,3)  &
      + A(1,2)*A(2,3)*A(3,1)  &
      + A(1,3)*A(2,1)*A(3,2)  &
      - A(1,3)*A(2,2)*A(3,1)

    IF (ABS(DET) .LE. EPS) THEN
      AINV = 0.0D0
      OK_FLAG = .FALSE.
      RETURN
    END IF

    COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    AINV = TRANSPOSE(COFACTOR) / DET

    OK_FLAG = .TRUE.

    RETURN
  end subroutine M33INV

  subroutine init_random_seed()
    ! Code entnommen von: https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
      form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
          + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
          + dt(3) * 24_int64 * 60 * 60 * 1000 &
          + dt(5) * 60 * 60 * 1000 &
          + dt(6) * 60 * 1000 + dt(7) * 1000 &
          + dt(8)
      end if
      pid = getpid()
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
        seed(i) = lcg(t)
      end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
        s = 104729
      else
        s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed
end module mod_Functions
