module mod_Constants
  implicit none
  double precision, parameter :: PI = 4.0d+00 * (4.0d+00 * atan (1.0d+00/5.0d+00) - atan (1.0d+00/239.0d+00)) !nach John Machin
  double precision, parameter :: G  = 9.81d+00 !m/s²
  double precision, parameter :: RHO_W = 1025d+00 !saltwater
end module mod_Constants
