module physical_constants

  implicit none

  double precision :: PI = 3.141592653589793238462643383279502884197d0
  double precision :: CM_PER_FM = 1.0d-13

  ! From PDG 2014
  double precision :: CLIGHT  = 299792458d2 ! cm/s
  double precision :: KB_MeV_PER_K = 8.6173324d-11 ! MeV/K
  double precision :: GNEWTON  = 6.67384d-8 ! cm^3/g/s^2
  double precision :: ERG_PER_MeV = 1.602176565d-6
  double precision :: N_AVO = 6.02214129d23

  ! From AME 2012 Atomic Mass Evaluation
  double precision :: GRAM_PER_AMU = 1660538.921d-30 ! g

end module physical_constants
