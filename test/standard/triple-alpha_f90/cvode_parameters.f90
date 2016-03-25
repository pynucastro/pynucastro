module cvode_parameters

  implicit none

  type :: cv_pars_t
     double precision :: dens, temp, ye
     double precision :: DT_SAVE ! Save solution each DT_SAVE (s)
     integer :: NDT_SAVE ! (Such that NDT_SAVE*DT_SAVE=T_TARGET)
     double precision :: R_TOL = 1.0d-10 ! Relative tolerance (was using 1d-12)
     double precision :: A_TOL = 1.0d-10 ! Absolute tolerance (was using 1d-12)
     integer :: MAX_NSTEPS = 1500000 ! Maximum number of steps to solution
     integer :: SET_JAC_FLAG = 1 ! Do use my supplied Jacobian
     integer :: ITASK = 1 ! Normal mode (overshoot and interpolate)
  end type cv_pars_t
  
  type :: cv_data_t
     ! CVODE
     integer :: KEY = 1 ! Use CVODE
     ! NOTE: Below, integer*8 necessary to match the C int type expected
     integer*8 :: NEQ = 3 ! Size of ODE system
     integer ::   IER     ! Error flag
     integer ::   METH = 2 ! Use stiff BDF integration
     integer :: ITMETH = 2 ! Use Newton iteration
     integer ::  IA_TOL = 1 ! Use scalar absolute tolerance
     integer, dimension(21) ::         IOUT ! Optional integer outputs
     double precision, dimension(6) :: ROUT ! Optional real outputs
     integer :: NRTFN = 1 ! Number of root functions to solve during CVODE
     ! User integer parameters
     integer, dimension(1) :: IPAR = (/ 0 /)
     ! User real parameters
     double precision, dimension(1) :: RPAR = (/ 0.0d0 /)

     ! Integration Data Structures and Parameters
     double precision :: T0 = 0.0d0
     double precision, dimension(:), allocatable :: Y0
     double precision :: T_TARGET ! Try to integrate to TFIN = NDTV*DTSV
     double precision :: TK ! Stores the time at the current interval
     integer :: KFIN ! Stores how many intervals we went (K = 1 to KFIN)
     double precision               :: T ! Holds the Time integrated to
     double precision, dimension(:), allocatable :: Y ! Holds the solution Y integrated to
     double precision, dimension(:,:), allocatable :: net_history
  end type cv_data_t

  type(cv_data_t), save :: cv_data
  type(cv_pars_t), save :: cv_pars

contains
  subroutine cvode_init(pfile_unit)
    integer, intent(in) :: pfile_unit
    ! Read runtime parameters
    namelist /intparams/ cv_pars
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=intparams)
  end subroutine cvode_init

end module cvode_parameters
