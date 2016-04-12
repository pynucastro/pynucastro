program integrator
  use cvode_parameters, only: cv_data, cv_pars, net_initial_abundances
  use network
  use net_rates, only: init_reaclib, net_screening_init
  use data_wrangler, only: store_solution, write_solution
  use parameters, only: init_parameters

  implicit none

  integer KOUNTER
  integer jk

  ! Initialization
  call init_parameters()
  call init_network()
  call init_reaclib()
  call net_screening_init()

  ! Allocate data storage
  allocate( cv_data%net_history(cv_data%NEQ+1, cv_pars%NDT_SAVE) )
  allocate( cv_data%Y0(cv_data%NEQ) )
  allocate( cv_data%Y(cv_data%NEQ) )

  ! Print Physical Parameters
  write(*,*) 'Integrating starting with:'
  write(*,'(A,ES25.14)') 'T = ', cv_pars%temp
  write(*,'(A,ES25.14)') 'Density = ', cv_pars%dens

  ! Initialize the Integration
  cv_data%Y0(jn)   = net_initial_abundances%yn
  cv_data%Y0(jp)   = net_initial_abundances%yp
  cv_data%Y0(jhe4)   = net_initial_abundances%yhe4
  cv_data%Y0(jc12)   = net_initial_abundances%yc12
  cv_data%Y0(jne20)   = net_initial_abundances%yne20
  cv_data%Y0(jne23)   = net_initial_abundances%yne23
  cv_data%Y0(jna23)   = net_initial_abundances%yna23
  cv_data%Y0(jmg23)   = net_initial_abundances%ymg23
  cv_data%Y0(jenuc)   = net_initial_abundances%yenuc

  write(*,*) 'Initialization at time 0: '
  do jk = 1, cv_data%NEQ
     write(*,*) 'index: ', jk, ' : ', cv_data%Y0(jk)
  end do

  ! Integration target T_TARGET = NDT_SAVE*DT_SAVE
  cv_data%T_TARGET = DBLE(cv_pars%NDT_SAVE)*cv_pars%DT_SAVE

  ! Setup CVODE
  call FNVINITS(cv_data%KEY, cv_data%NEQ, cv_data%IER)
  call FCVMALLOC(cv_data%T0, cv_data%Y0, cv_data%METH, cv_data%ITMETH, &
       cv_data%IA_TOL, cv_pars%R_TOL, cv_pars%A_TOL, &
       cv_data%IOUT, cv_data%ROUT, cv_data%IPAR, cv_data%RPAR, cv_data%IER)
  if (cv_data%IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVMALLOC!'
     stop
  end if
  call FCVROOTINIT(cv_data%NRTFN, cv_data%IER)
  call FCVSETIIN('MAX_NSTEPS', cv_pars%MAX_NSTEPS, cv_data%IER)
  if (cv_data%IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVSETIN!'
     stop
  end if
  call FCVLAPACKDENSE(cv_data%NEQ, cv_data%IER)
  if (cv_data%IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVLAPACKDENSE!'
     stop
  end if
  call FCVLAPACKDENSESETJAC(cv_pars%SET_JAC_FLAG, cv_data%IER)
  if (cv_data%IER /= 0) then
     write(*,*) 'ERROR: occurred in FCVLAPACKDENSESETJAC!'
     stop
  end if

  call store_solution(cv_data%Y0, cv_data%net_history, 1, 0.0d0)
  ! Do the integration
  do KOUNTER = 1, cv_pars%NDT_SAVE
     cv_data%TK = min(cv_data%T0 + DBLE(KOUNTER)*cv_pars%DT_SAVE, cv_data%T_TARGET)
     write(*,*) 'Attempting integration to: ', cv_data%TK
     call FCVODE(cv_data%TK, cv_data%T, cv_data%Y, cv_pars%ITASK, cv_data%IER)
     if (cv_data%IER == 2) then
        ! The root was found
        ! Y = Y(TK) not Y(T) so we have to find Y(T)
        call FCVODE(cv_data%T, cv_data%T, cv_data%Y, cv_pars%ITASK, cv_data%IER)
        if (cv_data%IER == 0) then
           call store_solution(cv_data%Y, cv_data%net_history, KOUNTER, cv_data%T)
           cv_data%KFIN = KOUNTER
           write(*,*) 'Root Found! Integration halting.'
        else
           write(*,*) 'Error: Root find was unsuccessful!'
        end if
        exit
     else
        ! Store solution values if this was a successful integration
        call store_solution(cv_data%Y, cv_data%net_history, KOUNTER, cv_data%T)
        cv_data%KFIN = KOUNTER
     end if
  end do

  if (cv_data%KFIN == cv_pars%NDT_SAVE) then
     write(*,*) 'Desired integration time was reached!'
  end if

  ! Deallocate CVODE Memory
  call FCVROOTFREE
  call FCVFREE

  ! Write profile to output file
  call write_solution(cv_data%net_history, cv_data%KFIN)

  ! Print output to console
  write(*,'(A,ES25.14)') 'Integrated to time: ', cv_data%T
  write(*,'(A,ES25.14)') 'Final enuc: ', cv_data%Y(jenuc)
  write(*,'(A,ES25.14)') 'Average enuc_dot: ', cv_data%Y(jenuc)/cv_data%T
  write(*,'(A,ES25.14)') 'n: ', cv_data%Y(jn)
  write(*,'(A,ES25.14)') 'p: ', cv_data%Y(jp)
  write(*,'(A,ES25.14)') 'he4: ', cv_data%Y(jhe4)
  write(*,'(A,ES25.14)') 'c12: ', cv_data%Y(jc12)
  write(*,'(A,ES25.14)') 'ne20: ', cv_data%Y(jne20)
  write(*,'(A,ES25.14)') 'ne23: ', cv_data%Y(jne23)
  write(*,'(A,ES25.14)') 'na23: ', cv_data%Y(jna23)
  write(*,'(A,ES25.14)') 'mg23: ', cv_data%Y(jmg23)

  deallocate( cv_data%net_history )

end program integrator
