program integrator
  use cvode_parameters
  use network
  use net_rates
  use data_wrangler
  use parameters

  implicit none

  integer KOUNTER
  integer jk

  ! Initialize parameters
  call init_parameters()
  call net_meta%initialize()

  ! Allocate data storage
  allocate( cv_data%net_history(net_meta%neqs+1, cv_pars%NDT_SAVE) )
  allocate( cv_data%Y0(net_meta%neqs) )
  allocate( cv_data%Y(net_meta%neqs) )

  ! Print Physical Parameters
  write(*,*) 'Integrating starting with:'
  write(*,'(A,ES25.14)') 'T = ', cv_pars%temp
  write(*,'(A,ES25.14)') 'Density = ', cv_pars%dens

  ! Initialize the Integration
  cv_data%Y0(net_meta%ip)   = net_initial_abundances%yp
  cv_data%Y0(net_meta%ihe4)   = net_initial_abundances%yhe4
  cv_data%Y0(net_meta%ic12)   = net_initial_abundances%yc12
  cv_data%Y0(net_meta%ic13)   = net_initial_abundances%yc13
  cv_data%Y0(net_meta%in13)   = net_initial_abundances%yn13
  cv_data%Y0(net_meta%in14)   = net_initial_abundances%yn14
  cv_data%Y0(net_meta%in15)   = net_initial_abundances%yn15
  cv_data%Y0(net_meta%io14)   = net_initial_abundances%yo14
  cv_data%Y0(net_meta%io15)   = net_initial_abundances%yo15
  cv_data%Y0(net_meta%ienuc)   = net_initial_abundances%yenuc

  write(*,*) 'Initialization at time 0: '
  do jk = 1, number_equations
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
  write(*,'(A,ES25.14)') 'Final enuc: ', cv_data%Y(net_meta%ienuc)
  write(*,'(A,ES25.14)') 'Average enuc_dot: ', cv_data%Y(net_meta%ienuc)/cv_data%T
  write(*,'(A,ES25.14)') 'p: ', cv_data%Y(net_meta%ip)
  write(*,'(A,ES25.14)') 'he4: ', cv_data%Y(net_meta%ihe4)
  write(*,'(A,ES25.14)') 'c12: ', cv_data%Y(net_meta%ic12)
  write(*,'(A,ES25.14)') 'c13: ', cv_data%Y(net_meta%ic13)
  write(*,'(A,ES25.14)') 'n13: ', cv_data%Y(net_meta%in13)
  write(*,'(A,ES25.14)') 'n14: ', cv_data%Y(net_meta%in14)
  write(*,'(A,ES25.14)') 'n15: ', cv_data%Y(net_meta%in15)
  write(*,'(A,ES25.14)') 'o14: ', cv_data%Y(net_meta%io14)
  write(*,'(A,ES25.14)') 'o15: ', cv_data%Y(net_meta%io15)

  deallocate( cv_data%net_history )

end program integrator
