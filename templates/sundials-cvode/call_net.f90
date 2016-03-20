program call_net
  ! This is a short program to call the network and evaluate rates
  ! Useful for debugging or closer inspection of the rate calculations

  use net_rates

  implicit none

  ! Declare variables
  character(len=*), parameter :: parameter_file_name = 'net.par'
  integer, parameter :: parameter_file_unit = 1
  double precision :: temp
  double precision, parameter :: temp_ini = 1.0d7
  double precision, parameter :: temp_fin = 1.0d10
  integer, parameter :: temp_num = 1000
  integer :: itemp, i
  double precision, dimension(number_reactions+1, temp_num) :: rxn_rate_data
  double precision :: ri
  integer :: NPROFILE
  character(len=50) :: rfmt
  character(len=50) :: hfmt
  character(len=*), parameter  :: profile_file_name = 'net_rates.dat'

  ! Initialize the network
  open(unit=parameter_file_unit, file=parameter_file_name, recl=1024, delim='APOSTROPHE')
  call init_net_pars(parameter_file_unit)
  close(unit=parameter_file_unit)
  call net_meta%initialize()
  
  ! Evaluate the temperature dependence of the rates
  do itemp = 1, temp_num
     temp = temp_ini + dble(itemp-1)*(temp_fin - temp_ini)/dble(temp_num-1)
     ri = 0.0d0
     do i = 1, net_meta%nrxn
        call net_meta%evaluate(temp, i, ri)
        rxn_rate_data(i, itemp) = ri
     end do
     rxn_rate_data(net_meta%nrxn+1, itemp) = temp
  end do

  ! Write the rates to a file for further analysis
  ! Set format for writing column entries
  NPROFILE = net_meta%nrxn + 1
  write(hfmt,'(A,I5,A)') '(', NPROFILE, '(A18))'
  write(rfmt,'(A,I5,A)') '(', NPROFILE, '(ES25.14))'
  open(unit=2, file=profile_file_name, recl=(25*NPROFILE+10), form='formatted')
  ! Write the header
  write(2, fmt=hfmt) 'c12(p,g)n13', 'n13( ,e+)c13', 'c13(p,g)n14', 'n14(p,g)o15', 'o15( ,e+)n15', 'n15(p,a)c12', 'n13(p,g)o14', 'o14( ,e+)n14', 'temp'
  ! Write the rates
  do itemp = 1, temp_num
     write(2, fmt=rfmt) (rxn_rate_data(i, itemp), i=1, NPROFILE)
  end do
  close(unit=2)
end program call_net
