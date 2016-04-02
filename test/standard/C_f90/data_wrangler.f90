module data_wrangler
  use net_rates

  implicit none

contains

  subroutine store_solution(ysol, net_history, k, t)
    double precision, dimension(:), intent(in)   :: ysol 
    double precision, dimension(:,:), intent(inout) :: net_history
    integer, intent(in) :: k
    double precision, intent(in) :: t
    integer :: i

    do i = 1, number_equations
       net_history(i, k) = ysol(i)
    end do
    net_history(number_equations+1, k) = t
  end subroutine store_solution
  
  subroutine write_solution(net_history, KFIN)
    double precision, dimension(:,:), intent(in) :: net_history
    integer, intent(in) :: KFIN
    integer :: K, J, NPROFILE

    ! Output Profile Misc
    character(len=50) :: rfmt
    character(len=50) :: hfmt
    character(len=*), parameter  :: profile_file_name = 'net_history.dat'

    NPROFILE = number_equations + 1
    write(*,*) 'NPROFILE: ', NPROFILE
    write(*,*) 'size(net_history, 1): ', size(net_history,1)
    write(*,*) 'size(net_history, 2): ', size(net_history,2)
    !! Save profile to file
    ! Set format for writing column entries
    write(hfmt,'(A,I5,A)') '(', NPROFILE, '(A18))'
    write(rfmt,'(A,I5,A)') '(', NPROFILE, '(ES25.14))'
    open(unit=2, file=profile_file_name, recl=(25*NPROFILE+10), form='formatted')
    write(2, fmt=hfmt) 'Y_n', 'Y_p', 'Y_he4', 'Y_c12', 'Y_ne20', 'Y_na23', 'Y_mg23', 'E_nuc', 'Time'
    do K = 1, KFIN
       write(2, fmt=rfmt) (net_history(J, K), J=1, NPROFILE)
    end do
    close(unit=2)
  end subroutine write_solution

end module data_wrangler
