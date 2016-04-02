module table_rates

  ! Table is expected to be in terms of log10(dens*ye) and log10(temp)
  
  implicit none

  public table_meta
  

  integer, parameter :: num_tables   = 0
  integer, parameter :: jtab_mu      = 1
  integer, parameter :: jtab_dq      = 2
  integer, parameter :: jtab_vs      = 3
  integer, parameter :: jtab_rate    = 4
  integer, parameter :: jtab_nuloss  = 5
  integer, parameter :: jtab_gamma   = 6


  type :: table_info
     double precision, dimension(:,:,:), allocatable :: rate_table
     double precision, dimension(:), allocatable :: rhoy_table
     double precision, dimension(:), allocatable :: temp_table
     character(len=50) :: rate_table_file
     integer :: num_header 
     integer :: num_rhoy
     integer :: num_temp
     integer :: num_vars
   contains
     procedure :: initialize => init_tab_info
     procedure :: terminate => term_tab_info
     procedure :: bl_lookup => bl_lu_tab_info
  end type table_info

  type(table_info), target, dimension(num_tables) :: table_meta

contains

  subroutine init_table_meta()
    integer :: n
    do n = 1, num_tables
       call table_meta(n)%initialize()
    end do
  end subroutine init_table_meta

  subroutine term_table_meta()
    integer :: n
    do n = 1, num_tables
       call table_meta(n)%terminate()
    end do
  end subroutine term_table_meta
  
  subroutine init_tab_info(self)
    class(table_info) :: self
    double precision, target, dimension(:,:,:), allocatable :: rate_table_scratch
    integer :: i, j, k

    allocate( self%rate_table( self%num_temp, self%num_rhoy, self%num_vars ) )
    allocate( self%rhoy_table( self%num_rhoy ) )
    allocate( self%temp_table( self%num_temp ) )
    allocate( rate_table_scratch( self%num_temp, self%num_rhoy, self%num_vars+2 ) )

    open(unit=11, file=self%rate_table_file)
    do i = 1, self%num_header
       read(11,*)
    end do
    do j = 1, self%num_rhoy
       do i = 1, self%num_temp
          read(11,*) ( rate_table_scratch(i, j, k), k=1, self%num_vars+2 )
       end do
       if (j/=self%num_rhoy) then
          read(11,*)
       end if
    end do
    close(11)

    self%rate_table(:,:,:) = rate_table_scratch(:,:,3:self%num_vars+2)
    do i = 1, self%num_rhoy
       self%rhoy_table(i) = rate_table_scratch( 1, i, 1 )
    end do
    do i = 1, self%num_temp
       self%temp_table(i) = rate_table_scratch( i, 1, 2 )
    end do
    deallocate( rate_table_scratch )
  end subroutine init_tab_info

  subroutine term_tab_info(self)
    class(table_info) :: self

    deallocate( self%rate_table )
    deallocate( self%rhoy_table )
    deallocate( self%temp_table )
  end subroutine term_tab_info

  subroutine vector_index_lu(vector, fvar, index)
    ! Returns the greatest index of vector for which vector(index) < fvar.
    ! Return 1 if fvar < vector(1)
    ! Return size(vector)-1 if fvar > vector(size(vector))
    ! The interval [index, index+1] brackets fvar for fvar within the range of vector.
    double precision, intent(in) :: vector(:)
    double precision, intent(in) :: fvar
    integer, intent(out) :: index
    integer :: n, i, j, nup, ndn

    n = size(vector)
    if ( fvar .lt. vector(1) ) then
       index = 1
       return
    else if ( fvar .gt. vector(n) ) then
       index = n
       return
    else
       nup = n
       ndn = 1
       do i = 1, n
          j = ndn + (nup - ndn)/2
          if ( fvar .lt. vector(j) ) then
             nup = j
          else
             ndn = j
          end if
          if ( ((nup - ndn) .eq. 1) ) then
             index = ndn
             return
          end if
       end do
    end if
  end subroutine vector_index_lu
  
  subroutine bl_lu_tab_info(self, rhoy, temp, ivar, fvar)
    class(table_info) :: self
    double precision, intent(in) :: rhoy, temp
    integer, intent(in) :: ivar            ! variable index
    double precision, intent(out) :: fvar  ! return variable value
    double precision :: logrhoy, logtemp
    double precision :: fab, fcd, fa, fb, fc, fd
    double precision :: temp_lo, temp_hi, rhoy_lo, rhoy_hi
    integer :: irhoy_lo, irhoy_hi, itemp_lo, itemp_hi

    ! The table is in terms of log(rho*ye) and log(temp) so convert
    logrhoy = log10(rhoy)
    logtemp = log10(temp)

    ! Get box-corner points for interpolation
    ! This deals with out-of-range inputs via linear extrapolation
    call vector_index_lu(self%rhoy_table, logrhoy, irhoy_lo)
    call vector_index_lu(self%temp_table, logtemp, itemp_lo)
    irhoy_hi = irhoy_lo + 1
    itemp_hi = itemp_lo + 1
    
    ! Bilinear interpolation within the box (log=log10)
    ! The desired point is denoted by ABCD, within the box.
    ! The value of ivar at ABCD is denoted by fvar.
    ! T ^   B .      . C
    !   |
    ! g |  AB   ABCD   CD
    ! o |     .      . 
    ! l |   A          D
    !   |___________________> log rho*Ye
    temp_lo = self%temp_table( itemp_lo )
    temp_hi = self%temp_table( itemp_hi )
    rhoy_lo = self%rhoy_table( irhoy_lo )
    rhoy_hi = self%rhoy_table( irhoy_hi )
    fa = self%rate_table( itemp_lo, irhoy_lo, ivar )
    fb = self%rate_table( itemp_hi, irhoy_lo, ivar )
    fc = self%rate_table( itemp_hi, irhoy_hi, ivar )
    fd = self%rate_table( itemp_lo, irhoy_hi, ivar )
    fab = ( fa * ( temp_hi - logtemp ) + fb * ( logtemp - temp_lo ) )/( temp_hi - temp_lo )
    fcd = ( fd * ( temp_hi - logtemp ) + fc * ( logtemp - temp_lo ) )/( temp_hi - temp_lo )
    fvar = ( fab * ( rhoy_hi - logrhoy ) + fcd * ( logrhoy - rhoy_lo ) )/( rhoy_hi - rhoy_lo )
  end subroutine bl_lu_tab_info

end module table_rates
