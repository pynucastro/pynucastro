module table_rates
  ! Table is expected to be in terms of dens*ye and temp (non-logarithmic, cgs units)
  ! Table energy units are expected in terms of ergs
  
  implicit none

  public table_meta
  
  integer, parameter :: num_tables   = 2
  integer, parameter :: jtab_mu      = 1
  integer, parameter :: jtab_dq      = 2
  integer, parameter :: jtab_vs      = 3
  integer, parameter :: jtab_rate    = 4
  integer, parameter :: jtab_nuloss  = 5
  integer, parameter :: jtab_gamma   = 6

  ! k_drate_dt is used only for calculating the derivative
  ! of rate with temperature from the table, it isn't an index
  ! into the table but into the 'entries' array in, eg. get_entries.
  integer, parameter :: k_drate_dt   = 7
  integer, parameter :: add_vars     = 1 ! 1 Additional Var in entries
  
  integer, parameter :: j_na23_ne23   = 1
  integer, parameter :: j_ne23_na23   = 2

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
     procedure :: reaction => get_reaction
     procedure :: entries => get_entries
  end type table_info

  type(table_info), target, dimension(num_tables) :: table_meta

contains

  subroutine init_table_meta()
    integer :: n
    table_meta(j_na23_ne23)%rate_table_file = '23Na-23Ne_electroncapture.dat'
    table_meta(j_na23_ne23)%num_header = 7
    table_meta(j_na23_ne23)%num_rhoy = 152
    table_meta(j_na23_ne23)%num_temp = 39
    table_meta(j_na23_ne23)%num_vars = 6

    table_meta(j_ne23_na23)%rate_table_file = '23Ne-23Na_betadecay.dat'
    table_meta(j_ne23_na23)%num_header = 6
    table_meta(j_ne23_na23)%num_rhoy = 152
    table_meta(j_ne23_na23)%num_temp = 39
    table_meta(j_ne23_na23)%num_vars = 6

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

  subroutine bl_clamp(xlo, xhi, flo, fhi, x, f)
    ! Perform bilinear interpolation within the interval [xlo, xhi]
    ! where the function values at the endpoints are defined by:
    ! flo = f(xlo)
    ! fhi = f(xhi)
    ! Returns f(x), the values flo and fhi interpolated at x
    ! f(x) = flo if x <= xlo
    ! f(x) = fhi if x >= xhi
    double precision, intent(in)  :: xlo, xhi, flo, fhi, x
    double precision, intent(out) :: f
    if ( x .le. xlo ) then
       f = flo
    else if ( x .ge. xhi ) then
       f = fhi
    else
       f = ( flo * ( xhi - x ) + fhi * ( x - xlo ) ) / ( xhi - xlo )
    end if
  end subroutine bl_clamp

  subroutine bl_extrap(xlo, xhi, flo, fhi, x, f)
    ! Perform bilinear interpolation within the interval [xlo, xhi]
    ! where the function values at the endpoints are defined by:
    ! flo = f(xlo)
    ! fhi = f(xhi)
    ! Returns f(x), the values flo and fhi interpolated at x
    ! If x <= xlo or x >= xhi, f(x) is extrapolated at x
    double precision, intent(in)  :: xlo, xhi, flo, fhi, x
    double precision, intent(out) :: f
    f = ( flo * ( xhi - x ) + fhi * ( x - xlo ) ) / ( xhi - xlo )
  end subroutine bl_extrap
  
  subroutine get_entries(self, rhoy, temp, entries)
    class(table_info) :: self
    double precision, intent(in) :: rhoy, temp

    double precision, dimension(self%num_vars+1), intent(out) :: entries
    ! The last element of entries is the derivative of rate with temperature
    ! drate_dt, evaluated by central differencing at the box corners
    ! and then performing a bilinear interpolation on those central differences.

    double precision :: f_im1, f_i, f_ip1, f_ip2
    double precision :: t_im1, t_i, t_ip1, t_ip2
    double precision :: drdt_i, drdt_ip1
    double precision :: temp_lo, temp_hi, rhoy_lo, rhoy_hi
    integer :: irhoy_lo, irhoy_hi, itemp_lo, itemp_hi
    integer :: ivar

    ! Get box-corner points for interpolation
    ! This deals with out-of-range inputs via linear extrapolation
    call vector_index_lu(self%rhoy_table, rhoy, irhoy_lo)
    call vector_index_lu(self%temp_table, temp, itemp_lo)
    irhoy_hi = irhoy_lo + 1
    itemp_hi = itemp_lo + 1

    ! Bilinear interpolation within the box
    ! The desired point is denoted by ABCD, within the box.
    ! The value of ivar at ABCD is denoted by fvar.
    ! T ^   B .      . C
    !   |
    !   |  AB   ABCD   CD
    !   |     .      . 
    !   |   A          D
    !   |___________________> rho*Ye
    temp_lo = self%temp_table( itemp_lo )
    temp_hi = self%temp_table( itemp_hi )
    rhoy_lo = self%rhoy_table( irhoy_lo )
    rhoy_hi = self%rhoy_table( irhoy_hi )

    ! Interpolate for each table entry
    do ivar = 1, self%num_vars
       call bl_extrap(temp_lo, temp_hi, &
            self%rate_table( itemp_lo, irhoy_lo, ivar ), &
            self%rate_table( itemp_hi, irhoy_lo, ivar ), &
            temp, f_i)
       call bl_extrap(temp_lo, temp_hi, &
            self%rate_table( itemp_lo, irhoy_hi, ivar ), &
            self%rate_table( itemp_hi, irhoy_hi, ivar ), &
            temp, f_ip1)
       call bl_extrap(rhoy_lo, rhoy_hi, f_i, f_ip1, rhoy, entries(ivar))
    end do

    ! Calculate the derivative of rate with temperature, d(rate)/d(t)
    ! (Clamp interpolations in rhoy to avoid unphysical temperature derivatives)
    if (( itemp_lo .eq. 1 ) .or. ( itemp_lo .eq. self%num_temp-1 )) then
       ! We're at the first or last table cell (in temperature)
       ! First do bilinear interpolation in rhoy for the table at tlo and thi
       call bl_clamp(rhoy_lo, rhoy_hi, &
            self%rate_table( itemp_lo, irhoy_lo, jtab_rate ), &
            self%rate_table( itemp_lo, irhoy_hi, jtab_rate ), &
            rhoy, f_i)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            self%rate_table( itemp_hi, irhoy_lo, ivar ), &
            self%rate_table( itemp_hi, irhoy_hi, ivar ), &
            rhoy, f_ip1)
       ! Approximate d(rate)/d(t) via forward differencing
       entries(k_drate_dt) = (f_ip1 - f_i) / (temp_hi - temp_lo)
    else
       ! Approximate d(rate)/d(t) via bilinear interpolation on central differences
       t_im1 = self%temp_table( itemp_lo-1 )
       t_i   = self%temp_table( itemp_lo )
       t_ip1 = self%temp_table( itemp_hi )
       t_ip2 = self%temp_table( itemp_lo+2 )
       call bl_clamp(rhoy_lo, rhoy_hi, &
            self%rate_table( itemp_lo-1, irhoy_lo, jtab_rate ), &
            self%rate_table( itemp_lo-1, irhoy_hi, jtab_rate ), &
            rhoy, f_im1)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            self%rate_table( itemp_lo, irhoy_lo, jtab_rate ), &
            self%rate_table( itemp_lo, irhoy_hi, jtab_rate ), &
            rhoy, f_i)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            self%rate_table( itemp_hi, irhoy_lo, jtab_rate ), &
            self%rate_table( itemp_hi, irhoy_hi, jtab_rate ), &
            rhoy, f_ip1)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            self%rate_table( itemp_lo+2, irhoy_lo, jtab_rate ), &
            self%rate_table( itemp_lo+2, irhoy_hi, jtab_rate ), &
            rhoy, f_ip2)
       ! Get central difference derivatives at the box corners
       drdt_i   = (f_ip1 - f_im1) / (t_ip1 - t_im1)
       drdt_ip1 = (f_ip2 - f_i)   / (t_ip2 - t_i)
       ! Interpolate in temperature
       ! (Since we're inside the table in temp, use bl_extrap, it's faster)
       call bl_extrap(t_i, t_ip1, drdt_i, drdt_ip1, temp, entries(k_drate_dt))
    end if
  end subroutine get_entries

  subroutine get_reaction(self, rhoy, temp, reactvec)
    use network, only: nreactvec
    class(table_info) :: self
    double precision, intent(in) :: rhoy, temp
    double precision, dimension(nreactvec), intent(out) :: reactvec
    double precision, dimension(self%num_vars+add_vars) :: entries
    
    ! Get the table entries at this rhoy, temp
    call self%entries(rhoy, temp, entries)

    ! Recast entries into reactvec
    reactvec(1) = entries(jtab_rate)
    reactvec(2) = entries(k_drate_dt)
    reactvec(3) = 1.0d0 
    reactvec(4) = 0.0d0
    reactvec(5) = entries(jtab_dq) 
    reactvec(6) = entries(jtab_gamma) - entries(jtab_nuloss)
  end subroutine get_reaction
  
end module table_rates
