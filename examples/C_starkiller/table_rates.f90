module table_rates
  ! Table is expected to be in terms of dens*ye and temp (non-logarithmic, cgs units)
  ! Table energy units are expected in terms of ergs
  
  implicit none

  public table_meta, tabular_evaluate

  private num_tables, jtab_mu, jtab_dq, jtab_vs, jtab_rate, jtab_nuloss, jtab_gamma
  private k_drate_dt, add_vars
  private table_read_meta
  
  integer, parameter :: num_tables   = 0
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


  type :: table_info
     double precision, dimension(:,:,:), allocatable :: rate_table
     double precision, dimension(:), allocatable :: rhoy_table
     double precision, dimension(:), allocatable :: temp_table
     integer :: num_rhoy
     integer :: num_temp
     integer :: num_vars
  end type table_info

  type :: table_read_info
     character(len=50) :: rate_table_file
     integer :: num_header 
  end type table_read_info

  type(table_info), dimension(num_tables) :: table_meta
  type(table_read_info), dimension(num_tables) :: table_read_meta

  ! Create the device pointers for this array of derived type.
  !$acc declare create(table_meta)    
  
contains

  subroutine init_tabular()
    integer :: n

    
    do n = 1, num_tables
       call init_tab_info(table_meta(n), table_read_meta(n))
       ! For scalars or arrays with size known at compile-time, do update device
       ! to move them to the device and point the derived type pointers at them.
       !$acc update device(table_meta(n)%num_rhoy)
       !$acc update device(table_meta(n)%num_temp)
       !$acc update device(table_meta(n)%num_vars)

       ! For dynamic arrays, do enter data copyin to move their data to the device
       ! and then point the derived type pointers to these arrays on the device.
       ! If you do update device instead, the device gets the host memory addresses
       ! for these dynamic arrays instead of device memory addresses.
       !$acc enter data copyin(table_meta(n)%rate_table)
       !$acc enter data copyin(table_meta(n)%rhoy_table)
       !$acc enter data copyin(table_meta(n)%temp_table)
    end do
  end subroutine init_tabular

  subroutine term_table_meta()
    integer :: n
    do n = 1, num_tables
       call term_tab_info(table_meta(n))
    end do
  end subroutine term_table_meta

  subroutine init_tab_info(self, self_read)
    type(table_info) :: self
    type(table_read_info) :: self_read
    double precision, target, dimension(:,:,:), allocatable :: rate_table_scratch
    integer :: i, j, k

    allocate( self%rate_table( self%num_temp, self%num_rhoy, self%num_vars ) )
    allocate( self%rhoy_table( self%num_rhoy ) )
    allocate( self%temp_table( self%num_temp ) )
    allocate( rate_table_scratch( self%num_temp, self%num_rhoy, self%num_vars+2 ) )

    open(unit=11, file=self_read%rate_table_file)
    do i = 1, self_read%num_header
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
    type(table_info) :: self

    deallocate( self%rate_table )
    deallocate( self%rhoy_table )
    deallocate( self%temp_table )
  end subroutine term_tab_info

  subroutine vector_index_lu(vector, fvar, index)
    !$acc routine seq

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
    else if ( fvar .gt. vector(n) ) then
       index = n - 1
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
    !$acc routine seq
    
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
    !$acc routine seq
    
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
    !$acc routine seq
    
    type(table_info) :: self
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
    ! write(*,*) 'upper self temp table: ', self%temp_table(39)
    ! write(*,*) 'temp: ', temp
    ! write(*,*) 'itemp_lo: ', itemp_lo
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
            self%rate_table( itemp_hi, irhoy_lo, jtab_rate ), &
            self%rate_table( itemp_hi, irhoy_hi, jtab_rate ), &
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

  subroutine tabular_evaluate(self, rhoy, temp, reactvec)
    !$acc routine seq
    
    use actual_network, only: num_rate_groups
    implicit none
    
    type(table_info) :: self
    double precision, intent(in) :: rhoy, temp
    double precision, dimension(num_rate_groups+2), intent(inout) :: reactvec
    double precision, dimension(self%num_vars+add_vars) :: entries
    
    ! Get the table entries at this rhoy, temp
    call get_entries(self, rhoy, temp, entries)

    ! Recast entries into reactvec
    reactvec(1) = entries(jtab_rate)
    reactvec(2) = entries(k_drate_dt)
    reactvec(3) = 1.0d0 
    reactvec(4) = 0.0d0
    reactvec(5) = entries(jtab_dq) 
    reactvec(6) = entries(jtab_gamma) - entries(jtab_nuloss)
  end subroutine tabular_evaluate
  
end module table_rates
