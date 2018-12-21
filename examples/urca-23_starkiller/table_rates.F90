module table_rates
  ! Table is expected to be in terms of dens*ye and temp (non-logarithmic, cgs units)
  ! Table energy units are expected in terms of ergs

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public tabular_evaluate
  public j_na23_ne23
  public j_ne23_na23

  public jtab_mu, jtab_dq, jtab_vs, jtab_rate, jtab_nuloss, jtab_gamma

  private num_tables
  private k_drate_dt, add_vars

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

  real(rt), allocatable :: rate_table_j_na23_ne23(:,:,:), rhoy_table_j_na23_ne23(:), temp_table_j_na23_ne23(:)
  integer, allocatable  :: num_rhoy_j_na23_ne23, num_temp_j_na23_ne23, num_vars_j_na23_ne23
  character(len=50)     :: rate_table_file_j_na23_ne23
  integer               :: num_header_j_na23_ne23
  logical, parameter    :: invert_chemical_potential_j_na23_ne23 = .false.

  real(rt), allocatable :: rate_table_j_ne23_na23(:,:,:), rhoy_table_j_ne23_na23(:), temp_table_j_ne23_na23(:)
  integer, allocatable  :: num_rhoy_j_ne23_na23, num_temp_j_ne23_na23, num_vars_j_ne23_na23
  character(len=50)     :: rate_table_file_j_ne23_na23
  integer               :: num_header_j_ne23_na23
  logical, parameter    :: invert_chemical_potential_j_ne23_na23 = .true.


#ifdef AMREX_USE_CUDA

  attributes(managed) :: rate_table_j_na23_ne23, rhoy_table_j_na23_ne23, temp_table_j_na23_ne23
  attributes(managed) :: num_rhoy_j_na23_ne23, num_temp_j_na23_ne23, num_vars_j_na23_ne23

  attributes(managed) :: rate_table_j_ne23_na23, rhoy_table_j_ne23_na23, temp_table_j_ne23_na23
  attributes(managed) :: num_rhoy_j_ne23_na23, num_temp_j_ne23_na23, num_vars_j_ne23_na23

#endif

contains

  subroutine init_tabular()
    integer :: n

    allocate(num_temp_j_na23_ne23)
    allocate(num_rhoy_j_na23_ne23)
    allocate(num_vars_j_na23_ne23)
    num_temp_j_na23_ne23 = 39
    num_rhoy_j_na23_ne23 = 152
    num_vars_j_na23_ne23 = 6
    num_header_j_na23_ne23 = 7
    rate_table_file_j_na23_ne23 = trim("23Na-23Ne_electroncapture.dat")
    allocate(rate_table_j_na23_ne23(num_temp_j_na23_ne23, num_rhoy_j_na23_ne23, num_vars_j_na23_ne23))
    allocate(rhoy_table_j_na23_ne23(num_rhoy_j_na23_ne23))
    allocate(temp_table_j_na23_ne23(num_temp_j_na23_ne23))
    call init_tab_info(rate_table_j_na23_ne23, rhoy_table_j_na23_ne23, temp_table_j_na23_ne23, num_rhoy_j_na23_ne23, num_temp_j_na23_ne23, num_vars_j_na23_ne23, rate_table_file_j_na23_ne23, num_header_j_na23_ne23, invert_chemical_potential_j_na23_ne23)

    allocate(num_temp_j_ne23_na23)
    allocate(num_rhoy_j_ne23_na23)
    allocate(num_vars_j_ne23_na23)
    num_temp_j_ne23_na23 = 39
    num_rhoy_j_ne23_na23 = 152
    num_vars_j_ne23_na23 = 6
    num_header_j_ne23_na23 = 6
    rate_table_file_j_ne23_na23 = trim("23Ne-23Na_betadecay.dat")
    allocate(rate_table_j_ne23_na23(num_temp_j_ne23_na23, num_rhoy_j_ne23_na23, num_vars_j_ne23_na23))
    allocate(rhoy_table_j_ne23_na23(num_rhoy_j_ne23_na23))
    allocate(temp_table_j_ne23_na23(num_temp_j_ne23_na23))
    call init_tab_info(rate_table_j_ne23_na23, rhoy_table_j_ne23_na23, temp_table_j_ne23_na23, num_rhoy_j_ne23_na23, num_temp_j_ne23_na23, num_vars_j_ne23_na23, rate_table_file_j_ne23_na23, num_header_j_ne23_na23, invert_chemical_potential_j_ne23_na23)


  end subroutine init_tabular


  subroutine term_table_meta()

    deallocate(num_temp_j_na23_ne23)
    deallocate(num_rhoy_j_na23_ne23)
    deallocate(num_vars_j_na23_ne23)
    deallocate(rate_table_j_na23_ne23)
    deallocate(rhoy_table_j_na23_ne23)
    deallocate(temp_table_j_na23_ne23)

    deallocate(num_temp_j_ne23_na23)
    deallocate(num_rhoy_j_ne23_na23)
    deallocate(num_vars_j_ne23_na23)
    deallocate(rate_table_j_ne23_na23)
    deallocate(rhoy_table_j_ne23_na23)
    deallocate(temp_table_j_ne23_na23)


  end subroutine term_table_meta


  subroutine init_tab_info(rate_table, rhoy_table, temp_table, &
                           num_rhoy, num_temp, num_vars, &
                           rate_table_file, num_header, invert_chemical_potential)
    integer  :: num_rhoy, num_temp, num_vars, num_header
    real(rt) :: rate_table(num_temp, num_rhoy, num_vars), rhoy_table(num_rhoy), temp_table(num_temp)
    character(len=50) :: rate_table_file
    logical :: invert_chemical_potential

    real(rt), allocatable :: rate_table_scratch(:,:,:)
    integer :: i, j, k

    allocate(rate_table_scratch(num_temp, num_rhoy, num_vars+2))

    open(unit=11, file=rate_table_file)
    do i = 1, num_header
       read(11,*)
    end do
    do j = 1, num_rhoy
       do i = 1, num_temp
          read(11,*) ( rate_table_scratch(i, j, k), k=1, num_vars+2 )
       end do
       if (j/=num_rhoy) then
          read(11,*)
       end if
    end do
    close(11)

    rate_table(:,:,:) = rate_table_scratch(:,:,3:num_vars+2)

    ! Set sign for chemical potential contribution to energy generation for
    ! electron capture vs beta decays.
    if (invert_chemical_potential) then
       rate_table(:,:,jtab_mu) = -rate_table(:,:,jtab_mu)
       rate_table(:,:,jtab_vs) = -rate_table(:,:,jtab_vs)
    end if

    do i = 1, num_rhoy
       rhoy_table(i) = rate_table_scratch(1, i, 1)
    end do
    do i = 1, num_temp
       temp_table(i) = rate_table_scratch(i, 1, 2)
    end do

    deallocate(rate_table_scratch)

  end subroutine init_tab_info


  subroutine vector_index_lu(vector, fvar, index)
    !$acc routine seq

    ! Returns the greatest index of vector for which vector(index) < fvar.
    ! Return 1 if fvar < vector(1)
    ! Return size(vector)-1 if fvar > vector(size(vector))
    ! The interval [index, index+1] brackets fvar for fvar within the range of vector.
    real(rt), intent(in) :: vector(:)
    real(rt), intent(in) :: fvar
    integer, intent(out) :: index
    integer :: n, i, j, nup, ndn

    !$gpu

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
    real(rt), intent(in)  :: xlo, xhi, flo, fhi, x
    real(rt), intent(out) :: f

    !$gpu

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
    real(rt), intent(in)  :: xlo, xhi, flo, fhi, x
    real(rt), intent(out) :: f

    !$gpu

    f = ( flo * ( xhi - x ) + fhi * ( x - xlo ) ) / ( xhi - xlo )
  end subroutine bl_extrap


  subroutine get_entries(rate_table, rhoy_table, temp_table, &
                         num_rhoy, num_temp, num_vars, &
                         rhoy, temp, entries)

    integer  :: num_rhoy, num_temp, num_vars
    real(rt) :: rate_table(num_temp, num_rhoy, num_vars), rhoy_table(num_rhoy), temp_table(num_temp)
    real(rt), intent(in) :: rhoy, temp
    real(rt), dimension(num_vars+1), intent(out) :: entries

    ! The last element of entries is the derivative of rate with temperature
    ! drate_dt, evaluated by central differencing at the box corners
    ! and then performing a bilinear interpolation on those central differences.

    real(rt) :: f_im1, f_i, f_ip1, f_ip2
    real(rt) :: t_im1, t_i, t_ip1, t_ip2
    real(rt) :: drdt_i, drdt_ip1
    real(rt) :: temp_lo, temp_hi, rhoy_lo, rhoy_hi
    integer :: irhoy_lo, irhoy_hi, itemp_lo, itemp_hi
    integer :: ivar

    !$gpu

    ! Get box-corner points for interpolation
    ! This deals with out-of-range inputs via linear extrapolation
    call vector_index_lu(rhoy_table, rhoy, irhoy_lo)
    call vector_index_lu(temp_table, temp, itemp_lo)

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
    temp_lo = temp_table( itemp_lo )
    temp_hi = temp_table( itemp_hi )
    rhoy_lo = rhoy_table( irhoy_lo )
    rhoy_hi = rhoy_table( irhoy_hi )

    ! Interpolate for each table entry
    do ivar = 1, num_vars
       call bl_extrap(temp_lo, temp_hi, &
            rate_table( itemp_lo, irhoy_lo, ivar ), &
            rate_table( itemp_hi, irhoy_lo, ivar ), &
            temp, f_i)
       call bl_extrap(temp_lo, temp_hi, &
            rate_table( itemp_lo, irhoy_hi, ivar ), &
            rate_table( itemp_hi, irhoy_hi, ivar ), &
            temp, f_ip1)
       call bl_extrap(rhoy_lo, rhoy_hi, f_i, f_ip1, rhoy, entries(ivar))
    end do

    ! Calculate the derivative of rate with temperature, d(rate)/d(t)
    ! (Clamp interpolations in rhoy to avoid unphysical temperature derivatives)
    if (( itemp_lo .eq. 1 ) .or. ( itemp_lo .eq. num_temp-1 )) then
       ! We're at the first or last table cell (in temperature)
       ! First do bilinear interpolation in rhoy for the table at tlo and thi
       call bl_clamp(rhoy_lo, rhoy_hi, &
            rate_table( itemp_lo, irhoy_lo, jtab_rate ), &
            rate_table( itemp_lo, irhoy_hi, jtab_rate ), &
            rhoy, f_i)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            rate_table( itemp_hi, irhoy_lo, jtab_rate ), &
            rate_table( itemp_hi, irhoy_hi, jtab_rate ), &
            rhoy, f_ip1)
       ! Approximate d(rate)/d(t) via forward differencing
       entries(k_drate_dt) = (f_ip1 - f_i) / (temp_hi - temp_lo)
    else
       ! Approximate d(rate)/d(t) via bilinear interpolation on central differences
       t_im1 = temp_table( itemp_lo-1 )
       t_i   = temp_table( itemp_lo )
       t_ip1 = temp_table( itemp_hi )
       t_ip2 = temp_table( itemp_lo+2 )
       call bl_clamp(rhoy_lo, rhoy_hi, &
            rate_table( itemp_lo-1, irhoy_lo, jtab_rate ), &
            rate_table( itemp_lo-1, irhoy_hi, jtab_rate ), &
            rhoy, f_im1)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            rate_table( itemp_lo, irhoy_lo, jtab_rate ), &
            rate_table( itemp_lo, irhoy_hi, jtab_rate ), &
            rhoy, f_i)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            rate_table( itemp_hi, irhoy_lo, jtab_rate ), &
            rate_table( itemp_hi, irhoy_hi, jtab_rate ), &
            rhoy, f_ip1)
       call bl_clamp(rhoy_lo, rhoy_hi, &
            rate_table( itemp_lo+2, irhoy_lo, jtab_rate ), &
            rate_table( itemp_lo+2, irhoy_hi, jtab_rate ), &
            rhoy, f_ip2)
       ! Get central difference derivatives at the box corners
       drdt_i   = (f_ip1 - f_im1) / (t_ip1 - t_im1)
       drdt_ip1 = (f_ip2 - f_i)   / (t_ip2 - t_i)
       ! Interpolate in temperature
       ! (Since we're inside the table in temp, use bl_extrap, it's faster)
       call bl_extrap(t_i, t_ip1, drdt_i, drdt_ip1, temp, entries(k_drate_dt))
    end if
  end subroutine get_entries


  subroutine tabular_evaluate(rate_table, rhoy_table, temp_table, &
                              num_rhoy, num_temp, num_vars, &
                              rhoy, temp, reactvec)

    use actual_network, only: num_rate_groups
    use extern_probin_module, only: disable_fermi_heating
    implicit none

    integer  :: num_rhoy, num_temp, num_vars, num_header
    real(rt) :: rate_table(num_temp, num_rhoy, num_vars), rhoy_table(num_rhoy), temp_table(num_temp)

    real(rt), intent(in)    :: rhoy, temp
    real(rt), intent(inout) :: reactvec(num_rate_groups+2)
    real(rt) :: entries(num_vars+add_vars)

    !$gpu

    ! Get the table entries at this rhoy, temp
    call get_entries(rate_table, rhoy_table, temp_table, &
                     num_rhoy, num_temp, num_vars, &
                     rhoy, temp, entries)

    ! Recast entries into reactvec
    reactvec(1) = entries(jtab_rate)
    reactvec(2) = entries(k_drate_dt)
    reactvec(3) = 1.0d0
    reactvec(4) = 0.0d0
    reactvec(5) = entries(jtab_dq)
    if (.not. disable_fermi_heating) then
       reactvec(5) = reactvec(5) + entries(jtab_mu) - entries(jtab_vs)
    end if
    reactvec(6) = entries(jtab_gamma) - entries(jtab_nuloss)

  end subroutine tabular_evaluate

end module table_rates
