module partition_function_table

   implicit none
   Real, allocatable :: Temperature_PF(:)
   Real(kind = 16), allocatable :: Partition_Functions(:,:)
   Integer, private :: N_nuclei
   Integer, private :: Num_Points

contains

   subroutine read_partition_file()
      !A subroutine to read in partition function data into a matrix, Partition_Functions(i,j)
      !Where i corresponds to the ith element and j is the partition function value at Temperature_PF(j)
      Implicit None

      Character(len=100) :: File
      Integer :: ierror
      Integer :: luninp
      Character(len=100) :: io_emsg
      Integer :: i,j
      Real, allocatable :: Data(:)

      Real :: size_real

      File = "partition_function_table.dat"

      Open(newunit=luninp, file=File, Status='old', Action='Read', iostat=ierror, FORM="FORMATTED", iomsg=io_emsg)

      if (ierror /= 0) then
         Write(*,*) "Failed to open"
         Stop 1
      Endif

       !First index is the number of nuclei
       Read(luninp,*,iostat=ierror) size_real
       N_nuclei = INT(size_real)

       !Next is the number of PF data points
       Read(luninp,*,iostat=ierror) size_real
       Num_Points = INT(size_real)

       Allocate(Partition_Functions(N_nuclei, Num_Points))
       Allocate(Temperature_PF(Num_Points))
       Allocate(Data(Num_Points))

       !Next is Temperature
       Read(luninp,*,iostat=ierror) Temperature_PF

       !Lastly is Partition Function Data
       Do i= 1, N_nuclei
         Read(luninp, *, iostat=ierror) Data
         Do j = 1, Num_Points
            Partition_Functions(i,j) = Data(j)
         EndDo
       EndDo

     close(unit=luninp)

   end subroutine read_partition_file

   function interpolate_partition_functions (Element, Temp)  result(PF_Val)
      !Interpolates existing Partition Function given an Element (index) and a Temp value.
      implicit none

      Integer, intent(in) :: Element
      Real(8), intent(in) :: Temp
      Integer :: closet_index
      Real, dimension(3) :: PF_Val

      Real :: dx, dxm, dxp
      Real :: Fim, Fi, Fip
      Real, dimension(Num_Points) :: Distance
      Real :: A, B, C ! Solution

      !Find closes grid point to temp provided.
      Distance = ABS(Temperature_PF-Temp)
      closet_index = MINLOC(Distance, 1)

      !Edge grid cases.
      If (closet_index == 1) then
         dx  = Temperature_PF(closet_index+1)
         dxp = Temperature_PF(closet_index+2)
         dxm = Temperature_PF(closet_index)
         Fi  = Partition_Functions(element, closet_index+1)
         Fip = Partition_Functions(element, closet_index+2)
         Fim = Partition_Functions(element, closet_index)
      Else if (closet_index == Num_Points) then
         dx  = Temperature_PF(closet_index-1)
         dxp = Temperature_PF(closet_index)
         dxm = Temperature_PF(closet_index-2)
         Fi  = Partition_Functions(element, closet_index-1)
         Fip = Partition_Functions(element, closet_index)
         Fim = Partition_Functions(element, closet_index-2)
      Else
         dx  = Temperature_PF(closet_index)
         dxp = Temperature_PF(closet_index+1)
         dxm = Temperature_PF(closet_index-1)
         Fi  = Partition_Functions(element, closet_index)
         Fip = Partition_Functions(element, closet_index+1)
         Fim = Partition_Functions(element, closet_index-1)
      Endif

      !Quadratic Interpolation
      A = ((Fim*(dx - dxm) - dxm*(Fi - Fim))*((dx - dxm)*(-1.0d0/2.0d0*dxm** &
         2 + (1.0d0/2.0d0)*dxp**2) - ((1.0d0/2.0d0)*dx**2 - 1.0d0/2.0d0* &
         dxm**2)*(-dxm + dxp)) - ((1.0d0/2.0d0)*dxm**2*(dx - dxm) - dxm*(( &
         1.0d0/2.0d0)*dx**2 - 1.0d0/2.0d0*dxm**2))*(-(Fi - Fim)*(-dxm + &
         dxp) + (-Fim + Fip)*(dx - dxm)))/((dx - dxm)*((dx - dxm)*(-1.0d0/ &
         2.0d0*dxm**2 + (1.0d0/2.0d0)*dxp**2) - ((1.0d0/2.0d0)*dx**2 - &
         1.0d0/2.0d0*dxm**2)*(-dxm + dxp)))

      B = ((Fi - Fim)*((dx - dxm)*(-1.0d0/2.0d0*dxm**2 + (1.0d0/2.0d0)*dxp** &
          2) - ((1.0d0/2.0d0)*dx**2 - 1.0d0/2.0d0*dxm**2)*(-dxm + dxp)) - &
          ((1.0d0/2.0d0)*dx**2 - 1.0d0/2.0d0*dxm**2)*(-(Fi - Fim)*(-dxm + &
          dxp) + (-Fim + Fip)*(dx - dxm)))/((dx - dxm)*((dx - dxm)*(-1.0d0/ &
          2.0d0*dxm**2 + (1.0d0/2.0d0)*dxp**2) - ((1.0d0/2.0d0)*dx**2 - &
          1.0d0/2.0d0*dxm**2)*(-dxm + dxp)))

      C = (-(Fi - Fim)*(-dxm + dxp) + (-Fim + Fip)*(dx - dxm))/((dx - dxm)*( &
          -1.0d0/2.0d0*dxm**2 + (1.0d0/2.0d0)*dxp**2) - ((1.0d0/2.0d0)*dx** &
          2 - 1.0d0/2.0d0*dxm**2)*(-dxm + dxp))

      PF_Val(1) = A + B*Temp + 0.5*C*Temp**2 !F(Temp)
      PF_Val(2) = B + C*Temp !F'(Temp)
      PF_Val(3) = C !F''(Temp)

   end function interpolate_partition_functions

end module partition_function_table
