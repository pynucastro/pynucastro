      double precision function vddot (n,dx,incx,dy,incy)

!$acc routine seq

c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      vddot = 0.0d0
      dtemp = 0.0d0
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      enddo
      vddot = dtemp
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if ( m .eq. 0 ) goto 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      enddo
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      enddo
   60 vddot = dtemp
      end


