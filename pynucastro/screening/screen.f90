module screening_module

  implicit none

  private
  public :: screen5, screenz, add_screening_factor, screening_init, &
            plasma_state, fill_plasma_state

  ! Constants
  double precision, parameter :: ZERO    =  0.0d0
  double precision, parameter :: ONE     =  1.0d0
  double precision, parameter :: TWO     =  2.0d0
  double precision, parameter :: THREE   =  3.0d0
  double precision, parameter :: FOUR    =  4.0d0
  double precision, parameter :: FIVE    =  5.0d0
  double precision, parameter :: SIX     =  6.0d0
  double precision, parameter :: SEVEN   =  7.0d0
  double precision, parameter :: EIGHT   =  8.0d0
  double precision, parameter :: NINE    =  9.0d0
  double precision, parameter :: TEN     = 10.0d0

  double precision, parameter :: ELEVEN  = 11.0d0
  double precision, parameter :: TWELVE  = 12.0d0
  double precision, parameter :: FIFTEEN = 15.0d0
  double precision, parameter :: SIXTEEN = 16.0d0

  double precision, parameter :: HALF    = 0.5d0
  double precision, parameter :: THIRD   = ONE/THREE
  double precision, parameter :: FOURTH  = 0.25d0
  double precision, parameter :: FIFTH   = ONE/FIVE
  double precision, parameter :: SIXTH   = ONE/SIX
  double precision, parameter :: SEVENTH = ONE/SEVEN
  double precision, parameter :: EIGHTH  = 0.125d0
  double precision, parameter :: NINETH  = ONE/NINE
  double precision, parameter :: TENTH   = 0.10d0
  double precision, parameter :: TWELFTH = ONE/TWELVE

  double precision, parameter :: TWO3RD    = TWO/THREE
  double precision, parameter :: FOUR3RD   = FOUR/THREE
  double precision, parameter :: FIVE3RD   = FIVE/THREE
  double precision, parameter :: FIVE6TH   = FIVE/SIX

  double precision, parameter :: THREE4TH  = 0.75d0

  double precision, parameter :: FIVE12TH = FIVE/TWELVE
  double precision, parameter :: SEVEN12TH = SEVEN/TWELVE

  double precision, parameter :: FIVE32ND = FIVE/32.0d0

  !! Pi
  double precision, parameter :: M_PI    = &
       3.141592653589793238462643383279502884197d0
  double precision, parameter :: M_SQRT_PI  = &
       1.772453850905516027298167483341145182798d0

  !! Roots
  double precision, parameter :: M_SQRT_2  = &
       1.414213562373095048801688724209698078570d0
  
  integer, parameter :: nscreen_max = 500
  integer            :: nscreen = 0
  
  double precision, parameter :: fact       = 1.25992104989487d0
  double precision, parameter :: co2        = THIRD * 4.248719d3
  double precision, parameter :: gamefx     = 0.3d0
  double precision, parameter :: gamefs     = 0.8d0
  double precision, parameter :: blend_frac = 0.05d0

  double precision :: z1scr(nscreen_max)
  double precision :: z2scr(nscreen_max)
  double precision :: a1scr(nscreen_max)
  double precision :: a2scr(nscreen_max)
  
  ! zs13    = (z1+z2)**(1./3.)
  ! zhat    = combination of z1 and z2 raised to the 5/3 power
  ! zhat2   = combination of z1 and z2 raised to the 5/12 power
  ! lzav    = log of effective charge
  ! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power

  double precision :: zs13(nscreen_max)
  double precision :: zs13inv(nscreen_max)
  double precision :: zhat(nscreen_max)
  double precision :: zhat2(nscreen_max)
  double precision :: lzav(nscreen_max)
  double precision :: aznut(nscreen_max)

  type :: plasma_state

     double precision :: qlam0z
     double precision :: qlam0zdt
     double precision :: qlam0zdd

     double precision :: taufac
     double precision :: taufacdt

     double precision :: aa
     double precision :: daadt
     double precision :: daadd

  end type plasma_state

contains

  subroutine screening_init()
    ! This routine assumes that we have already filled z1scr, z2scr,
    ! a1scr, and a2scr.
    
    integer :: i
    
    do i = 1, nscreen

       zs13(i)    = (z1scr(i) + z2scr(i))**THIRD
       zs13inv(i) = ONE/zs13(i)
       zhat(i)    = (z1scr(i) + z2scr(i))**FIVE3RD  - z1scr(i)**FIVE3RD - z2scr(i)**FIVE3RD
       zhat2(i)   = (z1scr(i) + z2scr(i))**FIVE12TH - z1scr(i)**FIVE12TH -z2scr(i)**FIVE12TH
       lzav(i)    = FIVE3RD * log(z1scr(i)*z2scr(i)/(z1scr(i) + z2scr(i)))
       aznut(i)   = (z1scr(i)**2 * z2scr(i)**2 * a1scr(i)*a2scr(i) / (a1scr(i) + a2scr(i)))**THIRD
       
    enddo

  end subroutine screening_init


  subroutine add_screening_factor(z1, a1, z2, a2)

    ! this is only called at initialization

    double precision :: z1, a1, z2, a2

    nscreen = nscreen + 1
    
    z1scr(nscreen) = z1
    a1scr(nscreen) = a1
    z2scr(nscreen) = z2
    a2scr(nscreen) = a2

  end subroutine add_screening_factor


  subroutine fill_plasma_state(state, temp, dens, y)

    use network, only: nspec, zion

    ! Input variables

    type (plasma_state) :: state
    double precision :: temp, dens, y(nspec)

    ! Local variables

    double precision :: abar, zbar, z2bar
    double precision :: ytot, rr, tempi, dtempi, deni
    double precision :: pp, qq, dppdt, xni
!    double precision :: dppdd

    abar   = ONE / sum(y)
    zbar   = sum(zion * y) * abar
    z2bar  = sum(zion**2 * y) * abar    
    
    ytot             = ONE / abar
    rr               = dens * ytot
    tempi            = ONE / temp
    dtempi           = -tempi * tempi
    deni             = ONE / dens

    pp               = sqrt(rr*tempi*(z2bar + zbar))
    qq               = HALF/pp *(z2bar + zbar)
    dppdt            = qq*rr*dtempi
    !dppdd            = qq * ytot * tempi

    state % qlam0z   = 1.88d8 * tempi * pp
    state % qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
    !state % qlam0zdd = 1.88d8 * tempi * dppdd

    state % taufac   = co2 * tempi**THIRD
    state % taufacdt = -THIRD * state % taufac * tempi

    qq               = rr * zbar
    xni              = qq**THIRD
    !dxnidd           = THIRD * xni * deni

    state % aa       = 2.27493d5 * tempi * xni
    state % daadt    = 2.27493d5 * dtempi * xni
    !state % daadd    = 2.27493d5 * tempi * dxnidd

  end subroutine fill_plasma_state



  subroutine screen5(state,jscreen,scor,scordt,scordd)

    implicit none

    ! this subroutine calculates screening factors and their derivatives
    ! for nuclear reaction rates in the weak, intermediate and strong regimes.
    ! based on graboske, dewit, grossman and cooper apj 181 457 1973 for
    ! weak screening. based on alastuey and jancovici apj 226 1034 1978,
    ! with plasma parameters from itoh et al apj 234 1079 1979, for strong
    ! screening.

    ! input:
    ! state   = plasma state (T, rho, abar, zbar, etc.)
    ! jscreen = counter of which reaction is being calculated

    ! output:
    ! scor    = screening correction
    ! scordt  = derivative of screening correction with temperature
    ! scordd  = derivative of screening correction with density


    ! declare the pass        
    integer             :: jscreen
    type (plasma_state) :: state
    double precision    :: scor, scordt, scordd


    ! local variables
    double precision :: z1, a1, z2, a2
    
    double precision :: bb,cc,dccdt, &
                        qq,dqqdt,rr,drrdt, &
                        ss,dssdt,tt,dttdt,uu,duudt, &
                        vv,dvvdt,a3,da3, &
                        h12w,dh12wdt,h12,dh12dt, &
                        h12x,dh12xdt,alfa,beta, &
                        gamp,gampdt, &
                        gamef,gamefdt, &
                        tau12,tau12dt,alph12,alph12dt, &
                        xlgfac,dxlgfacdt, &
                        gamp14,gamp14dt
!    double precision :: dccdd,dqqdd,dvvdd,drrdd,dssdd,dttdd,duudd
!    double precision :: dh12dd,dh12wdd,dh12xdd,alph12dd
!    double precision :: gampdd,gamefdd,dxlgcfacdd,gamp14dd

    ! Get the ion data based on the input index
    
    z1 = z1scr(jscreen)
    a1 = a1scr(jscreen)
    z2 = z2scr(jscreen)
    a2 = a2scr(jscreen)

    ! calculate individual screening factors
    bb       = z1 * z2
    gamp     = state % aa
    gampdt   = state % daadt
    !gampdd   = state % daadd

    qq       = fact * bb * zs13inv(jscreen)
    gamef    = qq * gamp
    gamefdt  = qq * gampdt
    !gamefdd  = qq * gampdd

    tau12    = state % taufac * aznut(jscreen)
    tau12dt  = state % taufacdt * aznut(jscreen)

    qq       = ONE/tau12
    alph12   = gamef * qq
    alph12dt = (gamefdt - alph12*tau12dt) * qq
    !alph12dd = gamefdd * qq



    ! limit alph12 to 1.6 to prevent unphysical behavior.
    ! this should really be replaced by a pycnonuclear reaction rate formula
    if (alph12 .gt. 1.6) then
       alph12   = 1.6d0
       alph12dt = ZERO
       !alph12dd = ZERO

       gamef    = 1.6d0 * tau12
       gamefdt  = 1.6d0 * tau12dt
       !gamefdd  = ZERO

       qq       = zs13(jscreen)/(fact * bb)
       gamp     = gamef * qq
       gampdt   = gamefdt * qq
       !gampdd   = ZERO
    end if



    ! weak screening regime
    h12w    = bb * state % qlam0z
    dh12wdt = bb * state % qlam0zdt
    !dh12wdd = bb * qlam0zdd

    h12     = h12w
    dh12dt  = dh12wdt
    !dh12dd  = dh12wdd



    ! intermediate and strong sceening regime
    if (gamef .gt. gamefx) then

       gamp14   = gamp**FOURTH
       rr       = ONE/gamp
       qq       = 0.25d0*gamp14*rr
       gamp14dt = qq * gampdt
       !gamp14dd = qq * gampdd

       cc       =   0.896434d0 * gamp * zhat(jscreen) &
            - 3.44740d0  * gamp14 * zhat2(jscreen) &
            - 0.5551d0   * (log(gamp) + lzav(jscreen)) &
            - 2.996d0

       dccdt    =   0.896434d0 * gampdt * zhat(jscreen) &
            - 3.44740d0  * gamp14dt * zhat2(jscreen) &
            - 0.5551d0*rr*gampdt

       !dccdd    =   0.896434d0 * gampdd * zhat(jscreen) &
       !     - 3.44740d0  * gamp14dd * zhat2(jscreen) &
       !     - 0.5551d0*rr*gampdd

       a3     = alph12 * alph12 * alph12
       da3    = 3.0d0 * alph12 * alph12

       qq     = 0.014d0 + 0.0128d0*alph12
       dqqdt  = 0.0128d0*alph12dt
       !dqqdd  = 0.0128d0*alph12dd

       rr     = FIVE32ND - alph12*qq
       drrdt  = -(alph12dt*qq + alph12*dqqdt)
       !drrdd  = -(alph12dd*qq + alph12*dqqdd)

       ss     = tau12*rr
       dssdt  = tau12dt*rr + tau12*drrdt
       !dssdd  = tau12*drrdd

       tt     =  -0.0098d0 + 0.0048d0*alph12
       dttdt  = 0.0048d0*alph12dt
       !dttdd  = 0.0048d0*alph12dd

       uu     =  0.0055d0 + alph12*tt
       duudt  = alph12dt*tt + alph12*dttdt
       !duudd  = alph12dd*tt + alph12*dttdd

       vv   = gamef * alph12 * uu
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt
       !dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd

       h12     = cc - a3 * (ss + vv)
       rr      = da3 * (ss + vv)
       dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
       !dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

       rr     =  ONE - 0.0562d0*a3
       ss     =  -0.0562d0*da3
       drrdt  = ss*alph12dt
       !drrdd  = ss*alph12dd

       if (rr .ge. 0.77d0) then
          xlgfac    = rr
          dxlgfacdt = drrdt
          !dxlgfacdd = drrdd
       else
          xlgfac    = 0.77d0
          dxlgfacdt = ZERO
          !dxlgfacdd = ZERO
       end if


       h12    = log(xlgfac) + h12
       rr     = ONE/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt
       !dh12dd = rr*dxlgfacdd + dh12dd


       if (gamef .le. gamefs) then
          rr     =  2.0d0*(gamefs - gamef)
          drrdt  = -2.0d0*gamefdt
          !drrdd  = -2.0d0*gamefdd

          ss     = 2.0d0*(gamef - gamefx)
          dssdt  = 2.0d0*gamefdt
          !dssdd  = 2.0d0*gamefdd


          ! store current values for possible blending
          h12x    = h12
          dh12xdt = dh12dt
          !dh12xdd = dh12dd

          vv     = h12
          h12    = h12w*rr + vv*ss
          dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
          !dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd

          ! blend the transition region - from bill paxton
          if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
             alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
             alfa   = HALF * (ONE - cos(M_PI*alfa))
             beta   = ONE - alfa
             h12    = alfa * h12 + beta * h12x
             dh12dt = alfa * dh12dt + beta * dh12xdt
             !dh12dd = alfa * dh12dd + beta * dh12xdd
          end if
       end if


       ! end of intermediate and strong screening if
    end if


    ! machine limit the output
    ! further limit to avoid the pycnonuclear regime
    h12    = max(min(h12, 30.0d0), ZERO)
    scor   = exp(h12)
    if (h12 .eq. 30.0d0) then
       scordt = ZERO
       !scordd = ZERO
    else
       scordt = scor * dh12dt
       !scordd = scor * dh12dd
    end if

  end subroutine screen5

  subroutine screenz (t,d,z1,z2,a1,a2,ymass,scfac,dscfacdt)

    use network, only: aion, zion, nspec

    implicit none

    double precision :: t, d, z1, z2, a1, a2
    double precision :: ymass(nspec)
    double precision :: scfac
    double precision :: dscfacdt

    ! this subroutine calculates screening factors for nuclear reaction
    ! rates in the weak, intermediate , and strong regimes given the
    ! temperature (t--degk), the density (d--g/cc), the atomic numbers
    ! and weights of the elements in the reaction channel with the
    ! largest coulomb barrier (z1,z2,a1,a2), and the mean plasma
    ! parameters calculated in main and passed over in common aver:
    ! (mean atomic number--zbar, mean square of the atomic
    ! number--z2bar, mean atomic weight--abar, and total number of moles
    ! of nuclei per gram--ytot1).  the unscreened rate is to be
    ! multiplied by the dimensionless the treatment is based on
    ! graboske, dewit, grossman, and cooper ap j. 181,457 (1973) for
    ! weak screening and on alastuey and jancovici, ap.j. 226, 1034,
    ! 1978, with plasma parameters from itoh, totsuji, setsuo, and
    ! dewitt, ap.j. 234, 1079,1979, for strong screening (rkw
    ! modification).

    !.... last revision 15 nov 1982

    double precision abar, zbar, ytot1, z2bar, theta

    integer iy

    double precision qlam0, ztilda, qlam0z, gamp, taufac
    double precision dqlam0dt, dqlam0zdt, dgampdt, dtaufacdt
    double precision zhat, zhat2, gamef, tau12, alph12
    double precision dgamefdt, dtau12dt, dalph12dt 
    double precision h12w, h12, c, h12fac
    double precision dh12wdt, dh12dt, dcdt


    ! calculate averages for screening routine
    ! nb  y = x/a with x the mass fraction
    ! zi and ai are the nuclear chage and atomic mass number
    ! respectively

    ! this part came in through a common block in Kepler -- do it
    ! directly here
    abar=0.d0
    zbar=0.d0
    ytot1=0.d0
    z2bar=0.d0

    do iy = 1, nspec
       ytot1 = ytot1 + ymass(iy)
       z2bar = z2bar + zion(iy)**2*ymass(iy)
       abar = abar + aion(iy)*ymass(iy)
       zbar = zbar + zion(iy)*ymass(iy)
    enddo

    z2bar=z2bar/ytot1
    abar=abar/ytot1
    zbar=zbar/ytot1

    ! resume original Kepler screen...
    theta=1.d0
    ytot1=1.d0/abar

    !.... calculate average plasma parameters
    !....
    if ((z1*z2) > 0.d0) then

       qlam0=1.88d+8*sqrt(d/(abar*t**3))
       dqlam0dt=-1.5d+0*qlam0 / t

       ztilda=sqrt(z2bar+zbar*theta)

       qlam0z=qlam0*ztilda
       dqlam0zdt=ztilda*dqlam0dt

       gamp=2.27493d+5*(d*zbar*ytot1)**THIRD/t
       dgampdt=-gamp/t

       taufac=4.248719d+3/t**THIRD
       dtaufacdt=-THIRD*taufac/t

       !.... calculate screening factor
       !.... approx. for strong screening only good for alpha .lt. 1.6

       zhat=(z1+z2)**FIVE3RD-z1**FIVE3RD-z2**FIVE3RD
       zhat2=(z1+z2)**FIVE12TH-z1**FIVE12TH-z2**FIVE12TH

       gamef=2.d0**THIRD*gamp*z1*z2/(z1+z2)**THIRD
       dgamefdt=gamef*dgampdt/gamp

       tau12=taufac*(z1**2*z2**2*a1*a2/(a1+a2))**THIRD
       dtau12dt=tau12*dtaufacdt/taufac

       alph12=3.d0*gamef/tau12
       dalph12dt=alph12*(dgamefdt/gamef - dtau12dt/tau12)

       !....
       !.... limit alph12 to 1.6 to prevent unphysical behavior
       !.... (h dec. as rho inc.) at high rho.  this should really
       !.... be replaced by a pycnonuclear reaction rate formula.
       !....
       if (alph12 > 1.6d0) then

          alph12=1.6d0
          dalph12dt=0.0d0

          gamef=1.6d0*tau12/3.d0
          dgamefdt=gamef*dtau12dt/tau12

          gamp=gamef*(z1+z2)**THIRD/(2.d0**THIRD*z1*z2)
          dgampdt=gamp*dgamefdt/gamef

       endif

       h12w=z1*z2*qlam0z
       dh12wdt=h12w*dqlam0zdt/qlam0z

       h12=h12w
       dh12dt=dh12wdt

       if (gamef > 0.3d0) then 

          c=0.896434d0*gamp*zhat-3.44740d0*gamp**FOURTH*zhat2- &
               0.5551d0*(log(gamp)+FIVE3RD*log(z1*z2/(z1+z2)))-2.996d0

          dcdt=0.896434d0*dgampdt*zhat- &
               3.44740d0*FOURTH*gamp**(FOURTH-1.0d0)*zhat2*dgampdt- &
               0.5551d0*dgampdt/gamp

          h12=c-(tau12/3.d0)*(5.d0*alph12**3/32.d0-0.014d0*alph12**4 &
               -0.0128d0*alph12**5)-gamef*(0.0055d0*alph12**4 &
               -0.0098d0*alph12**5+0.0048d0*alph12**6)

          dh12dt=dcdt - ((dtau12dt*alph12**3 + 3.0d0*tau12*alph12**2* &
               dalph12dt)*(5.d0/32.d0 - 0.014d0*alph12 - &
               0.0128d0*alph12**2) + tau12*alph12**3*dalph12dt*(-0.014d0 &
               - 2.d0*0.0128d0*alph12))/3.d0 -(dgamefdt*alph12**4 + 4.d0 &
               *gamef*alph12**3*dalph12dt)*(0.0055d0 - 0.0098d0*alph12 - &
               0.0048d0*alph12**2) - gamef*alph12**4*dalph12dt*(-0.0098d0 &
               + 2.d0*0.0048d0*alph12)

          h12fac=0.77d0

          h12=log(max(1.d+0-0.0562d+0*alph12**3,h12fac))+h12
          if (1.d+0-0.0562d+0*alph12**3 .gt. h12fac) then
             dh12dt=(-3.d0*0.0562d0*alph12**2*dalph12dt)/ &
                  (1.d0-0.0562d0*alph12**3) + dh12dt
          endif

          if(gamef <= 0.8d0) then

             h12=h12w*((0.8d0-gamef)/0.5d0)+h12*((gamef-0.3d0)/0.5d0)
             dh12dt=((dh12wdt*(0.8d0-gamef) - h12w*dgamefdt + dh12dt* &
                  (gamef-0.3d0) + h12*dgamefdt)/0.5d0)
          endif
       endif

       if (h12.gt.300.d0) then
          h12=300.d0
          dh12dt=0.d0
       endif

       if (h12.lt.0.d0) then
          h12=0.d0
          dh12dt=0.d0
       endif

       scfac=exp(h12)
       dscfacdt=scfac*dh12dt

    else
       scfac=1.d0
       dscfacdt=0.d0
    endif

    return

  end subroutine screenz
  

end module screening_module
