  
!  ************** testsymp.for **************************************************
  ! Program Suite testsymp (fixed form Fortran)
  !
  ! Simple driver and test program for symplectic integration routines
  ! of E. Teloy. See paper of A. Seiter and Ch. Schlier
  ! "High-Order Symplectic Integration: An Assessment"
  ! Computer Physics Communications, accepted
  !
  ! Written in standard Fortran 90 (fixed form). To change into free form
  ! module coefficients must be changed as follows: blancs in the numbers
  ! must be deleted, and the continuation codes (&) moved.
  !
  ! Select double or quadruple precision (mk = 8 or 16) in module
  ! coordinates. These kinds may be processor or compiler dependent.
  !
  ! Set dimension of problem (qdim) also in module coordinates.
  !
  ! One out of the documentd five integrators can be selected by
  ! uncommenting them in subroutine algo12.
  !
  ! The main program is placed below, since some compilers want to compile
  ! the modules first!
  !
  ! ---------------------------------------------------------------------
  !
  module coordinates

    !     Defines: mk, qdim, xq, fv, xtime, xdt used by the other program units

    implicit none

    !     Select double or quadruple precision here:
    !     integer, parameter :: mk = 16
    integer, parameter :: mk = 8

    !     Dimension of problem is defined here:
    integer, parameter :: d = 1
    integer, parameter :: nop = 2
    integer, parameter :: dof = nop*d
    integer, parameter :: qdim = dof*2

    !     Coordinates (q:1-6, p:7-12) and forces ("right sides")
    real (kind = mk), dimension (1:2*qdim) :: xq
    real (kind = mk), dimension (1:2*qdim) :: fv = 0.0_mk
    ! Time and time step
    real (kind = mk) :: xtime, xdt

    ! Masses
    real (kind = mk), dimension (1:nop) :: mi = 1.0_mk
    real (kind = mk), dimension (1:dof) :: mass = 0.0_mk

    !Fast Lyapunov Indicator
    real (kind = mk) :: FLI 

    ! Constants for potentials
    real (kind = mk) :: De = 1.0_mk
    real (kind = mk) :: a = 1.0_mk
    real (kind = mk) :: re = 1.0_mk
    real (kind = mk) :: kMeyer = 1.0_mk !Meyer potential mixing coeffecient
    real (kind = mk) :: kext = 1.0_mk !Spring coeffecient for external harm osc trap

    !Define Initial conditions
    real(kind = mk), dimension(1:dof) :: qi = 0.0_mk
    real(kind = mk), dimension(1:dof) :: vqi = 0.0_mk
    real(kind = mk), dimension(1:dof) :: momentumi = 0.0_mk
    real(kind = mk), dimension(1:dof) :: vmomentumi = 0.0_mk

    contains
            subroutine masses
                integer i
                do i =1,nop
                        mass((i-1)*d+1:i*d) = mi(i)
                enddo
            end subroutine masses

  end module coordinates

  !----------------------------------------------------------------------

  module coefficients

    ! Coefficients for five symplectic integration routines of E. Teloy
    ! See paper of A. Seiter and Ch. Schlier
    ! "High-Order Symplectic Integration: An Assessment"
    ! Computer Physics Communications, to be published
    !
    ! Listed is only the first half of the coefficients, ai(0) to ai(nc),
    ! where 2*nc is the number of substages of the integrator.
    ! The full set is computed on first start in subroutine algo12.
    !
    ! This module is correct only in fixed form Fortran 90!
    ! To change it to free form the blanks in the coefficients must be
    ! deleted, and the "&"s must be moved from before the continuation lines
    ! to the end of the lines which are to be continued.

    use coordinates, only: mk

    implicit none

    integer, parameter ::     &
     typ6a = 1, nc6a = 9,     &
     typ6b = 2, nc6b = 8,     &
     typ8a = 1, nc8a = 17,    &
     typ8b = 2, nc8b = 17,    &
     typ8c = 2, nc8c = 17

    real (kind = mk), dimension(0:9), parameter ::    &
    ai6a = (/                                         &     
      0.09517625454177405267746114335519342_mk, &
      0.66629689399770780134207498907168068_mk, &
     -0.12795028552368677941219191621429411_mk, &
      0.02461890095210508713078430308713062_mk, &
      0.10597295345325113143793587608716998_mk, &
     -0.41072553361795113231992873918199025_mk, &
      0.44822227660082748416851634186561201_mk, &
      0.65772926205091317768935130009339042_mk, &
     -0.02142119907216588887172144509368130_mk, &
     -0.87583904676554986768456370614042295_mk /)

    real (kind = mk), dimension(0:8), parameter ::    &
     ai6b = (/                                        &
      0.06942944346252987735848865824703402_mk, &
      0.28487837717280084052745346456657828_mk, &
     -0.13315519831598209409961309951373512_mk, &
      0.32783975759612945412054678367325547_mk, &
      0.00129038917981078974230481746443284_mk, &
     -0.38122104271932629475622784374211274_mk, &
      0.42243536567364142699881962380226825_mk, &
      0.26850290795039600010822759550227899_mk, &
      0.28000000000000000000000000000000000_mk /)

    real (kind = mk), parameter ::                     &
     ai8a(0:17) = (/                                   &
      0.04020757626295627296653921454892367_mk,  &
      0.10968252140081995880852111452131455_mk,  &
      0.17023759564885894706257453906663563_mk,  &
      0.36756158806337006433149757369026277_mk,  &
      0.24370233998503432353195633486895307_mk,  &
     -0.04544131419758065661437375963088864_mk,  &
      0.56601963795366046019899599701939548_mk,  &
      0.00022167162169864039643822185570309_mk,  &
     -0.58169695762497039518529999797620005_mk,  &
      0.05519927098092328759679762829526377_mk,  &
     -0.24138639830477987453171482029238617_mk,  &
     -0.12513929981618023524050370745321727_mk,  &
      0.36115097569793127373014000321599616_mk,  &
     -0.04284389352937610255914308734324331_mk,  &
     -0.53225450460377165284025446933453953_mk,  &
     -0.00393367299329157410510456094858013_mk,  &
      0.47401973498508064506706319888322175_mk,  &
      0.36938625693923323477174115402677030_mk /)

    real (kind = mk), parameter ::                     &
     ai8b(0:17) = (/                                   &
      0.03676680389912337302666154929429291_mk,  &
      0.11072655003739784175754797312279745_mk,  &
      0.16040429374255560219395381214509780_mk,  &
      0.61101267825171523627962718607785428_mk,  &
     -0.00472877643941287918639412436088645_mk,  &
     -0.19202809069032535396838334049379558_mk,  &
      0.02983098489335056954884440558763334_mk,  &
     -0.25979073929811660257162833544861286_mk,  &
      0.19135844311091097984885756175207225_mk,  &
      0.38384564066882093754274499421236298_mk,  &
     -0.03781968145745128677723635761417376_mk,  &
      0.32661664886778120135972921761872954_mk,  &
      0.00351845996378093605518443870229385_mk,  &
     -0.53463443374897025678663398242742174_mk,  &
      0.13067013867271618676514580608303276_mk,  &
     -0.39935632081078281354806842349635698_mk,  &
     -0.01000066638557348147501709158936269_mk,  &
      0.90721613344495961987012942166888585_mk /)

    real (kind = mk), parameter ::                     &
     ai8c(0:17) = (/                                   &
      0.04463795052359022755913999625733590_mk,  &
      0.13593258071690959145543264213495574_mk,  &
      0.21988440427147072254445535069606167_mk,  &
      0.13024946780523828601621193778196846_mk,  &
      0.10250365693975069608261241007779814_mk,  &
      0.43234521869358547487983257884877035_mk,  &
     -0.00477482916916881658022489063962934_mk,  &
     -0.58253476904040845493112837930861212_mk,  &
     -0.03886264282111817697737420875189743_mk,  &
      0.31548728537940479698273603797274199_mk,  &
      0.18681583743297155471526153503972746_mk,  &
      0.26500275499062083398346002963079872_mk,  &
     -0.02405084735747361993573587982407554_mk,  &
     -0.45040492499772251180922896712151891_mk,  &
     -0.05897433015592386914575323926766330_mk,  &
     -0.02168476171861335324934388684707580_mk,  &
      0.07282080033590128173761892641234244_mk,  &
      0.55121429634197067334405601381594315_mk /)

  end module coefficients
  !
  ! ---------------------------------------------------------------------
  !
  program testsymp

    ! This driver computes one cycle of a harmonic oscillator and returns
    ! some test data. The "right sides" are computed by subroutines pots1
    ! and potsMeyer2 below. The precision of the calculation (mk) and the
    ! dimension of the problem (qdim) are set in module coordinates.
    !
    ! Valid in fixed and free form Fortran 90
    !
    use coordinates

    implicit none

    integer :: vals(1:8)
    INTEGER(KIND=4) :: nn, ncall
    real(kind = mk) :: oscs = 0.0_mk
    integer :: i,j,m,k
    
    real(kind=mk), parameter ::                                       &
         & pi = 3.1415926535897932384626433832795029_mk
    real(kind = mk) :: qf, pf, phif, tf, ef, ei, normv, normv0, T, Ti
    real(kind = mk) :: Vtot, Vtoti = 0.0_mk

    external pots1, potsMeyer2
    integer, external :: vectpos
    real (kind = mk), external :: r, V, dVdr, d2Vdr2
    call masses

    ! Input of stepsize:
    write (*,*) "Total number of steps: "
    read  (*,*) ncall
    write (*,*) "Step size: "
    read (*,*) xdt

    !Redefine initial conditions
    qi(1) = -1.0_mk !x(1,1)
    qi(2) = 1.0_mk !x(1,2)
    vqi(1) = 0.1_mk !v_x
    vqi(2) = 0.0_mk !v_y
    momentumi(1) = 0.0_mk !p(1,1)
    momentumi(2) = 0.5_mk !p(1,2)
    vmomentumi(1) = -0.001_mk !v_px
    vmomentumi(2) = 0.001_mk !v_py

    ! Assign initial conditions:
    xq(1:dof) = qi !q's
    xq(dof+1:qdim) = vqi !v_q's
    xq(qdim+1:qdim+dof) = momentumi !p's
    xq(qdim+dof+1:2*qdim) = vmomentumi !v_p's

    !data file for results
    open(1, file = 'Morse.dat')
    write(1,*) xtime, qi(1), qi(2)

    ! Calculate initial energy
    Ti = sum(momentumi**2/(2*mass))
    do j = 1,nop
        do i = j+1,nop
            Vtoti = Vtoti + V(r(i,j))
        enddo
    enddo

    ! Integration loop:
    call algini                      ! First call
    intloop: do nn = 2, ncall
       call algrun ! Other calls
     enddo intloop

    ! End check :
    !tf = ncall*xdt
    !qf = xq(1);
    !pf = xq(qdim+1)
    !phif = atan2(pf,qf)
    !Check Meyer energy conservation
    T = sum(xq(qdim+1:qdim+dof)**2/(2*mass))
    do j = 1,nop
        do i = j+1,nop
            Vtot = Vtot + V(r(i,j))
        enddo
    enddo
    ef = T + Vtot
    !write(*,*) "T: ", T, "Vtot: ", Vtot

    !write(*,*) r(1,2)

    !normv0 = sqrt(sum(vqi**2)+sum(vmomentumi**2))
    !write (*,*) "tf-2pi:    ", tf - 2.0_mk*pi
    !write(*,*) "fv: ",fv
    !write (*,*) "Ef:   ", ef, "Ei:    ", ei
    !write (*,*) "Ef-E0:     ", ef - ei
    !write (*,*) "qf-1:      ", qf - 1.0_mk
    !write (*,*) "pf:        ", pf
    !write (*,*) "phif(deg): ", phif*180.0_mk/pi

    !write (9,*) "tf-2pi:    ", tf - 2.0_mk*pi
    !write (9,*) "Ef-E0:     ", ef - 1.0_mk
    !write (9,*) "qf-1:      ", qf - 1.0_mk
    !write (9,*) "pf:        ", pf
    !write (9,*) "phif(deg): ", phif*180.0_mk/pi
    call date_and_time(values=vals)
    !write (9,*) vals(1:3), "  ", vals(5:8) 
    !write (9,*)

    !CLOSE (9)
    close (1)

  end program testsymp

  !----------------------------------------------------------------------

  subroutine algo12

    ! This subroutine performs one full integration step.
    ! The increments fvh are added to fv only after the final substep to
    ! preserve as much accuracy as possisble.
    ! The use of the FSAL property has not been implemented.

    ! Uncomment one of the five following lines to determine which
    ! integrator you want to use:

    !      use coefficients, only : typ =>typ6a, nc => nc6a, ai => ai6a
    !      use coefficients, only : typ =>typ6b, nc => nc6b, ai => ai6b
    !      use coefficients, only : typ =>typ8a, nc => nc8a, ai => ai8a
    !      use coefficients, only : typ =>typ8b, nc => nc8b, ai => ai8b
    use coefficients, only : typ =>typ8c, nc => nc8c, ai => ai8c

    ! Exchange of parameters and coordinates with other program units:
    use coordinates

    implicit none

    integer :: nn, ii
    real (kind = mk) normv, normv0

    real (kind = mk), dimension (0:2*nc) :: hhc
    real (kind = mk), dimension (1:2*qdim) :: xqh, fvh

    external pots1, potsMorse2

    entry algini

    write (*,*) "typ ", typ, " nc: ", nc
    write (9,*) "typ ", typ, " nc: ", nc
    ! Compute full set of coefficients:
    hhc(0:nc) = ai(0:nc)
    hhc(nc+1:2*nc) = ai(nc-1:0:-1)
    ! Check integrity of coefficients (sum must be 2.0):
          write (*,*)"sum of hh ", sum(hhc)
          write (*,*)

    hhc(0:2*nc) = xdt*hhc(0:2*nc)

    entry algrun

    xqh(1:2*qdim) = xq(1:2*qdim)
    fvh(1:2*qdim) = 0.0e0_mk

    choice: if (typ == 1) then      ! Choose algorithm: typ = 1
       do ii = 0, 2*nc-2, 2
          call pots1      ! dT/dp
          fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(ii)
          xq(1:qdim) = xqh(1:qdim) + fvh(1:qdim)
          call potsMorse2     ! -dV/dq
          fvh((qdim+1):2*qdim)=fvh((qdim+1):2*qdim)+fv(1:qdim)*hhc(ii+1)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + fvh( qdim+1:2*qdim)
       enddo
       call pots1      ! dT/dp
       fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(2*nc)

    else  choice                    ! Choose algorithm: typ = 2

       do ii = 0, 2*nc-2, 2
          call potsMorse2  ! -dV/dq
          fvh( qdim+1:2*qdim) = fvh( qdim+1:2*qdim) + fv(1:qdim)*hhc(ii)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + fvh( qdim+1:2*qdim)
          call pots1   ! dT/dp
          fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(ii+1)
          xq(1:qdim) = xqh(1:qdim) + fvh(1:qdim)
       enddo
       call potsMorse2  ! -dV/dq
       fvh(qdim+1:2*qdim) = fvh(qdim+1:2*qdim) + fv(1:qdim)*hhc(2*nc)

    endif choice

    xq(1:2*qdim) = xqh(1:2*qdim) + fvh(1:2*qdim)
    xtime = xtime + xdt

    ! Calculate Lyapunov Indicator
    !normv0 = sqrt(vqi(1)**2+vqi(2)**2+vmomentumi(1)**2+vmomentumi(2)**2)
    !normv = sqrt(xq(3)**2+xq(4)**2+xq(7)**2+xq(8)**2)
    !FLI = log(normv/normv0)/xtime
    write(1,*) xtime, xq(1), xq(2) 
    return
  end subroutine algo12

  !----------------------------------------------------------------------

  ! Subroutine for the p derivatives of a well behaved Hamiltonian
  subroutine pots1        ! qdots and v_qdots
    use coordinates
    call masses
    fv(qdim+1:qdim+dof) = xq(qdim+1:qdim+dof)/mass !qdot = dH/dp = p/m
    fv(qdim+dof+1:2*qdim) = xq(qdim+dof+1:2*qdim)/mass !v_qdot = dK/dv_p = v_p*d2H/dp2 = v_p/m
    return
  end subroutine pots1

  !Subroutine for the Morse potential for n particles
    subroutine potsMorse2
    use coordinates
    implicit none
    real (kind = mk), external :: r, V, dVdr, d2Vdr2, dVextdx, d2Vextdx2
    integer, external :: vectpos
    real (kind = mk), dimension(1:nop) :: jsum, lsum, intsum = 0.0_mk
    integer :: n,j,m,k,o,l !j,k and l are particle numbers, m and o are coordinate numbers and n is  the position in the phase space-vector
    call masses

    ! Calculate all the pdots
    do n = 1,dof
       call indicies(n,m,k)
       do j = 1,nop !consider interaction with each other particle
          intsum(j) = -dVdr(r(k,j))*(xq(vectpos(m,k))-xq(vectpos(m,j)))/r(k,j)
          intsum(k) = 0.0_mk !except the kth particle
       enddo
       fv(n) = sum(intsum) - dVextdx(xq(n)) !sum over all particles (execpt kth)
    enddo

    ! Calculate all the v_pdots
    do n = dof+1,qdim
        call indicies(n-dof,m,k)
        do o = 1,d
           do j = 1, nop
              jsum(j) = d2Vdr2(r(k,j))*(xq(vectpos(m,k))-xq(vectpos(m,j)))*(xq(vectpos(o,k))-xq(vectpos(o,j)))/r(k,j)**2
              jsum(k) = 0.0_mk
           enddo
           do l = 1,nop
              lsum(l) = d2Vdr2(r(l,k))*(xq(vectpos(m,l))-xq(vectpos(m,k)))*(xq(vectpos(o,k))-xq(vectpos(o,l)))/r(l,k)**2
              lsum(k) = 0.0_mk
           enddo
           intsum(o) = -xq(vectpos(o,k)+dof)*sum(jsum)-sum(lsum)
        enddo
        fv(n) = sum(intsum) - d2Vextdx2(xq(n-dof))
    enddo
    return
  end subroutine potsMorse2

  !Subroutine for the Harmonic oscillator potential
  subroutine pots2       ! pdots and v_pdots
    use coordinates
    call masses
    fv(1:dof) = -xq(1:dof)*mass !pdot = -dH/dq = -q*m
    fv(dof+1:qdim) = -xq(dof+1:qdim)*mass !v_pdot = -dK/dv_q = -v_q*d2H/dq2 = -v_q*m
    return
  end subroutine pots2

  !Subroutines for the Meyer Hamiltonian
  subroutine potsMeyer2
    use coordinates     ! pdots and v_pdots 
    fv(1) = -(xq(1)+8*kMeyer*xq(1)*xq(2)**2) !p1dot = -(q1 + 8*k*q1*q2^2)
    fv(2) = -(xq(2)+8*kMeyer*xq(2)*xq(1)**2) !p2dot = -(q2 + 8*k*q2*q1^2)
    fv(3) = -xq(3)*(1+8*kMeyer*xq(2)**2)-xq(4)*16*kMeyer*xq(1)*xq(2) !v_p1dot = -vq1(16*k*q1*q2)-vq2(1+8*kq2^2)
    fv(4) = -xq(4)*(1+8*kMeyer*xq(1)**2)-xq(3)*16*kMeyer*xq(1)*xq(2) !v_p2dot = -vq1(1+8*kq1^2)-vq2(16*k*q1*q2)
  end subroutine potsMeyer2

  ! Function to compute separation distance between the ith and jth particles
  function r(i,j)
    use coordinates
    implicit none
    integer, intent(in) :: i, j
    real (kind = mk) :: r
    r = sqrt(sum((xq((i-1)*d+1:i*d)-xq((j-1)*d+1:j*d))**2))
  end function

  ! Functions to return the values of the Morse Potential and its derivatives
  function V(sep)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: sep
    real (kind = mk) :: V
    V = 0.5_mk*sep**2 
    !V = De*(exp(-2*a*(sep-re))-2*exp(-a*(sep-re)))
  end function

  function dVdr(sep)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: sep
    real (kind = mk) :: dVdr
    dVdr = sep
    !dVdr = De*(-2*a*exp(-2*a*(sep-re))+2*a*exp(-a*(sep-re)))
  end function

  function d2Vdr2(sep)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: sep
    real (kind = mk) :: d2Vdr2
    d2Vdr2 = 1.0_mk
    !d2Vdr2 = De*(4*a**2*exp(-2*a*(sep-re))-2*a**2*exp(-a*(sep-re)))
  end function

  ! Functions to return values of external potential (like a harm osc trap)
  function Vext(x)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: x
    real (kind = mk) :: Vext
    Vext = 0.5_mk*kext*x**2
  end function

  function dVextdx(x)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: x
    real (kind = mk) :: dVextdx
    dVextdx = x*kext
  end function

  function d2Vextdx2(x)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: x
    real (kind = mk) :: d2Vextdx2
    d2Vextdx2 = kext*1.0_mk
  end function

  ! Function which returns the coordinate number and particle number based on the position in the 
  ! phase-space vector
  subroutine indicies(n, coordnum, partnum)
    use coordinates
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: partnum, coordnum
    partnum = n/d + 1
    coordnum = mod(n,d)
    if (coordnum == 0) then
            coordnum = d
            partnum = partnum-1
    end if
  end subroutine indicies

  ! Function which returns the position in the phase-space vector as a function of particle and 
  ! coordinate numbers
  function vectpos(coordnum,partnum)
    use coordinates
    implicit none
    integer, intent(in) :: coordnum, partnum
    integer :: vectpos
    vectpos = d*(partnum-1)+coordnum
  end function 
