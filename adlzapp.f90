  
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

    !     Degrees of freedom is defined here (qdim is defined to accomadate for xq containing v's):
    integer, parameter :: dof = 1
    integer, parameter :: qdim = 2*dof
    integer, parameter :: d = 1 !Spatial degrees of freedom
    integer, parameter :: nop = dof/d !Determine number of particles from given information
    

    !     Coordinates (q:1-6, p:7-12) and forces/velocities ("right sides" defined as fv)
    !The first qdim/2 elements of xq are x, second are v_x, third are p, fourth are v_p
    real (kind = mk), dimension (1:2*qdim) :: xq, fv
    ! Time and time step
    real (kind = mk) :: xtime, xdt

    ! Create array of masses (each element corresponds to each particle)
    real (kind = mk), dimension(1:nop) :: mi = 1.0_mk
        
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
    ! and pots2 below. The precision of the calculation (mk) and the
    ! dimension of the problem (qdim) are set in module coordinates.
    !
    ! Valid in fixed and free form Fortran 90
    !
    use coordinates,  only : mk, qdim, dof, d, xq, fv, xtime, xdt

    implicit none

    integer :: vals(1:8)
    INTEGER(KIND=4) :: nn, ncall
    real(kind=mk) :: ttotal

    real(kind=mk), parameter ::                                       &
         & pi = 3.1415926535897932384626433832795029_mk
    real(kind = mk) :: qf, pf, phif, tf, ef

    external pots1, pots2

    ! Input of stepsize:
    write (*,*) "Total number of steps: "
    read  (*,*) ncall
    write (*,*) "Total time:"
    read  (*,*) ttotal
    xdt = ttotal/real(ncall,kind=mk)
    xtime = 0.0_mk
    write(6,*) xtime, xdt

    ! File for reporting results:
    open(9, file="testsymp.res", position="append")
    call date_and_time(values=vals)
    write (9,*) vals(1:3), "  ", vals(5:8)
    write (9,*) "Dimension: ", qdim, "    Calls per cycle ", ncall

    ! Initial conditions:
    xq(1:dof) =1.0_mk !x's
    xq(dof+1:qdim) =0.0_mk !v_x's
    xq(qdim+1:qdim+dof) = 0.0_mk !p's
    xq(qdim+dof+1:2*qdim) =0.0_mk !v_p's

    ! Redefinition of time step to get ncall steps for one cycle:
    !xdt = xdt*2.0_mk*pi !use only if you just want one period
    write(102,*) xtime, xq(1:dof), xq(qdim+1:qdim+dof) !write initial time, x, p
    ! Integration loop:
    call algini                      ! First call
    intloop: do nn = 2, ncall
       call algrun ! Other calls
       write(102,*) xtime !write all other times
     enddo intloop

    ! End check (only first coordinate pair) :
    tf = ncall*xdt
    qf = xq(1)
    pf = xq(qdim+1)
    phif = atan2(pf,qf)
    ef = qf*qf+pf*pf
    write (*,*) "tf-2pi:    ", tf - 2.0_mk*pi
    write (*,*) "Ef-E0:     ", ef - 1.0_mk
    write (*,*) "qf-1:      ", qf - 1.0_mk
    write (*,*) "pf:        ", pf
    write (*,*) "phif(deg): ", phif*180.0_mk/pi

    write (9,*) "tf-2pi:    ", tf - 2.0_mk*pi
    write (9,*) "Ef-E0:     ", ef - 1.0_mk
    write (9,*) "qf-1:      ", qf - 1.0_mk
    write (9,*) "pf:        ", pf
    write (9,*) "phif(deg): ", phif*180.0_mk/pi
    call date_and_time(values=vals)
    write (9,*) vals(1:3), "  ", vals(5:8) 
    write (9,*)

    CLOSE (9)

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
    use coordinates,  only : mk, qdim, dof, d, xq, fv, xtime, xdt

    implicit none

    integer :: nn, ii

    real (kind = mk), dimension (0:2*nc) :: hhc
    real (kind = mk), dimension (1:2*qdim) :: xqh, fvh

    external pots1, pots2

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
          call pots2     ! -dV/dq
          fvh((qdim+1):2*qdim)=fvh((qdim+1):2*qdim)+fv(1:qdim)*hhc(ii+1)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + fvh( qdim+1:2*qdim)
       enddo
       call pots1      ! dT/dp
       fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(2*nc)

    else  choice                    ! Choose algorithm: typ = 2

       do ii = 0, 2*nc-2, 2
          call pots2  ! -dV/dq
          fvh( qdim+1:2*qdim) = fvh( qdim+1:2*qdim) + fv(1:qdim)*hhc(ii)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + fvh( qdim+1:2*qdim)
          call pots1   ! dT/dp
          fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(ii+1)
          xq(1:qdim) = xqh(1:qdim) + fvh(1:qdim)
       enddo
       call pots2  ! -dV/dq
       fvh(qdim+1:2*qdim) = fvh(qdim+1:2*qdim) + fv(1:qdim)*hhc(2*nc)

    endif choice

    xq(1:2*qdim) = xqh(1:2*qdim) + fvh(1:2*qdim)
    xtime = xtime + xdt
    ! Write out the trajectory
    write(102,*) xq(1:dof), xq(qdim+1:qdim+dof) 

    return
  end subroutine algo12

  !----------------------------------------------------------------------

  ! Subroutines for the "right sides" of the harmonic oscillator with
  ! unit frequency:

  subroutine pots1        ! dH/dp and dK/dvp
    use coordinates, only : mk, dof, nop, d, qdim, xq, fv, mi, d
    implicit none
    integer i
    real (kind = mk), dimension(1:dof) :: m = 0.0_mk
    do i = 1,nop
        m((i-1)*d+1:i*d) = mi(i) 
    enddo
    !Define the xdots = dH/dp = p/m
    fv(qdim+1:qdim+dof) = xq(qdim+1:qdim+dof)/m
    !Define the v_xdots = dK/dvp = v_p/m
    fv(qdim+dof+1:2*qdim) = xq(qdim+dof+1:2*qdim)/m    
    return
  end subroutine pots1

  subroutine pots2       ! -dH/dx and -dK/dvx
    use coordinates,  only : mk, dof, nop, d, qdim, xq, fv, mi, d
    implicit none
    integer i
    real (kind = mk), dimension(1:dof) :: m = 0.0_mk
    do i = 1,nop
        m((i-1)*d+1:i*d) = mi(i)
    enddo
    fv(1:dof) = -xq(1:dof)*m !Define the pdots = -dH/dx = -x*m
    fv(dof+1:qdim) = -xq(dof+1:qdim)*m !Define the vpdots = -dK/dvx = -v_x*d^2H/dx^2 = -v_x*m
    return
  end subroutine pots2
!  ****
