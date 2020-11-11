!! This version has no FLI !!  
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
         !integer, parameter :: mk = 16
    
         integer, parameter :: mk = 8

    !     Dimension of problem is defined here:
    integer, parameter :: d = 1
    integer, parameter :: nop = 3
    integer, parameter :: dof = nop*d
    integer, parameter :: qdim = dof

    !     Coordinates (q:1-6, p:7-12) and forces ("right sides")
    real (kind = mk), dimension (1:2*qdim) :: xq ! The first half will be x's and the second p's
    real (kind = mk), dimension (1:2*qdim) :: fv = 0.0_mk
    ! Time and time step
    INTEGER(KIND=4) :: nn, ncall
    real (kind = mk) :: xtime, xdt, start, finish, armtime
    real (kind = mk), dimension(:), allocatable :: rtimestemp, rtimes, coltimes, Vpoints

    ! Conversion Factors
    real(kind=mk), parameter ::                                       &
         & pi = 3.1415926535897932384626433832795029_mk
    real (kind = mk) :: c = 299792458_mk
    real (kind = mk) :: meperamu = 1822.89_mk
    real (kind = mk) :: hartperinvcm = 4.55633525291e-6_mk
    real (kind = mk) :: Hzperinvcm = 2.99792458e6_mk
    real (kind = mk) :: angsperbohr = 0.529177249_mk
    real (kind = mk) :: nanoKperEh = 3.15775024804e14_mk
    real (kind = mk) :: secperatu = 2.41888432658e-17_mk
    real (kind = mk) :: JperEh = 4.3597447222071e-18_mk
    real (kind = mk) :: kgperamu = 1.66054e-27_mk

    ! Paramters of the Potentials
    ! Dissociation Energies in cm^-1 (converted to E_h)
    real (kind = mk) :: DeK2 = 4405.389_mk*4.55633525291e-6_mk
    real (kind = mk) :: DeRb2 = 3965.8_mk*4.55633525291e-6_mk
    real (kind = mk) :: DeKRb = 4180.417_mk*4.55633525291e-6_mk
    ! Equilibrium Distances in Angstroms (converted to bohr)
    real (kind = mk) :: reK2 = 3.956_mk/0.529177249_mk
    real (kind = mk) :: reRb2 = 4.233_mk/0.529177249_mk
    real (kind = mk) :: reKRb = 4.160_mk/0.529177249_mk
    real (kind = mk) :: aK2 = 1.0_mk !in Bohr^(-1)
    real (kind = mk) :: aKRb = 1.0_mk
    real (kind = mk) :: aRb2 = 1.0_mk
    real (kind = mk) :: omegaK2cm = 92.021_mk !in cm^-1
    real (kind = mk) :: omegaKRBcm = 75.5_mk
    real (kind = mk) :: omegaRb2cm = 57.31_mk
    real (kind = mk) :: omegaK2, omegaKRb, omegaRb2

    ! Jacobi Polar Coordinates
    !real (kind = mk), dimension(1:4000000) :: Vpoints
    real (kind = mk) :: cutoffrad = 150.0_mk
    real (kind = mk) :: hyprad = 100.0_mk
    real (kind = mk) :: hypang, hypangint = 0.0_mk
    real (kind = mk) :: hypraddot = -0.001_mk
    real (kind = mk) :: hypangdot = 0.0_mk
    real (kind = mk) :: KEy1 = 0.0_mk !Enter energies in nK
    real (kind = mk) :: KEy2 = 1e10_mk
    real (kind = mk) :: eps = 0.0_mk
    real (kind = mk), dimension(:), allocatable :: radvec, radmax
    logical, dimension(:), allocatable :: rchoose

    ! Masses (in amu)
    real (kind = mk), dimension (1:nop) :: miamu = 1.0_mk
    real (kind = mk), dimension (1:nop) :: mi = 1.0_mk
    real (kind = mk), dimension (1:dof) :: mass = 0.0_mk
    real (kind = mk) :: Mtot, mu, mu12, mu123, mu23, muKRb
    real (kind = mk) :: mPotas = 39.0983_mk
    real (kind = mk) :: mRb = 85.4678_mk

    ! Fast Lyapunov Indicator
    real (kind = mk) :: FLI 

    ! Constants for potentials
    real (kind = mk) :: kMeyer = 1.0_mk !Meyer potential mixing coeffecient
    real (kind = mk) :: omegaext = 0.0_mk !Spring coeffecient for external harm osc trap (in Kelvin/Angs^2)

    ! Define which particle is which
    integer, dimension(1:nop) :: partkind 
    integer :: coltype
    integer :: count12, count13, count23, countind, numcols, eps1, angstep, laststep

    ! Define Initial conditions
    real(kind = mk), dimension(1:dof) :: qi = 0.0_mk
    real(kind = mk), dimension(1:dof) :: vqi = 0.0_mk
    real(kind = mk), dimension(1:dof) :: momentumi = 0.0_mk
    real(kind = mk), dimension(1:dof) :: vmomentumi = 0.0_mk

    ! Define coordinates and energies to be defined later
    real (kind = mk) :: tf, ef, normv, xcm, y1, y2, normv0, y1dot, y2dot, &
         Ptoti, Ptot, E, T, Vtot, theta13, theta12, theta23, ang1, ang2 = 0.0_mk

    ! For RNG
    integer, parameter :: rngsize = 1000
    double precision :: rngvec(rngsize*2)
    real (kind = mk) :: rvec(size(rngvec))
    real (kind = mk), dimension(:), allocatable :: rvecselect

    ! For gaussian
    real (kind = mk) :: agauss = 1.0_mk, bgauss = 0.0_mk, cgauss = 1.0_mk

    contains
            subroutine masses
                integer i
                partkind = (/1, 0, 0/) !1 for Rubidium, 0 for Potassium
                
                ! Assign masses based on particle type
                do i = 1,nop
                  if (partkind(i) == 0) then
                    mi(i) = mPotas*meperamu
                  else
                    mi(i) = mRb*meperamu
                  endif
                  mass((i-1)*d+1:i*d) = mi(i)
                enddo
                Mtot = mi(1) + mi(2) + mi(3)
                mu12 = (mi(1)*mi(2))/(mi(1) + mi(2))
                mu123 = (mi(1) + mi(2))*mi(3)/Mtot
                mu = sqrt(mu12*mu123)
                mu23 = (mi(2) + mi(3))/(mi(2)*mi(3))
                muKRb = mPotas*mRb/(mPotas + mRb) 
    
                theta13 = atan(-mu/mi(1))
                theta23 = atan(mu/mi(2))

                !omega's converted from cm^-1 to Hz to 1/atu
                !omegaK2 = 2.0*pi*100.0*c*omegaK2cm*secperatu
                !omegaKRb = 2.0*pi*100.0*c*omegaKRbcm*secperatu
                !omegaRb2 = 2.0*pi*100.0*c*omegaRb2cm*secperatu
                omegaK2 = omegaK2cm*Hzperinvcm*secperatu
                omegaKRb = omegaKRbcm*Hzperinvcm*secperatu
                omegaRb2 = omegaRb2cm*Hzperinvcm*secperatu

                !a's calculated and converted to Bohr^(-1)
                aK2 = 0.3639!omegaK2*sqrt((mPotas*meperamu/2.0_mk)/(2.0_mk*DeK2))
                aKRb = 0.3898!omegaKRb*sqrt((muKRb*meperamu)/(2.0_mk*DeKRb))
                aRb2 = omegaRb2*sqrt((mRb*meperamu/2.0_mk)/(2.0_mk*DeRb2))
 
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
    integer :: i,j,m,k,l,jj
    
    real(kind = mk) :: Vtoti, Ei, Ti = 0.0_mk

    external pots1, potsMeyer2
    integer, external :: vectpos
    real (kind = mk), external :: r, V, dVdr, d2Vdr2, Vext, dVextdx, d2Vextdx2, De, re, apar
    call CPU_time(start)
    call masses

    ! Input of stepsize:
    !write (*,*) "Total time length: "
    !read  (*,*) tf
    !write (*,*) "Step size: "
    !read (*,*) xdt
    tf = 400000000_mk
    !do jj = 0,2
    xdt = 100!**real(jj,mk)
    write(*,*)"xdt:",  xdt
    ncall = int(tf/xdt)

    write(*,*) "weK2", omegaK2,"weRb2", omegaRb2, "weKRB", omegaKRb
    write(*,*) "aKRb: ", aKRb, "aK2: ", aK2, "aRb2", aRb2
    !write(*,*) theta13, theta23
    ! Create array of ran nums for hypang int condit
    !call rgnf_lux(rngvec,rngsize*2)
    !rvec = real(rngvec,mk)*0.001_mk 
    !do i = 1, size(rvec)
    !eps = rvec(i)
    write(*,*) "The final time: ", tf*secperatu*1e9_mk
    ang1 = 1.4858642565_mk
    ang2 = 1.5015699085_mk
    !ang1 = 0.5566701165_mk
    !ang2 = 0.5697164945_mk
    numcols = 1!00
    call rgnf_lux(rngvec,numcols)
    rngvec = rngvec*(ang2-ang1)+ang1
    !allocate(Vpoints(ncall)) 
    !numcols = 5697-5567+1
    angstep = 15705652
    count12 = 0
    do j = 0, 0 ! Loop for order of magnitude of energy
    k = 1
    write(*,*) "j: ", j
    do m = 1,1
    count13 = 0 
    count23 = 0
    countind = 0
    allocate(coltimes(numcols))
    do i = 1, numcols
    !Vpoints = 0.0_mk
    hypang = rngvec(i)
    !hypang = real(i,mk)*1e-4_mk
    write(*,*) "hypang: ", hypang
    !eps  = real(i,mk)/2000_mk
    !eps = 0.0055_mk
    armtime = 0.0_mk
    xtime = 0.0_mk
    hyprad = 100.0_mk
    KEy2 = 1e6_mk!(1.0_mk + real(j,mk)*10.0_mk)*1e7_mk
   ! KEy2 = (10**(real(j,mk)))/real(m,mk)
   ! allocate (radvec(ncall),rtimestemp(ncall),rchoose(ncall))
   ! radvec(1) = hyprad

    ! Redefine hyperangle based on hyperradius
    !hypang = atan(mu/mi(1)) + re(2,3)/hyprad !particles 2,3
    !hypang = pi/2.0_mk - re(1,2)/hyprad + eps !particles 1,2
    hypangint = hypang

    ! ---------------------------- Set initial conditions here -------------------------------------
    ! (in Bohr)
    ! Initial conditions determined from hyper-radial coordinates
    qi(1) = sqrt(mu)*hyprad*(sqrt(mu12)*cos(hypang)/mi(1) + sqrt(mu123)*sin(hypang)/(mi(1)+mi(2)))
    qi(2) = sqrt(mu)*hyprad*(-sqrt(mu12)*cos(hypang)/mi(2) + sqrt(mu123)*sin(hypang)/(mi(1)+mi(2)))
    qi(3) = -sqrt(mu*mu123)*hyprad*sin(hypang)/mi(3)
    xq(1:qdim) = qi !q's

    y1dot = sqrt(2*KEy1/nanoKperEh/mu)
    y2dot = -sqrt(2*KEy2/nanoKperEh/mu)
    momentumi(1) = mi(1)*sqrt(mu)*(sqrt(mu12)*y1dot/mi(1) + sqrt(mu123)*y2dot/(mi(1)+mi(2)))
    momentumi(2) = mi(2)*sqrt(mu)*(-sqrt(mu12)*y1dot/mi(2) + sqrt(mu123)*y2dot/(mi(1)+mi(2)))
    momentumi(3) = -sqrt(mu*mu123)*y2dot
    !Ptoti = momentumi(1)+momentumi(2)+momentumi(3)
    !write(*,*) 'Ptot: ', Ptoti
    !momentumi(1) = momentumi(1) - (mi(1)/Mtot)*Ptot
    !momentumi(2) = momentumi(2) - (mi(2)/Mtot)*Ptot
    !momentumi(3) = momentumi(3) - (mi(3)/Mtot)*Ptot
    !Ptot = momentumi(1)+momentumi(2)+momentumi(3)
    !write(*,*) 'Ptot new: ', Ptot
    xq(qdim+1:2*qdim) = momentumi !p's
    call energy
    ei = E
    Ti = T
    Vtoti = Vtot
    write(102,*) abs(xq(1)-xq(2)), Vtoti
    write(*,*) 'Ei1: ', Ei, 'Ti1: ', Ti, 'Vtoti1: ', Vtoti

!    call fstfwrd
    ! (p's in m_e*Bohr/atm. time units)
!    y1dot = sqrt(2*mu*(E12-V(abs(qi(1)-qi(2)), De(1,2), re(1,2), re(1,2))))/mu
!    y2dot = sqrt(2*mu*(E123-V(abs(qi(1)-qi(3)), De(1,3), re(1,3), re(1,3))- &
!            V(abs(qi(2)-qi(3)), De(2,3), re(2,3), re(2,3))))/mu
    !y1dot = sqrt(2*KEy1/nanoKperEh/mu)
    !y2dot = -sqrt(2*KEy2/nanoKperEh/mu)
    !momentumi(1) = mi(1)*sqrt(mu)*(sqrt(mu12)*y1dot/mi(1) + sqrt(mu123)*y2dot/(mi(1)+mi(2)))
    !momentumi(2) = mi(2)*sqrt(mu)*(-sqrt(mu12)*y1dot/mi(2) + sqrt(mu123)*y2dot/(mi(1)+mi(2)))
    !momentumi(3) = -sqrt(mu*mu123)*y2dot
    !momentumi = (/-1000.0_mk, 1000.0_mk, 1000.0_mk/)

    ! Define p's as functions of the hyper-coordinates
    !momentumi(1) = mi(1)*sqrt(mu)*hyprad*hypangdot*(-sin(hypang)*sqrt(mu12)/mi(1) + & 
    !        sqrt(mu123)*cos(hypang)/(mi(1) + mi(2))) + mi(1)*hypraddot/hyprad*qi(1)
    !momentumi(2) = mi(2)*sqrt(mu)*hyprad*hypangdot*(sin(hypang)*sqrt(mu12)/mi(1) + & 
    !        sqrt(mu123)*cos(hypang)/(mi(1) + mi(2))) + mi(2)*hypraddot/hyprad*qi(2)
    !momentumi(3) = -mi(3)*sqrt(mu)*hyprad*sqrt(mu123)*cos(hypang)*hypangdot/mi(3) - &
    !        mi(3)*hypraddot/hyprad*qi(3)

    ! Move to Xcm frame
    !Ptoti = momentumi(1)+momentumi(2)+momentumi(3)
    !write(*,*) 'Ptot: ', Ptoti
    !momentumi(1) = momentumi(1) - (mi(1)/Mtot)*Ptot
    !momentumi(2) = momentumi(2) - (mi(2)/Mtot)*Ptot
    !momentumi(3) = momentumi(3) - (mi(3)/Mtot)*Ptot
    !Ptot = momentumi(1)+momentumi(2)+momentumi(3)
    !write(*,*) 'Ptot new: ', Ptot

    ! Load initial conditions:
    !xq(qdim+1:2*qdim) = momentumi !p's

    ! Define Jacobi Coordinates
    call Jacobi

    ! Files for results
    open(1, file = 'traj.dat')
    !write(1,*) xtime, qi
    open(2, file = 'energy.dat')
    open(3, file = 'Jacobi.dat')
    open(4, file = 'radmax.dat')
    write(3,*) y1, y2, 0
!    write(*,*) 'Xcm: ', xcm

    ! Calculate initial energy
    !call energy
    !ei = E
    !Ti = T
    !Vtoti = Vtot
!    write(*,*) "Ti: ", T, "Vtoti: ", Vtot, 'Ei: ', E
!    write(*,*) 'Ei1 - Ei: ', Ei - E

    ! Integration loop:
    call algini                      ! First call
    intloop: do nn = 2, ncall
       call algrun ! Other calls
       if ((sqrt(y1**2 + y2**2) > cutoffrad)) then
               exit
       !else if (sqrt(y1**2 + y2**2) > 75) then
               
       endif
   !    radvec(nn) = sqrt(y1**2 + y2**2)
     enddo intloop

    ! Energy conservation check :
!    call energy
!    write (*,*) "Ef:   ", E, "Ei:    ", ei
!    write (*,*) "Ef-Ei:     ", e - ei

    ! Calculate final hyper-radial coordinates
    hypang = atan(y2/y1)
    hyprad = sqrt(y1**2 + y2**2)

    ! Determine exit channel
!    if (abs(y1) < 30.0_mk) then 
!            write(*,*) '1-2 dimer was formed.'
!            coltype = 1
!             count12 = count12 + 1
!    else if ( (hypang<(theta13+15.0/hyprad)).AND.(hypang>(theta13-15.0_mk/hyprad)) ) then
!            write(*,*) '1-3 dimer was formed.'
!            coltype = 2
!             count13 = count13 + 1
!    else if ( (hypang<(theta23+15.0/hyprad)).AND.(hypang>(theta23-15.0_mk/hyprad)) ) then
!            write(*,*) '2-3 dimer was formed.'
!            coltype = 3
!             count23 = count23 + 1
!    else
!            write(*,*) 'Result was indeterminate.'
!            coltype = 4
!             countind = countind + 1
!    endif

    if (hyprad < 150_mk) then
            count12 = count12 + 1
    endif
    !call rcheck

    !do j = 1,size(radvec)
    !   write(4,*) rtimestemp(j), radvec(j)
    !enddo
    !write(4,*) 
    !deallocate(radvec,rtimes,radmax,rtimestemp,rchoose)

    !write(1,*) eps, xtime*secperatu*1e9_mk - 0.0274_mk, coltype
    coltimes(i) = xtime*secperatu*1e9_mk
    !k = k + 1

    !do l = 1,ncall
    !   if ( (Vpoints(ncall-l)>Vpoints(ncall-l-1)) .AND. (Vpoints(ncall-l)>Vpoints(ncall-l+1)) ) then
    !      write(102,*) Vpoints(ncall-l), real(ncall-l,mk)*real(xdt,mk)*secperatu*1e9_mk
    !      exit
    !   endif
    !enddo
    write(2,*) hypangint, Vtoti, xtime*secperatu*1e9_mk
    call energy
    write(*,*) 'coltimes(i):  ', coltimes(i)
    write(*,*) "Ei - Ef: ", ei - E 
    enddo
!    write(*,*) count12, count13, count23, countind
!    write(1,*) KEy2, real(count12,mk)/101.0_mk, real(count13,mk)/101.0_mk, &
!          &  real(count23,mk)/101.0_mk, real(countind,mk)/101.0_mk
    write(*,*) "The collision energy in K is: ", KEy2*1e-9_mk
    write(*,*) "The average collision time is: ", sum(coltimes)/real(size(coltimes),mk)
    !write(*,*) count12, count13, count23, countind
    deallocate(coltimes)
    !enddo
    enddo
    enddo
    !call potmap
    !do i = 1,ncall
    !   if ( (Vpoints(ncall-i)>Vpoints(ncall-i-1)) .AND. (Vpoints(ncall-i)>Vpoints(ncall-i+1)) ) then
    !      write(102,*) Vpoints(ncall-i), real(ncall-i,mk)*real(xdt,mk)*secperatu*1e9_mk
    !      exit
    !   endif
    !enddo
!    write(*,*) maxval(Vpoints),maxloc(Vpoints)

    ! Create plots of the potentials
    do k = 0,1000
       !write(102,*) real(k,mk)/10_mk, V(real(k,mk)/10_mk,DeK2,reK2,aK2)
       write(103,*) real(k,mk)/10_mk, V(real(k,mk)/10_mk,DeKRb,reKRb,aKRb)
    enddo

    write(*,*) "Number of traj's that maxed out: ", count12
    close (1)
    close (2)
    close (3)
    close (4)
    call CPU_time(finish)
    write(*,*) 'I took ', finish-start, 'sec'
  end program testsymp

  !----------------------------------------------------------------------

  subroutine algo12

    ! This subroutine performs one full integration step.
    ! The increments fvh are added to fv only after the final substep to
    ! preserve as much accuracy as possisble.
    ! The use of the FSAL property has not been implemented.
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

    integer :: ii, i, j
    real (kind = mk), external :: r, V, dVdr, d2Vdr2, Vext, dVextdx, d2Vextdx2, De, re, apar

    real (kind = mk), dimension (0:2*nc) :: hhc
    real (kind = mk), dimension (1:2*qdim) :: xqh, fvh

    external pots1 ,potsMeyer2, potsnbody2, pots2

    entry algini

!    write (*,*) "typ ", typ, " nc: ", nc
!    write (9,*) "typ ", typ, " nc: ", nc
    ! Compute full set of coefficients:
    hhc(0:nc) = ai(0:nc)
    hhc(nc+1:2*nc) = ai(nc-1:0:-1)
    ! Check integrity of coefficients (sum must be 2.0):
!          write (*,*)"sum of hh ", sum(hhc)
!          write (*,*)

    hhc(0:2*nc) = xdt*hhc(0:2*nc)

    entry algrun

    xqh(1:2*qdim) = xq(1:2*qdim)
    fvh(1:2*qdim) = 0.0e0_mk

    choice: if (typ == 1) then      ! Choose algorithm: typ = 1
       do ii = 0, 2*nc-2, 2
          call pots1      ! dT/dp
          fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(ii)
          xq(1:qdim) = xqh(1:qdim) + fvh(1:qdim)
          call potsnbody2     ! -dV/dq
          fvh((qdim+1):2*qdim)=fvh((qdim+1):2*qdim)+fv(1:qdim)*hhc(ii+1)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + fvh( qdim+1:2*qdim)
       enddo
       call pots1      ! dT/dp
       fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(2*nc)

    else  choice                    ! Choose algorithm: typ = 2

       do ii = 0, 2*nc-2, 2
          call potsnbody2  ! -dV/dq
          fvh( qdim+1:2*qdim) = fvh( qdim+1:2*qdim) + fv(1:qdim)*hhc(ii)
          xq( qdim+1:2*qdim) = xqh( qdim+1:2*qdim) + fvh( qdim+1:2*qdim)
          call pots1   ! dT/dp
          fvh(1:qdim) = fvh(1:qdim) + fv(qdim+1:2*qdim)*hhc(ii+1)
          xq(1:qdim) = xqh(1:qdim) + fvh(1:qdim)
       enddo
       call potsnbody2  ! -dV/dq
       fvh(qdim+1:2*qdim) = fvh(qdim+1:2*qdim) + fv(1:qdim)*hhc(2*nc)

    endif choice

    xq(1:2*qdim) = xqh(1:2*qdim) + fvh(1:2*qdim)
    xtime = xtime + xdt

    ! Caclulate Jacobi Coordinates
    call Jacobi
    call Energy
    !Vpoints(nn) = Vtot

    ! Output results
!    write(101,*) xtime*secperatu*1e9_mk, xq(1:dof)
    !write(101,*) xq(1), 0, 0, 'i'
    !write(102,*) xq(2), 0, 0, 'i'
    !write(103,*) xq(3), 0, 0, 'i'
 
    !write(1,*) xtime*secperatu*1e9_mk, sqrt(y1**2 + y2**2)
    !--------------!!!!! Turn this one on for using pot.gnu !!!!!------------------
!    write(3,*) y1, y2, 0
    !write(102,*) xtime*secperatu*1e9_mk, Vpoints(nn)
    !write(4,*) xtime, y2, y1*cos(theta23)+y2*sin(theta23), sqrt(y1**2 + y2**2)
    return
  end subroutine algo12

  !----------------------------------------------------------------------

  ! Subroutine for the p derivatives of a well behaved Hamiltonian
  subroutine pots1        ! qdots and v_qdots
    use coordinates
    call masses
    fv(qdim+1:2*qdim) = xq(qdim+1:2*qdim)/mass !qdot = dH/dp = p/m
    !fv(qdim+dof+1:2*qdim) = xq(qdim+dof+1:2*qdim)/mass !v_qdot = dK/dv_p = v_p*d2H/dp2 = v_p/m
    return
  end subroutine pots1

  !Subroutine for arbitrary pairwise interaction for n particles
    subroutine potsnbody2
    use coordinates
    implicit none
    real (kind = mk), external :: r, V, dVdr, d2Vdr2, dVextdx, d2Vextdx2, De, re, apar
    integer, external :: vectpos
    real (kind = mk), dimension(1:nop) :: jsum, lsum, intsum = 0.0_mk
    integer :: n,j,m,k,o,l !j,k and l are particle numbers, m and o are coordinate numbers and n is  the position in the phase space-vector
    call masses

    ! Calculate all the pdots
    do n = 1,qdim
       call indicies(n,m,k)
       do j = 1,nop !consider interaction with each other particle
          intsum(j) = -dVdr(r(k,j),De(k,j),re(k,j),apar(k,j))*(xq(vectpos(m,k))-xq(vectpos(m,j)))/r(k,j)
          intsum(k) = 0.0_mk !except the kth particle
       enddo
       fv(n) = sum(intsum) - dVextdx(n) !sum over all particles (execpt kth)
    enddo

    ! Calculate all the v_pdots
    !do n = dof+1,qdim
    !    call indicies(n-dof,m,k)
    !    do o = 1,d
    !       do j = 1, nop
    !          jsum(j) = d2Vdr2(r(k,j),De(k,j),re(k,j),re(k,j))*(xq(vectpos(m,k))-xq(vectpos(m,j)))* &
    !                  (xq(vectpos(o,k))-xq(vectpos(o,j)))/r(k,j)**2
    !          jsum(k) = 0.0_mk
    !       enddo
    !       do l = 1,nop
    !          lsum(l) = d2Vdr2(r(l,k),De(l,k),re(l,k),re(l,k))*(xq(vectpos(m,l))-xq(vectpos(m,k)))* &
    !                 (xq(vectpos(o,k))-xq(vectpos(o,l)))/r(l,k)**2
    !          lsum(k) = 0.0_mk
    !       enddo
    !       intsum(o) = -xq(vectpos(o,k)+dof)*sum(jsum)-sum(lsum)
    !    enddo
    !    fv(n) = sum(intsum) - d2Vextdx2(n-dof)
    !enddo
    return
  end subroutine potsnbody2

  ! Subroutine for the Harmonic oscillator potential
  subroutine pots2       ! pdots and v_pdots
    use coordinates
    call masses
    fv(1:dof) = -xq(1:dof)*mass !pdot = -dH/dq = -q*m
    !fv(dof+1:qdim) = -xq(dof+1:qdim)*mass !v_pdot = -dK/dv_q = -v_q*d2H/dq2 = -v_q*m
    return
  end subroutine pots2

  ! Subroutines for the Meyer Hamiltonian
  subroutine potsMeyer2
    use coordinates     ! pdots and v_pdots 
    fv(1) = -(xq(1)+8*kMeyer*xq(1)*xq(2)**2) !p1dot = -(q1 + 8*k*q1*q2^2)
    fv(2) = -(xq(2)+8*kMeyer*xq(2)*xq(1)**2) !p2dot = -(q2 + 8*k*q2*q1^2)
    fv(3) = -xq(3)*(1+8*kMeyer*xq(2)**2)-xq(4)*16*kMeyer*xq(1)*xq(2) !v_p1dot = -vq1(16*k*q1*q2)-vq2(1+8*kq2^2)
    fv(4) = -xq(4)*(1+8*kMeyer*xq(1)**2)-xq(3)*16*kMeyer*xq(1)*xq(2) !v_p2dot = -vq1(1+8*kq1^2)-vq2(16*k*q1*q2)
    return
  end subroutine potsMeyer2

  ! Function to compute separation distance between the ith and jth particles
  function r(i,j)
    use coordinates
    implicit none
    integer, intent(in) :: i, j
    real (kind = mk) :: r
    r = sqrt(sum((xq((i-1)*d+1:i*d)-xq((j-1)*d+1:j*d))**2))
  end function

  ! Functions to return the values of the chosen potential and its derivatives
  ! ------------------ Choose which potential the particles interact by here ---------------------
  function V(sep, De, re, a)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: sep, De, re, a
    real (kind = mk) :: V
    !V = 0.5_mk*sep**2 !harm osc potential
    !V = 0.0_mk !free particle
    V = De*(1.0_mk-exp(-(sep-re)*a))**2.0_mk - De
    !V =  De*(exp(-2*a*(sep-re))-2*exp(-a*(sep-re))) !Morse potential
  end function

  function dVdr(sep, De, re, a)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: sep, De, re, a
    real (kind = mk) :: dVdr
    !dVdr = sep !harm osc potential
    !dVdr = 0.0_mk !free particle
    dVdr = (2.0_mk*De*a)*(1.0_mk-exp(-(sep-re)*a))*exp(-(sep-re)*a)
    !dVdr = De*(-2*a*exp(-2*a*(sep-re))+2*a*exp(-a*(sep-re))) !Morse potential
  end function

  function d2Vdr2(sep, De, re, a)
    use coordinates
    implicit none
    real (kind = mk), intent(in) :: sep, De, re, a
    real (kind = mk) :: d2Vdr2
    !d2Vdr2 = 1.0_mk !harm osc potential
    !d2Vdr2 = 0.0_mk !free partilce
    d2Vdr2 = (2.0_mk*De*a**2.0_mk)*exp(-2.0_mk*(sep-re)*a)*(-exp((sep-re)*a)+2.0_mk)
    !d2Vdr2 = De*(4*a**2*exp(-2*a*(sep-re))-2*a**2*exp(-a*(sep-re))) !Morse potential
  end function
  ! ------------------------------------------------------------------------------------------------ 

  ! Functions to return values of external potential (like a harm osc trap)
  function Vext(n)
    use coordinates
    implicit none
    integer, intent(in) :: n
    real (kind = mk) :: Vext
    call masses
    Vext = 0.5_mk*mass(n)*omegaext**2*xq(n)**2
  end function

  function dVextdx(n)
    use coordinates
    implicit none
    integer, intent(in) :: n
    real (kind = mk) :: dVextdx
    call masses
    dVextdx = mass(n)*omegaext**2*xq(n)
  end function

  function d2Vextdx2(n)
    use coordinates
    implicit none
    integer, intent(in) :: n
    real (kind = mk) :: d2Vextdx2
    call masses
    d2Vextdx2 = mass(n)*omegaext**2
  end function
  ! ------------------------------------------------------------------------------------------------

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

  ! Function to choose the correct dissociation energy in Hartree
  function De(part1, part2)
    use coordinates
    implicit none
    integer, intent(in) :: part1, part2
    real (kind = mk) :: De
    if (partkind(part1) == 0 .AND. partkind(part2) == 0) then 
         De = DeK2
    else if ((partkind(part1)==0 .AND.partkind(part2)==1).OR.(partkind(part1)==1 .AND.partkind(part2)==0)) then
         De = DeKRb
    else !(part1 == 1 .AND. part2 == 1) 
         De = DeRb2
    endif
   end function

   ! Function to choose the correct a_Morse parameter
   function apar(part1, part2)
     use coordinates
     implicit none
     integer, intent(in) :: part1, part2
     real (kind = mk) :: apar
     call masses
     if (partkind(part1) == 0 .AND. partkind(part2) == 0) then
          apar = aK2
     else if ((partkind(part1)==0 .AND.partkind(part2)==1).OR.(partkind(part1)==1 .AND.partkind(part2)==0)) then
          apar = aKRb
     else !(part1 == 1 .AND. part2 == 1) 
          apar = aRb2
     endif
   end function


   ! Function to choose the correct equilibrium distance in Bohr
  function re(part1, part2)
    use coordinates
    implicit none
    integer, intent(in) :: part1, part2
    real (kind = mk) :: re
    if (partkind(part1) == 0 .AND. partkind(part2) == 0) then
        re = reK2
    else if ((partkind(part1)==0 .AND.partkind(part2)==1).OR.(partkind(part1)==1 .AND.partkind(part2)==0)) then 
        re = reKRb
    else !(part1 == 1 .AND. part2 == 1) 
        re = reRb2
    endif
  end function

   ! Subroutine to compute the energy of the system (in Hartree)
   subroutine energy
     use coordinates
     real (kind = mk), external :: r, De, re, V, Vext, apar
     call masses
     Vtot = 0.0_mk
     T = sum((xq(qdim+1:qdim+dof))**2/(2*mass))
     do i = 1,nop
         do j = i+1,nop
             Vtot = Vtot + V(r(i,j), De(i,j), re(i,j), apar(i,j))
         enddo
     enddo
     do j = 1,dof
         Vtot = Vtot + Vext(j)
     enddo
     E = T + Vtot
     return
   end subroutine energy

   ! Function to create a plot of the potential surface
   subroutine potmap
     use coordinates
     real (kind = mk) x1, x2, x3, r12, r13, r23, V12, V23, V13
     real (kind = mk), external :: r, De, re, V, Vext, apar
     real (kind = mk), dimension(:), allocatable :: potvec
     call masses
         ! Create map of potential
    do i = -100,100 !y1
        do j = -100,100 !y2
            Vtot = 0.0_mk
            x1 = sqrt(mu)*(sqrt(mu12)*real(i,mk)/(mi(1)*10.0_mk) + sqrt(mu123)*real(j,mk)/ &
                    ((mi(1)+mi(2))*10.0_mk))
            x2 = sqrt(mu)*(-sqrt(mu12)*real(i,mk)/(mi(2)*10.0_mk) + sqrt(mu123)*real(j,mk)/ &
                    ((mi(1)+mi(2))*10.0_mk))
            x3 = sqrt(mu)*(-sqrt(mu123)*real(j,mk)/(mi(3)*10.0_mk))
            r12 = abs(x1 - x2)
            r13 = abs(x1 - x3)
            r23 = abs(x2 - x3)
            V12 = V(r12,De(1,2),re(1,2),apar(1,2))
            V13 = V(r13,De(1,3),re(1,3),apar(1,3))
            V23 = V(r23,De(2,3),re(2,3),apar(2,3))
            Vtot = V12 + V13 + V23 + &
                    0.5_mk*mu*omegaext**2*((real(i,mk)/10.0_mk)**2+(real(j,mk)/10.0_mk)**2)
            write(4,*) real(i,mk)/10.0_mk, real(j,mk)/10.0_mk, Vtot
        enddo
        write(4,*)
    enddo
    do i = -5*1571,5*1571
       x1 = sqrt(mu)*100_mk*(sqrt(mu12)*cos(real(i,mk)/5e3_mk)/mi(1) + &
               sqrt(mu123)*sin(real(i,mk)/5e3_mk)/(mi(1)+mi(2)))
       x2 = sqrt(mu)*100_mk*(-sqrt(mu12)*cos(real(i,mk)/5e3_mk)/mi(2) + &
                    sqrt(mu123)*sin(real(i,mk)/5e3_mk)/(mi(1)+mi(2)))
       x3 = sqrt(mu)*100_mk*(-sqrt(mu123)*sin(real(i,mk)/5e3_mk)/mi(3))
       r12 = abs(x1 - x2)
       r13 = abs(x1 - x3)
       r23 = abs(x2 - x3)
       V12 = V(r12,De(1,2),re(1,2),apar(1,2))
       V13 = V(r13,De(1,3),re(1,3),apar(1,3))
       V23 = V(r23,De(2,3),re(2,3),apar(2,3))
       Vtot = V12 + V13 + V23
       write(2,*) real(i,mk)/5e3_mk,  Vtot
    enddo
    return
    end subroutine potmap

    ! Function to calculate jacobi coordinates
    subroutine Jacobi
      use coordinates
      implicit none
      call masses
      xcm = sqrt(Mtot/mu)*(mi(1)*xq(1) + mi(2)*xq(2) + mi(3)*xq(3))/Mtot
      y1 = sqrt(mu12/mu)*(xq(1)-xq(2))
      y2 = sqrt(mu123/mu)*((mi(1)*xq(1) + mi(2)*xq(2))/(mi(1) + mi(2)) - xq(3))
    end subroutine Jacobi

    ! Subroutine to skip through slow initial oscillations
    subroutine fstfwrd
      use coordinates
      implicit none
      real (kind = mk), external :: re
      real (kind = mk) Vtemp, Vfwrd, deltaV
      call energy
      Vtemp = Vtot !Calculate initial potential energy
      hyprad = 50.0_mk ! Redefine position to skip slow beginning oscs
      hypang = pi/2.0_mk - re(1,2)/hyprad + eps*hyprad/100.0_mk !particles 1,2
      qi(1) = sqrt(mu)*hyprad*(sqrt(mu12)*cos(hypang)/mi(1) + sqrt(mu123)*sin(hypang)/(mi(1)+mi(2)))
      qi(2) = sqrt(mu)*hyprad*(-sqrt(mu12)*cos(hypang)/mi(2) + sqrt(mu123)*sin(hypang)/(mi(1)+mi(2)))
      qi(3) = -sqrt(mu*mu123)*hyprad*sin(hypang)/mi(3)
      !write(*,*) y1, y2
      !y1 = cos(hypang)*hyprad
      !y2 = 50.0_mk
      !qi(1) = sqrt(mu)*(sqrt(mu12)*y1/mi(1) + sqrt(mu123)*y2/(mi(1)+mi(2)))
      !qi(2) = sqrt(mu)*(-sqrt(mu12)*y1/mi(2) + sqrt(mu123)*y2/(mi(1)+mi(2)))
      !qi(3) = -sqrt(mu*mu123)*y2/mi(3)
      xq(1:dof) = qi
      call energy
      Vfwrd = Vtot !Calculate final potential energy
      deltaV = Vfwrd - Vtemp
      KEy2 = KEy2 - deltaV*nanoKperEh
      !write(*,*) 'Vint: ', Vtemp, 'Vfin: ', Vfwrd, 'deltaV: ', deltaV
    end subroutine fstfwrd

    subroutine rng
      use coordinates
      implicit none
      real (kind = mk), external :: gauss
      integer :: i
      logical :: choose(rngsize*2)
      ! Create vector of random numbers, first rngsize elements are x's, sencond rngsize are y's
      call rgnf_lux(rngvec,rngsize*2)
      rvec = real(rngvec,mk)
      rvec(1:size(rvec)/2) = 8.0_mk*cgauss*rvec(1:size(rvec)/2) - 4.0_mk*cgauss + bgauss
      !write(*,*) rvec
      do i = 1, rngsize
         if ( gauss(rvec(i)) < rvec(i+rngsize) ) then 
            choose(i) = .FALSE.
            choose(i+rngsize) = choose(i)
         else 
            choose(i) = .TRUE.
            choose(i+rngsize) = choose(i)
         endif
      enddo
      rvecselect = pack(rvec, choose)
      !write(*,*) choose
      write(*,*) 'the number of remaining points is: ', size(rvecselect)/2
      write(*,*) rngsize
      write(*,*) 'the fraction of points that survived is: ', real(size(rvecselect))/real(rngsize*2)
      open (1, file = 'rng.dat')
      open (2, file = 'rngunfilterd.dat')
      do i = 1, size(rvecselect)/2
         write(1,*) rvecselect(i), rvecselect(size(rvecselect)/2+i)
      enddo
      do i = 1, rngsize
         write(2,*) rvec(i), rvec(rngsize + i)
      enddo
      !write(102,*) rvecselect(size(rvecselect)/2+1:size(rvecselect))
      close(1)
      close(2)
    end subroutine rng

    function gauss(x)
      use coordinates
      implicit none
      real (kind = mk), intent(in) :: x
      real (kind = mk) :: gauss
      gauss = agauss*exp(-((x-bgauss)**2)/(2*cgauss**2))
    end function

    subroutine rcheck
      use coordinates
      implicit none
      integer :: i
      do i = 2, ncall
         if ((radvec(i) > radvec(i-1)).AND.(radvec(i) > radvec(i+1))) then 
            rtimestemp(i) = real(i,mk)*xdt
            rchoose(i) = .TRUE.
            if (radvec(i) > (cutoffrad-1.0_mk)) then
               rchoose(i) = .FALSE.
            endif
         else
            rchoose(i) = .FALSE.
         endif
      enddo
      rtimes = pack(rtimestemp, rchoose)
      radmax = pack(radvec, rchoose)
      deallocate(rchoose,rtimestemp,radvec)
      allocate(rchoose(size(radmax)))
      rchoose = .TRUE.
      do i = 1, size(radmax)
        if (radmax(i) < 50) then
                exit
        endif
        rchoose(i) = .FALSE.
      enddo
      radvec = pack(radmax,rchoose)
      rtimestemp = pack(rtimes,rchoose)
      return
    end subroutine rcheck
