MODULE constants

  IMPLICIT NONE
  PRIVATE dp
  INTEGER, PARAMETER :: dp=8
  INTEGER, PARAMETER :: mnsp=250, mre=2000
  INTEGER i
  ! photolysis rate indices (taken from photolysis.txt)
  INTEGER, dimension(35) :: jidx = (/1,2,3,4,5,6,7,8,11,12,13,14,15,16,&
       17,18,19,21,22,23,24,31,32,33,34,35,41,51,52,53,54,55,56,57,61/)
  ! generic reaction rate variables
  REAL(dp) kro2no, kro2ho2, kapho2, kapno, kro2no3, kno3al, kdec, &
    krosec, kalkoxy, kalkpxy, kroprim, &
    KCH3O2,K298CH3O2,KMT18,NCD, &
    NC,NC1,NC2,NC3,NC4,NC7,NC8,NC9,NC10,K120,K12I,KR12,FC12,NC12,F12, &
    NC13,NC14,K150,K15I,KR15,FC15,NC15,F15,NC16,K170,K17I,KR17,FC17, &
    NC17,F17,KPPN0,KPPNI,KRPPN,FCPPN,NCPPN,FPPN,KBPPN
  ! variables for calculation of kfpan and kbpan
  REAL(dp) kfpan, kbpan
  REAL(dp) kc0, kci, krc, fcc, fc
  REAL(dp) kd0, kdi, krd, fcd, fd
  ! variables for calculation of kmt01
  REAL(dp) kmt01
  REAL(dp) k10, k1i, kr1, fc1, f1
  ! variables for calculation of kmt02
  REAL(dp) kmt02
  REAL(dp) k20, k2i, kr2, fc2, f2
  ! variables for calculation of kmt03
  REAL(dp) kmt03
  REAL(dp) k30, k3i, kr3, fc3, f3
  ! variables for calculation of kmt04
  REAL(dp) kmt04
  REAL(dp) k40, k4i, kr4, fc4, f4
  ! variables for calculation of kmt05
  REAL(dp) kmt05
  ! variables for calculation of kmt06
  REAL(dp) kmt06
  ! variables for calculation of kmt07
  REAL(dp) kmt07
  REAL(dp) k70, k7i, kr7, fc7, f7
  ! variables for calculation of kmt08
  REAL(dp) kmt08
  REAL(dp) k80, k8i, kr8, fc8, f8
  ! variables for calculation of kmt09
  REAL(dp) kmt09
  REAL(dp) k90, k9i, kr9, fc9, f9
  ! variables for calculation of kmt10
  REAL(dp) kmt10
  REAL(dp) k100, k10i, kr10, fc10, f10
  ! variables for calculation of kmt11
  REAL(dp) kmt11
  REAL(dp) k1,k2,k3,k4
  ! variables for calculation of kmt12
  REAL(dp) kmt12
  REAL(dp) k0, ki, x, ssign,f
  ! variables for calculation of kmt13
  REAL(dp) kmt13
  REAL(dp) k130, k13i, kr13, fc13, f13
  ! variables for calculation of kmt14
  REAL(dp) kmt14
  REAL(dp) k140, k14i, kr14, fc14, f14
  ! variables for calculation of kmt15
  REAL(dp) kmt15
  ! variables for calculation of kmt16
  REAL(dp) kmt16
  REAL(dp) k160, k16i, kr16, fc16, f16
  ! variables for calculation of kmt17
  REAL(dp) kmt17
  ! variables for calculation of photolysis reaction rates
  REAL(dp)  l(61), mm(61), nn(61), j(61)
  INTEGER k

CONTAINS

  !***************************************************************************

  SUBROUTINE mcm_constants(time, temp, M, N2, O2, RO2, H2O)
    ! calculates rate constants from arrhenius informtion

    REAL(dp) time, temp, M, N2, O2, RO2, H2O

    ! ************************************************************************
    ! define generic reaction rates.
    ! ************************************************************************

    ! constants used in calculation of reaction rates
    M  = 2.55E+19
    N2 = 0.79*M
    O2 = 0.2095*M

    ! kro2no : ro2      + no      = ro      + no2
    !        : ro2      + no      = rono2
    ! iupac 1992
    kro2no    = 2.40d-12*EXP(360.0/temp)

    ! kro2ho2: ro2      + ho2     = rooh    + o2
    ! mcm protocol [1997]
    kro2ho2   = 2.91d-13*EXP(1300.0/temp)

    ! kapho2 : rcoo2    + ho2     = products
    ! mcm protocol [1997]
    kapho2    = 4.30d-13*EXP(1040.0/temp)

    ! kapno  : rcoo2    + no      = products
    ! mej [1998]
    kapno = 8.10d-12*EXP(270.0/temp)

    ! kro2no3: ro2      + no3     = products
    ! mcm protocol [1997]
    kro2no3   = 2.50d-12

    ! kno3al : no3      + rcho    = rcoo2   + hno3
    ! mcm protocol [1997]
    kno3al    = 1.44d-12*EXP(-1862.0/temp)

    ! kdec   : ro                 = products
    ! mcm protocol [1997]
    kdec      = 1.00d+06
    krosec = 1.50e-14*EXP(-200.0/temp)

    kalkoxy=3.70d-14*EXP(-460.0/temp)*o2
    kalkpxy=1.80d-14*EXP(-260.0/temp)*o2

    ! -------------------------------------------------------------------
    ! complex reactions
    ! -------------------------------------------------------------------

    ! kfpan kbpan
    ! formation and decomposition of pan
    ! iupac 1997
    kc0     = 2.70d-28*m*(temp/298.0)**(-7.1)
    kci     = 1.21d-11*(temp/298.0)**(-0.9)
    krc     = kc0/kci
    fcc     = 0.30
    fc      = 10**(dlog10(fcc)/(1+(dlog10(krc))**2))
    kfpan   = (kc0*kci)*fc/(kc0+kci)

    kd0     = 4.90d-03*m*EXP(-12100.0/temp)
    kdi     = 5.40d+16*EXP(-13830.0/temp)
    krd     = kd0/kdi
    fcd     = 0.30
    fd      = 10**(dlog10(fcd)/(1+(dlog10(krd))**2))
    kbpan   = (kd0*kdi)*fd/(kd0+kdi)

    ! kmt01  : o        + no      = no2
    ! iupac 1997
    k10     = 1.00d-31*m*(temp/300.0)**(-1.6)

    k1i     = 3.00d-11*(temp/300.0)**(0.3)
    kr1     = k10/k1i
    fc1     = EXP(-temp/1850.0)
    f1      = 10**(dlog10(fc1)/(1+(dlog10(kr1))**2))
    kmt01   = (k10*k1i)*f1/(k10+k1i)

    ! kmt02  : o        + no2     = no3
    ! iupac 1997
    k20     = 9.00d-32*m*(temp/300.0)**(-2.0)
    k2i     = 2.20d-11
    kr2     = k20/k2i
    fc2     = EXP(-temp/1300.0)
    f2      = 10**(dlog10(fc2)/(1+(dlog10(kr2))**2))
    kmt02   = (k20*k2i)*f2/(k20+k2i)

    ! kmt03  : no2      + no3     = n2o5
    ! iupac 1997
    k30     = 2.80d-30*m*(temp/300.0)**(-3.5)
    k3i     = 2.00d-12*(temp/300.0)**(0.2)
    kr3     = k30/k3i
    fc3     = 2.5*EXP(-1950.0/temp) + 0.9*EXP(-temp/430.0)
    f3      = 10**(dlog10(fc3)/(1+(dlog10(kr3))**2))
    kmt03   = (k30*k3i)*f3/(k30+k3i)

    ! kmt04  : n2o5               = no2     + no3
    ! iupac 1997
    k40     = 1.00d-03*m*(temp/300.0)**(-3.5)*EXP(-11000.0/temp)
    k4i     = 9.70d+14*(temp/300.0)**(0.1)*EXP(-11080.0/temp)
    kr4     = k40/k4i
    fc4     = 2.5*EXP(-1950.0/temp) + 0.9*EXP(-temp/430.0)
    f4      = 10**(dlog10(fc4)/(1+(dlog10(kr4))**2))
    kmt04   = (k40*k4i)*f4/(k40+k4i)

    ! kmt05  : oh       + co(+o2) = ho2     + co2
    ! iupac 1997
    kmt05  = 1 + ((0.6*m)/(2.652d+19*(273.0/temp)))

    ! kmt06  : ho2      + ho2     = h2o2    + o2
    ! water enhancement factor
    ! iupac 1992

    kmt06  = 1 + (1.40d-21*EXP(2200.0/temp)*h2o)

    ! kmt07  : oh       + no      = hono

    ! iupac 1992
    k70     = 7.40d-31*m*(temp/300.0)**(-2.4)
    k7i     = 3.30d-11
    kr7     = k70/k7i
    fc7     = EXP(-temp/1420.0)
    f7      = 10**(dlog10(fc7)/(1+(dlog10(kr7))**2))
    kmt07   = (k70*k7i)*f7/(k70+k7i)

    ! kmt08  : oh       + no2     = hno3

    ! nasa 2000
    k80     = 2.40d-30*m*(temp/300.0)**(-3.1)
    k8i     = 1.70d-11*(temp/300.0)**(-2.1)
    kr8     = k80/k8i
    fc8     = 0.6
    f8      = 10**(dlog10(fc8)/(1+(dlog10(kr8))**2))
    kmt08   = (k80*k8i)*f8/(k80+k8i)

    ! kmt09  : ho2      + no2     = ho2no2
    ! iupac 1997

    k90     = 1.80d-31*m*(temp/300.0)**(-3.2)
    k9i     = 4.70d-12
    kr9     = k90/k9i
    fc9     = 0.6
    f9      = 10**(dlog10(fc9)/(1+(dlog10(kr9))**2))
    kmt09   = (k90*k9i)*f9/(k90+k9i)

    ! kmt10  : ho2no2             = ho2     + no2
    ! iupac 1997

    k100     = 4.10d-05*m*EXP(-10650.0/temp)
    k10i     = 5.70d+15*EXP(-11170.0/temp)
    kr10     = k100/k10i
    fc10     = 0.5
    f10      = 10**(dlog10(fc10)/(1+(dlog10(kr10))**2))
    kmt10    = (k100*k10i)*f10/(k100+k10i)

    ! kmt11  : oh       + hno3    = h2o     + no3
    ! iupac 1992

    k1     = 7.20d-15*EXP(785.0/temp)
    k3     = 1.90d-33*EXP(725.0/temp)
    k4     = 4.10d-16*EXP(1440.0/temp)
    k2     = (k3*m)/(1+(k3*m/k4))
    kmt11  = k1 + k2

    ! kmt12

    k0 = 3.00d-31*((temp/300.0)**(-3.3))*m
    ki = 1.50d-12
    fc = 0.6
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    kmt12=(k0*ki*f)/(k0+ki)

    ! kmt13  : ch3o2    + no2     = ch3o2no2
    ! iupac 1997

    k130     = 2.50d-30*((temp/298.0)**(-5.5))*m
    k13i     = 7.50d-12
    kr13     = k130/k13i
    fc13     = 0.36
    f13      = 10**(dlog10(fc13)/(1+(dlog10(kr13))**2))
    kmt13    = (k130*k13i)*f13/(k130+k13i)

    ! kmt14  : ch3o2no2           = ch3o2   + no2
    ! iupac 1997

    k140     = 9.00d-05*EXP(-9690.0/temp)*m
    k14i     = 1.10d+16*EXP(-10560.0/temp)
    kr14     = k140/k14i
    fc14     = 0.36
    f14      = 10**(dlog10(fc14)/(1+(dlog10(kr14))**2))
    kmt14    = (k140*k14i)*f14/(k140+k14i)

    ! kmt15

    k0=6.00d-29*((temp/298.0)**(-4.0))*m
    ki=9.00d-12*((temp/298.0)**(-1.1))
    fc=0.7
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    kmt15=(k0*ki*f)/(k0+ki)

    ! kmt16  :  oh  +  c3h6
    ! atkinson 1994

    k160     = 3.00d-27*((temp/298.0)**(-3.0))*m
    k16i     = 2.80d-11*((temp/298.0)**(-1.3))
    kr16     = k160/k16i
    fc16     = 0.5
    f16      = 10**(dlog10(fc16)/(1+(dlog10(kr16))**2))
    kmt16    = (k160*k16i)*f16/(k160+k16i)

    ! kmt17

    k0 = 5.00d-30*((temp/298.0)**(-1.5))*m
    ki = 9.40d-12*EXP(-700.0/temp)
    fc = (EXP(-temp/580.0) + EXP(-2320.0/temp))
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    kmt17=(k0*ki*f)/(k0+ki)

    !       mcm 2001

    kroprim  = 3.70d-14*EXP(-460.0/temp)
    krosec   = 1.80d-14*EXP(-260.0/temp)

    ! ************************************************************************
    ! define fixed photolysis reaction rates
    ! ************************************************************************

    DO i = 1, size(jidx)
        j(jidx(i)) = 1.e-30
    END DO

  END SUBROUTINE mcm_constants

END MODULE constants
