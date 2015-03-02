module simpol_module

  implicit none

  integer, parameter :: ngroups = 31
  integer, parameter :: precision = selected_real_kind(14,300)
  real(precision), dimension(ngroups) :: bk1, bk2, bk3, bk4
  
  parameter &
       (bk1= (/-4.26938D+02,-4.11248D+02,-1.46442D+02,&
       3.50262D+01,-8.72770D+01,&
       5.73335D+00,-2.61268D+02,-7.25373D+02,-7.29501D+02,-1.37456D+01,&
       -7.98796D+02,-3.93345D+02,-1.44334D+02,4.05265D+01,-7.07406D+01,&
       -7.83648D+02,-5.63872D+02,-4.53961D+02,3.71375D+01,-5.03710D+02,&
       -3.59763D+01,-6.09432D+02,-1.02367D+02,-1.93802D+03,-5.26919D+00,&
       -2.84042D+02,1.50093D+02,-2.03387D+01,-8.38064D+02,-5.27934D+01,&
       -1.61520D+03/),&
       bk2=(/2.89223D-01,8.96919D-01,1.54528D+00,-9.20839D-01,1.78059D+00,&
       1.69764D-02,-7.63282D-01,8.26326D-01,9.86017D-01,5.23486D-01,&
       -1.09436D+00,-9.51778D-01,-1.85617D+00,-2.43780D+00,-1.06674D+00,&
       -1.03439D+00,-7.18416D-01,-3.26105D-01,-2.66753D+00,1.04092D+00,&
       -4.08458D-01,1.50436D+00,-7.16253D-01,6.48262D-01,3.06435D-01,&
       -6.25424D-01,2.39875D-02,-5.48718D+00,-1.09600D+00,-4.63689D-01,&
       9.01669D-01/),&
       bk3=(/4.42057D-03,-2.48607D-03,1.71021D-03,2.24399D-03,-3.07187D-03,&
       -6.28957D-04,-1.68213D-03,2.50957D-03,-2.92664D-03,5.50298D-04,&
       5.24132D-03,-2.19071D-03,-2.37491D-05,3.60133D-03,3.73104D-03,&
       -1.07148D-03,2.63016D-03,-1.39780D-04,1.01483D-03,-4.12746D-03,&
       1.67264D-03,-9.09024D-04,-2.90670D-04,1.73245D-03,3.25397D-03,&
       -8.22474D-04,-3.37969D-03,8.39075D-03,-4.24385D-04,-5.11647D-03,&
       1.44536D-03/),&
       bk4=(/2.92846D-01,1.40312D-01,-2.78291D-01,-9.36300D-02,-1.04341D-01,&
       7.55434D-03,2.89038D-01,-2.32304D-01,1.78077D-01,-2.76950D-01,&
       -2.28040D-01,3.05843D-01,2.88290D-01,9.86422D-02,-1.44003D-01,&
       3.15535D-01,-4.99470D-02,-3.93916D-02,2.14233D-01,1.82790D-01,&
       -9.98919D-02,-1.35495D-01,-5.88556D-01,3.47940D-02,-6.81506D-01,&
       -8.80549D-02,1.52789D-02,1.07884D-01,2.81812D-01,3.84965D-01,&
       2.66889D-01/))
!!$  data bk1 /ngroups*0.0/,&
!!$       bk2 /ngroups*0.0/,&
!!$       bk3 /ngroups*0.0/,&
!!$       bk4 /ngroups*0.0/

contains

  function simpolvp(temp, nuk) result(vp_atm)
! arguments
    real(precision), intent(in) :: temp
    integer, dimension(ngroups), intent(in) :: nuk
! local variables
    real(precision), dimension(ngroups) :: bk
    real(precision) :: vp_atm
! body    
    bk = bk1/temp + bk2 + bk3*temp + bk4*log(temp)
    vp_atm = 10**sum(nuk*bk)
  end function simpolvp

  function press2conc(temp,press) result(conc)
! arguments
    real(precision), intent(in) :: press, temp        ! atm, K
! local variables
    real(precision)             :: conc               ! molecules/cm^3
    real(precision), parameter  :: Av = 6.0221413d+23 ! molecules/mole
    real(precision), parameter  :: R = 82.05746       ! cm^3 atm/ mole/ K
! body
    conc = press/temp * Av/R
  end function press2conc

end module simpol_module
