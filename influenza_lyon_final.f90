!===============================================================================
! INFLUENZA R_eff: PERIODIC CALIBRATION + MCMC DIAGNOSTICS
! Compile:
!   gfortran -O3 -std=f2008 influenza_periodic.f90 -o influenza_periodic
!
! Run:
!   ./influenza_periodic
!===============================================================================

module precision_mod
  implicit none

  integer, parameter :: dp=selected_real_kind(15,307)
  integer, parameter :: i8=selected_int_kind(18)

  real(dp), parameter :: PI=3.1415926535897932384626433832795_dp

end module precision_mod


module settings_mod

  use precision_mod

  implicit none

  integer, parameter :: NG=4
  integer, parameter :: NY=8
  integer, parameter :: NDIM=4
  integer, parameter :: NDATA=13

  integer, parameter :: DAYS=7
  integer, parameter :: SHIFTS=3

  real(dp), parameter :: SHIFT_DURATION=1.0_dp/3.0_dp
  real(dp), parameter :: DAYS_PER_WEEK=7.0_dp
  real(dp), parameter :: SHIFT_DURATION_WEEKS=SHIFT_DURATION/7.0_dp

  real(dp), parameter :: SEASONAL_PERIOD_DAYS=364.0_dp

  real(dp), parameter :: BURN_IN_WEEKS=15.0_dp
  real(dp), parameter :: MODEL_MAX_WEEK=80.0_dp

  real(dp), parameter :: BURN_IN_DAYS=BURN_IN_WEEKS*7.0_dp
  real(dp), parameter :: MODEL_DAYS=MODEL_MAX_WEEK*7.0_dp
  real(dp), parameter :: SIM_END_DAYS=BURN_IN_DAYS+MODEL_DAYS

  integer, parameter :: TOTAL_SHIFTS= &
       int(SIM_END_DAYS/SHIFT_DURATION+0.5_dp)

  integer, parameter :: N_EVAL=TOTAL_SHIFTS+1


  !=====================================================================
  ! MCMC SETTINGS
  !=====================================================================

  integer, parameter :: NWALKERS=64
  integer, parameter :: NSTEPS=220000
  integer, parameter :: BURN_MCMC=65000
  integer, parameter :: THIN=20

  real(dp), parameter :: A_STRETCH=2.0_dp


  !=====================================================================
  ! POSTERIOR TRAJECTORY SETTINGS
  !=====================================================================

  integer, parameter :: N_POSTERIOR_MODELS=2000


  !=====================================================================
  ! MLE / PROFILE SETTINGS
  !=====================================================================

  integer, parameter :: N_MLE_STARTS=12
  integer, parameter :: N_PROFILE=31

  integer, parameter :: MIN_VALID_POSTERIOR=100


  !=====================================================================
  ! PARAMETER BOUNDS
  !=====================================================================

  real(dp), parameter :: ALPHA_LO=0.2_dp
  real(dp), parameter :: ALPHA_HI=1.0_dp

  real(dp), parameter :: SEASON_LO=0.05_dp
  real(dp), parameter :: SEASON_HI=0.8_dp

  real(dp), parameter :: SIGMA_LO=0.001_dp
  real(dp), parameter :: SIGMA_HI=1.0_dp

  real(dp), parameter :: SHIFT_LO=-8.0_dp
  real(dp), parameter :: SHIFT_HI=8.0_dp

  real(dp), parameter :: SHIFT_PRIOR_SIGMA=3.0_dp

  real(dp), parameter :: NEG_TOL=1.0e-9_dp


  !=====================================================================
  ! POPULATION
  !=====================================================================

  real(dp), parameter :: POP(NG)= &
       [11.0_dp,27.0_dp,8.0_dp,29.0_dp]


  !=====================================================================
  ! TRANSMISSION MATRICES
  !=====================================================================

  real(dp), parameter :: B1(NG,NG)= &
       0.9_dp*reshape([ &
       1.5_dp,0.0_dp,0.0_dp,0.0_dp, &
       0.0_dp,1.4_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.0_dp,1.2_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.0_dp,1.0_dp],[NG,NG])

  real(dp), parameter :: B2(NG,NG)= &
       1.1_dp*reshape([ &
       1.7_dp,0.0_dp,0.0_dp,0.0_dp, &
       0.0_dp,1.5_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.0_dp,1.3_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.0_dp,1.1_dp],[NG,NG])


  !=====================================================================
  ! RECOVERY RATES
  !=====================================================================

  real(dp), parameter :: GAM1(NG)= &
       1.3_dp*[0.15_dp,0.14_dp,0.12_dp,0.10_dp]

  real(dp), parameter :: GAM2(NG)= &
       [0.14_dp,0.13_dp,0.11_dp,0.09_dp]

end module settings_mod


module data_mod

  use precision_mod
  use settings_mod

  implicit none

  real(dp), parameter :: WEEKS_DATA(NDATA)=[ &
   50.02_dp-34.0_dp, &
   51.00_dp-34.0_dp, &
   52.01_dp-34.0_dp, &
   53.01_dp-34.0_dp, &
   54.01_dp-34.0_dp, &
   55.00_dp-34.0_dp, &
   56.01_dp-34.0_dp, &
   57.03_dp-34.0_dp, &
   58.01_dp-34.0_dp, &
   59.00_dp-34.0_dp, &
   60.01_dp-34.0_dp, &
   61.00_dp-34.0_dp, &
   62.02_dp-34.0_dp]

  real(dp), parameter :: R_DATA(NDATA)=[ &
   1.4355_dp, &
   1.3805_dp, &
   1.1171_dp, &
   1.0857_dp, &
   0.9717_dp, &
   1.1152_dp, &
   1.2154_dp, &
   1.0601_dp, &
   0.8557_dp, &
   0.9579_dp, &
   0.8204_dp, &
   0.7948_dp, &
   1.0503_dp]

end module data_mod


module rng_mod

  use precision_mod

  implicit none

contains


  function urand(seed) result(u)

    integer,intent(inout)::seed

    real(dp)::u

    integer(i8)::a,m,s

    a=1103515245_i8
    m=2147483647_i8

    s=mod(a*int(seed,i8)+12345_i8,m)

    if(s<=0_i8)s=1_i8

    seed=int(s)

    u=real(seed,dp)/real(m,dp)

    u=min(max(u,1.0e-15_dp),1.0_dp-1.0e-15_dp)

  end function urand


  function nrand(seed) result(z)

    integer,intent(inout)::seed

    real(dp)::z
    real(dp)::u1,u2

    u1=urand(seed)
    u2=urand(seed)

    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*PI*u2)

  end function nrand

end module rng_mod


module stats_mod

  use precision_mod

  implicit none

contains


  function pearson(x,y) result(r)

    real(dp),intent(in)::x(:),y(:)

    real(dp)::r
    real(dp)::mx,my,sx,sy

    mx=sum(x)/real(size(x),dp)
    my=sum(y)/real(size(y),dp)

    sx=sqrt(sum((x-mx)**2))
    sy=sqrt(sum((y-my)**2))

    if(sx<=0.0_dp .or. sy<=0.0_dp)then

      r=0.0_dp

    else

      r=sum((x-mx)*(y-my))/(sx*sy)

    end if

  end function pearson


  function percentile(v,p) result(q)

    real(dp),intent(in)::v(:),p

    real(dp)::q,a

    real(dp),allocatable::s(:)

    integer::n,k

    n=size(v)

    if(n==0)then
      q=0.0_dp
      return
    end if

    allocate(s(n))

    s=v

    call sort(s)

    if(n==1)then
      q=s(1)
      deallocate(s)
      return
    end if

    a=p*real(n-1,dp)+1.0_dp

    k=max(1,min(n-1,int(floor(a))))

    q=s(k)+(a-real(k,dp))*(s(k+1)-s(k))

    deallocate(s)

  end function percentile


  recursive subroutine sort(a)

    real(dp),intent(inout)::a(:)

    call qsort(a,1,size(a))

  end subroutine sort


  recursive subroutine qsort(a,l,r)

    real(dp),intent(inout)::a(:)

    integer,intent(in)::l,r

    integer::i,j

    real(dp)::p,tmp

    if(l>=r)return

    i=l
    j=r

    p=a((l+r)/2)

    do

      do while(a(i)<p)
        i=i+1
      end do

      do while(a(j)>p)
        j=j-1
      end do

      if(i<=j)then

        tmp=a(i)
        a(i)=a(j)
        a(j)=tmp

        i=i+1
        j=j-1

      end if

      if(i>j)exit

    end do

    if(l<j)call qsort(a,l,j)

    if(i<r)call qsort(a,i,r)

  end subroutine qsort

end module stats_mod


module contact_mod

  use precision_mod
  use settings_mod

  implicit none

  real(dp),save::C_CONTACT(DAYS,SHIFTS,NG,NG)
  real(dp),save::L_SHIFTS(SHIFTS,NG,NG)

contains


  subroutine initialize_contacts()

    C_CONTACT=0.0_dp
    L_SHIFTS=0.0_dp


    !===================================================================
    ! SHIFT-SPECIFIC INTERVENTION MATRICES
    !===================================================================

    L_SHIFTS(1,:,:)=1.3_dp*reshape([ &
     0.50_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.45_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.40_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.20_dp],[4,4])

    L_SHIFTS(2,:,:)=1.3_dp*reshape([ &
     0.40_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.40_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.35_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.18_dp],[4,4])

    L_SHIFTS(3,:,:)=1.3_dp*reshape([ &
     0.30_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.25_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.20_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.15_dp],[4,4])


    !===================================================================
    ! DAY 1
    !===================================================================

    C_CONTACT(1,1,:,:)=reshape([ &
     0.089334_dp,0.030682_dp,0.168875_dp,0.226144_dp, &
     0.256768_dp,0.044045_dp,0.384233_dp,0.068801_dp, &
     0.033724_dp,0.0_dp,0.148653_dp,0.042188_dp, &
     0.003302_dp,0.009303_dp,0.239060_dp,0.033885_dp],[4,4])

    C_CONTACT(1,2,:,:)=reshape([ &
     0.044020_dp,0.014489_dp,0.010500_dp,0.276993_dp, &
     0.046295_dp,0.046848_dp,0.107251_dp,0.004278_dp, &
     0.020234_dp,0.003369_dp,0.158113_dp,0.019922_dp, &
     0.0_dp,0.005582_dp,0.043102_dp,0.016697_dp],[4,4])

    C_CONTACT(1,3,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.062522_dp,0.0_dp,0.055101_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.000825_dp,0.0_dp,0.058210_dp,0.0_dp],[4,4])


    !===================================================================
    ! DAY 2
    !===================================================================

    C_CONTACT(2,1,:,:)=reshape([ &
     0.104870_dp,0.045170_dp,0.173250_dp,0.402778_dp, &
     0.266790_dp,0.036838_dp,0.443762_dp,0.070583_dp, &
     0.028328_dp,0.001685_dp,0.124328_dp,0.062109_dp, &
     0.007429_dp,0.007815_dp,0.248391_dp,0.039778_dp],[4,4])

    C_CONTACT(2,2,:,:)=reshape([ &
     0.102281_dp,0.040057_dp,0.063000_dp,0.388058_dp, &
     0.096884_dp,0.062464_dp,0.321752_dp,0.025667_dp, &
     0.040469_dp,0.005054_dp,0.210817_dp,0.055078_dp, &
     0.0_dp,0.011164_dp,0.090203_dp,0.038796_dp],[4,4])

    C_CONTACT(2,3,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.024818_dp,0.0_dp,0.049690_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.004127_dp,0.0_dp,0.023106_dp,0.0_dp],[4,4])


    !===================================================================
    ! DAY 3
    !===================================================================

    C_CONTACT(3,1,:,:)=reshape([ &
     0.091923_dp,0.043466_dp,0.180250_dp,0.453627_dp, &
     0.275858_dp,0.127331_dp,0.411291_dp,0.073435_dp, &
     0.060703_dp,0.015162_dp,0.429742_dp,0.059766_dp, &
     0.007429_dp,0.016746_dp,0.256834_dp,0.034868_dp],[4,4])

    C_CONTACT(3,2,:,:)=reshape([ &
     0.111344_dp,0.011080_dp,0.063875_dp,0.492432_dp, &
     0.068249_dp,0.058460_dp,0.141197_dp,0.026023_dp, &
     0.040469_dp,0.006739_dp,0.197303_dp,0.015234_dp, &
     0.002889_dp,0.011164_dp,0.063542_dp,0.042234_dp],[4,4])

    C_CONTACT(3,3,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.056794_dp,0.002402_dp,0.085112_dp,0.0_dp, &
     0.004047_dp,0.0_dp,0.008108_dp,0.0_dp, &
     0.0_dp,0.001116_dp,0.052877_dp,0.0_dp],[4,4])


    !===================================================================
    ! DAY 4
    !===================================================================

    C_CONTACT(4,1,:,:)=reshape([ &
     0.146301_dp,0.031534_dp,0.171500_dp,0.355943_dp, &
     0.214291_dp,0.102505_dp,0.343399_dp,0.069870_dp, &
     0.056656_dp,0.020217_dp,0.345956_dp,0.043359_dp, &
     0.006604_dp,0.015629_dp,0.199513_dp,0.055493_dp],[4,4])

    C_CONTACT(4,2,:,:)=reshape([ &
     0.066030_dp,0.019602_dp,0.043750_dp,0.289036_dp, &
     0.050590_dp,0.025226_dp,0.098395_dp,0.017824_dp, &
     0.021583_dp,0.010108_dp,0.085138_dp,0.026953_dp, &
     0.005366_dp,0.005954_dp,0.047101_dp,0.025046_dp],[4,4])

    C_CONTACT(4,3,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.028159_dp,0.003203_dp,0.023615_dp,0.0_dp, &
     0.001349_dp,0.0_dp,0.010811_dp,0.0_dp, &
     0.002476_dp,0.000372_dp,0.026217_dp,0.0_dp],[4,4])


    !===================================================================
    ! DAY 5
    !===================================================================

    C_CONTACT(5,1,:,:)=reshape([ &
     0.075092_dp,0.025568_dp,0.187250_dp,0.188677_dp, &
     0.215723_dp,0.067670_dp,0.529857_dp,0.076287_dp, &
     0.032375_dp,0.131408_dp,0.228385_dp,0.035156_dp, &
     0.002889_dp,0.008931_dp,0.200846_dp,0.028483_dp],[4,4])

    C_CONTACT(5,2,:,:)=reshape([ &
     0.001295_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.001432_dp,0.001201_dp,0.001968_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.004054_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.001333_dp,0.000491_dp],[4,4])

    C_CONTACT(5,3,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.000477_dp,0.000400_dp,0.000492_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.001351_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.000444_dp,0.0_dp],[4,4])


    !===================================================================
    ! DAY 6
    !===================================================================

    C_CONTACT(6,1,:,:)=reshape([ &
     0.029778_dp,0.010227_dp,0.074375_dp,0.074935_dp, &
     0.085907_dp,0.026828_dp,0.211549_dp,0.030301_dp, &
     0.012141_dp,0.052226_dp,0.090543_dp,0.014062_dp, &
     0.000825_dp,0.003349_dp,0.079983_dp,0.011295_dp],[4,4])

    C_CONTACT(6,2,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.000477_dp,0.000400_dp,0.000492_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.001351_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.000444_dp,0.0_dp],[4,4])


    !===================================================================
    ! DAY 7
    !===================================================================

    C_CONTACT(7,1,:,:)=reshape([ &
     0.022010_dp,0.007670_dp,0.056000_dp,0.056202_dp, &
     0.064431_dp,0.020021_dp,0.158908_dp,0.022815_dp, &
     0.009443_dp,0.038749_dp,0.067569_dp,0.010547_dp, &
     0.000825_dp,0.002605_dp,0.059987_dp,0.008349_dp],[4,4])

    C_CONTACT(7,2,:,:)=reshape([ &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.000492_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp],[4,4])

  end subroutine initialize_contacts


  subroutine get_CL_periodic(t,C,L)

    real(dp),intent(in)::t

    real(dp),intent(out)::C(NG,NG)
    real(dp),intent(out)::L(NG,NG)

    real(dp)::tm

    integer::day,sh

    tm=modulo(t,7.0_dp)

    day=int(floor(tm))+1

    sh=int(floor( &
         (tm-floor(tm))/SHIFT_DURATION+1.0e-10_dp))+1

    day=max(1,min(DAYS,day))
    sh=max(1,min(SHIFTS,sh))

    C=C_CONTACT(day,sh,:,:)
    L=L_SHIFTS(sh,:,:)

  end subroutine get_CL_periodic

end module contact_mod


module model_mod

  use, intrinsic :: ieee_arithmetic

  use precision_mod
  use settings_mod
  use contact_mod

  implicit none

contains


  function seasonal_factor(t,a) result(s)

    real(dp),intent(in)::t,a

    real(dp)::s

    s=1+a*sin( &
         2.0_dp*PI*t/SEASONAL_PERIOD_DAYS)

  end function seasonal_factor


  subroutine rhs(t,y,alpha,a,C,L,dy,ok)

    real(dp),intent(in)::t
    real(dp),intent(in)::y(NY)
    real(dp),intent(in)::alpha,a
    real(dp),intent(in)::C(NG,NG)
    real(dp),intent(in)::L(NG,NG)

    real(dp),intent(out)::dy(NY)

    logical,intent(out)::ok

    real(dp)::I1(NG)
    real(dp)::I2(NG)
    real(dp)::Sus(NG)

    real(dp)::lam1(NG)
    real(dp)::lam2(NG)

    real(dp)::q1(NG)
    real(dp)::q2(NG)

    real(dp)::seas

    integer::i,j


    I1=y(1:NG)
    I2=y(NG+1:NY)

    Sus=POP-I1-I2

    ok=.false.
    dy=0.0_dp


    if(minval(I1)<-NEG_TOL .or. &
       minval(I2)<-NEG_TOL .or. &
       minval(Sus)<-NEG_TOL) return


    seas=seasonal_factor(t,a)


    !===================================================================
    ! B1 and B2 are diagonal.
    !
    ! These explicit loops are mathematically identical to
    !
    !   (I-L) C B1
    !
    ! and
    !
    !   (I-alpha L) C B2
    !
    ! but avoid repeated small MATMUL operations.
    !===================================================================

    do j=1,NG

      q1(j)=B1(j,j)*I1(j)/POP(j)

      q2(j)=B2(j,j)*I2(j)/POP(j)

    end do


    do i=1,NG

      lam1(i)=0.0_dp
      lam2(i)=0.0_dp

      do j=1,NG

        lam1(i)=lam1(i)+C(i,j)*q1(j)

        lam2(i)=lam2(i)+C(i,j)*q2(j)

      end do

      lam1(i)=seas*(1.0_dp-L(i,i))*lam1(i)

      lam2(i)=seas*(1.0_dp-alpha*L(i,i))*lam2(i)

    end do


    do i=1,NG

      dy(i)=Sus(i)*lam1(i)-GAM1(i)*I1(i)

      dy(NG+i)=Sus(i)*lam2(i)-GAM2(i)*I2(i)

    end do


    if(any(.not.ieee_is_finite(dy))) return

    ok=.true.

  end subroutine rhs


  subroutine integrate(alpha,a,wm,rm,n,ok)

    real(dp),intent(in)::alpha,a

    real(dp),intent(out)::wm(:)
    real(dp),intent(out)::rm(:)

    integer,intent(out)::n

    logical,intent(out)::ok


    real(dp)::y(NY)

    real(dp)::k1(NY)
    real(dp)::k2(NY)
    real(dp)::k3(NY)
    real(dp)::k4(NY)

    real(dp)::yt(NY)

    real(dp)::C(NG,NG)
    real(dp)::L(NG,NG)

    real(dp)::t
    real(dp)::dt

    real(dp)::Itot
    real(dp)::gw

    real(dp)::dy(NY)

    logical::good

    integer::sh


    ok=.false.
    n=0

    y=0.0_dp


    ! Initial infected fractions
    y(1)=5.0e-4_dp*POP(1)
    y(NG+1)=2.0e-4_dp*POP(1)


    t=0.0_dp


    do sh=1,TOTAL_SHIFTS+1


      !=================================================================
      ! Calculate R_eff on the model grid after burn-in
      !=================================================================

      if(t>=BURN_IN_DAYS-1.0e-10_dp)then

        Itot=sum(y)

        if(Itot<=0.0_dp)return


        gw=( &
          sum(GAM1*y(1:NG))+ &
          sum(GAM2*y(NG+1:NY)) &
          )/Itot


        if(gw<=0.0_dp)return


        call get_CL_periodic(t,C,L)

        call rhs(t,y,alpha,a,C,L,dy,good)

        if(.not.good)return


        n=n+1

        wm(n)=(t-BURN_IN_DAYS)/7.0_dp

        rm(n)=1.0_dp+sum(dy)/(Itot*gw)

      end if


      if(sh>TOTAL_SHIFTS)exit


      !=================================================================
      ! Contact/intervention matrices remain fixed over this shift.
      !=================================================================

      call get_CL_periodic(t,C,L)

      dt=SHIFT_DURATION


      ! RK4 stage 1
      call rhs(t,y,alpha,a,C,L,k1,good)

      if(.not.good)return


      ! RK4 stage 2
      yt=y+0.5_dp*dt*k1

      call rhs( &
           t+0.5_dp*dt,yt,alpha,a,C,L,k2,good)

      if(.not.good)return


      ! RK4 stage 3
      yt=y+0.5_dp*dt*k2

      call rhs( &
           t+0.5_dp*dt,yt,alpha,a,C,L,k3,good)

      if(.not.good)return


      ! RK4 stage 4
      yt=y+dt*k3

      call rhs( &
           t+dt,yt,alpha,a,C,L,k4,good)

      if(.not.good)return


      y=y+(dt/6.0_dp)*( &
           k1+2.0_dp*k2+2.0_dp*k3+k4)

      y=max(y,0.0_dp)

      t=t+dt

    end do


    ok=(n>=NDATA)

  end subroutine integrate


  subroutine predict(theta,weeks,pred,ok)

    real(dp),intent(in)::theta(NDIM)
    real(dp),intent(in)::weeks(:)

    real(dp),intent(out)::pred(:)

    logical,intent(out)::ok

    real(dp)::wm(N_EVAL)
    real(dp)::rm(N_EVAL)

    real(dp)::x
    real(dp)::pos
    real(dp)::w

    integer::n
    integer::i
    integer::j


    pred=0.0_dp
    ok=.false.


    call integrate( &
         theta(1),theta(2),wm,rm,n,ok)

    if(.not.ok)return


    !===================================================================
    ! The integration grid is uniform:
    !
    !       Delta t = 1/21 week.
    !
    ! Therefore direct indexing replaces the previous linear search.
    ! This changes only the implementation, not the interpolation.
    !===================================================================

    do i=1,size(weeks)

      x=weeks(i)-theta(4)


      if(x<wm(1) .or. x>wm(n))then

        ok=.false.

        return

      end if


      pos=x/SHIFT_DURATION_WEEKS

      j=int(floor(pos))+1


      if(j<1)j=1


      if(j>=n)then

        j=n-1

        w=1.0_dp

      else

        w=pos-floor(pos)

      end if


      pred(i)=rm(j)+ &
           w*(rm(j+1)-rm(j))

    end do


    ok=.true.

  end subroutine predict

end module model_mod


module fit_mod

  use precision_mod
  use settings_mod
  use data_mod
  use model_mod
  use rng_mod
  use stats_mod

  implicit none

contains


  logical function bounds(x)

    real(dp),intent(in)::x(NDIM)

    bounds= &
      x(1)>=ALPHA_LO .and. x(1)<=ALPHA_HI .and. &
      x(2)>=SEASON_LO .and. x(2)<=SEASON_HI .and. &
      x(3)>=SIGMA_LO .and. x(3)<=SIGMA_HI .and. &
      x(4)>=SHIFT_LO .and. x(4)<=SHIFT_HI

  end function bounds


  function nll(x,y) result(f)

    real(dp),intent(in)::x(NDIM)
    real(dp),intent(in)::y(NDATA)

    real(dp)::f
    real(dp)::p(NDATA)

    real(dp)::res
    real(dp)::ll

    logical::ok

    integer::i


    f=1.0e100_dp


    if(.not.bounds(x))return

    if(x(3)<=0.0_dp)return


    call predict(x,WEEKS_DATA,p,ok)

    if(.not.ok)return


    ll=0.0_dp


    do i=1,NDATA

      res=y(i)-p(i)

      ll=ll-0.5_dp*( &
        (res/x(3))**2+ &
        log(2.0_dp*PI*x(3)**2))

    end do


    f=-ll

  end function nll


  subroutine nelder(x0,y,best,fbest,maxit)

    real(dp),intent(in)::x0(NDIM)
    real(dp),intent(in)::y(NDATA)

    real(dp),intent(out)::best(NDIM)
    real(dp),intent(out)::fbest

    integer,intent(in)::maxit


    real(dp)::x(NDIM,NDIM+1)
    real(dp)::f(NDIM+1)

    real(dp)::c(NDIM)
    real(dp)::xr(NDIM)
    real(dp)::xe(NDIM)
    real(dp)::xc(NDIM)

    real(dp)::fr
    real(dp)::fe
    real(dp)::fc

    real(dp)::scale

    integer::i
    integer::j
    integer::it

    integer::lo
    integer::hi
    integer::nhi


    x(:,1)=x0

    x(:,2)=x0+[0.02_dp,0.0_dp,0.0_dp,0.0_dp]

    x(:,3)=x0+[0.0_dp,0.02_dp,0.0_dp,0.0_dp]

    x(:,4)=x0+[0.0_dp,0.0_dp,0.005_dp,0.0_dp]

    x(:,5)=x0+[0.0_dp,0.0_dp,0.0_dp,0.1_dp]


    do i=1,5

      f(i)=nll(x(:,i),y)

    end do


    do it=1,maxit


      lo=minloc(f,1)

      hi=maxloc(f,1)

      nhi=1


      do j=1,5

        if(j/=hi .and. &
           (nhi==hi .or. f(j)>f(nhi)))then

          nhi=j

        end if

      end do


      scale=0.0_dp


      do j=1,5

        scale=max( &
          scale, &
          sum(abs(x(:,j)-x(:,lo))))

      end do


      if(scale<1.0e-8_dp)exit


      c=0.0_dp


      do j=1,5

        if(j/=hi)c=c+x(:,j)

      end do


      c=c/4.0_dp


      ! Reflection
      xr=c+(c-x(:,hi))

      fr=nll(xr,y)


      if(fr<f(lo))then


        ! Expansion
        xe=c+2.0_dp*(xr-c)

        fe=nll(xe,y)


        if(fe<fr)then

          x(:,hi)=xe
          f(hi)=fe

        else

          x(:,hi)=xr
          f(hi)=fr

        end if


      else if(fr<f(nhi))then


        x(:,hi)=xr
        f(hi)=fr


      else


        if(fr<f(hi))then

          xc=c+0.5_dp*(xr-c)

        else

          xc=c+0.5_dp*(x(:,hi)-c)

        end if


        fc=nll(xc,y)


        if(fc<f(hi))then

          x(:,hi)=xc
          f(hi)=fc

        else


          do j=1,5

            if(j/=lo)then

              x(:,j)=x(:,lo)+ &
                   0.5_dp*(x(:,j)-x(:,lo))

              f(j)=nll(x(:,j),y)

            end if

          end do

        end if

      end if

    end do


    lo=minloc(f,1)

    best=x(:,lo)

    fbest=f(lo)

  end subroutine nelder


  subroutine mle(start,best,fbest,seed)

    real(dp),intent(in)::start(NDIM)

    real(dp),intent(out)::best(NDIM)
    real(dp),intent(out)::fbest

    integer,intent(inout)::seed


    real(dp)::x0(NDIM)
    real(dp)::xb(NDIM)
    real(dp)::fb

    integer::s


    best=start
    fbest=1.0e100_dp


    do s=1,N_MLE_STARTS


      if(s==1)then

        x0=start

      else


        x0(1)= &
          ALPHA_LO+ &
          0.05_dp*(ALPHA_HI-ALPHA_LO)+ &
          0.90_dp*(ALPHA_HI-ALPHA_LO)*urand(seed)


        x0(2)= &
          SEASON_LO+ &
          0.05_dp*(SEASON_HI-SEASON_LO)+ &
          0.90_dp*(SEASON_HI-SEASON_LO)*urand(seed)


        x0(3)= &
          exp( &
          log(SIGMA_LO*2.0_dp)+ &
          (log(SIGMA_HI)-log(SIGMA_LO*2.0_dp))* &
          urand(seed))


        x0(4)= &
          SHIFT_LO+ &
          0.02_dp*(SHIFT_HI-SHIFT_LO)+ &
          0.96_dp*(SHIFT_HI-SHIFT_LO)*urand(seed)

      end if


      call nelder( &
           x0,R_DATA,xb,fb,800)


      if(fb<fbest)then

        fbest=fb
        best=xb

      end if

    end do

  end subroutine mle


  function logpost(x,y) result(lp)

    real(dp),intent(in)::x(NDIM)
    real(dp),intent(in)::y(NDATA)

    real(dp)::lp


    if(.not.bounds(x))then

      lp=-1.0e100_dp

      return

    end if


    ! Gaussian likelihood + Gaussian prior on shift.
    lp= &
      -0.5_dp*(x(4)/SHIFT_PRIOR_SIGMA)**2 &
      -log(SHIFT_PRIOR_SIGMA*sqrt(2.0_dp*PI)) &
      -nll(x,y)

  end function logpost


  subroutine mcmc( &
       center, &
       samples, &
       accept, &
       seed, &
       walker_acceptance, &
       rhat)

    real(dp),intent(in)::center(NDIM)

    real(dp),allocatable,intent(out)::samples(:,:)

    real(dp),intent(out)::accept

    real(dp),intent(out)::walker_acceptance(NWALKERS)

    real(dp),intent(out)::rhat(NDIM)

    integer,intent(inout)::seed


    !===================================================================
    ! Walker-major state was changed to walker-contiguous columns:
    !
    !       ch(parameter,walker)
    !
    ! This avoids repeated non-contiguous array sections such as ch(i,:)
    ! in calls to logpost.
    !===================================================================

    real(dp)::ch(NDIM,NWALKERS)

    real(dp)::lp(NWALKERS)

    real(dp)::p(NDIM)

    real(dp)::newlp

    real(dp)::z
    real(dp)::u


    real(dp),allocatable::saved_chain(:,:,:)

    integer,allocatable::walker_seed(:)


    integer::i
    integer::k
    integer::g
    integer::partner

    integer::nkeep
    integer::idx

    integer::nprop
    integer::nacc


    integer::walker_prop(NWALKERS)
    integer::walker_acc(NWALKERS)


    !===================================================================
    ! Independent RNG stream for each walker.
    !===================================================================

    allocate(walker_seed(NWALKERS))


    do i=1,NWALKERS

      walker_seed(i)= &
        abs(int( &
        2147483000.0_dp*urand(seed)))+1

    end do


    !===================================================================
    ! Initial walker ensemble
    !===================================================================

    do i=1,NWALKERS

      ch(:,i)=center


      ch(1,i)=min( &
        max( &
          center(1)+0.01_dp*nrand(seed), &
          ALPHA_LO+1.0e-6_dp), &
        ALPHA_HI-1.0e-6_dp)


      ch(2,i)=min( &
        max( &
          center(2)+0.01_dp*nrand(seed), &
          SEASON_LO+1.0e-6_dp), &
        SEASON_HI-1.0e-6_dp)


      ch(3,i)=min( &
        max( &
          center(3)*exp(0.05_dp*nrand(seed)), &
          SIGMA_LO*1.01_dp), &
        SIGMA_HI*0.99_dp)


      ch(4,i)=min( &
        max( &
          center(4)+0.05_dp*nrand(seed), &
          SHIFT_LO+1.0e-6_dp), &
        SHIFT_HI-1.0e-6_dp)


      lp(i)=logpost(ch(:,i),R_DATA)

    end do


    !===================================================================
    ! Storage
    !===================================================================

    nkeep=NWALKERS*((NSTEPS-BURN_MCMC)/THIN)


    allocate(samples(nkeep,NDIM))


    allocate(saved_chain( &
         (NSTEPS-BURN_MCMC)/THIN, &
         NWALKERS, &
         NDIM))


    idx=0

    nprop=0
    nacc=0

    walker_prop=0
    walker_acc=0


    !===================================================================
    ! AFFINE-INVARIANT STRETCH-MOVE MCMC
    !
    ! The Bayesian sampler itself is unchanged.
    !
    ! Two complementary groups are updated sequentially.
    !===================================================================

    do k=1,NSTEPS


      do g=1,2


        do i=g,NWALKERS,2


          !-------------------------------------------------------------
          ! Select a walker from the complementary group.
          !-------------------------------------------------------------

          partner=1+int( &
             real(NWALKERS,dp)* &
             urand(walker_seed(i)))


          do while( &
             partner==i .or. &
             mod(partner,2)==g)

            partner=1+int( &
               real(NWALKERS,dp)* &
               urand(walker_seed(i)))

          end do


          !-------------------------------------------------------------
          ! Goodman-Weare stretch variable
          !-------------------------------------------------------------

          z=((A_STRETCH-1.0_dp)* &
             urand(walker_seed(i))+1.0_dp)**2/A_STRETCH


          !-------------------------------------------------------------
          ! Proposed position
          !-------------------------------------------------------------

          p=ch(:,partner)+ &
             z*(ch(:,i)-ch(:,partner))


          !-------------------------------------------------------------
          ! Posterior
          !-------------------------------------------------------------

          newlp=logpost(p,R_DATA)


          nprop=nprop+1
          walker_prop(i)=walker_prop(i)+1


          if(newlp>-1.0e99_dp)then


            u=urand(walker_seed(i))


            if( &
              log(u) < &
              real(NDIM-1,dp)*log(z)+ &
              newlp-lp(i))then


              ch(:,i)=p

              lp(i)=newlp

              nacc=nacc+1

              walker_acc(i)=walker_acc(i)+1

            end if

          end if


        end do

      end do


      !=================================================================
      ! Save thinned chains
      !=================================================================

      if( &
        k>BURN_MCMC .and. &
        mod(k-BURN_MCMC,THIN)==0)then


        idx=idx+1


        do i=1,NWALKERS

          samples( &
            (idx-1)*NWALKERS+i,:)=ch(:,i)


          saved_chain(idx,i,:)=ch(:,i)

        end do

      end if


    end do


    !===================================================================
    ! Return a valid scalar seed to preserve caller interface.
    !===================================================================

    seed=max( &
      1, &
      int(mod( &
      real(sum(walker_seed),dp), &
      2147483000.0_dp)))


    !===================================================================
    ! Acceptance
    !===================================================================

    accept= &
      real(nacc,dp)/real(max(1,nprop),dp)


    do i=1,NWALKERS

      walker_acceptance(i)= &
        real(walker_acc(i),dp)/ &
        real(max(1,walker_prop(i)),dp)

    end do


    !===================================================================
    ! Split-R-hat
    !===================================================================

    call compute_split_rhat(saved_chain,rhat)


    deallocate(saved_chain)
    deallocate(walker_seed)

  end subroutine mcmc


  subroutine compute_split_rhat(samples,rhat)

    real(dp),intent(in)::samples(:,:,:)

    real(dp),intent(out)::rhat(NDIM)

    integer::ns
    integer::nw
    integer::p

    integer::m
    integer::ii
    integer::jj
    integer::kk

    real(dp)::grand_mean
    real(dp)::w
    real(dp)::b
    real(dp)::var_plus

    real(dp),allocatable::chain_mean(:)
    real(dp),allocatable::chain_var(:)


    ns=size(samples,1)
    nw=size(samples,2)
    p=size(samples,3)


    if(ns<4 .or. nw<2)then

      rhat=huge(1.0_dp)

      return

    end if


    m=2*nw


    allocate(chain_mean(m))
    allocate(chain_var(m))


    do kk=1,min(p,NDIM)


      do ii=1,nw


        !---------------------------------------------------------------
        ! First half
        !---------------------------------------------------------------

        jj=2*ii-1


        chain_mean(jj)= &
          sum(samples(1:ns/2,ii,kk))/ &
          real(ns/2,dp)


        chain_var(jj)= &
          sum( &
          (samples(1:ns/2,ii,kk)-chain_mean(jj))**2)/ &
          real(ns/2-1,dp)


        !---------------------------------------------------------------
        ! Second half
        !---------------------------------------------------------------

        jj=2*ii


        chain_mean(jj)= &
          sum(samples(ns/2+1:2*(ns/2),ii,kk))/ &
          real(ns/2,dp)


        chain_var(jj)= &
          sum( &
          (samples(ns/2+1:2*(ns/2),ii,kk)-chain_mean(jj))**2)/ &
          real(ns/2-1,dp)


      end do


      grand_mean= &
        sum(chain_mean)/real(m,dp)


      w= &
        sum(chain_var)/real(m,dp)


      b= &
        real(ns/2,dp)* &
        sum((chain_mean-grand_mean)**2)/ &
        real(m-1,dp)


      if(w>0.0_dp)then


        var_plus= &
          (real(ns/2-1,dp)/real(ns/2,dp))*w + &
          b/real(ns/2,dp)


        rhat(kk)= &
          sqrt(max(var_plus/w,0.0_dp))


      else


        if(b<=0.0_dp)then

          rhat(kk)=1.0_dp

        else

          rhat(kk)=huge(1.0_dp)

        end if

      end if


    end do


    if(p<NDIM)then

      rhat(p+1:NDIM)=huge(1.0_dp)

    end if


    deallocate(chain_mean)
    deallocate(chain_var)

  end subroutine compute_split_rhat


  subroutine write_case(seed)

    integer,intent(inout)::seed


    real(dp)::start(NDIM)
    real(dp)::best(NDIM)

    real(dp)::fbest

    real(dp)::pred(NDATA)
    real(dp)::mle_curve(NDATA)

    real(dp)::rmse
    real(dp)::corr

    real(dp)::acc

    real(dp)::rhat(NDIM)
    real(dp)::max_rhat

    real(dp)::walker_acceptance(NWALKERS)


    logical::convergence_ok
    logical::ok


    real(dp),allocatable::samp(:,:)
    real(dp),allocatable::curves(:,:)

    real(dp),allocatable::pc(:)
    real(dp),allocatable::pe(:)


    real(dp)::q025
    real(dp)::q50
    real(dp)::q975
    real(dp)::meanv


    integer::i
    integer::j
    integer::nv
    integer::unit

    character(len=80)::fn

    integer::posterior_seed


    !===================================================================
    ! Initial MLE guess
    !===================================================================

    start=[ &
      0.8_dp, &
      0.35_dp, &
      0.15_dp, &
      1.8_dp]


    !===================================================================
    ! MLE
    !===================================================================

    call mle( &
         start, &
         best, &
         fbest, &
         seed)


    call predict( &
         best, &
         WEEKS_DATA, &
         pred, &
         ok)


    if(.not.ok)stop 'MLE prediction failed'


    ! IMPORTANT:
    ! Preserve the MLE curve before posterior predictions overwrite pred.

    mle_curve=pred


    rmse=sqrt( &
      sum((R_DATA-mle_curve)**2)/ &
      real(NDATA,dp))


    corr=pearson( &
      R_DATA, &
      mle_curve)


    write(*,'(A)')' '
    write(*,'(A)')'==============================================================='
    write(*,'(A)')' MLE RESULTS: PERIODIC MODEL'
    write(*,'(A)')'==============================================================='


    write(*,'(A,F14.7)') &
      ' alpha         = ',best(1)


    write(*,'(A,F14.7)') &
      ' season_amp    = ',best(2)


    write(*,'(A,F14.7)') &
      ' sigma         = ',best(3)


    write(*,'(A,F14.7)') &
      ' shift_weeks   = ',best(4)


    write(*,'(A,F14.7)') &
      ' logLikelihood = ',-fbest


    write(*,'(A,F14.7)') &
      ' RMSE          = ',rmse


    write(*,'(A,F14.7)') &
      ' Pearson r     = ',corr


    write(*,'(A)')'==============================================================='
    write(*,'(A)')' Starting MCMC...'
    write(*,'(A)')'==============================================================='


    !===================================================================
    ! MCMC
    !===================================================================

    call mcmc( &
         best, &
         samp, &
         acc, &
         seed, &
         walker_acceptance, &
         rhat)


    max_rhat=maxval(rhat)

    convergence_ok=(max_rhat<=1.05_dp)


    write(*,'(A,F10.6)') &
      'MCMC acceptance = ',acc


    write(*,'(A,4(1X,F10.6))') &
      'R-hat = ',rhat


    write(*,'(A,F10.6)') &
      'Max R-hat = ',max_rhat


    if(convergence_ok)then

      write(*,'(A)') &
        'MCMC convergence = GOOD'

    else

      write(*,'(A)') &
        'MCMC convergence = NOT CONVERGED'

    end if


    !===================================================================
    ! Posterior trajectory evaluation
    !===================================================================

    allocate( &
      curves( &
        NDATA, &
        N_POSTERIOR_MODELS))


    allocate(pc(N_POSTERIOR_MODELS))
    allocate(pe(N_POSTERIOR_MODELS))


    nv=0

    posterior_seed=314159


    do i=1,N_POSTERIOR_MODELS


      j=1+int( &
        real(size(samp,1),dp)* &
        urand(posterior_seed))


      j=max(1,min(size(samp,1),j))


      call predict( &
        samp(j,:), &
        WEEKS_DATA, &
        pred, &
        ok)


      if(ok)then


        nv=nv+1


        curves(:,nv)=pred


        pc(nv)=pearson( &
          R_DATA, &
          pred)


        pe(nv)=sqrt( &
          sum((R_DATA-pred)**2)/ &
          real(NDATA,dp))

      end if


    end do


    if(nv<MIN_VALID_POSTERIOR)then

      stop 'Too few valid posterior trajectories'

    end if


    !===================================================================
    ! MLE summary
    !===================================================================

    fn='Influenza_MLE_Summary.dat'


    open( &
      newunit=unit, &
      file=trim(fn), &
      status='replace')


    write(unit,'(A)') &
      '# alpha season_amp sigma shift_weeks logLikelihood RMSE Pearson_r'


    write(unit,'(7ES24.14)') &
      best(1), &
      best(2), &
      best(3), &
      best(4), &
      -fbest, &
      rmse, &
      corr


    close(unit)


    !===================================================================
    ! Fit results
    !
    ! MLE curve is explicitly stored in mle_curve.
    !===================================================================

    fn='Influenza_Fit_Results.csv'


    open( &
      newunit=unit, &
      file=trim(fn), &
      status='replace')


    write(unit,'(A)') &
      'week,observed,mle,posterior_median,posterior_lower,posterior_upper'


    do i=1,NDATA


      q50=percentile( &
        curves(i,1:nv), &
        0.5_dp)


      q025=percentile( &
        curves(i,1:nv), &
        0.025_dp)


      q975=percentile( &
        curves(i,1:nv), &
        0.975_dp)


      write(unit,'(6ES24.14)') &
        WEEKS_DATA(i), &
        R_DATA(i), &
        mle_curve(i), &
        q50, &
        q025, &
        q975

    end do


    close(unit)


    !===================================================================
    ! Posterior samples
    !===================================================================

    fn='Influenza_Posterior_Samples.dat'


    open( &
      newunit=unit, &
      file=trim(fn), &
      status='replace')


    write(unit,'(A)') &
      '# alpha season_amp sigma shift_weeks'


    do i=1,size(samp,1)

      write(unit,'(4ES24.14)') &
        samp(i,:)

    end do


    close(unit)


    !===================================================================
    ! Posterior performance
    !===================================================================

    fn='Influenza_Posterior_Performance.dat'


    open( &
      newunit=unit, &
      file=trim(fn), &
      status='replace')


    write(unit,'(A)') &
      '# Pearson_r RMSE'


    do i=1,nv

      write(unit,'(2ES24.14)') &
        pc(i), &
        pe(i)

    end do


    close(unit)


    !===================================================================
    ! MCMC diagnostics
    !===================================================================

    fn='Influenza_MCMC_Diagnostics.dat'


    open( &
      newunit=unit, &
      file=trim(fn), &
      status='replace')


    write(unit,'(A)') &
      '# overall_acceptance posterior_samples max_rhat convergence'


    write(unit,'(3ES24.14,1X,L1)') &
      acc, &
      real(size(samp,1),dp), &
      max_rhat, &
      convergence_ok


    write(unit,'(A)') &
      '# walker acceptance'


    do i=1,NWALKERS

      write(unit,'(I6,1X,ES24.14)') &
        i, &
        walker_acceptance(i)

    end do


    write(unit,'(A)') &
      '# split_Rhat alpha season_amp sigma shift_weeks'


    write(unit,'(4ES24.14)') &
      rhat


    close(unit)


    !===================================================================
    ! Profile likelihoods
    !===================================================================

    call write_profile( &
      'Influenza_Profile_alpha.dat', &
      best, &
      1)


    call write_profile( &
      'Influenza_Profile_season_amp.dat', &
      best, &
      2)


    call write_profile( &
      'Influenza_Profile_sigma.dat', &
      best, &
      3)


    call write_profile( &
      'Influenza_Profile_shift_weeks.dat', &
      best, &
      4)


    !===================================================================
    ! Bayesian summary
    !===================================================================

    fn='Influenza_Bayesian_Results_Summary.txt'


    open( &
      newunit=unit, &
      file=trim(fn), &
      status='replace')


    write(unit,'(A)') &
      'INFLUENZA R_eff CALIBRATION - PERIODIC MODEL'


    write(unit,'(A)') &
      '21 shift-specific contact matrices and 3 shift-specific intervention matrices'


    write(unit,'(A)') &
      'Force of infection: (I-L(t)) C(t) B1 and (I-alpha L(t)) C(t) B2'


    write(unit,'(A,F14.7)') &
      'MLE alpha = ',best(1)


    write(unit,'(A,F14.7)') &
      'MLE season_amp = ',best(2)


    write(unit,'(A,F14.7)') &
      'MLE sigma = ',best(3)


    write(unit,'(A,F14.7)') &
      'MLE shift_weeks = ',best(4)


    write(unit,'(A,F14.7)') &
      'MLE logLikelihood = ',-fbest


    write(unit,'(A,F14.7)') &
      'MLE RMSE = ',rmse


    write(unit,'(A,F14.7)') &
      'MLE Pearson_r = ',corr


    write(unit,'(A,F14.7)') &
      'MCMC acceptance = ',acc


    write(unit,'(A,F14.7)') &
      'MCMC max R-hat = ',max_rhat


    if(convergence_ok)then

      write(unit,'(A)') &
        'MCMC convergence = GOOD'

    else

      write(unit,'(A)') &
        'MCMC convergence = NOT_CONVERGED'

    end if


    !===================================================================
    ! Posterior parameter summaries
    !===================================================================

    do i=1,NDIM


      meanv= &
        sum(samp(:,i))/ &
        real(size(samp,1),dp)


      q025=percentile( &
        samp(:,i), &
        0.025_dp)


      q50=percentile( &
        samp(:,i), &
        0.5_dp)


      q975=percentile( &
        samp(:,i), &
        0.975_dp)


      write(unit, &
        '(A,I2,A,F14.7,A,F14.7,A,F14.7,A,F14.7)') &
        'parameter ', &
        i, &
        ' mean=', &
        meanv, &
        ' median=', &
        q50, &
        ' CrI_low=', &
        q025, &
        ' CrI_high=', &
        q975

    end do


    write(unit,'(A,F14.7,A,F14.7)') &
      'Posterior RMSE CrI = ', &
      percentile(pe(1:nv),0.025_dp), &
      ' to ', &
      percentile(pe(1:nv),0.975_dp)


    write(unit,'(A,F14.7,A,F14.7)') &
      'Posterior Pearson r CrI = ', &
      percentile(pc(1:nv),0.025_dp), &
      ' to ', &
      percentile(pc(1:nv),0.975_dp)


    close(unit)


    deallocate(samp)
    deallocate(curves)
    deallocate(pc)
    deallocate(pe)

  end subroutine write_case


  subroutine write_profile(base,best,which)

    character(len=*),intent(in)::base

    real(dp),intent(in)::best(NDIM)

    integer,intent(in)::which


    character(len=120)::out

    real(dp)::x(NDIM)
    real(dp)::f
    real(dp)::fmin

    real(dp)::lo
    real(dp)::hi
    real(dp)::v

    integer::i
    integer::u


    select case(which)

    case(1)

      lo=ALPHA_LO
      hi=ALPHA_HI

    case(2)

      lo=SEASON_LO
      hi=SEASON_HI

    case(3)

      lo=SIGMA_LO
      hi=SIGMA_HI

    case default

      lo=SHIFT_LO
      hi=SHIFT_HI

    end select


    out=trim(base)


    fmin=nll(best,R_DATA)


    open( &
      newunit=u, &
      file=trim(out), &
      status='replace')


    write(u,'(A)') &
      '# parameter value negative_loglikelihood delta_nll'


    do i=1,N_PROFILE


      v=lo+(hi-lo)* &
        real(i-1,dp)/ &
        real(N_PROFILE-1,dp)


      x=best

      x(which)=v


      f=nll(x,R_DATA)


      write(u,'(3ES24.14)') &
        v, &
        f, &
        f-fmin

    end do


    close(u)

  end subroutine write_profile

end module fit_mod


program influenza_periodic

  use precision_mod
  use settings_mod
  use contact_mod
  use fit_mod

  implicit none

  integer::seed


  call initialize_contacts()


  seed=42


  write(*,'(A)') &
    '==============================================================='


  write(*,'(A)') &
    ' INFLUENZA R_eff: PERIODIC CALIBRATION + MCMC DIAGNOSTICS'


  write(*,'(A)') &
    '==============================================================='


  write(*,'(A,I6)') &
    ' Data points = ',NDATA


  write(*,'(A,I6)') &
    ' Walkers    = ',NWALKERS


  write(*,'(A,I7)') &
    ' MCMC steps = ',NSTEPS


  write(*,'(A,I7)') &
    ' Burn-in    = ',BURN_MCMC


  write(*,'(A,I6)') &
    ' Thin       = ',THIN


  write(*,'(A)') &
    ' Bootstrap  = REMOVED'


  write(*,'(A)') &
    ' Model      = PERIODIC 21-contact / 3-intervention shifts'


  write(*,'(A)') &
    '==============================================================='


  call write_case(seed)


  write(*,'(A)') &
    'Completed: periodic Influenza calibration and diagnostics.'

end program influenza_periodic
