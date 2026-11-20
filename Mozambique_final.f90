!===============================================================================
! COVID-19 RURAL vs URBAN - COMPLETE BAYESIAN CALIBRATION FRAMEWORK
! Corrected, self-contained Fortran 2008 implementation of the supplied model.
!
! Preserved model/data/settings from the supplied COVID source:
!   NDIM=5: alpha, season, A_week, rho, sigma
!   N_AGE=10, NDATA=30, NWALKERS=64, NSTEPS=30000,
!   BURN_MCMC=35000, THIN=20, N_POSTERIOR_MODELS=2000
!
! Corrections:
!   * Correct parameter-bound construction
!   * MLE optimized using the same Gaussian likelihood as MCMC
!   * sigma profiled analytically during MLE optimization
!   * correct parameter-specific R-hat assignment
!   * corrected chain-printing variable
!
! Other model/data/settings/output filenames are unchanged.
!
! Output files:
!   RURAL_COVID_Fit_Results.csv
!   RURAL_COVID_Posterior_Samples.dat
!   RURAL_COVID_MCMC_Diagnostics.dat
!   RURAL_COVID_Parameter_Summary.dat
!   RURAL_COVID_Summary.txt
!
!   URBAN_COVID_Fit_Results.csv
!   URBAN_COVID_Posterior_Samples.dat
!   URBAN_COVID_MCMC_Diagnostics.dat
!   URBAN_COVID_Parameter_Summary.dat
!   URBAN_COVID_Summary.txt
!
!   COVID_Rural_Urban_Combined_Results.csv
!
! No external libraries are required.
!===============================================================================


module precision_mod
  implicit none

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: i8 = selected_int_kind(18)

  real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp

end module precision_mod


!===============================================================================
! SETTINGS
!===============================================================================

module settings_mod

  use precision_mod
  implicit none

  integer, parameter :: N_AGE = 10
  integer, parameter :: NY = 20

  real(dp), parameter :: DAYS_PER_WEEK = 7.0_dp
  real(dp), parameter :: SIM_END_DAYS = 600.0_dp

  integer, parameter :: N_T_EVAL = 600

  integer, parameter :: NDIM = 5
  integer, parameter :: NWALKERS = 64
  integer, parameter :: NCHAINS = 4

  integer, parameter :: NSTEPS = 120000
  integer, parameter :: BURN_MCMC = 35000
  integer, parameter :: THIN = 20

  integer, parameter :: N_POSTERIOR_MODELS = 2000

  integer, parameter :: NDATA = 30

  !---------------------------------------------------------------------------
  ! Parameter bounds
  !---------------------------------------------------------------------------

  real(dp), parameter :: ALPHA_LO = 0.2_dp
  real(dp), parameter :: ALPHA_HI = 1.0_dp

  real(dp), parameter :: SEASON_LO = 0.01_dp
  real(dp), parameter :: SEASON_HI = 0.9_dp

  real(dp), parameter :: A_WEEK_LO = 0.01_dp
  real(dp), parameter :: A_WEEK_HI = 0.9_dp

  real(dp), parameter :: RHO_LO = 10.0_dp
  real(dp), parameter :: RHO_HI = 2000.0_dp

  real(dp), parameter :: SIGMA_LO = 1.0e-5_dp
  real(dp), parameter :: SIGMA_HI = 1.0_dp

  integer, parameter :: MLE_MAXITER = 25
  integer, parameter :: MLE_POPSIZE = 10

  !---------------------------------------------------------------------------
  ! Population
  !---------------------------------------------------------------------------

  real(dp), parameter :: POP_RAW(10) = [ &
    103885.0_dp*0.02_dp, &
    103885.0_dp*0.02_dp, &
    103885.0_dp*0.15_dp, &
    103885.0_dp*0.25_dp, &
    103885.0_dp*0.56_dp, &
    103885.0_dp*0.30_dp, &
    103885.0_dp*0.70_dp, &
    226409.0_dp*0.70_dp, &
    (226409.0_dp*0.30_dp)+254504.0_dp, &
    173257.0_dp ]

  real(dp), parameter :: POP(10) = POP_RAW/sum(POP_RAW)

  !---------------------------------------------------------------------------
  ! Transmission / recovery parameters
  !---------------------------------------------------------------------------

  real(dp), parameter :: BETA1_BASE(10) = [ &
    0.2_dp, 0.3_dp, 0.23_dp, 0.38_dp, 0.33_dp, &
    0.36_dp, 0.35_dp, 0.25_dp, 0.36_dp, 0.35_dp ]

  real(dp), parameter :: GAMMA1(10) = [ &
    0.15_dp, 0.18_dp, 0.2_dp, 0.22_dp, 0.25_dp, &
    0.3_dp, 0.28_dp, 0.24_dp, 0.2_dp, 0.18_dp ]

  real(dp), parameter :: GAMMA2(10) = 0.85_dp*GAMMA1

  real(dp), parameter :: L_DIAG(10) = [ &
    0.5_dp, 0.4_dp, 0.6_dp, 0.8_dp, 0.3_dp, &
    0.4_dp, 0.5_dp, 0.38_dp, 0.52_dp, 0.62_dp ]

end module settings_mod


!===============================================================================
! DATA
!===============================================================================

module data_mod

  use precision_mod
  use settings_mod

  implicit none

  real(dp), parameter :: WEEKS_DATA(NDATA) = [ &
    50.49_dp,51.54_dp,52.59_dp,53.53_dp,54.52_dp,55.51_dp,56.56_dp,57.55_dp,58.54_dp,59.53_dp, &
    60.52_dp,61.51_dp,62.62_dp,63.55_dp,64.48_dp,65.59_dp,66.57_dp,67.49_dp,68.47_dp,69.52_dp, &
    70.50_dp,71.50_dp,72.50_dp,73.50_dp,74.61_dp,75.60_dp,76.54_dp,77.59_dp,78.58_dp,79.63_dp ]

  real(dp), parameter :: Y_DATA(NDATA) = [ &
    13688.0_dp,11987.0_dp,10797.0_dp,10797.0_dp,6376.0_dp,8331.0_dp,8502.0_dp,8587.0_dp,6631.0_dp,8842.0_dp, &
    8672.0_dp,7566.0_dp,7311.0_dp,7311.0_dp,8417.0_dp,12753.0_dp,17854.0_dp,27801.0_dp,30352.0_dp,32222.0_dp, &
    39959.0_dp,35538.0_dp,29076.0_dp,22615.0_dp,20064.0_dp,17089.0_dp,13688.0_dp,12242.0_dp,10797.0_dp,9437.0_dp ]

  real(dp), parameter :: Y_DATA_NORM(NDATA) = Y_DATA/maxval(Y_DATA)

end module data_mod


!===============================================================================
! RANDOM NUMBER GENERATOR
!===============================================================================

module rng_mod

  use precision_mod
  implicit none

contains

  subroutine seed_rng(seed)

    integer, intent(in) :: seed

    integer :: n, i
    integer, allocatable :: s(:)

    call random_seed(size=n)

    allocate(s(n))

    do i=1,n
      s(i)=mod(seed+104729*i,2147483647)
    end do

    call random_seed(put=s)

    deallocate(s)

  end subroutine seed_rng


  function urand(seed) result(u)

    integer, intent(inout) :: seed

    real(dp) :: u

    integer(i8) :: a, m, s

    a=1103515245_i8
    m=2147483647_i8

    s=mod(a*int(seed,i8)+12345_i8,m)

    seed=int(s)

    u=max( &
      min(real(seed,dp)/real(m,dp),1.0_dp-1.0e-15_dp), &
      1.0e-15_dp )

  end function urand


  function grand(seed) result(z)

    integer, intent(inout) :: seed

    real(dp) :: z
    real(dp) :: u1,u2

    u1=urand(seed)
    u2=urand(seed)

    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*PI*u2)

  end function grand

end module rng_mod


!===============================================================================
! CONTACT MATRICES
!===============================================================================

module contact_mod

  use precision_mod
  use settings_mod

  implicit none

  real(dp), save :: RURAL_MATRIX(N_AGE,N_AGE)
  real(dp), save :: URBAN_MATRIX(N_AGE,N_AGE)

contains

  subroutine initialize_contacts()

    RURAL_MATRIX=reshape([ &
      0.2_dp,0.2_dp,0.3_dp,0.3_dp,0.6_dp,0.6_dp,0.5_dp,0.9_dp,1.3_dp,0.9_dp, &
      0.5_dp,0.5_dp,0.6_dp,0.5_dp,1.3_dp,1.4_dp,1.1_dp,1.8_dp,3.0_dp,0.8_dp, &
      0.5_dp,0.6_dp,0.6_dp,0.5_dp,1.1_dp,1.1_dp,1.1_dp,2.1_dp,1.1_dp,0.6_dp, &
      0.5_dp,0.7_dp,0.6_dp,0.5_dp,1.6_dp,1.5_dp,1.7_dp,0.4_dp,0.5_dp,0.3_dp, &
      0.4_dp,0.4_dp,0.5_dp,0.8_dp,3.9_dp,1.6_dp,1.0_dp,0.7_dp,1.1_dp,0.4_dp, &
      1.1_dp,1.3_dp,1.1_dp,1.7_dp,5.5_dp,3.4_dp,1.3_dp,0.6_dp,1.0_dp,0.5_dp, &
      0.6_dp,0.8_dp,1.2_dp,3.7_dp,0.6_dp,0.4_dp,0.3_dp,0.2_dp,0.2_dp,0.1_dp, &
      0.4_dp,0.5_dp,1.3_dp,0.5_dp,0.1_dp,0.2_dp,0.1_dp,0.2_dp,0.1_dp,0.1_dp, &
      0.1_dp,0.0_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp, &
      0.0_dp,0.0_dp,0.1_dp,0.1_dp,0.0_dp,0.0_dp,0.1_dp,0.0_dp,0.1_dp,0.0_dp ], &
      [N_AGE,N_AGE])


    URBAN_MATRIX=reshape([ &
      0.2_dp,0.1_dp,0.2_dp,0.2_dp,0.3_dp,0.5_dp,0.5_dp,0.6_dp,0.7_dp,0.8_dp, &
      0.2_dp,0.2_dp,0.3_dp,0.3_dp,0.7_dp,1.0_dp,1.0_dp,1.2_dp,1.7_dp,0.5_dp, &
      0.4_dp,0.6_dp,0.5_dp,0.6_dp,0.9_dp,1.0_dp,0.9_dp,1.8_dp,0.8_dp,0.4_dp, &
      0.4_dp,0.6_dp,0.7_dp,0.7_dp,0.8_dp,1.3_dp,1.0_dp,0.4_dp,0.7_dp,0.4_dp, &
      0.3_dp,0.3_dp,0.4_dp,0.4_dp,2.4_dp,1.2_dp,0.7_dp,0.6_dp,0.7_dp,0.3_dp, &
      0.4_dp,0.5_dp,0.5_dp,0.9_dp,4.3_dp,2.0_dp,0.4_dp,0.5_dp,0.5_dp,0.2_dp, &
      0.1_dp,0.3_dp,0.8_dp,2.1_dp,0.3_dp,0.1_dp,0.2_dp,0.2_dp,0.1_dp,0.1_dp, &
      0.1_dp,0.3_dp,1.2_dp,0.2_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.2_dp, &
      0.1_dp,0.1_dp,0.1_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.1_dp,0.0_dp,0.0_dp, &
      0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.1_dp ], &
      [N_AGE,N_AGE])


    call spectral_normalize(RURAL_MATRIX)
    call spectral_normalize(URBAN_MATRIX)

  end subroutine initialize_contacts


  subroutine spectral_normalize(C)

    real(dp), intent(inout) :: C(N_AGE,N_AGE)

    real(dp) :: v(N_AGE),vn(N_AGE),eig
    integer :: iter

    v=1.0_dp/real(N_AGE,dp)
    eig=0.0_dp

    do iter=1,100

      vn=matmul(C,v)

      eig=maxval(vn)

      if(eig<=0.0_dp)return

      vn=vn/eig

      if(maxval(abs(vn-v))<1.0e-6_dp)exit

      v=vn

    end do

    if(eig>0.0_dp)C=C/eig

  end subroutine spectral_normalize

end module contact_mod


!===============================================================================
! MODEL
!===============================================================================

module model_mod

  use precision_mod
  use settings_mod
  use data_mod

  implicit none

contains


  subroutine rhs(t,y,alpha,season,A_week,rho,C,dy)

    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(NY)

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: season
    real(dp), intent(in) :: A_week
    real(dp), intent(in) :: rho

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp), intent(out) :: dy(NY)

    real(dp) :: I1(N_AGE)
    real(dp) :: I2(N_AGE)
    real(dp) :: S(N_AGE)

    real(dp) :: Ct(N_AGE,N_AGE)

    real(dp) :: lam1(N_AGE)
    real(dp) :: lam2(N_AGE)

    real(dp) :: seasonal
    real(dp) :: weekly

    real(dp) :: control1(N_AGE)
    real(dp) :: control2(N_AGE)

    real(dp) :: tmp1(N_AGE)
    real(dp) :: tmp2(N_AGE)

    integer :: j

    I1=max(y(1:N_AGE),0.0_dp)

    I2=max(y(N_AGE+1:NY),0.0_dp)

    S=max(POP-I1-I2,1.0e-9_dp)

    seasonal=1.0_dp+season*sin(2.0_dp*PI*t/365.0_dp)

    weekly=1.0_dp+A_week*cos(2.0_dp*PI*t/7.0_dp)

    Ct=C*weekly

    do j=1,N_AGE

      tmp1(j)=sum(Ct(j,:)*(I1/POP))

      tmp2(j)=sum(Ct(j,:)*(I2/POP))

    end do

    lam1=tmp1*BETA1_BASE

    lam2=1.1_dp*tmp2*BETA1_BASE

    control1=1.0_dp-seasonal*L_DIAG

    control2=1.0_dp-alpha*seasonal*L_DIAG

    do j=1,N_AGE

      dy(j)=S(j)*control1(j)*lam1(j)-GAMMA1(j)*I1(j)

      dy(N_AGE+j)= &
        S(j)*control2(j)*lam2(j)-GAMMA2(j)*I2(j)

    end do

  end subroutine rhs


!===============================================================================
! SIMULATION
!===============================================================================

  subroutine simulate(alpha,season,A_week,rho,C,t_eval,y_out,ok)

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: season
    real(dp), intent(in) :: A_week
    real(dp), intent(in) :: rho

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp), intent(out) :: t_eval(:)
    real(dp), intent(out) :: y_out(:)

    logical, intent(out) :: ok

    real(dp) :: y(NY)
    real(dp) :: k1(NY)
    real(dp) :: k2(NY)
    real(dp) :: k3(NY)
    real(dp) :: k4(NY)
    real(dp) :: yt(NY)

    real(dp) :: t,dt

    integer :: i,n

    n=size(t_eval)

    ok=.false.

    if(n<2)return

    y=0.0_dp

    y(5)=1.0e-4_dp

    y(N_AGE+5)=1.0e-5_dp

    dt=SIM_END_DAYS/real(n-1,dp)

    t=0.0_dp

    do i=1,n

      t_eval(i)=t/DAYS_PER_WEEK

      call rhs( &
        t,y,alpha,season,A_week,rho,C,k1)

      yt=max(y+0.5_dp*dt*k1,0.0_dp)

      call rhs( &
        t+0.5_dp*dt,yt,alpha,season,A_week,rho,C,k2)

      yt=max(y+0.5_dp*dt*k2,0.0_dp)

      call rhs( &
        t+0.5_dp*dt,yt,alpha,season,A_week,rho,C,k3)

      yt=max(y+dt*k3,0.0_dp)

      call rhs( &
        t+dt,yt,alpha,season,A_week,rho,C,k4)

      y=max( &
        y+(dt/6.0_dp)*( &
        k1+2.0_dp*k2+2.0_dp*k3+k4), &
        0.0_dp)

      y_out(i)=sum(y)*rho

      t=t+dt

      if(any(y>1.0e6_dp) .or. any(y/=y))then

        ok=.false.
        return

      end if

    end do

    ok=.true.

  end subroutine simulate


!===============================================================================
! LINEAR INTERPOLATION
!===============================================================================

  function interp1d(x,y,xq) result(yq)

    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: xq

    real(dp) :: yq

    integer :: i,n

    n=size(x)

    yq=0.0_dp

    if(xq<x(1).or.xq>x(n))return

    do i=1,n-1

      if(xq>=x(i).and.xq<=x(i+1))then

        yq=y(i)+(y(i+1)-y(i))* &
          (xq-x(i))/(x(i+1)-x(i))

        return

      end if

    end do

  end function interp1d


!===============================================================================
! TIME ALIGNMENT
!===============================================================================

  subroutine align_model( &
      alpha,season,A_week,rho,C,fit,rmse,shift,ok)

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: season
    real(dp), intent(in) :: A_week
    real(dp), intent(in) :: rho

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp), intent(out) :: fit(NDATA)
    real(dp), intent(out) :: rmse
    real(dp), intent(out) :: shift

    logical, intent(out) :: ok

    real(dp), allocatable :: te(:)
    real(dp), allocatable :: yo(:)

    real(dp) :: yf(NDATA)
    real(dp) :: res(NDATA)
    real(dp) :: best_fit(NDATA)

    real(dp) :: s
    real(dp) :: best_rmse

    integer :: i,idx

    logical :: good

    allocate(te(N_T_EVAL),yo(N_T_EVAL))

    best_rmse=huge(1.0_dp)

    ok=.false.

    shift=0.0_dp

    fit=0.0_dp

    call simulate( &
      alpha,season,A_week,rho,C,te,yo,good)

    if(.not.good)then

      deallocate(te,yo)

      rmse=huge(1.0_dp)

      return

    end if

    do idx=1,100

      s=-15.0_dp+ &
        30.0_dp*real(idx-1,dp)/99.0_dp

      do i=1,NDATA

        yf(i)=interp1d( &
          te,yo,WEEKS_DATA(i)+s)

      end do

      if(any(yf<0.0_dp))cycle

      res=Y_DATA_NORM-yf

      rmse=sqrt( &
        sum(res**2)/real(NDATA,dp))

      if(rmse<best_rmse)then

        best_rmse=rmse

        best_fit=yf

        shift=s

      end if

    end do

    if(best_rmse<1.0e8_dp)then

      fit=best_fit

      rmse=best_rmse

      ok=.true.

    end if

    deallocate(te,yo)

  end subroutine align_model

end module model_mod


!===============================================================================
! STATISTICS
!===============================================================================

module stats_mod

  use precision_mod

  implicit none

contains


  function pearson(x,y) result(r)

    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)

    real(dp) :: r
    real(dp) :: mx,my,sx,sy

    mx=sum(x)/real(size(x),dp)

    my=sum(y)/real(size(y),dp)

    sx=sqrt(sum((x-mx)**2))

    sy=sqrt(sum((y-my)**2))

    if(sx<=0.0_dp.or.sy<=0.0_dp)then

      r=0.0_dp

    else

      r=sum((x-mx)*(y-my))/(sx*sy)

    end if

  end function pearson


!===============================================================================
! QUICKSORT
!===============================================================================

  recursive subroutine quicksort(a,l,r)

    real(dp), intent(inout) :: a(:)

    integer, intent(in) :: l,r

    integer :: i,j

    real(dp) :: p,tmp

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

    if(l<j)call quicksort(a,l,j)

    if(i<r)call quicksort(a,i,r)

  end subroutine quicksort


!===============================================================================
! PERCENTILE
!===============================================================================

  subroutine percentile(v,p,q)

    real(dp), intent(in) :: v(:)
    real(dp), intent(in) :: p

    real(dp), intent(out) :: q

    real(dp), allocatable :: a(:)

    integer :: n,k

    real(dp) :: z

    n=size(v)

    if(n==0)then

      q=0.0_dp

      return

    end if

    allocate(a(n))

    a=v

    call quicksort(a,1,n)

    if(n==1)then

      q=a(1)

      deallocate(a)

      return

    end if

    z=p*real(n-1,dp)+1.0_dp

    k=max(1,int(floor(z)))

    if(k>=n)then

      q=a(n)

    else

      q=a(k)+(z-real(k,dp))* &
        (a(k+1)-a(k))

    end if

    deallocate(a)

  end subroutine percentile


!===============================================================================
! SUMMARY
!===============================================================================

  subroutine summary(v,meanv,sd,q025,q50,q975)

    real(dp), intent(in) :: v(:)

    real(dp), intent(out) :: meanv
    real(dp), intent(out) :: sd
    real(dp), intent(out) :: q025
    real(dp), intent(out) :: q50
    real(dp), intent(out) :: q975

    integer :: n

    n=size(v)

    if(n==0)then

      meanv=0.0_dp
      sd=0.0_dp
      q025=0.0_dp
      q50=0.0_dp
      q975=0.0_dp

      return

    end if

    meanv=sum(v)/real(n,dp)

    if(n>1)then

      sd=sqrt( &
        sum((v-meanv)**2)/real(n-1,dp))

    else

      sd=0.0_dp

    end if

    call percentile(v,0.025_dp,q025)

    call percentile(v,0.5_dp,q50)

    call percentile(v,0.975_dp,q975)

  end subroutine summary

end module stats_mod


!===============================================================================
! LIKELIHOOD
!===============================================================================

module likelihood_mod

  use precision_mod
  use settings_mod
  use data_mod
  use model_mod

  implicit none

contains


!===============================================================================
! LOG PRIOR
!===============================================================================

  function log_prior(theta) result(lp)

    real(dp), intent(in) :: theta(NDIM)

    real(dp) :: lp

    if( &
      theta(1)>=ALPHA_LO .and. &
      theta(1)<=ALPHA_HI .and. &
      theta(2)>=SEASON_LO .and. &
      theta(2)<=SEASON_HI .and. &
      theta(3)>=A_WEEK_LO .and. &
      theta(3)<=A_WEEK_HI .and. &
      theta(4)>=RHO_LO .and. &
      theta(4)<=RHO_HI .and. &
      theta(5)>=SIGMA_LO .and. &
      theta(5)<=SIGMA_HI )then

      lp=0.0_dp

    else

      lp=-huge(1.0_dp)

    end if

  end function log_prior


!===============================================================================
! GAUSSIAN LOG LIKELIHOOD
!===============================================================================

  function log_likelihood(theta,C) result(ll)

    real(dp), intent(in) :: theta(NDIM)

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp) :: ll

    real(dp) :: fit(NDATA)
    real(dp) :: rmse
    real(dp) :: shift
    real(dp) :: res(NDATA)

    logical :: ok

    call align_model( &
      theta(1),theta(2),theta(3),theta(4), &
      C,fit,rmse,shift,ok)

    if(.not.ok)then

      ll=-huge(1.0_dp)

      return

    end if

    if(theta(5)<=0.0_dp)then

      ll=-huge(1.0_dp)

      return

    end if

    res=Y_DATA_NORM-fit

    ll=-0.5_dp*sum( &
      (res/theta(5))**2 + &
      log(2.0_dp*PI*theta(5)**2) )

  end function log_likelihood


!===============================================================================
! LOG POSTERIOR
!===============================================================================

  function log_posterior(theta,C) result(lp)

    real(dp), intent(in) :: theta(NDIM)

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp) :: lp

    lp=log_prior(theta)

    if(lp<-1.0e100_dp)return

    lp=lp+log_likelihood(theta,C)

  end function log_posterior


!===============================================================================
! MLE OBJECTIVE
!
! IMPORTANT:
! The MLE now uses the same Gaussian likelihood as the Bayesian model.
!
! sigma is profiled analytically:
!
!   sigma_hat^2 = mean(residual^2)
!
! and therefore the four-dimensional optimizer minimizes
!
!   - log L(theta, sigma_hat).
!
!===============================================================================

  function objective4(x,C) result(f)

    real(dp), intent(in) :: x(4)

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp) :: f

    real(dp) :: fit(NDATA)
    real(dp) :: rmse
    real(dp) :: shift
    real(dp) :: sigma_hat
    real(dp) :: sigma2

    logical :: ok

    call align_model( &
      x(1),x(2),x(3),x(4), &
      C,fit,rmse,shift,ok)

    if(.not.ok)then

      f=1.0e9_dp

      return

    end if

    !---------------------------------------------------------------------------
    ! Profile sigma analytically from the residuals.
    !---------------------------------------------------------------------------

    sigma_hat=max( &
      sqrt(sum((Y_DATA_NORM-fit)**2)/real(NDATA,dp)), &
      SIGMA_LO*10.0_dp )

    sigma_hat=min( &
      sigma_hat, &
      SIGMA_HI*0.999999_dp )

    sigma2=sigma_hat**2

    !---------------------------------------------------------------------------
    ! Negative Gaussian log-likelihood.
    !---------------------------------------------------------------------------

    f=0.5_dp*real(NDATA,dp)* &
      (1.0_dp+log(2.0_dp*PI*sigma2))

  end function objective4

end module likelihood_mod


!===============================================================================
! DIFFERENTIAL EVOLUTION
!===============================================================================

module optimization_mod

  use precision_mod
  use settings_mod
  use rng_mod
  use likelihood_mod

  implicit none

contains


  subroutine differential_evolution_4d( &
      C,bounds,xbest,fbest,maxiter,popsize,seed)

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp), intent(in) :: bounds(4,2)

    real(dp), intent(out) :: xbest(4)
    real(dp), intent(out) :: fbest

    integer, intent(in) :: maxiter
    integer, intent(in) :: popsize
    integer, intent(in) :: seed

    real(dp) :: pop(popsize,4)
    real(dp) :: fcost(popsize)

    real(dp) :: trial(4)
    real(dp) :: mutant(4)

    real(dp) :: ff
    real(dp) :: CR

    integer :: i,j,it
    integer :: a,b,cc,r
    integer :: sl

    sl=seed

    ff=0.8_dp

    CR=0.9_dp

    !---------------------------------------------------------------------------
    ! Initial population
    !---------------------------------------------------------------------------

    do i=1,popsize

      do j=1,4

        pop(i,j)= &
          bounds(j,1)+ &
          (bounds(j,2)-bounds(j,1))* &
          urand(sl)

      end do

      fcost(i)=objective4(pop(i,:),C)

    end do

    !---------------------------------------------------------------------------
    ! Differential evolution
    !---------------------------------------------------------------------------

    do it=1,maxiter

      do i=1,popsize

        do

          a=1+int(real(popsize,dp)*urand(sl))

          if(a/=i)exit

        end do

        do

          b=1+int(real(popsize,dp)*urand(sl))

          if(b/=i.and.b/=a)exit

        end do

        do

          cc=1+int(real(popsize,dp)*urand(sl))

          if(cc/=i.and.cc/=a.and.cc/=b)exit

        end do

        mutant= &
          pop(a,:)+ff*(pop(b,:)-pop(cc,:))

        r=1+int(4.0_dp*urand(sl))

        do j=1,4

          if(urand(sl)<CR.or.j==r)then

            trial(j)=mutant(j)

          else

            trial(j)=pop(i,j)

          end if

          trial(j)=min( &
            max(trial(j),bounds(j,1)), &
            bounds(j,2))

        end do

        if(objective4(trial,C)<fcost(i))then

          pop(i,:)=trial

          fcost(i)=objective4(trial,C)

        end if

      end do

    end do

    i=minloc(fcost,1)

    xbest=pop(i,:)

    fbest=fcost(i)

  end subroutine differential_evolution_4d

end module optimization_mod


!===============================================================================
! MCMC
!===============================================================================

module mcmc_mod

  use precision_mod
  use settings_mod
  use rng_mod
  use likelihood_mod
  use stats_mod

  implicit none

contains


!===============================================================================
! INITIAL WALKER POSITIONS
!===============================================================================

  subroutine make_start(center,chain,seed)

    real(dp), intent(in) :: center(NDIM)

    real(dp), intent(out) :: chain(NWALKERS,NDIM)

    integer, intent(in) :: seed

    integer :: i,j,sl

    real(dp) :: scale(NDIM)

    sl=seed

    scale=[ &
      0.001_dp, &
      0.001_dp, &
      0.001_dp, &
      0.5_dp, &
      0.001_dp ]

    do i=1,NWALKERS

      do

        do j=1,NDIM

          chain(i,j)= &
            center(j)+scale(j)*grand(sl)

        end do

        if(log_prior(chain(i,:))>-1.0e100_dp)exit

      end do

    end do

  end subroutine make_start


!===============================================================================
! RUN ONE MCMC CHAIN
!===============================================================================

  subroutine run_chain( &
      center,C,chain_samples,acceptance,seed)

    real(dp), intent(in) :: center(NDIM)

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    real(dp), intent(out) :: chain_samples(:,:)

    real(dp), intent(out) :: acceptance

    integer, intent(in) :: seed

    real(dp) :: chain(NWALKERS,NDIM)

    real(dp) :: lp(NWALKERS)

    real(dp) :: proposal(NDIM)

    real(dp) :: newlp

    real(dp) :: z
    real(dp) :: u
    real(dp) :: logalpha

    integer :: i,j,k
    integer :: sl
    integer :: na,np
    integer :: idx,keep

    sl=seed

    na=0
    np=0

    call make_start(center,chain,sl)

    do i=1,NWALKERS

      lp(i)=log_posterior(chain(i,:),C)

    end do

    idx=0

    keep=(NSTEPS-BURN_MCMC)/THIN

    do k=1,NSTEPS

      do i=1,NWALKERS

        do

          j=1+int(real(NWALKERS,dp)*urand(sl))

          if(j/=i)exit

        end do

        !-----------------------------------------------------------------------
        ! Goodman-Weare stretch move
        !-----------------------------------------------------------------------

        z=((2.0_dp-1.0_dp)*urand(sl)+1.0_dp)**2/2.0_dp

        proposal= &
          chain(j,:)+ &
          z*(chain(i,:)-chain(j,:))

        newlp=log_posterior(proposal,C)

        np=np+1

        if(newlp>-1.0e100_dp)then

          logalpha= &
            real(NDIM-1,dp)*log(z)+ &
            newlp-lp(i)

          u=log(urand(sl))

          if(u<logalpha)then

            chain(i,:)=proposal

            lp(i)=newlp

            na=na+1

          end if

        end if

      end do

      !-----------------------------------------------------------------------
      ! Store posterior samples
      !-----------------------------------------------------------------------

      if(k>BURN_MCMC.and. &
         mod(k-BURN_MCMC,THIN)==0)then

        do i=1,NWALKERS

          idx=idx+1

          chain_samples(idx,:)=chain(i,:)

        end do

      end if

      if(mod(k,500)==0.or.k==NSTEPS)then

        write(*,'(A,I5,A,I5,A,F6.1,A)') &
          '  MCMC ', &
          k,'/',NSTEPS,' (', &
          100.0_dp*real(k,dp)/real(NSTEPS,dp),'%)'

      end if

    end do

    acceptance= &
      real(na,dp)/real(max(1,np),dp)

  end subroutine run_chain


!===============================================================================
! GELMAN-RUBIN R-HAT
!===============================================================================

  subroutine gelman_rhat(chains,rhat)

    real(dp), intent(in) :: chains(:,:,:)

    real(dp), intent(out) :: rhat(NDIM)

    integer :: m,n,d,c

    real(dp) :: mean_c
    real(dp) :: grandm

    real(dp) :: W
    real(dp) :: B
    real(dp) :: varhat

    m=size(chains,1)

    n=size(chains,2)

    do d=1,NDIM

      W=0.0_dp

      grandm=0.0_dp

      do c=1,m

        mean_c= &
          sum(chains(c,:,d))/real(n,dp)

        grandm=grandm+mean_c

        W=W+ &
          sum((chains(c,:,d)-mean_c)**2)/ &
          real(max(1,n-1),dp)

      end do

      grandm=grandm/real(m,dp)

      W=W/real(m,dp)

      B=0.0_dp

      do c=1,m

        B=B+( &
          sum(chains(c,:,d))/real(n,dp)- &
          grandm )**2

      end do

      B=real(n,dp)*B/real(max(1,m-1),dp)

      varhat=( &
        real(n-1,dp)/real(n,dp))*W + &
        B/real(n,dp)

      !-----------------------------------------------------------------------
      ! CORRECTED: assign only rhat(d)
      !-----------------------------------------------------------------------

      if(W>0.0_dp)then

        rhat(d)=sqrt(varhat/W)

      else

        rhat(d)=huge(1.0_dp)

      end if

    end do

  end subroutine gelman_rhat

end module mcmc_mod


!===============================================================================
! MAIN PROGRAM
!===============================================================================

program COVID_RURAL_URBAN_COMPLETE

  use precision_mod
  use settings_mod
  use data_mod
  use contact_mod
  use rng_mod
  use model_mod
  use stats_mod
  use likelihood_mod
  use optimization_mod
  use mcmc_mod

  implicit none


!===============================================================================
! RESULT TYPE
!===============================================================================

  type :: result_type

    character(len=16) :: name

    real(dp) :: mle4(4)
    real(dp) :: mle5(NDIM)

    real(dp) :: mle_objective
    real(dp) :: mle_loglik

    real(dp) :: mle_r
    real(dp) :: mle_rmse
    real(dp) :: mle_shift

    real(dp) :: acceptance(NCHAINS)

    real(dp) :: rhat(NDIM)

    real(dp) :: post_mean(NDIM)
    real(dp) :: post_median(NDIM)

    real(dp) :: post_lo(NDIM)
    real(dp) :: post_hi(NDIM)

    real(dp) :: post_r_median
    real(dp) :: post_r_lo
    real(dp) :: post_r_hi

    real(dp) :: post_rmse_median
    real(dp) :: post_rmse_lo
    real(dp) :: post_rmse_hi

    real(dp) :: coverage

    real(dp) :: mle_curve(NDATA)

    real(dp), allocatable :: samples(:,:)

    real(dp), allocatable :: traj(:,:)

    real(dp), allocatable :: traj_med(:)
    real(dp), allocatable :: traj_lo(:)
    real(dp), allocatable :: traj_hi(:)

    real(dp), allocatable :: chain_samples(:,:,:)

  end type result_type


!===============================================================================
! VARIABLES
!===============================================================================

  real(dp) :: bounds(4,2)

  real(dp) :: xbest(4)
  real(dp) :: fbest

  real(dp) :: theta0(NDIM)

  real(dp) :: fit(NDATA)

  real(dp) :: rmse
  real(dp) :: shift
  real(dp) :: sigmahat
  real(dp) :: ll

  real(dp) :: q025
  real(dp) :: q50
  real(dp) :: q975

  real(dp) :: meanv
  real(dp) :: sdv

  real(dp) :: coverage

  integer :: i
  integer :: j

  integer :: nkeep
  integer :: nflat

  integer :: seed

  logical :: ok

  type(result_type) :: R
  type(result_type) :: U


!===============================================================================
! INITIALIZATION
!===============================================================================

  call initialize_contacts()

  call seed_rng(42)


!===============================================================================
! CORRECT PARAMETER BOUNDS
!
! Explicit assignments avoid Fortran column-major RESHAPE ordering problems.
!===============================================================================

  bounds(1,1)=ALPHA_LO
  bounds(1,2)=ALPHA_HI

  bounds(2,1)=SEASON_LO
  bounds(2,2)=SEASON_HI

  bounds(3,1)=A_WEEK_LO
  bounds(3,2)=A_WEEK_HI

  bounds(4,1)=RHO_LO
  bounds(4,2)=RHO_HI


!===============================================================================
! NUMBER OF POSTERIOR SAMPLES
!===============================================================================

  nkeep=(NSTEPS-BURN_MCMC)/THIN*NWALKERS

  nflat=nkeep*NCHAINS


!===============================================================================
! RURAL
!===============================================================================

  call analyze( &
    'RURAL', &
    RURAL_MATRIX, &
    R, &
    42)


!===============================================================================
! URBAN
!===============================================================================

  call analyze( &
    'URBAN', &
    URBAN_MATRIX, &
    U, &
    4242)


!===============================================================================
! COMBINED OUTPUT
!===============================================================================

  call write_combined(R,U)


!===============================================================================
! PRINT RESULTS
!===============================================================================

  call print_result(R)

  call print_result(U)


  print *,'==============================================================='

  print *, &
    'COVID-19 RURAL vs URBAN CALIBRATION COMPLETED'

  print *,'==============================================================='


contains


!===============================================================================
! ANALYZE ONE SCENARIO
!===============================================================================

  subroutine analyze(name,C,res,base_seed)

    character(len=*), intent(in) :: name

    real(dp), intent(in) :: C(N_AGE,N_AGE)

    type(result_type), intent(out) :: res

    integer, intent(in) :: base_seed


!===============================================================================
! LOCAL VARIABLES
!===============================================================================

    real(dp) :: chainmean(NDIM)

    real(dp) :: chains( &
      NCHAINS,nkeep,NDIM)

    real(dp) :: samples_all( &
      nflat,NDIM)

    real(dp) :: theta(NDIM)

    real(dp) :: models( &
      NDATA,N_POSTERIOR_MODELS)

    real(dp) :: perr_r( &
      N_POSTERIOR_MODELS)

    real(dp) :: perr_rm( &
      N_POSTERIOR_MODELS)

    real(dp) :: fit0(NDATA)

    real(dp) :: rr

    real(dp) :: qq1
    real(dp) :: qq2
    real(dp) :: qq3

    real(dp) :: meanx
    real(dp) :: sdx

    real(dp) :: rm
    real(dp) :: rs

    real(dp) :: shift0

    real(dp) :: coverage0

    real(dp) :: paramtmp(nflat)

    real(dp) :: model(NDATA)

    real(dp) :: resid(NDATA)

    real(dp) :: best_sigma

    real(dp) :: posterior_r( &
      nkeep*NCHAINS)

    real(dp) :: posterior_rmse( &
      nkeep*NCHAINS)

    real(dp) :: sample_fit(NDATA)

    real(dp) :: tr_lo(NDATA)

    real(dp) :: tr_hi(NDATA)

    real(dp) :: tr_med(NDATA)

    real(dp) :: rv(N_POSTERIOR_MODELS)

    real(dp) :: mv(N_POSTERIOR_MODELS)

    real(dp) :: z

    integer :: ch
    integer :: i
    integer :: d
    integer :: idx

    integer :: valid
    integer :: seedloc
    integer :: bi

    logical :: ok0


!===============================================================================
! SCENARIO NAME
!===============================================================================

    res%name=name


    print *,'==============================================================='

    print *, &
      ' COVID-19 '//trim(name)// &
      ' - MLE + BAYESIAN CALIBRATION'

    print *,'==============================================================='


!===============================================================================
! MLE
!===============================================================================

    call differential_evolution_4d( &
      C, &
      bounds, &
      res%mle4, &
      res%mle_objective, &
      MLE_MAXITER, &
      MLE_POPSIZE, &
      base_seed)


!===============================================================================
! MLE MODEL
!===============================================================================

    call align_model( &
      res%mle4(1), &
      res%mle4(2), &
      res%mle4(3), &
      res%mle4(4), &
      C, &
      fit0, &
      rm, &
      shift0, &
      ok0)


    if(.not.ok0)then

      stop 'MLE model evaluation failed.'

    end if


!===============================================================================
! MLE SIGMA
!===============================================================================

    best_sigma=max( &
      sqrt( &
        sum((Y_DATA_NORM-fit0)**2)/ &
        real(NDATA,dp)), &
      SIGMA_LO*10.0_dp)

    best_sigma=min( &
      best_sigma, &
      SIGMA_HI*0.999999_dp)


    res%mle5=[ &
      res%mle4, &
      best_sigma ]


    res%mle_shift=shift0

    res%mle_rmse=rm

    res%mle_r=pearson( &
      Y_DATA_NORM,fit0)

    res%mle_loglik=log_likelihood( &
      res%mle5,C)

    res%mle_curve=fit0


    theta0=res%mle5


!===============================================================================
! PRINT MLE
!===============================================================================

    print *,'MLE parameters:',theta0

    print *, &
      'MLE objective =', &
      res%mle_objective, &
      ' MLE log-likelihood =', &
      res%mle_loglik


!===============================================================================
! MCMC CHAINS
!===============================================================================

    do ch=1,NCHAINS

      !-----------------------------------------------------------------------
      ! Corrected variable: ch, not c
      !-----------------------------------------------------------------------

      print *, &
        'Running chain ', &
        ch,'/',NCHAINS, &
        ' for ',trim(name)

      call run_chain( &
        theta0, &
        C, &
        chains(ch,:,:), &
        res%acceptance(ch), &
        base_seed+10000*ch)

    end do


!===============================================================================
! R-HAT
!===============================================================================

    call gelman_rhat( &
      chains, &
      res%rhat)


!===============================================================================
! SAVE CHAIN SAMPLES
!===============================================================================

    allocate( &
      res%chain_samples( &
      NCHAINS,nkeep,NDIM))

    res%chain_samples=chains


!===============================================================================
! FLATTEN CHAINS
!===============================================================================

    allocate( &
      res%samples( &
      nflat,NDIM))

    idx=0

    do ch=1,NCHAINS

      do i=1,nkeep

        idx=idx+1

        res%samples(idx,:)= &
          chains(ch,i,:)

      end do

    end do


!===============================================================================
! PARAMETER POSTERIOR SUMMARIES
!===============================================================================

    do d=1,NDIM

      paramtmp=res%samples(:,d)

      call summary( &
        paramtmp, &
        meanx, &
        sdx, &
        qq1, &
        qq2, &
        qq3)

      res%post_mean(d)=meanx

      res%post_median(d)=qq2

      res%post_lo(d)=qq1

      res%post_hi(d)=qq3

    end do


!===============================================================================
! POSTERIOR TRAJECTORIES
!===============================================================================

    allocate( &
      res%traj( &
      NDATA,N_POSTERIOR_MODELS))

    valid=0

    seedloc=base_seed+777


    do i=1,N_POSTERIOR_MODELS

      idx=1+int( &
        real(nflat,dp)*urand(seedloc))

      idx=min( &
        max(idx,1), &
        nflat)


      call align_model( &
        res%samples(idx,1), &
        res%samples(idx,2), &
        res%samples(idx,3), &
        res%samples(idx,4), &
        C, &
        sample_fit, &
        rm, &
        shift0, &
        ok0)


      if(ok0)then

        valid=valid+1

        res%traj(:,valid)=sample_fit

        rv(valid)=pearson( &
          Y_DATA_NORM, &
          sample_fit)

        mv(valid)=sqrt( &
          sum((Y_DATA_NORM-sample_fit)**2)/ &
          real(NDATA,dp))

      end if

    end do


    if(valid<2)then

      stop 'Too few valid posterior trajectories.'

    end if


!===============================================================================
! POINTWISE POSTERIOR TRAJECTORY CrI
!===============================================================================

    do i=1,NDATA

      call summary( &
        res%traj(i,1:valid), &
        meanx, &
        sdx, &
        tr_lo(i), &
        tr_med(i), &
        tr_hi(i))

    end do


!===============================================================================
! STORE TRAJECTORY SUMMARIES
!===============================================================================

    allocate(res%traj_med(NDATA))

    allocate(res%traj_lo(NDATA))

    allocate(res%traj_hi(NDATA))


    res%traj_med=tr_med

    res%traj_lo=tr_lo

    res%traj_hi=tr_hi


!===============================================================================
! POSTERIOR PEARSON r
!===============================================================================

    call summary( &
      rv(1:valid), &
      meanx, &
      sdx, &
      qq1, &
      res%post_r_median, &
      qq3)

    res%post_r_lo=qq1

    res%post_r_hi=qq3


!===============================================================================
! POSTERIOR RMSE
!===============================================================================

    call summary( &
      mv(1:valid), &
      meanx, &
      sdx, &
      qq1, &
      res%post_rmse_median, &
      qq3)

    res%post_rmse_lo=qq1

    res%post_rmse_hi=qq3


!===============================================================================
! POINTWISE CrI COVERAGE
!===============================================================================

    coverage0= &
      100.0_dp*real( &
      count( &
      Y_DATA_NORM>=tr_lo.and. &
      Y_DATA_NORM<=tr_hi),dp)/ &
      real(NDATA,dp)

    res%coverage=coverage0


!===============================================================================
! SAVE
!===============================================================================

    call save_result(res)

  end subroutine analyze


!===============================================================================
! SAVE RESULT
!===============================================================================

  subroutine save_result(res)

    type(result_type), intent(in) :: res

    integer :: u
    integer :: i
    integer :: d

    character(len=80) :: fn


!===============================================================================
! FIT RESULTS CSV
!===============================================================================

    fn=trim(res%name)//'_COVID_Fit_Results.csv'

    u=31

    open( &
      u, &
      file=fn, &
      status='replace')

    write(u,'(A)') &
      'week,r_data,mle_curve,posterior_median,posterior_lower,posterior_upper'


    do i=1,NDATA

      write(u,'(F12.5,5(ES18.10,","))') &
        WEEKS_DATA(i), &
        Y_DATA_NORM(i), &
        res%mle_curve(i), &
        res%traj_med(i), &
        res%traj_lo(i), &
        res%traj_hi(i)

    end do


    close(u)


!===============================================================================
! POSTERIOR SAMPLES
!===============================================================================

    fn=trim(res%name)//'_COVID_Posterior_Samples.dat'

    u=32

    open( &
      u, &
      file=fn, &
      status='replace')

    write(u,'(A)') &
      '# alpha season A_week rho sigma'


    do i=1,size(res%samples,1)

      write(u,'(5ES18.10)') &
        res%samples(i,:)

    end do


    close(u)


!===============================================================================
! MCMC DIAGNOSTICS
!===============================================================================

    fn=trim(res%name)//'_COVID_MCMC_Diagnostics.dat'

    u=33

    open( &
      u, &
      file=fn, &
      status='replace')

    write(u,'(A)') &
      '# chain acceptance'


    do i=1,NCHAINS

      write(u,'(I5,1X,F14.8)') &
        i,res%acceptance(i)

    end do


    write(u,'(A)') &
      '# parameter Rhat'


    do d=1,NDIM

      write(u,'(A,1X,F14.8)') &
        pname(d), &
        res%rhat(d)

    end do


    close(u)


!===============================================================================
! PARAMETER SUMMARY
!===============================================================================

    fn=trim(res%name)//'_COVID_Parameter_Summary.dat'

    u=34

    open( &
      u, &
      file=fn, &
      status='replace')

    write(u,'(A)') &
      '# parameter MLE posterior_mean posterior_median CrI_low CrI_high Rhat'


    do d=1,NDIM

      write(u,'(A,6(1X,ES18.10))') &
        pname(d), &
        res%mle5(d), &
        res%post_mean(d), &
        res%post_median(d), &
        res%post_lo(d), &
        res%post_hi(d), &
        res%rhat(d)

    end do


    close(u)


!===============================================================================
! SUMMARY TXT
!===============================================================================

    fn=trim(res%name)//'_COVID_Summary.txt'

    u=35

    open( &
      u, &
      file=fn, &
      status='replace')


    write(u,'(A)') &
      trim(res%name)//' COVID-19 RESULTS'


    write(u,'(A,F16.8)') &
      'MLE log-likelihood = ', &
      res%mle_loglik


    write(u,'(A,F16.8)') &
      'MLE Pearson r = ', &
      res%mle_r


    write(u,'(A,F16.8)') &
      'MLE RMSE = ', &
      res%mle_rmse


    write(u,'(A,F16.8)') &
      'MLE shift weeks = ', &
      res%mle_shift


    write(u,'(A,F16.8)') &
      'Posterior Pearson r median = ', &
      res%post_r_median


    write(u,'(A,2F16.8)') &
      'Posterior Pearson r 95% CrI = ', &
      res%post_r_lo, &
      res%post_r_hi


    write(u,'(A,F16.8)') &
      'Posterior RMSE median = ', &
      res%post_rmse_median


    write(u,'(A,2F16.8)') &
      'Posterior RMSE 95% CrI = ', &
      res%post_rmse_lo, &
      res%post_rmse_hi


    write(u,'(A,F10.3)') &
      'Latent 95% CrI coverage (%) = ', &
      res%coverage


    write(u,'(A)') &
      'Chain acceptance fractions:'


    do i=1,NCHAINS

      write(u,'(I4,1X,F12.6)') &
        i,res%acceptance(i)

    end do


    write(u,'(A)') &
      'Rhat:'


    do d=1,NDIM

      write(u,'(A,1X,F12.6)') &
        pname(d), &
        res%rhat(d)

    end do


    close(u)

  end subroutine save_result


!===============================================================================
! COMBINED RURAL / URBAN RESULTS
!===============================================================================

  subroutine write_combined(R,U)

    type(result_type), intent(in) :: R
    type(result_type), intent(in) :: U

    integer :: unitc
    integer :: i


    unitc=36

    open( &
      unitc, &
      file='COVID_Rural_Urban_Combined_Results.csv', &
      status='replace')


    write(unitc,'(A)') &
      'week,observed,rural_median,rural_lower,rural_upper,urban_median,urban_lower,urban_upper'


    do i=1,NDATA

      write(unitc,'(8(ES18.10,","))') &
        WEEKS_DATA(i), &
        Y_DATA_NORM(i), &
        R%traj_med(i), &
        R%traj_lo(i), &
        R%traj_hi(i), &
        U%traj_med(i), &
        U%traj_lo(i), &
        U%traj_hi(i)

    end do


    close(unitc)

  end subroutine write_combined


!===============================================================================
! PRINT RESULTS
!===============================================================================

  subroutine print_result(res)

    type(result_type), intent(in) :: res

    integer :: d


    print *,'---------------------------------------------------------------'

    print *,trim(res%name)


    print '(A,ES14.6)', &
      'MLE log-likelihood = ', &
      res%mle_loglik


    print '(A,F10.6)', &
      'MLE Pearson r = ', &
      res%mle_r


    print '(A,F10.6)', &
      'MLE RMSE = ', &
      res%mle_rmse


    print '(A,F10.6)', &
      'MLE shift = ', &
      res%mle_shift


    do d=1,NDIM

      print '(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,F8.4)', &
        'parameter ', &
        res%mle5(d), &
        '  median ', &
        res%post_median(d), &
        '  CrI[', &
        res%post_lo(d), &
        ', ', &
        res%post_hi(d), &
        '] Rhat=', &
        res%rhat(d)

    end do


    print '(A,5F10.5)', &
      'acceptance = ', &
      res%acceptance


    print '(A,F10.6,A,F10.6,A,F10.6)', &
      'Posterior r = ', &
      res%post_r_median, &
      ' CrI[', &
      res%post_r_lo, &
      ', ', &
      res%post_r_hi, &
      ']'


    print '(A,F10.6,A,F10.6,A,F10.6)', &
      'Posterior RMSE = ', &
      res%post_rmse_median, &
      ' CrI[', &
      res%post_rmse_lo, &
      ', ', &
      res%post_rmse_hi, &
      ']'


    print '(A,F8.2,A)', &
      'Latent CrI coverage = ', &
      res%coverage, &
      ' %'

  end subroutine print_result


!===============================================================================
! PARAMETER NAMES
!===============================================================================

  function pname(d) result(s)

    integer, intent(in) :: d

    character(len=12) :: s

    select case(d)

    case(1)

      s='alpha'

    case(2)

      s='season'

    case(3)

      s='A_week'

    case(4)

      s='rho'

    case(5)

      s='sigma'

    end select

  end function pname


end program COVID_RURAL_URBAN_COMPLETE
