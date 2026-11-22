!===============================================================================
! Sars-Cov2 LYON
!===============================================================================


module precision_mod
  implicit none

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: i8 = selected_int_kind(18)

  real(dp), parameter :: PI = &
       3.1415926535897932384626433832795_dp

end module precision_mod


!===============================================================================
! SETTINGS
!===============================================================================

module settings_mod

  use precision_mod
  implicit none

  integer, parameter :: NG = 4
  integer, parameter :: NY = 8

  !---------------------------------------------------------------------------
  ! Weekly periodic contact structure
  !---------------------------------------------------------------------------

  real(dp), parameter :: SHIFT_DURATION = 1.0_dp/3.0_dp

  real(dp), parameter :: SEASONAL_PERIOD_DAYS = 52.0_dp*7.0_dp

  !---------------------------------------------------------------------------
  ! Burn-in and model duration
  !---------------------------------------------------------------------------

  real(dp), parameter :: BURN_IN_WEEKS = 15.0_dp
  real(dp), parameter :: BURN_IN_DAYS  = BURN_IN_WEEKS*7.0_dp

  real(dp), parameter :: MODEL_WEEKS = 80.0_dp
  real(dp), parameter :: MODEL_DAYS  = MODEL_WEEKS*7.0_dp

  integer, parameter :: N_T_EVAL = 2000

  !---------------------------------------------------------------------------
  ! Bayesian dimension
  !---------------------------------------------------------------------------

  integer, parameter :: NDIM = 4

  integer, parameter :: NWALKERS = 64
  integer, parameter :: NSTEPS   = 220000
  integer, parameter :: BURN_MCMC = 65000
  integer, parameter :: THIN = 20

  integer, parameter :: N_POSTERIOR_MODELS = 2000

  !---------------------------------------------------------------------------
  ! Bootstrap
  !---------------------------------------------------------------------------

  integer, parameter :: N_BOOTSTRAP = 200
  integer, parameter :: BOOTSTRAP_MAXITER = 800

  !---------------------------------------------------------------------------
  ! Prior bounds
  !---------------------------------------------------------------------------

  real(dp), parameter :: ALPHA_LO  = 0.2_dp
  real(dp), parameter :: ALPHA_HI  = 1.0_dp

  real(dp), parameter :: SEASON_LO = 0.05_dp
  real(dp), parameter :: SEASON_HI = 0.8_dp

  real(dp), parameter :: SIGMA_LO = 0.001_dp
  real(dp), parameter :: SIGMA_HI = 1.0_dp

  real(dp), parameter :: SHIFT_LO = -8.0_dp
  real(dp), parameter :: SHIFT_HI = 8.0_dp

  !---------------------------------------------------------------------------
  ! Population
  !---------------------------------------------------------------------------

  real(dp), parameter :: POP(NG) = [ &
       11.0_dp, 27.0_dp, 8.0_dp, 29.0_dp ]

  !---------------------------------------------------------------------------
  ! Baseline transmission parameters
  !---------------------------------------------------------------------------

  real(dp), parameter :: B1(NG) = 0.9_dp*[ &
       1.5_dp, 1.4_dp, 1.2_dp, 1.0_dp ]

  real(dp), parameter :: B2(NG) = 1.1_dp*[ &
       1.7_dp, 1.5_dp, 1.3_dp, 1.1_dp ]

  !---------------------------------------------------------------------------
  ! Recovery rates
  !---------------------------------------------------------------------------

  real(dp), parameter :: GAM1(NG) = 1.3_dp*[ &
       0.15_dp, 0.14_dp, 0.12_dp, 0.10_dp ]

  real(dp), parameter :: GAM2(NG) = [ &
       0.14_dp, 0.13_dp, 0.11_dp, 0.09_dp ]

  !---------------------------------------------------------------------------
  ! Number of observations
  !---------------------------------------------------------------------------

  integer, parameter :: NDATA = 25

end module settings_mod


!===============================================================================
! Sars-Cov2 LYON DATA
!===============================================================================

module data_mod

  use precision_mod
  use settings_mod
  implicit none

  real(dp), parameter :: WEEKS_DATA(NDATA) = [ &
    2.9723_dp, 3.9487_dp, 4.9854_dp, 5.9785_dp, 7.0143_dp, &
    8.0054_dp, 9.0152_dp,10.0074_dp,11.0239_dp,12.0182_dp, &
   13.0212_dp,14.0127_dp,15.0322_dp,16.0233_dp,17.0168_dp, &
   18.0287_dp,19.0207_dp,20.0111_dp,21.0311_dp,22.0657_dp, &
   23.0117_dp,24.0028_dp,24.9943_dp,26.0041_dp,27.0228_dp ]

  real(dp), parameter :: R_DATA(NDATA) = [ &
    1.4461_dp,1.1592_dp,1.1208_dp,0.9550_dp,1.0086_dp, &
    1.0238_dp,0.9777_dp,0.8902_dp,1.0619_dp,0.7827_dp, &
    1.3608_dp,1.3453_dp,1.2302_dp,1.2500_dp,1.0460_dp, &
    0.8005_dp,0.7390_dp,0.8156_dp,0.6591_dp,0.8231_dp, &
    0.8460_dp,0.8612_dp,0.8427_dp,0.7889_dp,0.7582_dp ]

end module data_mod


!===============================================================================
! CONTACT MATRICES
!===============================================================================

module contact_mod

  use precision_mod
  use settings_mod
  implicit none

  real(dp), save :: C_CONTACT(7,3,NG,NG)
  real(dp), save :: L_SHIFTS(3,NG,NG)

contains


  subroutine initialize_contacts()

    !---------------------------------------------------------------------------
    ! Shift-specific intervention/contact reduction matrices
    !---------------------------------------------------------------------------

    L_SHIFTS(1,:,:) = 1.3_dp*reshape([ &
       0.50_dp,0.0_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.45_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.40_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.0_dp,0.20_dp ],[NG,NG])

    L_SHIFTS(2,:,:) = 1.3_dp*reshape([ &
       0.40_dp,0.0_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.40_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.35_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.0_dp,0.18_dp ],[NG,NG])

    L_SHIFTS(3,:,:) = 1.3_dp*reshape([ &
       0.30_dp,0.0_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.25_dp,0.0_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.20_dp,0.0_dp, &
       0.0_dp,0.0_dp,0.0_dp,0.15_dp ],[NG,NG])

    ! Monday
    C_CONTACT(1,1,:,:) = reshape([ &
      0.089334_dp, 0.030682_dp, 0.168875_dp, 0.226144_dp, &
      0.256768_dp, 0.044045_dp, 0.384233_dp, 0.068801_dp, &
      0.033724_dp, 0.0_dp, 0.148653_dp, 0.042188_dp, &
      0.003302_dp, 0.009303_dp, 0.239060_dp, 0.033885_dp], [NG, NG])
    
    C_CONTACT(1,2,:,:) = reshape([ &
      0.044020_dp, 0.014489_dp, 0.010500_dp, 0.276993_dp, &
      0.046295_dp, 0.046848_dp, 0.107251_dp, 0.004278_dp, &
      0.020234_dp, 0.003369_dp, 0.158113_dp, 0.019922_dp, &
      0.0_dp, 0.005582_dp, 0.043102_dp, 0.016697_dp], [NG, NG])
    
    C_CONTACT(1,3,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.062522_dp, 0.0_dp, 0.055101_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.000825_dp, 0.0_dp, 0.058210_dp, 0.0_dp], [NG, NG])
    
    ! Tuesday
    C_CONTACT(2,1,:,:) = reshape([ &
      0.104870_dp, 0.045170_dp, 0.173250_dp, 0.402778_dp, &
      0.266790_dp, 0.036838_dp, 0.443762_dp, 0.070583_dp, &
      0.028328_dp, 0.001685_dp, 0.124328_dp, 0.062109_dp, &
      0.007429_dp, 0.007815_dp, 0.248391_dp, 0.039778_dp], [NG, NG])
    
    C_CONTACT(2,2,:,:) = reshape([ &
      0.102281_dp, 0.040057_dp, 0.063000_dp, 0.388058_dp, &
      0.096884_dp, 0.062464_dp, 0.321752_dp, 0.025667_dp, &
      0.040469_dp, 0.005054_dp, 0.210817_dp, 0.055078_dp, &
      0.0_dp, 0.011164_dp, 0.090203_dp, 0.038796_dp], [NG, NG])
    
    C_CONTACT(2,3,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.024818_dp, 0.0_dp, 0.049690_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.004127_dp, 0.0_dp, 0.023106_dp, 0.0_dp], [NG, NG])
    
    ! Wednesday
    C_CONTACT(3,1,:,:) = reshape([ &
      0.091923_dp, 0.043466_dp, 0.180250_dp, 0.453627_dp, &
      0.275858_dp, 0.127331_dp, 0.411291_dp, 0.073435_dp, &
      0.060703_dp, 0.015162_dp, 0.429742_dp, 0.059766_dp, &
      0.007429_dp, 0.016746_dp, 0.256834_dp, 0.034868_dp], [NG, NG])
    
    C_CONTACT(3,2,:,:) = reshape([ &
      0.111344_dp, 0.011080_dp, 0.063875_dp, 0.492432_dp, &
      0.068249_dp, 0.058460_dp, 0.141197_dp, 0.026023_dp, &
      0.040469_dp, 0.006739_dp, 0.197303_dp, 0.015234_dp, &
      0.002889_dp, 0.011164_dp, 0.063542_dp, 0.042234_dp], [NG, NG])
    
    C_CONTACT(3,3,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.056794_dp, 0.002402_dp, 0.085112_dp, 0.0_dp, &
      0.004047_dp, 0.0_dp, 0.008108_dp, 0.0_dp, &
      0.0_dp, 0.001116_dp, 0.052877_dp, 0.0_dp], [NG, NG])
    
    ! Thursday
    C_CONTACT(4,1,:,:) = reshape([ &
      0.146301_dp, 0.031534_dp, 0.171500_dp, 0.355943_dp, &
      0.214291_dp, 0.102505_dp, 0.343399_dp, 0.069870_dp, &
      0.056656_dp, 0.020217_dp, 0.345956_dp, 0.043359_dp, &
      0.006604_dp, 0.015629_dp, 0.199513_dp, 0.055493_dp], [NG, NG])
    
    C_CONTACT(4,2,:,:) = reshape([ &
      0.066030_dp, 0.019602_dp, 0.043750_dp, 0.289036_dp, &
      0.050590_dp, 0.025226_dp, 0.098395_dp, 0.017824_dp, &
      0.021583_dp, 0.010108_dp, 0.085138_dp, 0.026953_dp, &
      0.005366_dp, 0.005954_dp, 0.047101_dp, 0.025046_dp], [NG, NG])
    
    C_CONTACT(4,3,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.028159_dp, 0.003203_dp, 0.023615_dp, 0.0_dp, &
      0.001349_dp, 0.0_dp, 0.010811_dp, 0.0_dp, &
      0.002476_dp, 0.000372_dp, 0.026217_dp, 0.0_dp], [NG, NG])
    
    ! Friday
    C_CONTACT(5,1,:,:) = reshape([ &
      0.075092_dp, 0.025568_dp, 0.187250_dp, 0.188677_dp, &
      0.215723_dp, 0.067670_dp, 0.529857_dp, 0.076287_dp, &
      0.032375_dp, 0.131408_dp, 0.228385_dp, 0.035156_dp, &
      0.002889_dp, 0.008931_dp, 0.200846_dp, 0.028483_dp], [NG, NG])
    
    C_CONTACT(5,2,:,:) = reshape([ &
      0.001295_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.001432_dp, 0.001201_dp, 0.001968_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.004054_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.001333_dp, 0.000491_dp], [NG, NG])
    
    C_CONTACT(5,3,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.000477_dp, 0.000400_dp, 0.000492_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.001351_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.000444_dp, 0.0_dp], [NG, NG])
    
    ! Saturday
    C_CONTACT(6,1,:,:) = reshape([ &
      0.029778_dp, 0.010227_dp, 0.074375_dp, 0.074935_dp, &
      0.085907_dp, 0.026828_dp, 0.211549_dp, 0.030301_dp, &
      0.012141_dp, 0.052226_dp, 0.090543_dp, 0.014062_dp, &
      0.000825_dp, 0.003349_dp, 0.079983_dp, 0.011295_dp], [NG, NG])
    
    C_CONTACT(6,2,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.000477_dp, 0.000400_dp, 0.000492_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.001351_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.000444_dp, 0.0_dp], [NG, NG])
    
    C_CONTACT(6,3,:,:) = 0.0_dp
    
    ! Sunday
    C_CONTACT(7,1,:,:) = reshape([ &
      0.022010_dp, 0.007670_dp, 0.056000_dp, 0.056202_dp, &
      0.064431_dp, 0.020021_dp, 0.158908_dp, 0.022815_dp, &
      0.009443_dp, 0.038749_dp, 0.067569_dp, 0.010547_dp, &
      0.000825_dp, 0.002605_dp, 0.059987_dp, 0.008349_dp], [NG, NG])
    
    C_CONTACT(7,2,:,:) = reshape([ &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.000492_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [NG, NG])
    
    C_CONTACT(7,3,:,:) = 0.0_dp
  end subroutine initialize_contacts


  subroutine get_contact_and_L(t,C_out,L_out)

    real(dp),intent(in) :: t
    real(dp),intent(out) :: C_out(NG,NG)
    real(dp),intent(out) :: L_out(NG,NG)

    real(dp) :: t_mod
    integer :: day,shift_idx

    t_mod = modulo(t,7.0_dp)

    day = int(floor(t_mod))+1
    if(day < 1) day = 1
    if(day > 7) day = 7

    shift_idx = int(floor( &
       (t_mod-floor(t_mod))/SHIFT_DURATION))+1

    if(shift_idx < 1) shift_idx = 1
    if(shift_idx > 3) shift_idx = 3

    C_out = C_CONTACT(day,shift_idx,:,:)
    L_out = L_SHIFTS(shift_idx,:,:)

  end subroutine get_contact_and_L

end module contact_mod


!===============================================================================
! RANDOM NUMBER GENERATOR
!===============================================================================

module rng_mod

  use precision_mod
  implicit none

contains


  subroutine seed_rng(seed)

    integer,intent(in) :: seed

    integer :: n,i
    integer,allocatable :: s(:)

    call random_seed(size=n)

    allocate(s(n))

    do i=1,n
      s(i)=mod(seed+104729*i,2147483647)
    end do

    call random_seed(put=s)

    deallocate(s)

  end subroutine seed_rng


  function uniform_rand_seeded(seed) result(u)

    integer,intent(inout) :: seed

    real(dp) :: u
    integer(i8) :: a,m,s

    a=1103515245_i8
    m=2147483647_i8

    s=mod(a*int(seed,i8)+12345_i8,m)

    if(s <= 0_i8) s=1_i8

    seed=int(s)

    u=real(seed,dp)/real(m,dp)

    u=min(max(u,1.0e-15_dp),1.0_dp-1.0e-15_dp)

  end function uniform_rand_seeded


  function normal_rand_seeded(seed) result(z)

    integer,intent(inout) :: seed

    real(dp) :: z,u1,u2

    u1=uniform_rand_seeded(seed)
    u2=uniform_rand_seeded(seed)

    u1=max(u1,1.0e-15_dp)

    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*PI*u2)

  end function normal_rand_seeded

end module rng_mod


!===============================================================================
! Sars-Cov2 MODEL
!===============================================================================

module model_mod

  use precision_mod
  use settings_mod
  use contact_mod
  implicit none

contains


  subroutine rhs(t,y,alpha,season_amp,dy)

    real(dp),intent(in) :: t
    real(dp),intent(in) :: y(NY)
    real(dp),intent(in) :: alpha
    real(dp),intent(in) :: season_amp

    real(dp),intent(out) :: dy(NY)

    real(dp) :: I1(NG),I2(NG),S(NG)
    real(dp) :: lam1(NG),lam2(NG)
    real(dp) :: seasonal
    real(dp) :: Cmat(NG,NG),Lmat(NG,NG)

    integer :: j,k

    I1=max(y(1:NG),0.0_dp)
    I2=max(y(NG+1:NY),0.0_dp)

    S=max(POP-I1-I2,1.0e-10_dp)

    call get_contact_and_L(t,Cmat,Lmat)

    seasonal=1.0_dp+season_amp* &
       sin(2.0_dp*PI*t/SEASONAL_PERIOD_DAYS)

        lam1 = 0.0_dp
    lam2 = 0.0_dp

    do j=1,NG

       !--------------------------------------------------------------------
       ! Strain 1:
       !
       ! lambda_1 = (I - L) C B_1 (I_1 / POP)
       !
       ! L_j acts on the receiving group j.
       ! B_k acts on the infectious/source group k.
       !--------------------------------------------------------------------

       lam1(j) = (1.0_dp - Lmat(j,j)) * &
                 sum(Cmat(j,:) * &
                     (B1(:) * I1(:) / POP(:)))


       !--------------------------------------------------------------------
       ! Strain 2:
       !
       ! lambda_2 = (I - alpha L) C B_2 (I_2 / POP)
       !
       !--------------------------------------------------------------------

       lam2(j) = (1.0_dp - alpha*Lmat(j,j)) * &
                 sum(Cmat(j,:) * &
                     (B2(:) * I2(:) / POP(:)))

    end do

    lam1 = seasonal * lam1
    lam2 = seasonal * lam2
    do j=1,NG

      dy(j)=S(j)*lam1(j)-GAM1(j)*I1(j)

      dy(NG+j)=S(j)*lam2(j)- &
                GAM2(j)*I2(j)

    end do

  end subroutine rhs


!-------------------------------------------------------------------------------
! RK4 integration
!-------------------------------------------------------------------------------

  subroutine integrate_model(alpha,season_amp,t_eval,R_out,npts,ok)

    real(dp),intent(in) :: alpha,season_amp

    real(dp),intent(out) :: t_eval(:)
    real(dp),intent(out) :: R_out(:)

    integer,intent(out) :: npts
    logical,intent(out) :: ok

    real(dp) :: y(NY)
    real(dp) :: k1(NY),k2(NY),k3(NY),k4(NY)
    real(dp) :: y_temp(NY)

    real(dp) :: t,dt,total_days
    real(dp) :: Itot,gamma_weighted

    real(dp),allocatable :: t_local(:)
    real(dp),allocatable :: R_local(:)

    integer :: i,n_keep

    ok=.false.
    npts=0

    !---------------------------------------------------------------------------
    ! Initial infection
    !---------------------------------------------------------------------------

    y=0.0_dp

    y(1)=5.0e-4_dp*POP(1)
    y(NG+1)=2.0e-4_dp*POP(1)

    total_days=BURN_IN_DAYS+MODEL_DAYS

    dt=total_days/real(N_T_EVAL-1,dp)

    allocate(t_local(N_T_EVAL))
    allocate(R_local(N_T_EVAL))

    t=0.0_dp

    do i=1,N_T_EVAL

      t_local(i)=t

      Itot=max(sum(y),1.0e-12_dp)

      gamma_weighted=( &
          sum(GAM1*y(1:NG))+ &
          sum(GAM2*y(NG+1:NY)) )/Itot

      gamma_weighted=max(gamma_weighted,1.0e-6_dp)

      call rhs(t,y,alpha,season_amp,k1)

      ! Effective growth-rate reproduction estimate
      R_local(i)=1.0_dp+ &
         sum(k1)/(Itot*gamma_weighted)

      !-----------------------------------------------------------------------
      ! RK4
      !-----------------------------------------------------------------------

      y_temp=max(y+0.5_dp*dt*k1,0.0_dp)

      call rhs(t+0.5_dp*dt,y_temp, &
               alpha,season_amp,k2)

      y_temp=max(y+0.5_dp*dt*k2,0.0_dp)

      call rhs(t+0.5_dp*dt,y_temp, &
               alpha,season_amp,k3)

      y_temp=max(y+dt*k3,0.0_dp)

      call rhs(t+dt,y_temp, &
               alpha,season_amp,k4)

      y=max(y+(dt/6.0_dp)* &
          (k1+2.0_dp*k2+2.0_dp*k3+k4),0.0_dp)

      t=t+dt

    end do

    !---------------------------------------------------------------------------
    ! Remove burn-in
    !---------------------------------------------------------------------------

    n_keep=0

    do i=1,N_T_EVAL

      if(t_local(i) >= BURN_IN_DAYS) then

        n_keep=n_keep+1

        t_eval(n_keep)= &
           (t_local(i)-BURN_IN_DAYS)/7.0_dp

        R_out(n_keep)=R_local(i)

      end if

    end do

    npts=n_keep

    ok=(n_keep >= 10)

    deallocate(t_local)
    deallocate(R_local)

  end subroutine integrate_model

end module model_mod


!===============================================================================
! STATISTICS
!===============================================================================

module stats_mod

  use precision_mod
  implicit none

contains


  function pearson(x,y) result(r)

    real(dp),intent(in) :: x(:),y(:)

    real(dp) :: r,mx,my,sx,sy

    mx=sum(x)/real(size(x),dp)
    my=sum(y)/real(size(y),dp)

    sx=sqrt(sum((x-mx)**2))
    sy=sqrt(sum((y-my)**2))

    if(sx <= 0.0_dp .or. sy <= 0.0_dp) then

      r=0.0_dp

    else

      r=sum((x-mx)*(y-my))/(sx*sy)

    end if

  end function pearson


  subroutine sort_array(a)

    real(dp),intent(inout) :: a(:)

    call quicksort(a,1,size(a))

  end subroutine sort_array


  recursive subroutine quicksort(a,l,r)

    real(dp),intent(inout) :: a(:)

    integer,intent(in) :: l,r

    integer :: i,j
    real(dp) :: p,tmp

    if(l >= r) return

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

      if(i <= j) then

        tmp=a(i)
        a(i)=a(j)
        a(j)=tmp

        i=i+1
        j=j-1

      end if

      if(i > j) exit

    end do

    if(l < j) call quicksort(a,l,j)
    if(i < r) call quicksort(a,i,r)

  end subroutine quicksort


  subroutine percentile(values,p,q)

    real(dp),intent(in) :: values(:)
    real(dp),intent(in) :: p
    real(dp),intent(out) :: q

    real(dp),allocatable :: a(:)

    integer :: n,k
    real(dp) :: pos

    n=size(values)

    if(n == 0) then
      q=0.0_dp
      return
    end if

    allocate(a(n))
    a=values

    call sort_array(a)

    if(n == 1) then

      q=a(1)

    else

      pos=p*real(n-1,dp)+1.0_dp

      k=int(floor(pos))

      if(k < 1) k=1

      if(k >= n) then

        q=a(n)

      else

        q=a(k)+(pos-real(k,dp))* &
          (a(k+1)-a(k))

      end if

    end if

    deallocate(a)

  end subroutine percentile


  subroutine summary_stats(values,mean_val,std_val,q025,q50,q975)

    real(dp),intent(in) :: values(:)

    real(dp),intent(out) :: mean_val
    real(dp),intent(out) :: std_val
    real(dp),intent(out) :: q025,q50,q975

    integer :: n

    n=size(values)

    if(n == 0) then

      mean_val=0.0_dp
      std_val=0.0_dp
      q025=0.0_dp
      q50=0.0_dp
      q975=0.0_dp

      return

    end if

    mean_val=sum(values)/real(n,dp)

    if(n > 1) then

      std_val=sqrt( &
        sum((values-mean_val)**2)/real(n-1,dp))

    else

      std_val=0.0_dp

    end if

    call percentile(values,0.025_dp,q025)
    call percentile(values,0.50_dp,q50)
    call percentile(values,0.975_dp,q975)

  end subroutine summary_stats

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


  subroutine evaluate_model_at_shift( &
       alpha,season_amp,shift, &
       weeks_obs,R_pred,ok)

    real(dp),intent(in) :: alpha
    real(dp),intent(in) :: season_amp
    real(dp),intent(in) :: shift
    real(dp),intent(in) :: weeks_obs(:)

    real(dp),intent(out) :: R_pred(:)

    logical,intent(out) :: ok

    real(dp),allocatable :: weeks_model(:)
    real(dp),allocatable :: R_model(:)

    integer :: i,j,n_valid
    real(dp) :: xq,w

    logical :: good

    ok=.false.
    R_pred=0.0_dp

    allocate(weeks_model(N_T_EVAL))
    allocate(R_model(N_T_EVAL))

    call integrate_model( &
       alpha,season_amp, &
       weeks_model,R_model,n_valid,good)

    if(.not.good .or. n_valid < 10) then

      deallocate(weeks_model,R_model)
      return

    end if

    do i=1,size(weeks_obs)

      xq=weeks_obs(i)-shift

      if(xq < weeks_model(1) .or. &
         xq > weeks_model(n_valid)) then

        deallocate(weeks_model,R_model)
        return

      end if

      do j=1,n_valid-1

        if(xq >= weeks_model(j) .and. &
           xq <= weeks_model(j+1)) then

          w=(xq-weeks_model(j))/ &
            (weeks_model(j+1)-weeks_model(j))

          R_pred(i)=R_model(j)+ &
            w*(R_model(j+1)-R_model(j))

          exit

        end if

      end do

    end do

    ok=.true.

    deallocate(weeks_model,R_model)

  end subroutine evaluate_model_at_shift


!-------------------------------------------------------------------------------
! 4-D Gaussian likelihood
!-------------------------------------------------------------------------------

  function log_likelihood_4d(theta,w_obs,r_obs) result(ll)

    real(dp),intent(in) :: theta(4)
    real(dp),intent(in) :: w_obs(:)
    real(dp),intent(in) :: r_obs(:)

    real(dp) :: ll
    real(dp) :: R_pred(size(r_obs))
    real(dp) :: res,sigma

    integer :: i

    logical :: ok

    sigma=theta(3)

    if(sigma <= 0.0_dp) then

      ll=-huge(1.0_dp)
      return

    end if

    call evaluate_model_at_shift( &
       theta(1),theta(2),theta(4), &
       w_obs,R_pred,ok)

    if(.not.ok) then

      ll=-huge(1.0_dp)
      return

    end if

    ll=0.0_dp

    do i=1,size(r_obs)

      res=r_obs(i)-R_pred(i)

      ll=ll-0.5_dp*( &
          (res/sigma)**2 + &
          log(2.0_dp*PI*sigma**2) )

    end do

  end function log_likelihood_4d


  function log_prior_4d(theta) result(lp)

    real(dp),intent(in) :: theta(4)

    real(dp) :: lp

    if(theta(1) >= ALPHA_LO .and. &
       theta(1) <= ALPHA_HI .and. &
       theta(2) >= SEASON_LO .and. &
       theta(2) <= SEASON_HI .and. &
       theta(3) >= SIGMA_LO .and. &
       theta(3) <= SIGMA_HI .and. &
       theta(4) >= SHIFT_LO .and. &
       theta(4) <= SHIFT_HI) then

      lp=0.0_dp

    else

      lp=-huge(1.0_dp)

    end if

  end function log_prior_4d


  function log_posterior_4d(theta) result(lp)

    real(dp),intent(in) :: theta(4)

    real(dp) :: lp,ll

    lp=log_prior_4d(theta)

    if(lp < -1.0e100_dp) return

    ll=log_likelihood_4d( &
       theta,WEEKS_DATA,R_DATA)

    if(ll < -1.0e100_dp) then

      lp=-huge(1.0_dp)

    else

      lp=lp+ll

    end if

  end function log_posterior_4d


  function negative_loglike_4d(theta) result(f)

    real(dp),intent(in) :: theta(4)

    real(dp) :: f,ll

    if(log_prior_4d(theta) < -1.0e100_dp) then

      f=1.0e100_dp
      return

    end if

    ll=log_likelihood_4d( &
       theta,WEEKS_DATA,R_DATA)

    if(ll < -1.0e100_dp) then

      f=1.0e100_dp

    else

      f=-ll

    end if

  end function negative_loglike_4d

end module likelihood_mod


!===============================================================================
! OPTIMIZATION
!===============================================================================

module optimization_mod

  use precision_mod
  use settings_mod
  use likelihood_mod
  use rng_mod
  implicit none

  abstract interface

    function cost_func(x) result(f)

      use precision_mod

      real(dp),intent(in) :: x(4)

      real(dp) :: f

    end function cost_func

  end interface

contains


  subroutine nelder_mead_4d( &
       x0,func,xbest,fbest,maxiter,ftol)

    real(dp),intent(in) :: x0(4)

    procedure(cost_func) :: func

    real(dp),intent(out) :: xbest(4)
    real(dp),intent(out) :: fbest

    integer,intent(in) :: maxiter

    real(dp),intent(in),optional :: ftol

    real(dp) :: x(4,5)
    real(dp) :: f(5)

    real(dp) :: cent(4),xr(4),xe(4),xc(4)
    real(dp) :: fr,fe,fc,tol,scale

    integer :: j,it,ilo,ihi,inhi

    tol=1.0e-6_dp

    if(present(ftol)) tol=ftol

    x(:,1)=x0

    x(:,2)=x0+[0.02_dp,0.0_dp,0.0_dp,0.0_dp]
    x(:,3)=x0+[0.0_dp,0.02_dp,0.0_dp,0.0_dp]
    x(:,4)=x0+[0.0_dp,0.0_dp,0.01_dp,0.0_dp]
    x(:,5)=x0+[0.0_dp,0.0_dp,0.0_dp,0.05_dp]

    do j=1,5

      x(:,j)=min(max(x(:,j), &
           [ALPHA_LO,SEASON_LO,SIGMA_LO,SHIFT_LO]), &
           [ALPHA_HI,SEASON_HI,SIGMA_HI,SHIFT_HI])

      f(j)=func(x(:,j))

    end do

    do it=1,maxiter

      ilo=1
      ihi=1

      do j=2,5

        if(f(j)<f(ilo)) ilo=j
        if(f(j)>f(ihi)) ihi=j

      end do

      inhi=1

      do j=1,5

        if(j /= ihi) then

          if(inhi == ihi .or. &
             f(j)>f(inhi)) inhi=j

        end if

      end do

      scale=0.0_dp

      do j=1,5

        scale=max(scale, &
          sum(abs(x(:,j)-x(:,ilo))))

      end do

      if(scale < 1.0e-5_dp .and. &
         abs(f(ihi)-f(ilo)) < tol) exit

      cent=0.0_dp

      do j=1,5

        if(j /= ihi) &
          cent=cent+x(:,j)

      end do

      cent=cent/4.0_dp

      xr=cent+(cent-x(:,ihi))

      xr=min(max(xr, &
         [ALPHA_LO,SEASON_LO,SIGMA_LO,SHIFT_LO]), &
         [ALPHA_HI,SEASON_HI,SIGMA_HI,SHIFT_HI])

      fr=func(xr)

      if(fr < f(ilo)) then

        xe=cent+2.0_dp*(xr-cent)

        xe=min(max(xe, &
           [ALPHA_LO,SEASON_LO,SIGMA_LO,SHIFT_LO]), &
           [ALPHA_HI,SEASON_HI,SIGMA_HI,SHIFT_HI])

        fe=func(xe)

        if(fe < fr) then

          x(:,ihi)=xe
          f(ihi)=fe

        else

          x(:,ihi)=xr
          f(ihi)=fr

        end if

      else if(fr < f(inhi)) then

        x(:,ihi)=xr
        f(ihi)=fr

      else

        if(fr < f(ihi)) then

          xc=cent+0.5_dp*(xr-cent)

        else

          xc=cent+0.5_dp*(x(:,ihi)-cent)

        end if

        xc=min(max(xc, &
           [ALPHA_LO,SEASON_LO,SIGMA_LO,SHIFT_LO]), &
           [ALPHA_HI,SEASON_HI,SIGMA_HI,SHIFT_HI])

        fc=func(xc)

        if(fc < f(ihi)) then

          x(:,ihi)=xc
          f(ihi)=fc

        else

          do j=1,5

            if(j /= ilo) then

              x(:,j)=x(:,ilo)+ &
                     0.5_dp*(x(:,j)-x(:,ilo))

              f(j)=func(x(:,j))

            end if

          end do

        end if

      end if

    end do

    ilo=minloc(f,1)

    xbest=x(:,ilo)
    fbest=f(ilo)

  end subroutine nelder_mead_4d


  subroutine make_mcmc_start_4d(center,pos)

    real(dp),intent(in) :: center(4)

    real(dp),intent(out) :: pos(NWALKERS,4)

    integer :: i
    integer :: seed_local

    real(dp) :: z(4)

    seed_local=12345

    do i=1,NWALKERS

      do

        z(1)=normal_rand_seeded(seed_local)
        z(2)=normal_rand_seeded(seed_local)
        z(3)=normal_rand_seeded(seed_local)
        z(4)=normal_rand_seeded(seed_local)

        pos(i,:)=center+ &
          [0.03_dp,0.03_dp,0.01_dp,0.15_dp]*z

        if(log_prior_4d(pos(i,:)) > -1.0e100_dp) exit

      end do

    end do

  end subroutine make_mcmc_start_4d

end module optimization_mod


!===============================================================================
! MCMC
!===============================================================================

module mcmc_mod

  use precision_mod
  use settings_mod
  use likelihood_mod
  use rng_mod
  use optimization_mod
  implicit none

contains


  subroutine run_mcmc_4d( &
       start,flat_samples,acceptance,rhat,walker_acceptance)

    real(dp),intent(in) :: start(4)

    real(dp),allocatable,intent(out) :: flat_samples(:,:)

    real(dp),intent(out) :: acceptance
    real(dp),intent(out) :: rhat(4)
    real(dp),intent(out) :: walker_acceptance(NWALKERS)

    real(dp) :: chain(NWALKERS,4)
    real(dp) :: lp(NWALKERS)

    real(dp) :: proposal(4)
    real(dp) :: newlp,u,z,a

    integer :: i,j,k
    integer :: naccept,nprop
    integer :: walker_accept(NWALKERS),walker_prop(NWALKERS)
    integer :: seed_local
    integer :: n_keep,idx
    real(dp),allocatable :: saved_chain(:,:,:)

    seed_local=9173

    naccept=0
    nprop=0
    walker_accept=0
    walker_prop=0
    rhat=1.0_dp
    walker_acceptance=0.0_dp

    call make_mcmc_start_4d(start,chain)

    do i=1,NWALKERS

      lp(i)=log_posterior_4d(chain(i,:))

    end do

    n_keep=NWALKERS*((NSTEPS-BURN_MCMC)/THIN)

    allocate(flat_samples(n_keep,4))
    allocate(saved_chain(n_keep/NWALKERS,NWALKERS,4))

    idx=0

    do k=1,NSTEPS

      do i=1,NWALKERS

        do

          j=1+int( &
             real(NWALKERS,dp)* &
             uniform_rand_seeded(seed_local))

          if(j<1) j=1
          if(j>NWALKERS) j=NWALKERS

          if(j /= i) exit

        end do

        a=2.0_dp

        z=((a-1.0_dp)* &
           uniform_rand_seeded(seed_local)+1.0_dp)**2/a

        proposal=chain(j,:)+ &
           z*(chain(i,:)-chain(j,:))

        newlp=log_posterior_4d(proposal)

        nprop=nprop+1
        walker_prop(i)=walker_prop(i)+1

        if(newlp > -1.0e100_dp) then

          u=log(max( &
             uniform_rand_seeded(seed_local), &
             1.0e-15_dp))

          if(u < 3.0_dp*log(z)+ &
             newlp-lp(i)) then

            chain(i,:)=proposal
            lp(i)=newlp

            naccept=naccept+1
            walker_accept(i)=walker_accept(i)+1

          end if

        end if

      end do

      if(k > BURN_MCMC .and. &
         mod(k-BURN_MCMC,THIN)==0) then

        do i=1,NWALKERS

          idx=idx+1

          flat_samples(idx,:)=chain(i,:)
          saved_chain((k-BURN_MCMC)/THIN,i,:)=chain(i,:)

        end do

      end if

      if(mod(k,2000)==0 .or. k==NSTEPS) then

        write(*,'(A,I7,A,I7,A,F6.2,A)') &
          'MCMC: ',k,'/',NSTEPS,' (', &
          100.0_dp*real(k,dp)/real(NSTEPS,dp),'%)'

      end if

    end do

    acceptance=real(naccept,dp)/ &
               real(max(1,nprop),dp)

    do i=1,NWALKERS
      walker_acceptance(i)=real(walker_accept(i),dp)/ &
                           real(max(1,walker_prop(i)),dp)
    end do

    call compute_split_rhat(saved_chain,rhat)

    deallocate(saved_chain)

  end subroutine run_mcmc_4d


  subroutine compute_split_rhat(samples,rhat)

    real(dp),intent(in) :: samples(:,:,:)
    real(dp),intent(out) :: rhat(4)

    integer :: ns,nw,p,nhalf,m,ii,jj,kk
    real(dp) :: chain_mean(2*NWALKERS)
    real(dp) :: chain_var(2*NWALKERS)
    real(dp) :: grand_mean,w,b,var_plus

    ns=size(samples,1)
    nw=size(samples,2)
    p=size(samples,3)

    rhat=1.0_dp

    if(nw /= NWALKERS .or. p /= 4) return
    if(ns < 4) return

    nhalf=ns/2
    if(nhalf < 2) return

    m=2*nw

    do kk=1,p

      do ii=1,nw

        jj=2*ii-1
        chain_mean(jj)=sum(samples(1:nhalf,ii,kk))/real(nhalf,dp)
        chain_var(jj)=sum((samples(1:nhalf,ii,kk)-chain_mean(jj))**2)/ &
                      real(nhalf-1,dp)

        jj=2*ii
        chain_mean(jj)=sum(samples(nhalf+1:2*nhalf,ii,kk))/real(nhalf,dp)
        chain_var(jj)=sum((samples(nhalf+1:2*nhalf,ii,kk)-chain_mean(jj))**2)/ &
                      real(nhalf-1,dp)

      end do

      grand_mean=sum(chain_mean)/real(m,dp)

      w=sum(chain_var)/real(m,dp)

      b=real(nhalf,dp)*sum((chain_mean-grand_mean)**2)/real(m-1,dp)

      if(w > 0.0_dp) then
        var_plus=(real(nhalf-1,dp)/real(nhalf,dp))*w + &
                 b/real(nhalf,dp)
        rhat(kk)=sqrt(max(var_plus/w,0.0_dp))
      else
        if(b <= 0.0_dp) then
          rhat(kk)=1.0_dp
        else
          rhat(kk)=huge(1.0_dp)
        end if
      end if

    end do

  end subroutine compute_split_rhat

end module mcmc_mod


!===============================================================================
! BOOTSTRAP
!===============================================================================

module bootstrap_mod

  use precision_mod
  use settings_mod
  use data_mod
  use likelihood_mod
  use optimization_mod
  use rng_mod
  use stats_mod
  implicit none

contains


  subroutine fit_single_bootstrap( &
       mle_curve,mle_start,seed,beta, &
       boot_corr,boot_rmse,boot_curve,success)

    real(dp),intent(in) :: mle_curve(NDATA)
    real(dp),intent(in) :: mle_start(4)

    integer,intent(in) :: seed

    real(dp),intent(out) :: beta(4)
    real(dp),intent(out) :: boot_corr
    real(dp),intent(out) :: boot_rmse
    real(dp),intent(out) :: boot_curve(NDATA)

    logical,intent(out) :: success

    real(dp) :: R_boot(NDATA)
    real(dp) :: residuals(NDATA)

    real(dp) :: theta0(4)
    real(dp) :: xbest(4)

    real(dp) :: fbest
    real(dp) :: best_R(NDATA)

    integer :: i,idx,sl

    logical :: ok

    sl=seed

    success=.false.

    beta=0.0_dp
    boot_corr=0.0_dp
    boot_rmse=0.0_dp
    boot_curve=0.0_dp

    residuals=R_DATA-mle_curve

    do i=1,NDATA

      idx=1+int( &
        real(NDATA,dp)* &
        uniform_rand_seeded(sl))

      if(idx<1) idx=1
      if(idx>NDATA) idx=NDATA

      R_boot(i)=mle_curve(i)+residuals(idx)

    end do

    theta0=mle_start

    theta0(1)=theta0(1)+ &
       0.01_dp*normal_rand_seeded(sl)

    theta0(2)=theta0(2)+ &
       0.01_dp*normal_rand_seeded(sl)

    theta0(3)=theta0(3)+ &
       0.005_dp*normal_rand_seeded(sl)

    theta0(4)=theta0(4)+ &
       0.05_dp*normal_rand_seeded(sl)

    theta0=min(max(theta0, &
       [ALPHA_LO,SEASON_LO,SIGMA_LO,SHIFT_LO]), &
       [ALPHA_HI,SEASON_HI,SIGMA_HI,SHIFT_HI])

    call nelder_mead_4d( &
       theta0,boot_cost,xbest,fbest, &
       BOOTSTRAP_MAXITER,1.0e-7_dp)

    if(fbest > 1.0e99_dp) return

    call evaluate_model_at_shift( &
       xbest(1),xbest(2),xbest(4), &
       WEEKS_DATA,best_R,ok)

    if(.not.ok) return

    beta=xbest
    boot_curve=best_R
    boot_corr=pearson(R_DATA,best_R)
    boot_rmse=sqrt(sum((R_DATA-best_R)**2)/real(NDATA,dp))

    success=.true.

  contains


    function boot_cost(theta) result(f)

      real(dp),intent(in) :: theta(4)

      real(dp) :: f,ll

      if(log_prior_4d(theta) < -1.0e100_dp) then

        f=1.0e100_dp

        return

      end if

      ll=log_likelihood_4d( &
         theta,WEEKS_DATA,R_boot)

      if(ll < -1.0e100_dp) then

        f=1.0e100_dp

      else

        f=-ll

      end if

    end function boot_cost

  end subroutine fit_single_bootstrap

end module bootstrap_mod


!===============================================================================
! MAIN PROGRAM
!===============================================================================

program Sars-Cov2_lyon

  use precision_mod
  use settings_mod
  use data_mod
  use contact_mod
  use rng_mod
  use likelihood_mod
  use optimization_mod
  use mcmc_mod
  use bootstrap_mod
  use stats_mod

  use omp_lib

  implicit none


  real(dp) :: initial_guess(4)
  real(dp) :: best_theta(4)
  real(dp) :: xbest(4)

  real(dp) :: fbest
  real(dp) :: best_ll

  real(dp) :: best_R(NDATA)
  real(dp) :: mle_curve(NDATA)

  real(dp) :: mle_rmse
  real(dp) :: mle_corr

  logical :: ok


  real(dp),allocatable :: flat_samples(:,:)

  real(dp),allocatable :: posterior_models(:,:)
  real(dp),allocatable :: posterior_median(:)
  real(dp),allocatable :: posterior_lower(:)
  real(dp),allocatable :: posterior_upper(:)

  real(dp),allocatable :: posterior_corr_dist(:)
  real(dp),allocatable :: posterior_rmse_dist(:)

  real(dp) :: corr_025,corr_50,corr_975
  real(dp) :: rmse_025,rmse_50,rmse_975

  real(dp),allocatable :: bootstrap_params(:,:)
  real(dp),allocatable :: bootstrap_corr(:)
  real(dp),allocatable :: bootstrap_rmse(:)
  real(dp),allocatable :: bootstrap_curves(:,:)
  logical,allocatable :: bootstrap_success(:)

  real(dp),allocatable :: success_params(:,:)
  real(dp),allocatable :: success_corr(:)
  real(dp),allocatable :: success_rmse(:)
  real(dp),allocatable :: success_curves(:,:)

  real(dp),allocatable :: temp_curve(:)

  real(dp) :: mean_val,std_val
  real(dp) :: q025,q50,q975

  real(dp) :: p_boot(5)
  real(dp) :: corr_boot(5)

  real(dp) :: boot_alpha_ci(2)
  real(dp) :: boot_season_ci(2)
  real(dp) :: boot_sigma_ci(2)
  real(dp) :: boot_shift_ci(2)

  real(dp) :: acceptance
  real(dp) :: rhat(4)
  real(dp) :: walker_acceptance(NWALKERS)
  real(dp) :: max_rhat
  logical :: convergence_ok

  integer :: i,j,idx
  integer :: n_valid_post
  integer :: n_success
  integer :: boot_n_success
  integer :: ios

  integer :: posterior_seed

  character(len=16) :: labels(4)

  labels(1)='alpha'
  labels(2)='season_amp'
  labels(3)='sigma'
  labels(4)='shift_weeks'


  !===========================================================================
  ! INITIALIZATION
  !===========================================================================

  call initialize_contacts()

  call seed_rng(42)

  posterior_seed=314159


  write(*,*)
  write(*,*) '==============================================================='
  write(*,*) ' Sars-Cov2 LYON - COMPLETE BAYESIAN CALIBRATION'
  write(*,*) '==============================================================='
  write(*,'(A,I6)') ' Data points       = ',NDATA
  write(*,'(A,I6)') ' Walkers           = ',NWALKERS
  write(*,'(A,I7)') ' MCMC steps        = ',NSTEPS
  write(*,'(A,I7)') ' Burn-in           = ',BURN_MCMC
  write(*,'(A,I6)') ' Thin              = ',THIN
  write(*,'(A,I6)') ' Bootstrap         = ',N_BOOTSTRAP
  write(*,'(A)') ' Parameters        = alpha, season_amp, sigma, shift'
  write(*,'(A)') ' Shift treatment   = sampled nuisance parameter'
  write(*,'(A,F8.2,A,F8.2,A)') &
      ' Shift prior       = [',SHIFT_LO,',',SHIFT_HI,'] weeks'
  write(*,*) '==============================================================='
  write(*,*)


  !===========================================================================
  ! MLE
  !===========================================================================

  write(*,*) '--- MAXIMUM-LIKELIHOOD FIT ---'

  initial_guess=[ &
       0.8_dp, &
       0.35_dp, &
       0.15_dp, &
       1.8_dp ]

  call nelder_mead_4d( &
       initial_guess, &
       negative_loglike_4d, &
       xbest,fbest, &
       1000,1.0e-8_dp)

  best_theta=xbest

  best_ll=-fbest

  call evaluate_model_at_shift( &
       best_theta(1), &
       best_theta(2), &
       best_theta(4), &
       WEEKS_DATA,best_R,ok)

  if(.not.ok) then

    write(*,*) 'ERROR: MLE evaluation failed.'
    stop

  end if

  mle_curve=best_R

  mle_rmse=sqrt( &
     sum((R_DATA-mle_curve)**2)/real(NDATA,dp))

  mle_corr=pearson(R_DATA,mle_curve)

  write(*,'(A,F12.6)') &
     'alpha         = ',best_theta(1)

  write(*,'(A,F12.6)') &
     'season_amp    = ',best_theta(2)

  write(*,'(A,F12.6)') &
     'sigma         = ',best_theta(3)

  write(*,'(A,F12.6)') &
     'shift (weeks) = ',best_theta(4)

  write(*,'(A,F12.6)') &
     'logLik        = ',best_ll

  write(*,'(A,F12.6)') &
     'RMSE          = ',mle_rmse

  write(*,'(A,F12.6)') &
     'Pearson r     = ',mle_corr

  write(*,*)


  !===========================================================================
  ! MCMC
  !===========================================================================

  write(*,*) '--- BAYESIAN 4-D MCMC ---'

  call run_mcmc_4d( &
       best_theta,flat_samples,acceptance,rhat,walker_acceptance)

  write(*,'(A,F8.4)') &
       'Acceptance fraction = ',acceptance

  write(*,'(A,I10)') &
       'Posterior samples  = ',size(flat_samples,1)

  !---------------------------------------------------------------------------
  ! ADDITION ONLY: MCMC diagnostics. The sampler and existing calculations
  ! above are unchanged; this only records diagnostics from the same chains.
  !---------------------------------------------------------------------------
  max_rhat=maxval(rhat)
  convergence_ok=(max_rhat <= 1.05_dp)

  write(*,*)
  write(*,*) '--- MCMC DIAGNOSTICS (ADDED) ---'
  write(*,'(A,F10.6)') 'Overall acceptance fraction = ',acceptance
  write(*,*) 'Per-walker acceptance fractions:'
  do i=1,NWALKERS
    write(*,'(A,I4,A,F10.6)') 'Walker ',i,' : ',walker_acceptance(i)
  end do
  write(*,*) 'Split-R-hat:'
  do i=1,4
    write(*,'(A14,F12.6)') labels(i),rhat(i)
  end do
  write(*,'(A,F12.6)') 'Maximum R-hat = ',max_rhat
  if(convergence_ok) then
    write(*,*) 'Convergence assessment = GOOD (max R-hat <= 1.05)'
  else
    write(*,*) 'Convergence assessment = NOT CONVERGED (max R-hat > 1.05)'
  end if

  write(*,*)


  !===========================================================================
  ! POSTERIOR PARAMETERS
  !===========================================================================

  write(*,*) '--- POSTERIOR PARAMETERS ---'

  do i=1,4

    call summary_stats( &
       flat_samples(:,i), &
       mean_val,std_val, &
       q025,q50,q975)

    write(*,'(A14,A,F10.5,A,F10.5,A,F10.5,A)') &
       labels(i),' median=',q50, &
       ' 95% CrI=[',q025,',',q975,']'

  end do

  write(*,*)


  !===========================================================================
  ! POSTERIOR CORRELATION MATRIX
  !===========================================================================

  write(*,*) '--- POSTERIOR PARAMETER CORRELATION MATRIX ---'
  
  do i=1,4

    write(*,'(4F12.5)') &
      (pearson(flat_samples(:,i), &
               flat_samples(:,j)),j=1,4)

  end do

  write(*,*)


  !===========================================================================
  ! POSTERIOR TRAJECTORIES
  !===========================================================================

  write(*,*) '--- POSTERIOR TRAJECTORY EVALUATION ---'

  allocate(temp_curve(NDATA))

  allocate(posterior_models( &
       NDATA,N_POSTERIOR_MODELS))

  allocate(posterior_corr_dist( &
       N_POSTERIOR_MODELS))

  allocate(posterior_rmse_dist( &
       N_POSTERIOR_MODELS))

  n_valid_post=0

  do i=1,N_POSTERIOR_MODELS

    !-----------------------------------------------------------------------
    ! IMPORTANT:
    ! The original code used "idx" itself as the random seed before idx
    ! was initialized. This version uses a dedicated initialized seed.
    !-----------------------------------------------------------------------

    idx=1+int( &
       real(size(flat_samples,1),dp)* &
       uniform_rand_seeded(posterior_seed))

    idx=max(1,min(idx,size(flat_samples,1)))

    call evaluate_model_at_shift( &
       flat_samples(idx,1), &
       flat_samples(idx,2), &
       flat_samples(idx,4), &
       WEEKS_DATA,temp_curve,ok)

    if(ok) then

      n_valid_post=n_valid_post+1

      posterior_models(:,n_valid_post)=temp_curve

      posterior_corr_dist(n_valid_post)= &
         pearson(R_DATA,temp_curve)

      posterior_rmse_dist(n_valid_post)= &
         sqrt(sum((R_DATA-temp_curve)**2)/ &
              real(NDATA,dp))

    end if

  end do


  write(*,'(A,I6)') &
     'Valid posterior trajectories = ',n_valid_post


  if(n_valid_post > 0) then

    allocate(posterior_median(NDATA))
    allocate(posterior_lower(NDATA))
    allocate(posterior_upper(NDATA))

    do i=1,NDATA

      call percentile( &
         posterior_models(i,1:n_valid_post), &
         0.50_dp,posterior_median(i))

      call percentile( &
         posterior_models(i,1:n_valid_post), &
         0.025_dp,posterior_lower(i))

      call percentile( &
         posterior_models(i,1:n_valid_post), &
         0.975_dp,posterior_upper(i))

    end do


    call percentile( &
       posterior_corr_dist(1:n_valid_post), &
       0.025_dp,corr_025)

    call percentile( &
       posterior_corr_dist(1:n_valid_post), &
       0.50_dp,corr_50)

    call percentile( &
       posterior_corr_dist(1:n_valid_post), &
       0.975_dp,corr_975)


    call percentile( &
       posterior_rmse_dist(1:n_valid_post), &
       0.025_dp,rmse_025)

    call percentile( &
       posterior_rmse_dist(1:n_valid_post), &
       0.50_dp,rmse_50)

    call percentile( &
       posterior_rmse_dist(1:n_valid_post), &
       0.975_dp,rmse_975)


    write(*,*)
    write(*,*) 'Posterior Pearson r:'

    write(*,'(A,F10.6)') &
       'Median = ',corr_50

    write(*,'(A,F10.6,A,F10.6,A)') &
       '95% CrI = [',corr_025,',',corr_975,']'


    write(*,*)
    write(*,*) 'Posterior RMSE:'

    write(*,'(A,F10.6)') &
       'Median = ',rmse_50

    write(*,'(A,F10.6,A,F10.6,A)') &
       '95% CrI = [',rmse_025,',',rmse_975,']'

  else

    write(*,*) 'ERROR: No valid posterior trajectories were generated.'
    write(*,*) '       Posterior CrI for r/RMSE cannot be constructed.'
    stop 2

  end if

  deallocate(temp_curve)


  !===========================================================================
  ! RESIDUAL BOOTSTRAP
  !===========================================================================

  write(*,*)
  write(*,*) '--- RESIDUAL BOOTSTRAP ---'

  allocate(bootstrap_params(N_BOOTSTRAP,4))
  allocate(bootstrap_corr(N_BOOTSTRAP))
  allocate(bootstrap_rmse(N_BOOTSTRAP))
  allocate(bootstrap_curves(NDATA,N_BOOTSTRAP))
  allocate(bootstrap_success(N_BOOTSTRAP))

  bootstrap_params=0.0_dp
  bootstrap_corr=0.0_dp
  bootstrap_rmse=0.0_dp
  bootstrap_curves=0.0_dp
  bootstrap_success=.false.


  !$omp parallel do default(shared) &
  !$omp private(i,xbest,ok,q50) schedule(dynamic)

  do i=1,N_BOOTSTRAP

    call fit_single_bootstrap( &
       best_R,best_theta, &
       100000+i*7919, &
       xbest,q50,bootstrap_rmse(i),bootstrap_curves(:,i),ok)

    if(ok) then

      bootstrap_params(i,:)=xbest
      bootstrap_corr(i)=q50
      bootstrap_success(i)=.true.

    end if

  end do

  !$omp end parallel do


  n_success=count(bootstrap_success)

  write(*,'(A,I5,A,I5,A,F7.2,A)') &
     'Bootstrap success = ', &
     n_success,'/',N_BOOTSTRAP, &
     ' (',100.0_dp*real(n_success,dp)/ &
     real(N_BOOTSTRAP,dp),'%)'


  if(n_success > 0) then

    allocate(success_params(n_success,4))
    allocate(success_corr(n_success))
    allocate(success_rmse(n_success))
    allocate(success_curves(NDATA,n_success))

    idx=0

    do i=1,N_BOOTSTRAP

      if(bootstrap_success(i)) then

        idx=idx+1

        success_params(idx,:)= &
             bootstrap_params(i,:)

        success_corr(idx)=bootstrap_corr(i)
        success_rmse(idx)=bootstrap_rmse(i)
        success_curves(:,idx)=bootstrap_curves(:,i)

      end if

    end do


    write(*,*)
    write(*,*) 'Bootstrap 95% confidence intervals:'


    do i=1,4

      call summary_stats( &
         success_params(:,i), &
         mean_val,std_val, &
         q025,q50,q975)

      write(*,'(A14,A,F10.5,A,F10.5,A)') &
         labels(i),' = [',q025,',',q975,']'

      if(i == 1) then

        boot_alpha_ci=[q025,q975]

      else if(i == 2) then

        boot_season_ci=[q025,q975]

      else if(i == 3) then

        boot_sigma_ci=[q025,q975]

      else

        boot_shift_ci=[q025,q975]

      end if

    end do


    call summary_stats( &
       success_corr, &
       mean_val,std_val, &
       corr_boot(3),corr_boot(4),corr_boot(5))

    call summary_stats(success_rmse,mean_val,std_val, &
         q025,q50,q975)

    write(*,'(A14,A,F10.5,A,F10.5,A)') 'RMSE',' = [',q025,',',q975,']'


    write(*,'(A14,A,F10.5,A,F10.5,A)') &
       'Pearson r',' = [', &
       corr_boot(3),',',corr_boot(5),']'


    !-----------------------------------------------------------------------
    ! Save bootstrap samples
    !-----------------------------------------------------------------------

    open(unit=12, &
         file='Sars-Cov2_lyon_MLE_Bootstrap_Samples.dat', &
         status='replace',iostat=ios)

    if(ios == 0) then

      write(12,'(A)') &
        '# alpha season_amp sigma shift Pearson_r RMSE'

      do i=1,n_success

        write(12,'(6ES20.10)') &
          success_params(i,:),success_corr(i),success_rmse(i)

      end do

      close(12)

    end if


    boot_n_success=n_success

    deallocate(success_params)
    deallocate(success_corr)

  else

    boot_n_success=0

    write(*,*) &
      'WARNING: No successful bootstrap replicates.'

  end if


  !===========================================================================
  ! SAVE MCMC DIAGNOSTICS -- NEW FILE; EXISTING OUTPUT FILES UNCHANGED
  !===========================================================================

  open(unit=41, &
       file='Sars-Cov2_lyon_MCMC_Diagnostics_ADDED.dat', &
       status='replace',iostat=ios)

  if(ios == 0) then

    write(41,'(A)') '# Sars-Cov2 LYON - MCMC diagnostics (added only)'
    write(41,'(A)') '# Existing model/MCMC settings and posterior samples are unchanged.'
    write(41,'(A,F20.10)') 'overall_acceptance ',acceptance
    write(41,'(A,I10)') 'n_walkers ',NWALKERS
    write(41,'(A,I10)') 'n_steps ',NSTEPS
    write(41,'(A,I10)') 'burn_in ',BURN_MCMC
    write(41,'(A,I10)') 'thin ',THIN
    write(41,'(A,I10)') 'retained_steps_per_walker ',(NSTEPS-BURN_MCMC)/THIN
    write(41,'(A)') ''
    write(41,'(A)') '# walker_id acceptance_fraction'
    do i=1,NWALKERS
      write(41,'(I6,1X,ES20.10)') i,walker_acceptance(i)
    end do
    write(41,'(A)') ''
    write(41,'(A)') '# parameter Rhat'
    do i=1,4
      write(41,'(A16,1X,ES20.10)') labels(i),rhat(i)
    end do
    write(41,'(A)') ''
    write(41,'(A,1X,ES20.10)') 'max_Rhat',max_rhat
    if(convergence_ok) then
      write(41,'(A)') 'convergence GOOD'
    else
      write(41,'(A)') 'convergence NOT_CONVERGED'
    end if
    close(41)

  else
    write(*,*) 'WARNING: Could not write MCMC diagnostics file.'
  end if


  !===========================================================================
  ! SAVE POSTERIOR SAMPLES
  !===========================================================================

  open(unit=11, &
       file='Sars-Cov2_lyon_MLE_Posterior_Samples.dat', &
       status='replace',iostat=ios)

  if(ios == 0) then

    write(11,'(A)') &
      '# alpha season_amp sigma shift_weeks'
    
    do i=1,size(flat_samples,1)

      write(11,'(4ES20.10)') &
        flat_samples(i,:)

    end do

    close(11)

  end if


  !===========================================================================
  ! SAVE MLE
  !===========================================================================

  open(unit=13, &
       file='Sars-Cov2_lyon_MLE_Results.txt', &
       status='replace',iostat=ios)

  if(ios == 0) then

    write(13,'(A,1X,ES20.10)') 'alpha',best_theta(1)
    write(13,'(A,1X,ES20.10)') 'season_amp',best_theta(2)
    write(13,'(A,1X,ES20.10)') 'sigma',best_theta(3)
    write(13,'(A,1X,ES20.10)') 'shift',best_theta(4)
    write(13,'(A,1X,ES20.10)') 'log_likelihood',best_ll
    write(13,'(A,1X,ES20.10)') 'MLE_RMSE',mle_rmse
    write(13,'(A,1X,ES20.10)') 'MLE_Pearson_r',mle_corr
    write(13,'(A,1X,ES20.10)') 'Posterior_r_median',corr_50
    write(13,'(A,1X,ES20.10)') 'Posterior_r_lower',corr_025
    write(13,'(A,1X,ES20.10)') 'Posterior_r_upper',corr_975
    write(13,'(A,1X,ES20.10)') 'Posterior_RMSE_median',rmse_50
    write(13,'(A,1X,ES20.10)') 'Posterior_RMSE_lower',rmse_025
    write(13,'(A,1X,ES20.10)') 'Posterior_RMSE_upper',rmse_975

    close(13)

  end if


  !===========================================================================
  ! SAVE FIT CURVE
  !===========================================================================

  open(unit=20, &
       file='Sars-Cov2_lyon_MLE_Fit_Results.csv', &
       status='replace',iostat=ios)

  if(ios == 0) then

    write(20,'(A)') &
      'week,r_data,mle_curve,posterior_median,'// &
      'posterior_lower,posterior_upper'

    do i=1,NDATA

      write(20,'(F8.4,5(",",F10.6))') &
        WEEKS_DATA(i), &
        R_DATA(i), &
        mle_curve(i), &
        posterior_median(i), &
        posterior_lower(i), &
        posterior_upper(i)

    end do

    close(20)

  end if


  !===========================================================================
  ! SAVE POSTERIOR CORRELATION DISTRIBUTION
  !===========================================================================

  if(n_valid_post > 0) then

    open(unit=30, &
         file='Sars-Cov2_lyon_MLE_Posterior_Performance.txt', &
         status='replace',iostat=ios)

    if(ios == 0) then

      write(30,'(A)') '# Pearson_r RMSE'

      do i=1,n_valid_post

        write(30,'(2ES20.10)') &
          posterior_corr_dist(i), &
          posterior_rmse_dist(i)

      end do

      close(30)

    end if

  end if


  !===========================================================================
  ! FINAL SUMMARY
  !===========================================================================

  write(*,*)
  write(*,*) '==============================================================='
  write(*,*) ' Sars-Cov2 LYON - FINAL RESULTS'
  write(*,*) '==============================================================='

  write(*,*)
  write(*,*) 'MLE:'

  write(*,'(A,F12.6)') &
      'alpha         = ',best_theta(1)

  write(*,'(A,F12.6)') &
      'season_amp    = ',best_theta(2)

  write(*,'(A,F12.6)') &
      'sigma         = ',best_theta(3)

  write(*,'(A,F12.6)') &
      'shift (weeks) = ',best_theta(4)

  write(*,'(A,F12.6)') &
      'Pearson r     = ',mle_corr

  write(*,'(A,F12.6)') &
      'RMSE          = ',mle_rmse


  write(*,*)
  write(*,*) 'POSTERIOR 95% CREDIBLE INTERVALS:'

  do i=1,4

    call percentile( &
       flat_samples(:,i), &
       0.025_dp,q025)

    call percentile( &
       flat_samples(:,i), &
       0.975_dp,q975)

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
       labels(i),' = [',q025,',',q975,']'

  end do


  write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'Pearson r',' = [',corr_025,',',corr_975,']'

  write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'RMSE',' = [',rmse_025,',',rmse_975,']'


  if(boot_n_success > 0) then

    write(*,*)
    write(*,*) 'BOOTSTRAP 95% CONFIDENCE INTERVALS:'

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'alpha',' = [', &
      boot_alpha_ci(1),',', &
      boot_alpha_ci(2),']'

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'season_amp',' = [', &
      boot_season_ci(1),',', &
      boot_season_ci(2),']'

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'sigma',' = [', &
      boot_sigma_ci(1),',', &
      boot_sigma_ci(2),']'

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'shift',' = [', &
      boot_shift_ci(1),',', &
      boot_shift_ci(2),']'

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
      'Pearson r',' = [', &
      corr_boot(3),',', &
      corr_boot(5),']'

  end if


  write(*,*)
  write(*,*) 'MCMC convergence:'
  write(*,'(A,F8.4)') &
      'Acceptance fraction = ',acceptance
  do i=1,4
    write(*,'(A14,A,F10.6)') labels(i),' R-hat = ',rhat(i)
  end do
  write(*,'(A,F10.6)') 'Maximum R-hat = ',max_rhat
  if(convergence_ok) then
    write(*,*) 'Convergence assessment = GOOD'
  else
    write(*,*) 'Convergence assessment = NOT CONVERGED'
  end if

  write(*,*)
  write(*,'(A,I8)') &
      'Posterior samples = ',size(flat_samples,1)

  write(*,'(A,I8)') &
      'Valid posterior trajectories = ',n_valid_post

  if(boot_n_success > 0) then

    write(*,'(A,I6,A,I6)') &
      'Bootstrap success = ', &
      boot_n_success,'/',N_BOOTSTRAP

  end if


  write(*,*)
  write(*,*) '==============================================================='
  write(*,*) ' Output files:'
  write(*,*) '   Sars-Cov2_lyon_MCMC_Diagnostics_ADDED.dat'
  write(*,*) '   Sars-Cov2_lyon_MLE_Posterior_Samples.dat'
  write(*,*) '   Sars-Cov2_lyon_MLE_Results.txt'
  write(*,*) '   Sars-Cov2_lyon_MLE_Fit_Results.csv'
  write(*,*) '   Sars-Cov2_lyon_MLE_Posterior_Performance.txt'
  write(*,*) '   Sars-Cov2_lyon_MLE_Bootstrap_Samples.dat'
  write(*,*) '   Sars-Cov2_lyon_MLE_Bootstrap_Trajectories.dat'
  write(*,*) '   Sars-Cov2_lyon_MLE_Posterior_Trajectories.dat'
  write(*,*) '==============================================================='


  !===========================================================================
  ! CLEANUP
  !===========================================================================

  deallocate(flat_samples)

  deallocate(posterior_median)
  deallocate(posterior_lower)
  deallocate(posterior_upper)

  if(allocated(posterior_models)) &
      deallocate(posterior_models)

  if(allocated(posterior_corr_dist)) &
      deallocate(posterior_corr_dist)

  if(allocated(posterior_rmse_dist)) &
      deallocate(posterior_rmse_dist)

  deallocate(bootstrap_params)
  deallocate(bootstrap_corr)
  deallocate(bootstrap_rmse)
  deallocate(bootstrap_curves)
  deallocate(bootstrap_success)


end program Sars-Cov2_lyon
