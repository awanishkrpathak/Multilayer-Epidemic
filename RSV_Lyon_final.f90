!===============================================================================
! RSV LYON
!===============================================================================

module precision_mod
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: i8 = selected_int_kind(18)
  real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
end module precision_mod


module settings_mod
  use precision_mod
  implicit none

  integer, parameter :: NG = 4
  integer, parameter :: NY = 8

  real(dp), parameter :: SHIFT_DURATION = 1.0_dp / 3.0_dp
  real(dp), parameter :: DAYS_PER_WEEK = 7.0_dp

  real(dp), parameter :: SEASONAL_PERIOD_DAYS = 52.0_dp * DAYS_PER_WEEK

  real(dp), parameter :: BURN_IN_WEEKS = 15.0_dp
  real(dp), parameter :: BURN_IN_DAYS  = BURN_IN_WEEKS * DAYS_PER_WEEK

  real(dp), parameter :: MODEL_MAX_WEEK = 80.0_dp
  real(dp), parameter :: MODEL_DAYS = MODEL_MAX_WEEK * DAYS_PER_WEEK

  real(dp), parameter :: SIM_END_DAYS = BURN_IN_DAYS + MODEL_DAYS

  integer, parameter :: N_T_EVAL = 2000
  real(dp), parameter :: MAX_STEP = SHIFT_DURATION / 2.0_dp

  integer, parameter :: NDIM = 4
  integer, parameter :: NWALKERS = 64
  integer, parameter :: NSTEPS = 120000
  integer, parameter :: BURN_MCMC = 35000
  integer, parameter :: THIN = 20
  integer, parameter :: N_POSTERIOR_MODELS = 2000

  integer, parameter :: N_MLE_STARTS = 12
  integer, parameter :: N_PROFILE = 31
  integer, parameter :: MIN_VALID_POSTERIOR = 100

  real(dp), parameter :: ALPHA_LO  = 0.2_dp,  ALPHA_HI  = 1.0_dp
  real(dp), parameter :: SEASON_LO = 0.1_dp, SEASON_HI = 0.99_dp
  real(dp), parameter :: SIGMA_LO  = 0.001_dp, SIGMA_HI = 1.0_dp
  real(dp), parameter :: SHIFT_LO  = -8.0_dp, SHIFT_HI = 8.0_dp

  real(dp), parameter :: SHIFT_PRIOR_SIGMA = 3.0_dp

  real(dp), parameter :: POP(NG) = [11.0_dp, 27.0_dp, 8.0_dp, 29.0_dp]

  real(dp), parameter :: B1(NG, NG) = 0.9_dp*reshape([ &
    1.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 1.4_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 1.2_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [NG, NG])

  real(dp), parameter :: B2(NG, NG) = 1.3_dp*reshape([ &
    1.7_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 1.3_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 1.1_dp], [NG, NG])

  real(dp), parameter :: GAM1(NG) = 1.3*[0.15_dp, 0.14_dp, 0.12_dp, 0.10_dp]
  real(dp), parameter :: GAM2(NG) = [0.14_dp, 0.13_dp, 0.11_dp, 0.09_dp]

  integer, parameter :: NDATA = 11

end module settings_mod


module data_mod
  use precision_mod
  use settings_mod
  implicit none

  real(dp), parameter :: WEEKS_DATA(NDATA) = [ &
    45.0013670539986_dp - 34.0_dp, &
    45.9760765550239_dp - 34.0_dp, &
    46.9931647300068_dp - 34.0_dp, &
    47.9890635680109_dp - 34.0_dp, &
    48.9849624060150_dp - 34.0_dp, &
    50.0232399179768_dp - 34.0_dp, &
    51.0191387559809_dp - 34.0_dp, &
    51.9938482570062_dp - 34.0_dp, &
    53.0109364319891_dp - 34.0_dp, &
    54.0280246069720_dp - 34.0_dp, &
    55.0027341079973_dp - 34.0_dp ]

  real(dp), parameter :: R_DATA(NDATA) = [ &
    1.38744588744589_dp, &
    1.65079365079365_dp, &
    2.32178932178932_dp, &
    1.65079365079365_dp, &
    1.12049062049062_dp, &
    1.02669552669553_dp, &
    0.904040404040404_dp, &
    0.806637806637807_dp, &
    0.730880230880231_dp, &
    0.795815295815296_dp, &
    0.788600288600289_dp ]

end module data_mod


module contact_mod
  use precision_mod
  use settings_mod
  implicit none

  real(dp), save :: C_CONTACT(7, 3, NG, NG)
  real(dp), save :: L_SHIFTS(3, NG, NG)

contains

  subroutine initialize_contacts()

    L_SHIFTS(1,:,:) = 1.2_dp*reshape([ &
      0.50_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.45_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.40_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.20_dp], [NG, NG])

    L_SHIFTS(2,:,:) = 1.2_dp*reshape([ &
      0.40_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.40_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.35_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.18_dp], [NG, NG])

    L_SHIFTS(3,:,:) = 1.2_dp*reshape([ &
      0.30_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.20_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.15_dp], [NG, NG])

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


  subroutine get_contact_and_L(t, C_out, L_out)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: C_out(NG,NG), L_out(NG,NG)
    real(dp) :: t_mod
    integer :: day, shift_idx

    t_mod = modulo(t, 7.0_dp)
    day = int(floor(t_mod))

    if (day < 0) day = 0
    if (day > 6) day = 6

    shift_idx = int(floor((t_mod - floor(t_mod)) / SHIFT_DURATION + 1.0e-12_dp))

    if (shift_idx < 0) shift_idx = 0
    if (shift_idx > 2) shift_idx = 2

    C_out = C_CONTACT(day+1, shift_idx+1, :, :)
    L_out = L_SHIFTS(shift_idx+1, :, :)

  end subroutine get_contact_and_L

end module contact_mod


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

    do i = 1, n
      s(i) = mod(seed + 104729*i, 2147483647)
    end do

    call random_seed(put=s)
    deallocate(s)

  end subroutine seed_rng


  function uniform_rand_seeded(seed) result(u)
    integer, intent(inout) :: seed
    real(dp) :: u
    integer(i8) :: a, m, s

    a = 1103515245_i8
    m = 2147483647_i8

    s = mod(a * int(seed, i8) + 12345_i8, m)

    seed = int(s)
    u = real(seed, dp) / real(m, dp)

  end function uniform_rand_seeded


  function normal_rand_seeded(seed) result(z)
    integer, intent(inout) :: seed
    real(dp) :: z, u1, u2

    u1 = uniform_rand_seeded(seed)
    u2 = uniform_rand_seeded(seed)

    u1 = max(u1, 1.0e-15_dp)

    z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * PI * u2)

  end function normal_rand_seeded

end module rng_mod


module model_mod
  use precision_mod
  use settings_mod
  use data_mod
  use contact_mod
  implicit none

contains

  function seasonal_factor(t, season_amp) result(s)
    real(dp), intent(in) :: t, season_amp
    real(dp) :: s

    s = 1+ season_amp * sin(2.0_dp * PI * (t / DAYS_PER_WEEK) / 52.0_dp)

  end function seasonal_factor


  subroutine rhs(t, y, alpha, season_amp, dy)

    real(dp), intent(in) :: t, y(NY), alpha, season_amp
    real(dp), intent(out) :: dy(NY)

    real(dp) :: I1(NG), I2(NG), S(NG)
    real(dp) :: lam1(NG), lam2(NG)
    real(dp) :: seasonal, Cmat(NG,NG), Lmat(NG,NG)
    real(dp) :: control1(NG,NG), control2(NG,NG)
    real(dp) :: temp1(NG), temp2(NG)

    integer :: j, k

    I1 = max(y(1:NG), 0.0_dp)
    I2 = max(y(NG+1:NY), 0.0_dp)
    S  = max(1.0_dp - I1 - I2, 1.0e-12_dp)

    call get_contact_and_L(t, Cmat, Lmat)

    seasonal = seasonal_factor(t, season_amp)

    control1 = 0.0_dp
    control2 = 0.0_dp

    do j = 1, NG
      do k = 1, NG

        if (j == k) then
          control1(j,k) = 1.0_dp - Lmat(j,k)
          control2(j,k) = 1.0_dp - alpha * Lmat(j,k)
        else
          control1(j,k) = -Lmat(j,k)
          control2(j,k) = -alpha * Lmat(j,k)
        end if

      end do
    end do

    temp1 = matmul(B1, I1)
    temp1 = matmul(Cmat, temp1)
    lam1 = seasonal * matmul(control1, temp1)

    temp2 = matmul(B2, I2)
    temp2 = matmul(Cmat, temp2)
    lam2 = seasonal * matmul(control2, temp2)

    do j = 1, NG
      dy(j)    = 1.2_dp * S(j) * lam1(j) - GAM1(j) * I1(j)
      dy(NG+j) = 1.2_dp * S(j) * lam2(j) - GAM2(j) * I2(j)
    end do

  end subroutine rhs


  subroutine integrate_model(alpha, season_amp, t_eval, R_out, npts, ok)

    real(dp), intent(in) :: alpha, season_amp
    real(dp), intent(out) :: t_eval(:), R_out(:)
    integer, intent(out) :: npts
    logical, intent(out) :: ok

    real(dp) :: y(NY)
    real(dp) :: k1(NY), k2(NY), k3(NY), k4(NY), y_temp(NY)
    real(dp) :: t, dt, total_days, Itot, gamma_weighted

    integer :: i, n_keep

    real(dp), allocatable :: t_local(:), R_local(:)

    ok = .false.
    npts = 0

    y = 0.0_dp
    y(1) = 5.0e-4_dp
    y(NG+1) = 2.0e-4_dp

    total_days = SIM_END_DAYS
    dt = total_days / real(N_T_EVAL-1, dp)

    allocate(t_local(N_T_EVAL), R_local(N_T_EVAL))

    t = 0.0_dp

    do i = 1, N_T_EVAL

      t_local(i) = t

      Itot = max(sum(y(1:NY)), 1.0e-12_dp)

      gamma_weighted = &
        (sum(GAM1 * y(1:NG)) + sum(GAM2 * y(NG+1:NY))) / Itot

      gamma_weighted = max(gamma_weighted, 1.0e-6_dp)

      call rhs(t, y, alpha, season_amp, k1)

      R_local(i) = 1.0_dp + sum(k1) / (Itot * gamma_weighted)

      y_temp = max(y + 0.5_dp * dt * k1, 0.0_dp)
      call rhs(t + 0.5_dp * dt, y_temp, alpha, season_amp, k2)

      y_temp = max(y + 0.5_dp * dt * k2, 0.0_dp)
      call rhs(t + 0.5_dp * dt, y_temp, alpha, season_amp, k3)

      y_temp = max(y + dt * k3, 0.0_dp)
      call rhs(t + dt, y_temp, alpha, season_amp, k4)

      y = max(y + (dt / 6.0_dp) * &
          (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4), 0.0_dp)

      t = t + dt

    end do

    n_keep = 0

    do i = 1, N_T_EVAL

      if (t_local(i) >= BURN_IN_DAYS) then

        n_keep = n_keep + 1

        t_eval(n_keep) = &
          (t_local(i) - BURN_IN_DAYS) / DAYS_PER_WEEK

        R_out(n_keep) = R_local(i)

      end if

    end do

    npts = n_keep
    ok = .true.

    deallocate(t_local, R_local)

  end subroutine integrate_model


  subroutine evaluate_model_at_shift(alpha, season_amp, shift, &
                                     weeks_obs, R_pred, ok)

    real(dp), intent(in) :: alpha, season_amp, shift
    real(dp), intent(in) :: weeks_obs(:)
    real(dp), intent(out) :: R_pred(:)

    logical, intent(out) :: ok

    real(dp), allocatable :: weeks_model(:), R_model(:)

    integer :: i, j, n_valid

    real(dp) :: xq, w

    logical :: good

    ok = .false.
    R_pred = 0.0_dp

    allocate(weeks_model(N_T_EVAL), R_model(N_T_EVAL))

    call integrate_model(alpha, season_amp, weeks_model, &
                         R_model, n_valid, good)

    if (.not. good .or. n_valid < 10) then
      deallocate(weeks_model, R_model)
      return
    end if

    do i = 1, size(weeks_obs)

      xq = weeks_obs(i) - shift

      if (xq < weeks_model(1) .or. xq > weeks_model(n_valid)) then
        deallocate(weeks_model, R_model)
        return
      end if

      do j = 1, n_valid-1

        if (xq >= weeks_model(j) .and. &
            xq <= weeks_model(j+1)) then

          w = (xq - weeks_model(j)) / &
              (weeks_model(j+1) - weeks_model(j))

          R_pred(i) = R_model(j) + &
                      w * (R_model(j+1) - R_model(j))

          exit

        end if

      end do

    end do

    ok = .true.

    deallocate(weeks_model, R_model)

  end subroutine evaluate_model_at_shift

end module model_mod


module stats_mod
  use precision_mod
  use settings_mod
  implicit none

contains

  function pearson(x, y) result(r)

    real(dp), intent(in) :: x(:), y(:)

    real(dp) :: r, mx, my, sx, sy

    mx = sum(x) / real(size(x), dp)
    my = sum(y) / real(size(y), dp)

    sx = sqrt(sum((x - mx)**2))
    sy = sqrt(sum((y - my)**2))

    if (sx <= 0.0_dp .or. sy <= 0.0_dp) then
      r = 0.0_dp
    else
      r = sum((x - mx) * (y - my)) / (sx * sy)
    end if

  end function pearson


  subroutine percentile(values, p, q)

    real(dp), intent(in) :: values(:)
    real(dp), intent(in) :: p
    real(dp), intent(out) :: q

    real(dp), allocatable :: sorted_vals(:)

    integer :: n, k
    real(dp) :: a

    n = size(values)

    if (n == 0) then
      q = 0.0_dp
      return
    end if

    allocate(sorted_vals(n))

    sorted_vals = values

    call sort_array(sorted_vals)

    if (n == 1) then
      q = sorted_vals(1)
      deallocate(sorted_vals)
      return
    end if

    a = p * real(n-1, dp) + 1.0_dp

    k = int(floor(a))

    if (k < 1) k = 1

    if (k >= n) then
      q = sorted_vals(n)
    else
      q = sorted_vals(k) + &
          (a - real(k, dp)) * &
          (sorted_vals(k+1) - sorted_vals(k))
    end if

    deallocate(sorted_vals)

  end subroutine percentile


  subroutine sort_array(a)

    real(dp), intent(inout) :: a(:)

    call quicksort(a, 1, size(a))

  end subroutine sort_array


  recursive subroutine quicksort(a, l, r)

    real(dp), intent(inout) :: a(:)

    integer, intent(in) :: l, r

    integer :: i, j

    real(dp) :: p, tmp

    if (l >= r) return

    i = l
    j = r

    p = a((l+r)/2)

    do

      do while (a(i) < p)
        i = i + 1
      end do

      do while (a(j) > p)
        j = j - 1
      end do

      if (i <= j) then
        tmp = a(i)
        a(i) = a(j)
        a(j) = tmp

        i = i + 1
        j = j - 1
      end if

      if (i > j) exit

    end do

    if (l < j) call quicksort(a, l, j)
    if (i < r) call quicksort(a, i, r)

  end subroutine quicksort


  subroutine summary_stats(values, mean_val, std_val, &
                           q025, q50, q975)

    real(dp), intent(in) :: values(:)

    real(dp), intent(out) :: mean_val, std_val
    real(dp), intent(out) :: q025, q50, q975

    real(dp), allocatable :: vals_copy(:)

    integer :: n

    n = size(values)

    if (n == 0) then
      mean_val = 0.0_dp
      std_val  = 0.0_dp
      q025     = 0.0_dp
      q50      = 0.0_dp
      q975     = 0.0_dp
      return
    end if

    allocate(vals_copy(n))

    vals_copy = values

    mean_val = sum(vals_copy) / real(n, dp)

    if (n > 1) then
      std_val = sqrt(sum((vals_copy - mean_val)**2) / &
                      real(n-1, dp))
    else
      std_val = 0.0_dp
    end if

    call percentile(vals_copy, 0.025_dp, q025)
    call percentile(vals_copy, 0.50_dp,  q50)
    call percentile(vals_copy, 0.975_dp, q975)

    deallocate(vals_copy)

  end subroutine summary_stats


  !===========================================================================
  ! SPLIT-RHAT DIAGNOSTIC
  !
  ! flat_samples are stored as:
  !
  !   retained iteration 1: walkers 1,...,NWALKERS
  !   retained iteration 2: walkers 1,...,NWALKERS
  !   ...
  !
  ! Each walker therefore forms an individual chain.
  !
  ! The chains are split into two equal halves and the maximum split-Rhat
  ! across the four parameters is returned.
  !===========================================================================

  subroutine calculate_max_rhat(flat_samples, max_rhat)

    real(dp), intent(in) :: flat_samples(:,:)
    real(dp), intent(out) :: max_rhat

    integer :: n_total
    integer :: n_per_walker
    integer :: half_n
    integer :: m_split
    integer :: p
    integer :: walker
    integer :: h
    integer :: i
    integer :: start_idx
    integer :: end_idx
    integer :: chain_idx

    real(dp) :: W_within
    real(dp) :: B_between
    real(dp) :: var_hat
    real(dp) :: rhat
    real(dp) :: chain_mean
    real(dp) :: overall_mean
    real(dp) :: chain_var

    real(dp), allocatable :: chain_means(:)
    real(dp), allocatable :: chain_vars(:)

    max_rhat = 0.0_dp

    n_total = size(flat_samples, 1)

    if (size(flat_samples, 2) /= NDIM) then
      max_rhat = huge(1.0_dp)
      return
    end if

    if (mod(n_total, NWALKERS) /= 0) then
      max_rhat = huge(1.0_dp)
      return
    end if

    n_per_walker = n_total / NWALKERS

    if (n_per_walker < 4) then
      max_rhat = huge(1.0_dp)
      return
    end if

    half_n = n_per_walker / 2

    if (half_n < 2) then
      max_rhat = huge(1.0_dp)
      return
    end if

    m_split = 2 * NWALKERS

    allocate(chain_means(m_split))
    allocate(chain_vars(m_split))

    do p = 1, NDIM

      chain_idx = 0

      do walker = 1, NWALKERS

        do h = 1, 2

          chain_idx = chain_idx + 1

          if (h == 1) then

            start_idx = (1-1)*NWALKERS + walker
            end_idx = (half_n-1)*NWALKERS + walker

          else

            start_idx = half_n*NWALKERS + walker
            end_idx = (n_per_walker-1)*NWALKERS + walker

          end if

          chain_mean = 0.0_dp

          i = start_idx

          do while (i <= end_idx)

            chain_mean = chain_mean + flat_samples(i,p)

            i = i + NWALKERS

          end do

          chain_mean = chain_mean / real(half_n, dp)

          chain_means(chain_idx) = chain_mean

          chain_var = 0.0_dp

          i = start_idx

          do while (i <= end_idx)

            chain_var = chain_var + &
                        (flat_samples(i,p) - chain_mean)**2

            i = i + NWALKERS

          end do

          chain_var = chain_var / real(half_n-1, dp)

          chain_vars(chain_idx) = chain_var

        end do

      end do

      overall_mean = sum(chain_means) / real(m_split, dp)

      W_within = sum(chain_vars) / real(m_split, dp)

      B_between = real(half_n, dp) / real(m_split-1, dp) * &
                  sum((chain_means - overall_mean)**2)

      if (W_within > 0.0_dp) then

        var_hat = &
          ((real(half_n-1, dp) / real(half_n, dp)) * W_within) + &
          (B_between / real(half_n, dp))

        rhat = sqrt(max(var_hat / W_within, 0.0_dp))

      else

        if (B_between <= 0.0_dp) then
          rhat = 1.0_dp
        else
          rhat = huge(1.0_dp)
        end if

      end if

      max_rhat = max(max_rhat, rhat)

    end do

    deallocate(chain_means)
    deallocate(chain_vars)

  end subroutine calculate_max_rhat

end module stats_mod


module likelihood_mod
  use precision_mod
  use settings_mod
  use data_mod
  use model_mod
  use stats_mod
  implicit none

contains

  function log_likelihood_4d(theta) result(ll)

    real(dp), intent(in) :: theta(NDIM)

    real(dp) :: ll
    real(dp) :: R_pred(NDATA)
    real(dp) :: res
    real(dp) :: sigma

    real(dp) :: alpha
    real(dp) :: season_amp
    real(dp) :: shift

    logical :: ok

    integer :: i

    alpha = theta(1)
    season_amp = theta(2)
    sigma = theta(3)
    shift = theta(4)

    if (sigma <= 0.0_dp) then
      ll = -huge(1.0_dp)
      return
    end if

    call evaluate_model_at_shift(alpha, season_amp, shift, &
                                 WEEKS_DATA, R_pred, ok)

    if (.not. ok) then
      ll = -huge(1.0_dp)
      return
    end if

    ll = 0.0_dp

    do i = 1, NDATA

      res = R_DATA(i) - R_pred(i)

      ll = ll - 0.5_dp * &
           ((res / sigma)**2 + &
            log(2.0_dp * PI * sigma**2))

    end do

  end function log_likelihood_4d


  function log_prior_4d(theta) result(lp)

    real(dp), intent(in) :: theta(NDIM)

    real(dp) :: lp
    real(dp) :: shift_prior

    if (theta(1) > ALPHA_LO  .and. theta(1) < ALPHA_HI  .and. &
        theta(2) > SEASON_LO .and. theta(2) < SEASON_HI .and. &
        theta(3) > SIGMA_LO  .and. theta(3) < SIGMA_HI  .and. &
        theta(4) > SHIFT_LO  .and. theta(4) < SHIFT_HI) then

      lp = 0.0_dp

    else

      lp = -huge(1.0_dp)
      return

    end if

    shift_prior = -0.5_dp * &
                  (theta(4) / SHIFT_PRIOR_SIGMA)**2

    lp = lp + shift_prior

  end function log_prior_4d


  function log_posterior_4d(theta) result(lp)

    real(dp), intent(in) :: theta(NDIM)

    real(dp) :: lp
    real(dp) :: ll

    lp = log_prior_4d(theta)

    if (lp < -1.0e100_dp) return

    ll = log_likelihood_4d(theta)

    if (ll < -1.0e100_dp) then
      lp = -huge(1.0_dp)
    else
      lp = lp + ll
    end if

  end function log_posterior_4d


  function negative_loglike_4d(theta) result(f)

    real(dp), intent(in) :: theta(NDIM)

    real(dp) :: f
    real(dp) :: ll

    if (log_prior_4d(theta) < -1.0e100_dp) then
      f = 1.0e100_dp
      return
    end if

    ll = log_likelihood_4d(theta)

    if (ll < -1.0e100_dp) then
      f = 1.0e100_dp
    else
      f = -ll
    end if

  end function negative_loglike_4d

end module likelihood_mod


module optimization_mod
  use precision_mod
  use settings_mod
  use likelihood_mod
  use rng_mod
  implicit none

  abstract interface

    function cost_func(x) result(f)

      use precision_mod
      import :: NDIM

      real(dp), intent(in) :: x(NDIM)
      real(dp) :: f

    end function cost_func

  end interface

contains

  subroutine nelder_mead_4d(x0, func, xbest, fbest, maxiter, ftol)

    real(dp), intent(in) :: x0(NDIM)
    procedure(cost_func) :: func

    real(dp), intent(out) :: xbest(NDIM)
    real(dp), intent(out) :: fbest

    integer, intent(in) :: maxiter
    real(dp), optional, intent(in) :: ftol

    real(dp) :: x(NDIM, NDIM+1)
    real(dp) :: f(NDIM+1)

    real(dp) :: cent(NDIM)
    real(dp) :: xr(NDIM)
    real(dp) :: xe(NDIM)
    real(dp) :: xc(NDIM)

    real(dp) :: fr, fe, fc
    real(dp) :: tol, scale

    integer :: j, it
    integer :: ilo, ihi, inhi

    tol = 1.0e-6_dp

    if (present(ftol)) tol = ftol

    x(:,1) = x0
    x(:,2) = x0 + [0.02_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    x(:,3) = x0 + [0.0_dp, 0.02_dp, 0.0_dp, 0.0_dp]
    x(:,4) = x0 + [0.0_dp, 0.0_dp, 0.01_dp, 0.0_dp]
    x(:,5) = x0 + [0.0_dp, 0.0_dp, 0.0_dp, 0.05_dp]

    do j = 1, NDIM+1
      f(j) = func(x(:,j))
    end do

    do it = 1, maxiter

      ilo = 1
      ihi = 1

      do j = 2, NDIM+1

        if (f(j) < f(ilo)) ilo = j
        if (f(j) > f(ihi)) ihi = j

      end do

      inhi = 1

      do j = 1, NDIM+1

        if (j /= ihi) then

          if (inhi == ihi .or. f(j) > f(inhi)) then
            inhi = j
          end if

        end if

      end do

      scale = 0.0_dp

      do j = 1, NDIM+1

        scale = max(scale, &
                    abs(x(1,j)-x(1,ilo)) + &
                    abs(x(2,j)-x(2,ilo)) + &
                    abs(x(3,j)-x(3,ilo)) + &
                    abs(x(4,j)-x(4,ilo)))

      end do

      if (scale < 1.0e-5_dp .and. &
          abs(f(ihi)-f(ilo)) < tol) exit

      cent = 0.0_dp

      do j = 1, NDIM+1

        if (j /= ihi) then
          cent = cent + x(:,j)
        end if

      end do

      cent = cent / real(NDIM, dp)

      xr = cent + 1.0_dp * (cent - x(:,ihi))

      fr = func(xr)

      if (fr < f(ilo)) then

        xe = cent + 2.0_dp * (xr - cent)
        fe = func(xe)

        if (fe < fr) then
          x(:,ihi) = xe
          f(ihi) = fe
        else
          x(:,ihi) = xr
          f(ihi) = fr
        end if

      else if (fr < f(inhi)) then

        x(:,ihi) = xr
        f(ihi) = fr

      else

        if (fr < f(ihi)) then
          xc = cent + 0.5_dp * (xr - cent)
        else
          xc = cent + 0.5_dp * (x(:,ihi) - cent)
        end if

        fc = func(xc)

        if (fc < f(ihi)) then

          x(:,ihi) = xc
          f(ihi) = fc

        else

          do j = 1, NDIM+1

            if (j /= ilo) then

              x(:,j) = x(:,ilo) + &
                       0.5_dp * (x(:,j) - x(:,ilo))

              f(j) = func(x(:,j))

            end if

          end do

        end if

      end if

    end do

    ilo = minloc(f, 1)

    xbest = x(:,ilo)
    fbest = f(ilo)

  end subroutine nelder_mead_4d


  subroutine make_mcmc_start_4d(center, pos)

    real(dp), intent(in) :: center(NDIM)
    real(dp), intent(out) :: pos(NWALKERS, NDIM)

    integer :: i
    integer :: seed_local

    seed_local = 12345

    do i = 1, NWALKERS

      do

        pos(i,:) = center + &
          [0.01_dp, 0.01_dp, 0.003_dp, 0.05_dp] * &
          normal_rand_seeded(seed_local)

        if (log_prior_4d(pos(i,:)) > -1.0e100_dp) exit

      end do

    end do

  end subroutine make_mcmc_start_4d

end module optimization_mod


module mcmc_mod
  use precision_mod
  use settings_mod
  use likelihood_mod
  use rng_mod
  use optimization_mod
  implicit none

contains

  subroutine run_mcmc_4d(start, flat_samples, acceptance, &
                         walker_acceptance)

    real(dp), intent(in) :: start(NDIM)

    real(dp), allocatable, intent(out) :: flat_samples(:,:)

    real(dp), intent(out) :: acceptance

    real(dp), intent(out) :: walker_acceptance(NWALKERS)

    real(dp) :: chain(NWALKERS, NDIM)
    real(dp) :: lp(NWALKERS)
    real(dp) :: proposal(NDIM)

    real(dp) :: newlp
    real(dp) :: u, z, a

    integer :: i, j, k
    integer :: naccept, nprop
    integer :: seed_local
    integer :: n_keep, idx

    integer :: walker_prop(NWALKERS)
    integer :: walker_acc(NWALKERS)

    seed_local = 9173

    naccept = 0
    nprop = 0

    walker_prop = 0
    walker_acc = 0

    call make_mcmc_start_4d(start, chain)

    do i = 1, NWALKERS
      lp(i) = log_posterior_4d(chain(i,:))
    end do

    n_keep = NWALKERS * ((NSTEPS-BURN_MCMC)/THIN)

    allocate(flat_samples(n_keep, NDIM))

    idx = 0

    do k = 1, NSTEPS

      do i = 1, NWALKERS

        do

          j = 1 + int(real(NWALKERS, dp) * &
              uniform_rand_seeded(seed_local))

          if (j < 1) j = 1
          if (j > NWALKERS) j = NWALKERS

          if (j /= i) exit

        end do

        a = 2.0_dp

        z = ((a-1.0_dp) * &
             uniform_rand_seeded(seed_local) + 1.0_dp)**2 / a

        proposal = chain(j,:) + &
                   z * (chain(i,:) - chain(j,:))

        newlp = log_posterior_4d(proposal)

        nprop = nprop + 1
        walker_prop(i) = walker_prop(i) + 1

        if (newlp > -1.0e100_dp) then

          u = log(max(uniform_rand_seeded(seed_local), &
                      1.0e-15_dp))

          if (u < real(NDIM-1, dp) * log(z) + &
              newlp - lp(i)) then

            chain(i,:) = proposal
            lp(i) = newlp

            naccept = naccept + 1
            walker_acc(i) = walker_acc(i) + 1

          end if

        end if

      end do

      if (k > BURN_MCMC .and. &
          mod(k-BURN_MCMC, THIN) == 0) then

        idx = idx + 1

        do i = 1, NWALKERS
          flat_samples((idx-1)*NWALKERS+i,:) = chain(i,:)
        end do

      end if

      if (mod(k,500) == 0 .or. k == NSTEPS) then

        write(*,'(A,I6,A,I6,A,F5.1,A)') &
          'MCMC: ', k, '/', NSTEPS, ' (', &
          real(k,dp)/real(NSTEPS,dp)*100.0_dp, '%)'

      end if

    end do

    acceptance = real(naccept, dp) / &
                 real(max(1,nprop), dp)

    do i = 1, NWALKERS

      walker_acceptance(i) = &
        real(walker_acc(i),dp) / &
        real(max(1,walker_prop(i)),dp)

    end do

  end subroutine run_mcmc_4d

end module mcmc_mod


!===============================================================================
! MAIN PROGRAM - RSV LYON BAYESIAN CALIBRATION
!===============================================================================

program rsv_lyon_bayesian

  use precision_mod
  use settings_mod
  use data_mod
  use contact_mod
  use rng_mod
  use likelihood_mod
  use optimization_mod
  use mcmc_mod
  use stats_mod

  implicit none

  real(dp) :: initial_guess(NDIM)
  real(dp) :: best_theta(NDIM)
  real(dp) :: xbest(NDIM)

  real(dp) :: fbest
  real(dp) :: best_ll

  real(dp) :: best_R(NDATA)
  real(dp) :: mle_curve(NDATA)

  real(dp) :: mle_rmse
  real(dp) :: mle_corr

  logical :: ok

  real(dp), allocatable :: flat_samples(:,:)
  real(dp), allocatable :: posterior_models(:,:)

  real(dp), allocatable :: posterior_median(:)
  real(dp), allocatable :: posterior_lower(:)
  real(dp), allocatable :: posterior_upper(:)

  real(dp), allocatable :: posterior_corr_dist(:)
  real(dp), allocatable :: posterior_rmse_dist(:)

  real(dp) :: corr_025
  real(dp) :: corr_50
  real(dp) :: corr_975

  real(dp) :: rmse_025
  real(dp) :: rmse_50
  real(dp) :: rmse_975

  real(dp) :: acceptance
  real(dp) :: walker_acceptance(NWALKERS)

  real(dp) :: posterior_corr
  real(dp) :: posterior_rmse

  real(dp), allocatable :: temp_curve(:)

  integer :: n_valid_post
  integer :: i
  integer :: ios

  real(dp) :: mean_val
  real(dp) :: std_val
  real(dp) :: q025
  real(dp) :: q50
  real(dp) :: q975

  ! Added only for maximum split-Rhat calculation
  real(dp) :: max_rhat

  integer, allocatable :: draw_idx(:)
  integer :: post_seed

  character(len=16) :: labels(NDIM) = &
    ['alpha        ', 'season_amp   ', &
     'sigma        ', 'shift (weeks)']

  call initialize_contacts()
  call seed_rng(42)

  write(*,*) '==============================================================='
  write(*,*) ' RSV LYON - CORRECTED BAYESIAN CALIBRATION'
  write(*,*) ' Growth-Rate-Derived Reproduction Number'
  write(*,*) '==============================================================='

  write(*,'(A,I6)') ' Walkers   = ', NWALKERS
  write(*,'(A,I6)') ' Steps     = ', NSTEPS
  write(*,'(A,I6)') ' Burn-in   = ', BURN_MCMC
  write(*,'(A,I6)') ' Thin      = ', THIN
  write(*,'(A,I6)') ' MLE starts = ', N_MLE_STARTS
  write(*,'(A,I6)') ' Data pts  = ', NDATA

  write(*,*) ''

  !===========================================================================
  ! MAXIMUM LIKELIHOOD ESTIMATION
  !===========================================================================

  write(*,*) '--- MAXIMUM-LIKELIHOOD FIT ---'

  initial_guess = [0.8_dp, 0.35_dp, 0.15_dp, 1.8_dp]

  call nelder_mead_4d(initial_guess, negative_loglike_4d, &
                      xbest, fbest, 500, 1.0e-8_dp)

  best_theta = xbest
  best_ll = -fbest

  call evaluate_model_at_shift(best_theta(1), best_theta(2), &
                               best_theta(4), WEEKS_DATA, &
                               best_R, ok)

  if (.not. ok) then
    write(*,*) 'ERROR: MLE evaluation failed.'
    stop
  end if

  mle_curve = best_R

  mle_rmse = sqrt(sum((R_DATA-mle_curve)**2) / &
                  real(NDATA,dp))

  mle_corr = pearson(R_DATA,mle_curve)

  write(*,'(A,F12.6)') '  alpha         = ', best_theta(1)
  write(*,'(A,F12.6)') '  season_amp    = ', best_theta(2)
  write(*,'(A,F12.6)') '  sigma         = ', best_theta(3)
  write(*,'(A,F12.6)') '  shift (weeks) = ', best_theta(4)
  write(*,'(A,F12.6)') '  Log-Likelihood= ', best_ll
  write(*,'(A,F12.6)') '  MLE RMSE      = ', mle_rmse
  write(*,'(A,F12.6)') '  MLE Pearson r = ', mle_corr

  write(*,*) ''

  !===========================================================================
  ! BAYESIAN MCMC
  !===========================================================================

  write(*,*) '--- BAYESIAN 4D MCMC SAMPLING ---'

  call run_mcmc_4d(best_theta, flat_samples, &
                   acceptance, walker_acceptance)

  write(*,'(A,F8.4)') &
    'Mean acceptance fraction = ', acceptance

  ! Added diagnostic only
  call calculate_max_rhat(flat_samples, max_rhat)

  write(*,'(A,F10.6)') &
    'Maximum split-Rhat = ', max_rhat

  write(*,*) ''

  !===========================================================================
  ! POSTERIOR PARAMETER ESTIMATES
  !===========================================================================

  write(*,*) '--- BAYESIAN POSTERIOR PARAMETERS ---'

  do i = 1, NDIM

    call summary_stats(flat_samples(:,i), &
                       mean_val, std_val, &
                       q025, q50, q975)

    write(*,'(A14,A,F10.5,A,F10.5,A,F10.5,A)') &
      labels(i), ' median = ', q50, &
      ' | 95% CrI = [', q025, ', ', q975, ']'

  end do

  write(*,*) ''

  !===========================================================================
  ! POSTERIOR ENSEMBLE EVALUATION
  !===========================================================================

  write(*,*) '--- POSTERIOR ENSEMBLE EVALUATION ---'

  allocate(temp_curve(NDATA))
  allocate(draw_idx(N_POSTERIOR_MODELS))

  post_seed = 271828

  do i = 1, N_POSTERIOR_MODELS

    draw_idx(i) = 1 + int( &
      real(size(flat_samples,1),dp) * &
      uniform_rand_seeded(post_seed))

    if (draw_idx(i) < 1) &
      draw_idx(i) = 1

    if (draw_idx(i) > size(flat_samples,1)) &
      draw_idx(i) = size(flat_samples,1)

  end do

  n_valid_post = 0

  do i = 1, N_POSTERIOR_MODELS

    call evaluate_model_at_shift( &
      flat_samples(draw_idx(i),1), &
      flat_samples(draw_idx(i),2), &
      flat_samples(draw_idx(i),4), &
      WEEKS_DATA, temp_curve, ok)

    if (ok) n_valid_post = n_valid_post + 1

  end do

  write(*,'(A,I6)') &
    'Valid posterior models found: ', n_valid_post

  if (n_valid_post > 50) then

    allocate(posterior_models(NDATA,n_valid_post))
    allocate(posterior_corr_dist(n_valid_post))
    allocate(posterior_rmse_dist(n_valid_post))

    n_valid_post = 0

    do i = 1, N_POSTERIOR_MODELS

      call evaluate_model_at_shift( &
        flat_samples(draw_idx(i),1), &
        flat_samples(draw_idx(i),2), &
        flat_samples(draw_idx(i),4), &
        WEEKS_DATA, temp_curve, ok)

      if (ok) then

        n_valid_post = n_valid_post + 1

        posterior_models(:,n_valid_post) = temp_curve

        posterior_corr_dist(n_valid_post) = &
          pearson(R_DATA,temp_curve)

        posterior_rmse_dist(n_valid_post) = &
          sqrt(sum((R_DATA-temp_curve)**2) / &
          real(NDATA,dp))

      end if

      if (mod(i,50) == 0 .or. &
          i == N_POSTERIOR_MODELS) then

        write(*,'(A,I4,A,I4,A,I4)') &
          '  Processing: ', i, '/', &
          N_POSTERIOR_MODELS, &
          ' (valid: ', n_valid_post, ')'

      end if

    end do

    allocate(posterior_median(NDATA))
    allocate(posterior_lower(NDATA))
    allocate(posterior_upper(NDATA))

    do i = 1, NDATA

      call percentile(posterior_models(i,:), &
                      0.50_dp, posterior_median(i))

      call percentile(posterior_models(i,:), &
                      0.025_dp, posterior_lower(i))

      call percentile(posterior_models(i,:), &
                      0.975_dp, posterior_upper(i))

    end do

    call percentile(posterior_corr_dist, &
                    0.025_dp, corr_025)

    call percentile(posterior_corr_dist, &
                    0.50_dp, corr_50)

    call percentile(posterior_corr_dist, &
                    0.975_dp, corr_975)

    call percentile(posterior_rmse_dist, &
                    0.025_dp, rmse_025)

    call percentile(posterior_rmse_dist, &
                    0.50_dp, rmse_50)

    call percentile(posterior_rmse_dist, &
                    0.975_dp, rmse_975)

    posterior_corr = corr_50
    posterior_rmse = rmse_50

    write(*,*) ''

    write(*,'(A)') &
      'POSTERIOR CORRELATION DISTRIBUTION:'

    write(*,'(A,I6)') &
      '  Number of draws: ', n_valid_post

    write(*,'(A,F8.4)') &
      '  Median: ', corr_50

    write(*,'(A,F8.4,A,F8.4,A)') &
      '  95% CrI: [', corr_025, ', ', corr_975, ']'

    write(*,'(A,F8.4)') &
      '  Mean:   ', &
      sum(posterior_corr_dist) / &
      real(n_valid_post,dp)

    write(*,'(A,F8.4)') &
      '  Std Dev:', &
      sqrt(sum((posterior_corr_dist-corr_50)**2) / &
      real(n_valid_post-1,dp))

    write(*,*) ''

    write(*,'(A)') &
      'POSTERIOR RMSE DISTRIBUTION:'

    write(*,'(A,F8.4)') &
      '  Median: ', rmse_50

    write(*,'(A,F8.4,A,F8.4,A)') &
      '  95% CrI: [', rmse_025, ', ', rmse_975, ']'

    write(*,*) ''

    open(unit=30, &
         file='RSV_Posterior_Corr_Dist.dat', &
         status='replace')

    write(30,'(A)') '# corr rmse'

    do i = 1, n_valid_post

      write(30,'(2ES20.10)') &
        posterior_corr_dist(i), &
        posterior_rmse_dist(i)

    end do

    close(30)

    write(*,'(A)') &
      '  Saved posterior correlation distribution to RSV_Posterior_Corr_Dist.dat'

  else

    write(*,'(A)') &
      'WARNING: Few valid posterior models found!'

    allocate(posterior_median(NDATA))
    allocate(posterior_lower(NDATA))
    allocate(posterior_upper(NDATA))

    posterior_median = mle_curve
    posterior_lower = mle_curve * 0.8_dp
    posterior_upper = mle_curve * 1.2_dp

    posterior_corr = mle_corr
    posterior_rmse = mle_rmse

  end if

  deallocate(draw_idx)

  write(*,*) ''

  !===========================================================================
  ! SAVE OUTPUT FILES
  !===========================================================================

  write(*,*) '--- SAVING OUTPUT FILES ---'

  open(unit=11, &
       file='RSV_Posterior_Samples.dat', &
       status='replace', iostat=ios)

  if (ios == 0) then

    write(11,'(A)') &
      '# alpha season_amp sigma shift'

    do i = 1, size(flat_samples,1)

      write(11,'(4ES20.10)') &
        flat_samples(i,:)

    end do

    close(11)

    write(*,*) &
      '  Saved posterior samples to RSV_Posterior_Samples.dat'

  end if


  open(unit=13, &
       file='RSV_MLE_Summary.dat', &
       status='replace', iostat=ios)

  if (ios == 0) then

    write(13,'(A)') &
      '# alpha season_amp sigma shift logLik RMSE corr'

    write(13,'(8ES20.10)') &
      best_theta(1), best_theta(2), &
      best_theta(3), best_theta(4), &
      best_ll, mle_rmse, mle_corr

    close(13)

    write(*,*) &
      '  Saved MLE summary to RSV_MLE_Summary.dat'

  end if


  open(unit=20, &
       file='RSV_Fit_Results.csv', &
       status='replace', iostat=ios)

  if (ios == 0) then

    write(20,'(A)') &
      'week,r_data,mle_curve,post_median,post_lower,post_upper'

    do i = 1, NDATA

      write(20,'(F8.4,5(",",F8.4))') &
        WEEKS_DATA(i), &
        R_DATA(i), &
        mle_curve(i), &
        posterior_median(i), &
        posterior_lower(i), &
        posterior_upper(i)

    end do

    close(20)

    write(*,*) &
      '  Saved curve data to RSV_Fit_Results.csv'

  end if

  !===========================================================================
  ! FINAL SUMMARY
  !===========================================================================

  write(*,*) ''
  write(*,*) '==============================================================='
  write(*,*) ' FINAL RESULTS'
  write(*,*) '==============================================================='
  write(*,*) ''

  write(*,*) 'MLE:'

  write(*,'(A,F12.6)') &
    '  alpha         = ', best_theta(1)

  write(*,'(A,F12.6)') &
    '  season_amp    = ', best_theta(2)

  write(*,'(A,F12.6)') &
    '  sigma         = ', best_theta(3)

  write(*,'(A,F12.6)') &
    '  shift (weeks) = ', best_theta(4)

  write(*,'(A,F12.6)') &
    '  Pearson r     = ', mle_corr

  write(*,'(A,F12.6)') &
    '  RMSE          = ', mle_rmse

  write(*,*) ''

  write(*,*) 'POSTERIOR (95% Credible Interval):'

  do i = 1, NDIM

    call percentile(flat_samples(:,i), &
                    0.025_dp, q025)

    call percentile(flat_samples(:,i), &
                    0.975_dp, q975)

    write(*,'(A14,A,F10.5,A,F10.5,A)') &
      labels(i), ' = [', q025, ', ', q975, ']'

  end do

  write(*,'(A14,A,F10.5,A,F10.5,A)') &
    'Pearson r     ', &
    ' = [', corr_025, ', ', corr_975, ']'

  write(*,'(A14,A,F10.5)') &
    'RMSE         ', &
    ' = ', rmse_50

  write(*,*) ''

  write(*,'(A,F8.4)') &
    'MCMC acceptance fraction = ', acceptance

  !=======================================================================
  ! ONLY NEW RESULT ADDED
  !=======================================================================

  write(*,'(A,F10.6)') &
    'Maximum split-Rhat        = ', max_rhat

  write(*,*) ''

  write(*,*) 'Files written:'

  write(*,*) &
    '  1. RSV_Posterior_Samples.dat'

  write(*,*) &
    '  2. RSV_MLE_Summary.dat'

  write(*,*) &
    '  3. RSV_Fit_Results.csv'

  if (allocated(posterior_corr_dist)) then
    write(*,*) &
      '  4. RSV_Posterior_Corr_Dist.dat'
  end if

  write(*,*) ''

  write(*,*) '==============================================================='

  deallocate(flat_samples)

  deallocate(posterior_median)
  deallocate(posterior_lower)
  deallocate(posterior_upper)

  if (allocated(posterior_models)) &
    deallocate(posterior_models)

  if (allocated(posterior_corr_dist)) &
    deallocate(posterior_corr_dist)

  if (allocated(posterior_rmse_dist)) &
    deallocate(posterior_rmse_dist)

  deallocate(temp_curve)

end program rsv_lyon_bayesian
