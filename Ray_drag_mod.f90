module Ray_drag_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

contains

  subroutine Ray_drag (return_eng, cp, tau_ray, u_dt, ua, v_dt, va, dTdt)
    implicit none

    logical, intent(in) :: return_eng
    real(dp), intent(in) :: cp, tau_ray, ua, va

    real(dp), intent(inout) :: u_dt, v_dt, dTdt

    real(dp) :: dTdt_ray

    ! Add drag to the flow for both zonal and meridional component
    u_dt = -ua/tau_ray
    v_dt = -va/tau_ray

    ! Add dissipated energy to dTdt if optioned
    if (return_eng .eqv. .True.) then
      dTdt_ray = (ua**2 + va**2)/(tau_ray*cp)
      dTdt = dTdt + dTdt_ray
    end if

  end subroutine Ray_drag

end module Ray_drag_mod

module Mag_Ray_drag_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-23_dp
  real(dp), parameter :: h = 6.62607015e-34_dp
  real(dp), parameter :: me = 9.1093837015e-31_dp
  real(dp), parameter :: ev = 1.602176634e-19_dp
  real(dp), parameter :: amu = 1.66053906660e-27_dp

  real(dp), parameter :: eps(28) = (/ &
  &   12.00_dp, 10.914_dp, 0.96_dp, 1.38_dp, 2.70_dp, 8.46_dp, 7.83_dp, 8.69_dp, &
  &   4.40_dp, 8.06_dp, 6.22_dp, 7.55_dp, 6.43_dp, 7.51_dp, 5.41_dp, 7.12_dp, &
  &   5.31_dp, 6.38_dp, 5.07_dp, 6.30_dp, 3.14_dp, 4.97_dp, 3.90_dp, 5.62_dp, &
  &   5.42_dp, 7.46_dp, 4.94_dp, 6.20_dp /)

  real(dp), parameter :: Ij(28) = (/ &
  &   13.59844_dp, 24.58738_dp, 5.39171_dp, 9.32269_dp, 8.29803_dp, 11.26030_dp, &
  &   14.53414_dp, 13.61806_dp, 17.42282_dp, 21.5646_dp, 5.13908_dp, 7.64624_dp, &
  &   5.98577_dp, 8.15169_dp, 10.48669_dp, 10.36001_dp, 12.96764_dp, 15.75962_dp, &
  &   4.34066_dp, 6.11316_dp, 6.5615_dp, 6.8281_dp, 6.7462_dp, 6.7665_dp, 7.43402_dp, &
  &   7.9024_dp, 7.8810_dp, 7.6398_dp /)


contains

  subroutine Mag_Ray_drag (return_eng, saha, xe_in, lat, B, Tl, pl, mu, cp, u_dt, ua, dTdt)
    implicit none

    logical, intent(in) :: return_eng, saha
    real(dp), intent(in) :: cp, ua, lat, B, Tl, pl, mu, xe_in

    real(dp), intent(inout) :: u_dt, dTdt

    real(dp) :: nd, rho, xe, tau_mag, eta, dTdt_mag, nj_tot
    real(dp), dimension(28) :: xj, nj

    rho = (pl * mu * amu)/(kb * Tl)
    nd = rho/(mu * amu)

    nj(:) = 10.0_dp**(eps(:) - eps(1)) * nd
    nj_tot = sum(nj(:))
    print*, nj(:)/nj_tot

    if (saha .eqv. .True.) then
      ! Saha equation for element j
      xj(:) = (2.0_dp*pi*me)/(h**2)**(1.5_dp) * (kb*Tl)**(1.5_dp) * exp((Ij(:)*ev)/(kb*Tl))/nj(:)
      print*, exp((Ij(:)*ev)/(kb*Tl))/nj(:)
      ! Do algebra
      xj(:) = sqrt(xj(:)/(1.0_dp + xj(:)))
      print*, xj(:)
      ! Total electron number density is sum of VMR * xj
      xe = sum(nj(:)/nj_tot*xj(:))
    else
      xe = xe_in
    end if

    eta = (230.0_dp * sqrt(Tl)) / xe

    tau_mag = (4.0_dp * pi * rho * eta) / (B**2 * abs(sin(lat)))

    ! Add drag to the flow for only the zonal component
    u_dt = -ua/tau_mag

    ! Add dissipated energy to dTdt if optioned
    if (return_eng .eqv. .True.) then
      dTdt_mag = ua**2/tau_mag/cp
      dTdt = dTdt + dTdt_mag
    end if

    print*, 'rho: ', rho
    print*, 'nd: ', nd
    print*, 'xe: ', xe
    print*,' tau mag: ', tau_mag

  end subroutine Mag_Ray_drag

end module Mag_Ray_drag_mod


program test_Ray_drag_mods
  use Ray_drag_mod, only : Ray_drag
  use Mag_Ray_drag_mod, only : Mag_Ray_drag
  implicit none

  logical :: return_eng, saha
  double precision :: cp, tau_ray, u_dt, ua, v_dt, va, dTdt
  double precision :: lat, mu, Tl, pl, xe_in, B

  return_eng = .True.
  cp = 1.41d4
  tau_ray = 86400.0d0 !* 2.0d0
  u_dt = 0.0d0
  v_dt = 0.0d0
  ua = 1000.0d0
  va = 1000.0d0
  dTdt = 0.0d0

  call Ray_drag (return_eng, cp, tau_ray, u_dt, ua, v_dt, va, dTdt) 

  print*, 'tau_ray, cp, ua, va: ', tau_ray, cp, ua, va
  print*, 'u_dt: ', u_dt
  print*, 'v_dt: ', v_dt
  print*, 'dTdt: ', dTdt

  saha = .False.
  lat = 15.0d0 * (3.14d0/180.0d0)
  B = 10.0d0
  Tl = 1000.0d0
  pl = 1d5
  mu = 2.3d0/1000.0d0
  xe_in = 1.0d-6
  u_dt = 0.0d0
  v_dt = 0.0d0
  ua = 1000.0d0
  va = 1000.0d0
  dTdt = 0.0d0

  call Mag_Ray_Drag (return_eng, saha, xe_in, lat, B, Tl, pl, mu, cp, u_dt, ua, dTdt)

  print*, 'lat, B, Tl, pl, mu, cp, ua, va: ', abs(sin(lat)), B, Tl, pl, mu, cp, ua, va
  print*, 'u_dt: ', u_dt
  print*, 'v_dt: ', v_dt
  print*, 'dTdt: ', dTdt

end program test_Ray_drag_mods
