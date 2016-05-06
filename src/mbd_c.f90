module mbd_c

use mbd, only: Damping_t, Context_t, Geometry_t, MBDResults_t, Relay_t, bohr, &
    init_grid, destroy_grid
use iso_c_binding

implicit none

real(c_double), bind(c, name='bohr') :: bohr_c = bohr

type, bind(c) :: Damping_ct
    type(c_ptr) :: label
    real(c_double) :: d, s_R, a, beta
    type(c_ptr) :: alpha_0 = c_null_ptr, C6 = c_null_ptr, &
        R_vdw = c_null_ptr, overlap = c_null_ptr, custom = c_null_ptr, &
        custom_potential = c_null_ptr
    integer(c_int) :: n_atoms
end type

type, bind(c) :: Context_ct
    logical(c_bool) :: do_rpa, do_reciprocal, do_ewald
    real(c_double) :: ts_energy_accuracy, ts_shell_thickness, &
        dipole_low_dim_cutoff, mayer_scaling, ewald_real_cutoff_scaling, &
        ewald_rec_cutoff_scaling
    type(c_ptr) :: freq_grid, freq_grid_w
    integer(c_int) :: n_freq
end type

type, bind(c) :: Geometry_ct
    type(c_ptr) :: coords
    logical(c_bool) :: is_periodic
    real(c_double) :: vacuum_axes(3)
    real(c_double) :: lattice_vectors(3, 3)
    type(c_ptr) :: k_points
    integer(c_int) :: supercell(3)
    integer(c_int) :: n_atoms, n_kpoints
end type

type, bind(c) :: MBDResults_ct
    real(c_double) :: energy
    logical(c_bool) :: has_error
    type(c_ptr) :: omegas, modes, bandomegas, rpa_orders, modes_im, bandmodes
    integer(c_int) :: n_atoms, n_orders, n_kpoints
end type

type, bind(c) :: Relay_ct
    type(c_ptr) :: re, im
    integer(c_int) :: n_atoms
end type

contains

subroutine init_grid_c(ctx_c, n) bind(c, name='_init_grid')
    type(Context_ct), intent(inout) :: ctx_c
    integer(c_int), intent(in), value :: n

    type(Context_t) :: ctx

    ctx = context_c_to_f(ctx_c)
    call init_grid(ctx, n)
    ctx_c = context_f_to_c(ctx)
end subroutine

subroutine destroy_grid_c(ctx_c) bind(c, name='_destroy_grid')
    type(Context_ct), intent(inout) :: ctx_c

    type(Context_t) :: ctx

    ctx = context_c_to_f(ctx_c)
    call destroy_grid(ctx)
    ctx_c = context_f_to_c(ctx)
end subroutine

subroutine run_scs_c( &
    ctx_c, geom_c, alpha_c, damping_c, alpha_full_c, lambda, k_point_c &
) bind(c, name='_run_scs')
    type(Context_ct), intent(in) :: ctx_c
    type(Geometry_ct), intent(in) :: geom_c
    type(c_ptr), intent(in) :: alpha_c
    type(Damping_ct), intent(in) :: damping_c
    type(Relay_ct), intent(in) :: alpha_full_c
    real(c_double), intent(in) :: lambda
    type(c_ptr), intent(in) :: k_point_c
end subroutine

real(c_double) function get_C6_from_alpha_c( &
    ctx_c, alpha_c &
) bind(c, name='_get_C6_from_alpha')
    type(Context_ct), intent(in) :: ctx_c
    type(c_ptr), intent(in) :: alpha_c
end function

subroutine omega_eff_c(C6_c, alpha_c, omega_c) bind(c, name='_omega_eff')
    type(c_ptr), intent(in) :: C6_c
    type(c_ptr), intent(in) :: alpha_c
    type(c_ptr), intent(in) :: omega_c
end subroutine

subroutine run_mbd_c( &
    ctx_c, geom_c, alpha_0_c, omega_c, damping_c, results_c &
) bind(c, name='_run_mbd')
    type(Context_ct), intent(in) :: ctx_c
    type(Geometry_ct), intent(in) :: geom_c
    type(c_ptr), intent(in) :: alpha_0_c
    type(c_ptr), intent(in) :: omega_c
    type(Damping_ct), intent(in) :: damping_c
    type(MBDResults_ct), intent(in) :: results_c
end subroutine

type(Context_t) function context_c_to_f(ctx_c) result(ctx)
    type(Context_ct), intent(in) :: ctx_c

    ctx%do_rpa = ctx_c%do_rpa
    ctx%do_reciprocal = ctx_c%do_reciprocal
    ctx%do_ewald = ctx_c%do_ewald
    ctx%ts_energy_accuracy = ctx_c%ts_energy_accuracy
    ctx%ts_shell_thickness = ctx_c%ts_shell_thickness
    ctx%dipole_low_dim_cutoff = ctx_c%dipole_low_dim_cutoff
    ctx%mayer_scaling = ctx_c%mayer_scaling
    ctx%ewald_real_cutoff_scaling = ctx_c%ewald_real_cutoff_scaling
    ctx%ewald_rec_cutoff_scaling = ctx_c%ewald_rec_cutoff_scaling
    ctx%n_freq = ctx_c%n_freq
    call c_f_pointer(ctx_c%freq_grid, ctx%freq_grid, (/ ctx%n_freq /))
    call c_f_pointer(ctx_c%freq_grid_w, ctx%freq_grid_w, (/ ctx%n_freq /))
end function

type(Context_ct) function context_f_to_c(ctx) result(ctx_c)
    type(Context_t), intent(in) :: ctx

    ctx_c%do_rpa = ctx%do_rpa
    ctx_c%do_reciprocal = ctx%do_reciprocal
    ctx_c%do_ewald = ctx%do_ewald
    ctx_c%ts_energy_accuracy = ctx%ts_energy_accuracy
    ctx_c%ts_shell_thickness = ctx%ts_shell_thickness
    ctx_c%dipole_low_dim_cutoff = ctx%dipole_low_dim_cutoff
    ctx_c%mayer_scaling = ctx%mayer_scaling
    ctx_c%ewald_real_cutoff_scaling = ctx%ewald_real_cutoff_scaling
    ctx_c%ewald_rec_cutoff_scaling = ctx%ewald_rec_cutoff_scaling
    ctx_c%n_freq = ctx%n_freq
    ctx_c%freq_grid = c_loc(ctx%freq_grid)
    ctx_c%freq_grid_w = c_loc(ctx%freq_grid_w)
end function

end module
