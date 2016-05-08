module mbd_c

use mbd, only: Damping_t, Context_t, Geometry_t, MBDResults_t, Relay_t, bohr, &
    init_grid, destroy_grid, set_damping_parameters
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

subroutine c_init_ctx(c_ctx, n) bind(c, name='_init_ctx')
    type(Context_ct), intent(inout) :: c_ctx
    integer(c_int), intent(in), value :: n

    type(Context_t) :: f_ctx

    call init_grid(f_ctx, n)
    c_ctx = f_c_context(f_ctx)
end subroutine

subroutine c_destroy_grid(c_ctx) bind(c, name='_destroy_grid')
    type(Context_ct), intent(inout) :: c_ctx

    type(Context_t) :: f_ctx

    call c_f_pointer(c_ctx%freq_grid, f_ctx%freq_grid, (/ c_ctx%n_freq /))
    call c_f_pointer(c_ctx%freq_grid_w, f_ctx%freq_grid_w, (/ c_ctx%n_freq /))
    call destroy_grid(f_ctx)
end subroutine

subroutine c_run_scs( &
    c_ctx, geom_c, alpha_c, c_damping, alpha_full_c, lambda, k_point_c &
) bind(c, name='_run_scs')
    type(Context_ct), intent(in) :: c_ctx
    type(Geometry_ct), intent(in) :: geom_c
    type(c_ptr), intent(in) :: alpha_c
    type(Damping_ct), intent(in) :: c_damping
    type(Relay_ct), intent(in) :: alpha_full_c
    real(c_double), intent(in) :: lambda
    type(c_ptr), intent(in) :: k_point_c
end subroutine

real(c_double) function c_get_C6_from_alpha( &
    c_ctx, alpha_c &
) bind(c, name='_get_C6_from_alpha')
    type(Context_ct), intent(in) :: c_ctx
    type(c_ptr), intent(in) :: alpha_c
end function

subroutine c_omega_eff(C6_c, alpha_c, omega_c) bind(c, name='_omega_eff')
    type(c_ptr), intent(in) :: C6_c
    type(c_ptr), intent(in) :: alpha_c
    type(c_ptr), intent(in) :: omega_c
end subroutine

subroutine c_run_mbd( &
    c_ctx, geom_c, alpha_0_c, omega_c, c_damping, results_c &
) bind(c, name='_run_mbd')
    type(Context_ct), intent(in) :: c_ctx
    type(Geometry_ct), intent(in) :: geom_c
    type(c_ptr), intent(in) :: alpha_0_c
    type(c_ptr), intent(in) :: omega_c
    type(Damping_ct), intent(in) :: c_damping
    type(MBDResults_ct), intent(in) :: results_c
end subroutine

subroutine c_set_damping_parameters( &
    c_damping, c_method, c_xc &
) bind(c, name='_set_damping_parameters')
    type(Damping_ct), intent(inout) :: c_damping
    type(c_ptr), intent(in), value :: c_method, c_xc

    type(Damping_t) :: f_damping

    call set_damping_parameters(f_damping, c_f_string(c_method), c_f_string(c_xc))
    c_damping%d = f_damping%d
    c_damping%s_R = f_damping%s_R
    c_damping%a = f_damping%a
    c_damping%beta = f_damping%beta
end subroutine

type(Context_ct) function f_c_context(f_ctx) result(c_ctx)
    type(Context_t), intent(in) :: f_ctx

    c_ctx%do_rpa = f_ctx%do_rpa
    c_ctx%do_reciprocal = f_ctx%do_reciprocal
    c_ctx%do_ewald = f_ctx%do_ewald
    c_ctx%ts_energy_accuracy = f_ctx%ts_energy_accuracy
    c_ctx%ts_shell_thickness = f_ctx%ts_shell_thickness
    c_ctx%dipole_low_dim_cutoff = f_ctx%dipole_low_dim_cutoff
    c_ctx%mayer_scaling = f_ctx%mayer_scaling
    c_ctx%ewald_real_cutoff_scaling = f_ctx%ewald_real_cutoff_scaling
    c_ctx%ewald_rec_cutoff_scaling = f_ctx%ewald_rec_cutoff_scaling
    c_ctx%n_freq = f_ctx%n_freq
    c_ctx%freq_grid = c_loc(f_ctx%freq_grid)
    c_ctx%freq_grid_w = c_loc(f_ctx%freq_grid_w)
end function

type(Context_t) function c_f_context(c_ctx) result(f_ctx)
    type(Context_ct), intent(in) :: c_ctx

    f_ctx%do_rpa = c_ctx%do_rpa
    f_ctx%do_reciprocal = c_ctx%do_reciprocal
    f_ctx%do_ewald = c_ctx%do_ewald
    f_ctx%ts_energy_accuracy = c_ctx%ts_energy_accuracy
    f_ctx%ts_shell_thickness = c_ctx%ts_shell_thickness
    f_ctx%dipole_low_dim_cutoff = c_ctx%dipole_low_dim_cutoff
    f_ctx%mayer_scaling = c_ctx%mayer_scaling
    f_ctx%ewald_real_cutoff_scaling = c_ctx%ewald_real_cutoff_scaling
    f_ctx%ewald_rec_cutoff_scaling = c_ctx%ewald_rec_cutoff_scaling
    f_ctx%n_freq = c_ctx%n_freq
    call c_f_pointer(c_ctx%freq_grid, f_ctx%freq_grid, (/ c_ctx%n_freq /))
    call c_f_pointer(c_ctx%freq_grid_w, f_ctx%freq_grid_w, (/ c_ctx%n_freq /))
end function

type(Damping_t) function c_f_damping(c_damping) result(f_damping)
    type(Damping_ct), intent(in) :: c_damping

    integer :: n

    n = c_damping%n_atoms
    f_damping%d = c_damping%d
    f_damping%s_R = c_damping%s_R
    f_damping%a = c_damping%a
    f_damping%beta = c_damping%beta
    call c_f_pointer(c_damping%alpha_0, f_damping%alpha_0, (/ n /))
    call c_f_pointer(c_damping%C6, f_damping%C6, (/ n /))
    call c_f_pointer(c_damping%R_vdw, f_damping%R_vdw, (/ n /))
    call c_f_pointer(c_damping%overlap, f_damping%overlap, (/ n, n /))
    call c_f_pointer(c_damping%custom, f_damping%custom, (/ n, n /))
    call c_f_pointer( &
        c_damping%custom_potential, f_damping%custom_potential, (/ n, n, 3, 3 /) &
    )
end function

function c_f_string(c_str) result(f_str)
    type(c_ptr), intent(in) :: c_str
    character(len=:, kind=c_char), pointer :: f_str

    interface
        integer(c_size_t) function strlen(s) bind(c, name='strlen')
            import c_ptr, c_size_t
            type(c_ptr), intent(in), value :: s
        end function strlen
    end interface
    character(kind=c_char), pointer :: char_arr(:)

    call c_f_pointer(c_str, char_arr, [strlen(c_str)])
    call get_scalar(char_arr, size(char_arr), f_str)

    contains

    subroutine get_scalar(str, length, f_str)
        integer, intent(in) :: length
        character(kind=c_char, len=length), intent(in), target :: str(1)
        character(:, kind=c_char), intent(out), pointer :: f_str

        f_str => str(1)
    end subroutine
end function

end module
