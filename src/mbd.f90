module mbd

use mbd_interface, only: sync_sum, print_error, print_warning, &
    print_log, mute, unmute, my_task, n_tasks
use mbd_util, only: tostr, start_clock, stop_clock, measure_time, get_plock, &
    release_plock, have_plock, operator(.cprod.), diag, eye, invert, inverted, &
    diagonalize, sdiagonalize, diagonalized, sdiagonalized

implicit none

private

public :: init_grid, destroy_grid

real(8), parameter, public :: &
    bohr = 0.529177249d0, &
    pi = 3.1415926535897932385d0

type, public :: Damping_t
    character(len=30) :: label
    real(8) :: d, s_R, a, beta
    real(8), pointer :: &
        alpha_0(:) => null(), &
        C6(:) => null(), &
        R_vdw(:) => null(), &
        overlap(:, :) => null(), &
        custom(:, :) => null(), &
        custom_potential(:, :, :, :) => null()
end type

type, public :: Context_t
    logical :: &
        do_rpa = .false., &
        do_reciprocal = .true., &
        do_ewald = .true.
    real(8) :: &
        ts_energy_accuracy = 1.d-10, &
        ts_shell_thickness = 10.d0, &
        dipole_low_dim_cutoff = 100.d0/bohr, &
        mayer_scaling = 1.d0, &
        ewald_real_cutoff_scaling = 1.d0, &
        ewald_rec_cutoff_scaling = 1.d0
    integer :: n_freq = 0
    real(8), pointer :: &
        freq_grid(:) => null(), &
        freq_grid_w(:) => null()
end type

type, public :: Geometry_t
    real(8), pointer :: coords(:, :) => null()
    logical :: is_periodic = .false.
    logical :: vacuum_axes(3) = (/ .false., .false., .false. /)
    real(8) :: lattice_vectors(3, 3)
    real(8), pointer :: k_points(:, :) => null()
    integer :: supercell(3)
end type

type, public :: MBDResults_t
    real(8) :: energy
    logical :: has_error = .false.
    real(8), pointer :: &
        omegas(:) => null(), &
        modes(:, :) => null(), &
        bandomegas(:, :) => null(), &
        rpa_orders(:) => null()
    complex(8), pointer :: &
        modes_im(:, :) => null(), &  ! for modes of a single band
        bandmodes(:, :, :) => null()
end type

type, public :: Relay_t
    real(8), pointer :: re(:, :)
    complex(8), pointer :: cplx(:, :)
end type

contains

real(8) function get_ts_energy(ctx, geom, C6, alpha_0, damping) result(ene)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: C6(:), alpha_0(:)
    type(Damping_t), intent(in) :: damping

    real(8) :: C6_ij, r(3), r_norm, f_damp, R_vdw_ij, overlap_ij, &
        ene_shell, ene_pair, R_cell(3), ts_radius, ts_radius_inner, s_R, d
    integer :: i_shell, i_cell, i_atom, j_atom, range_cell(3), idx_cell(3)

    ene = 0.d0
    i_shell = 0
    s_R = damping%s_R
    d = damping%s_R
    call get_plock('ts')
    do
        i_shell = i_shell + 1
        ene_shell = 0.d0
        ts_radius = i_shell*ctx%ts_shell_thickness
        ts_radius_inner = ts_radius-ctx%ts_shell_thickness
        if (geom%is_periodic) then
            range_cell = supercell_circum(geom%lattice_vectors, ts_radius)
            where (geom%vacuum_axes) range_cell = 0
        else
            range_cell = (/ 0, 0, 0 /)
        end if
        idx_cell = (/ 0, 0, -1 /)
        do i_cell = 1, product(1+2*range_cell)
            call shift_cell(idx_cell, -range_cell, range_cell)
            if (have_plock('ts') .and. geom%is_periodic) then
                if (my_task /= modulo(i_cell, n_tasks)) cycle
            end if
            if (geom%is_periodic) then
                R_cell = matmul(idx_cell, geom%lattice_vectors)
            else
                R_cell = (/ 0.d0, 0.d0, 0.d0 /)
            end if
            do i_atom = 1, size(geom%coords, 1)
                if (have_plock('ts') .and. .not. geom%is_periodic) then
                    if (my_task /= modulo(i_atom, n_tasks)) cycle
                end if
                do j_atom = 1, i_atom
                    if (i_cell == 1) then
                        if (i_atom == j_atom) cycle
                    end if
                    r = geom%coords(i_atom, :)-geom%coords(j_atom, :)-R_cell
                    r_norm = sqrt(dot_product(r, r))
                    if (r_norm < ts_radius_inner) cycle
                    if (r_norm >= ts_radius) cycle
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha_0(i_atom), alpha_0(j_atom) &
                    )
                    if (associated(damping%R_vdw)) then
                        R_vdw_ij = damping%R_vdw(i_atom)+damping%R_vdw(j_atom)
                    end if
                    if (associated(damping%overlap)) then
                        overlap_ij = damping%overlap(i_atom, j_atom)
                    end if
                    select case (damping%label)
                    case ('fermi')
                        f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)
                    case ('fermi2')
                        f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)**2
                    case ('erf')
                        f_damp = damping_erf(r_norm, s_R*R_vdw_ij, d)
                    case ('1mexp')
                        f_damp = damping_1mexp(r_norm, s_R*R_vdw_ij, d)
                    case ('overlap')
                        f_damp = damping_overlap( &
                            r_norm, overlap_ij, C6_ij, s_R, d &
                        )
                    case ('custom')
                        f_damp = damping%custom(i_atom, j_atom)
                    end select
                    ene_pair = -C6_ij*f_damp/r_norm**6
                    if (i_atom == j_atom) then
                        ene_shell = ene_shell + ene_pair/2
                    else
                        ene_shell = ene_shell + ene_pair
                    endif
                end do ! j_atom
            end do ! i_atom
        end do ! i_cell
        if (have_plock('ts')) then
            call sync_sum(ene_shell)
        end if
        ene = ene + ene_shell
        if (.not. geom%is_periodic) exit
        if (i_shell > 1 .and. abs(ene_shell) < ctx%ts_energy_accuracy) then
            call print_log( &
                'Periodic TS converged in ' &
                // trim(tostr(i_shell)) // ' shells, ' &
                // trim(tostr(ts_radius*bohr)) // ' angstroms' &
            )
            exit
        endif
        call release_plock('ts')
    end do ! i_shell
end function get_ts_energy

subroutine init_grid(ctx, n)
    type(Context_t), intent(inout) :: ctx
    integer, intent(in) :: n

    if (associated(ctx%freq_grid)) deallocate (ctx%freq_grid)
    if (associated(ctx%freq_grid_w)) deallocate (ctx%freq_grid_w)
    ctx%n_freq = n
    allocate (ctx%freq_grid(0:n))
    allocate (ctx%freq_grid_w(0:n))
    ctx%freq_grid(0) = 0.d0
    ctx%freq_grid_w(0) = 0.d0
    call get_frequency_grid(n, 0.6d0, ctx%freq_grid(1:n), ctx%freq_grid_w(1:n))
    call print_log( &
        'Initialized a frequency grid of ' // trim(tostr(n)) // ' points.' &
    )
    call print_log( &
        'Relative quadrature error in C6 of carbon atom: '// &
        trim(tostr(test_frequency_grid(ctx))) &
    )
end subroutine

real(8) function test_frequency_grid(ctx) result(error)
    type(Context_t), intent(in) :: ctx

    real(8) :: alpha(0:ctx%n_freq)

    alpha = alpha_dynamic_ts(ctx%freq_grid, 21.d0, C6=99.5d0)
    error = abs(get_C6_from_alpha(ctx, alpha)/99.5d0-1.d0)
end function

subroutine destroy_grid(ctx)
    type(Context_t), intent(inout) :: ctx

    deallocate (ctx%freq_grid)
    deallocate (ctx%freq_grid_w)
end subroutine

subroutine add_dipole_matrix(ctx, geom, damping, relay, k_point)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    type(Damping_t), intent(in) :: damping
    type(Relay_t), intent(inout) :: relay
    real(8), intent(in), optional :: k_point(3)

    real(8) :: Tpp(3, 3), R_cell(3), r(3), r_norm, R_vdw_ij, C6_ij, &
        overlap_ij, sigma_ij, volume, ewald_alpha, real_space_cutoff, &
        beta, a
    complex(8) :: Tpp_c(3, 3)
    character(len=10) :: parallel_mode
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j
    logical :: is_lrange, want_ewald

    call get_plock('add_dipole')
    call start_clock('add_dipole_matrix')
    a = damping%a
    beta = damping%beta
    want_ewald = &
        geom%is_periodic .and. .not. any(geom%vacuum_axes) .and. ctx%do_ewald
    if (have_plock('add_dipole')) then
        if (geom%is_periodic .and. size(geom%coords, 1) < n_tasks) then
            parallel_mode = 'cells'
        else
            parallel_mode = 'atoms'
        end if
    else
        parallel_mode = ''
    end if
    if (have_plock('add_dipole')) then
        ! will be restored by syncing at the end
        if (present(k_point)) then
            relay%cplx = relay%cplx/n_tasks
        else
            relay%re = relay%re/n_tasks
        end if
    end if
    if (geom%is_periodic) then
        if (any(geom%vacuum_axes)) then
            real_space_cutoff = ctx%dipole_low_dim_cutoff
        else
            volume = max(dble(product(diagonalized(geom%lattice_vectors))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1.d0/3)
            real_space_cutoff = 6.d0/ewald_alpha*ctx%ewald_real_cutoff_scaling
            call print_log( &
                'Ewald: using alpha = ' // trim(tostr(ewald_alpha)) &
                // ', real cutoff = ' // trim(tostr(real_space_cutoff)) &
            )
        end if
        range_cell = supercell_circum(geom%lattice_vectors, real_space_cutoff)
        where (geom%vacuum_axes) range_cell = 0
    else
        range_cell(:) = 0
    end if
    if (geom%is_periodic) then
        call print_log( &
            'Ewald: summing real part in cell vector range of ' &
            // trim(tostr(1+2*range_cell(1))) // 'x' &
            // trim(tostr(1+2*range_cell(2))) // 'x' &
            // trim(tostr(1+2*range_cell(3))) &
        )
    end if
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        if (parallel_mode == 'cells') then
            if (my_task /= modulo(i_cell, n_tasks)) cycle
        end if
        if (geom%is_periodic) then
            R_cell = matmul(idx_cell, geom%lattice_vectors)
        else
            R_cell(:) = 0.d0
        end if
        do i_atom = 1, size(geom%coords, 1)
            if (parallel_mode == 'atoms') then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            do j_atom = 1, i_atom
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = geom%coords(i_atom, :)-geom%coords(j_atom, :)-R_cell
                r_norm = sqrt(dot_product(r, r))
                if (geom%is_periodic) then
                    if (r_norm > real_space_cutoff) cycle
                end if
                if (associated(damping%R_vdw)) then
                    R_vdw_ij = damping%R_vdw(i_atom)+damping%R_vdw(j_atom)
                end if
                if (associated(damping%alpha_0)) then
                    sigma_ij = sqrt(sum(get_sigma_selfint( &
                        damping%alpha_0((/ i_atom , j_atom /)))**2) &
                    )
                end if
                if (associated(damping%overlap)) then
                    overlap_ij = damping%overlap(i_atom, j_atom)
                end if
                if (associated(damping%C6)) then
                    C6_ij = combine_C6( &
                        damping%C6(i_atom), damping%C6(j_atom), &
                        damping%alpha_0(i_atom), damping%alpha_0(j_atom) &
                    )
                end if
                is_lrange = .true.
                select case (damping%label)
                case ('bare')
                    Tpp = T_bare(r)
                case ('dip,1mexp')
                    Tpp = T_1mexp_coulomb(r, beta*R_vdw_ij, a)
                case ('dip,erf')
                    Tpp = T_erf_coulomb(r, beta*R_vdw_ij, a)
                case ('dip,fermi')
                    Tpp = T_fermi_coulomb(r, beta*R_vdw_ij, a)
                case ('dip,overlap')
                    Tpp = T_overlap_coulomb(r, overlap_ij, C6_ij, beta, a)
                case ('1mexp,dip')
                    Tpp = damping_1mexp(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                case ('erf,dip')
                    Tpp = damping_erf(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                case ('fermi,dip')
                    Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                case ('fermi^2,dip')
                    Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)**2*T_bare(r)
                case ('overlap,dip')
                    Tpp = damping_overlap( &
                        r_norm, overlap_ij, C6_ij, beta, a &
                    )*T_bare(r)
                case ('custom,dip')
                    Tpp = damping%custom(i_atom, j_atom)*T_bare(r)
                case ('dip,custom')
                    Tpp = damping%custom_potential(i_atom, j_atom, :, :)
                case ('dip,gg')
                    Tpp = T_erf_coulomb(r, sigma_ij, 1.d0)
                case ('1mexp,dip,gg')
                    Tpp = (1.d0-damping_1mexp(r_norm, beta*R_vdw_ij, a)) &
                        * T_erf_coulomb(r, sigma_ij, 1.d0)
                    is_lrange = .false.
                case ('erf,dip,gg')
                    Tpp = (1.d0-damping_erf(r_norm, beta*R_vdw_ij, a)) &
                        * T_erf_coulomb(r, sigma_ij, 1.d0)
                    is_lrange = .false.
                case ('fermi,dip,gg')
                    Tpp = (1.d0-damping_fermi(r_norm, beta*R_vdw_ij, a)) &
                        * T_erf_coulomb(r, sigma_ij, 1.d0)
                    is_lrange = .false.
                case ('custom,dip,gg')
                    Tpp = (1.d0-damping%custom(i_atom, j_atom)) &
                        * T_erf_coulomb(r, sigma_ij, 1.d0)
                    is_lrange = .false.
                end select
                if (want_ewald .and. is_lrange) then
                    Tpp = Tpp + T_erfc(r, ewald_alpha)-T_bare(r)
                end if
                if (present(k_point)) then
                    Tpp_c = Tpp*exp(-cmplx(0.d0, 1.d0, 8)*( &
                        dot_product(k_point, r)))
                end if
                i = 3*(i_atom-1)
                j = 3*(j_atom-1)
                if (present(k_point)) then
                    relay%cplx(i+1:i+3, j+1:j+3) = &
                        relay%cplx(i+1:i+3, j+1:j+3) + Tpp_c
                    if (i_atom /= j_atom) then
                        relay%cplx(j+1:j+3, i+1:i+3) = &
                            relay%cplx(j+1:j+3, i+1:i+3) + transpose(Tpp_c)
                    end if
                else
                    relay%re(i+1:i+3, j+1:j+3) = &
                        relay%re(i+1:i+3, j+1:j+3) + Tpp
                    if (i_atom /= j_atom) then
                        relay%re(j+1:j+3, i+1:i+3) = &
                            relay%re(j+1:j+3, i+1:i+3) + transpose(Tpp)
                    end if
                end if
            end do ! j_atom
        end do ! i_atom
    end do ! i_cell
    if (have_plock('add_dipole')) then
        if (present(k_point)) then
            call sync_sum(relay%cplx)
        else
            call sync_sum(relay%re)
        end if
    end if
    call stop_clock('add_dipole_matrix')
    call release_plock('add_dipole')
    if (want_ewald .and. is_lrange) then
        call add_ewald_dipole_parts(ctx, geom, ewald_alpha, relay, k_point)
    end if
end subroutine

subroutine add_ewald_dipole_parts(ctx, geom, ewald_alpha, relay, k_point)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: ewald_alpha
    type(Relay_t), intent(inout) :: relay
    real(8), intent(in), optional :: k_point(3)

    logical :: do_surface
    real(8) :: rec_lattice_vectors(3, 3), volume, G_vector(3), r(3), k_total(3), &
        k_sq, rec_space_cutoff, Tpp(3, 3), k_prefactor(3, 3), elem
    complex(8) :: Tpp_c(3, 3)
    integer :: &
        i_atom, j_atom, i, j, i_xyz, j_xyz, idx_G_vector(3), i_G_vector, &
        range_G_vector(3)
    character(len=10) :: parallel_mode

    call get_plock('add_dipole_ewald')
    call start_clock('add_ewald_dipole_parts')
    if (have_plock('add_dipole_ewald')) then
        if (geom%is_periodic .and. size(geom%coords, 1) < n_tasks) then
            parallel_mode = 'gvectors'
        else
            parallel_mode = 'atoms'
        end if
    else
        parallel_mode = ''
    end if
    if (have_plock('add_dipole_ewald')) then
        ! will be restored by syncing at the end
        if (present(k_point)) then
            relay%cplx = relay%cplx/n_tasks
        else
            relay%re = relay%re/n_tasks
        end if
    end if
    rec_lattice_vectors = 2*pi*inverted(transpose(geom%lattice_vectors))
    volume = dble(product(diagonalized(geom%lattice_vectors)))
    rec_space_cutoff = 10.d0*ewald_alpha*ctx%ewald_real_cutoff_scaling
    range_G_vector = supercell_circum(rec_lattice_vectors, rec_space_cutoff)
    where (geom%vacuum_axes) range_G_vector = 0
    call print_log( &
        'Ewald: using reciprocal cutoff = ' &
        // trim(tostr(rec_space_cutoff)) &
    )
    call print_log( &
        'Ewald: summing reciprocal part in G vector range of ' &
        // trim(tostr(1+2*range_G_vector(1))) // 'x' &
        // trim(tostr(1+2*range_G_vector(2))) // 'x' &
        // trim(tostr(1+2*range_G_vector(3))) &
    )
    idx_G_vector = (/ 0, 0, -1 /)
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        if (parallel_mode == 'gvectors') then
            if (my_task /= modulo(i_G_vector, n_tasks)) cycle
        end if
        G_vector = matmul(idx_G_vector, rec_lattice_vectors)
        if (present(k_point)) then
            k_total = k_point+G_vector
        else
            k_total = G_vector
        end if
        k_sq = sum(k_total**2)
        if (sqrt(k_sq) > rec_space_cutoff) cycle
        k_prefactor(:, :) = 4*pi/volume*exp(-k_sq/(4*ewald_alpha**2))
        forall (i_xyz = 1:3, j_xyz = 1:3)
            k_prefactor(i_xyz, j_xyz) = &
                k_prefactor(i_xyz, j_xyz) * k_total(i_xyz)*k_total(j_xyz)/k_sq
        end forall
        do i_atom = 1, size(geom%coords, 1)
            if (parallel_mode == 'atoms') then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            do j_atom = 1, i_atom
                r = geom%coords(i_atom, :)-geom%coords(j_atom, :)
                if (present(k_point)) then
                    Tpp_c = k_prefactor*exp(cmplx(0.d0, 1.d0, 8) &
                        * dot_product(G_vector, r))
                else
                    Tpp = k_prefactor*cos(dot_product(G_vector, r))
                end if
                i = 3*(i_atom-1)
                j = 3*(j_atom-1)
                if (present(k_point)) then
                    relay%cplx(i+1:i+3, j+1:j+3) = &
                        relay%cplx(i+1:i+3, j+1:j+3) + Tpp_c
                    if (i_atom /= j_atom) then
                        relay%cplx(j+1:j+3, i+1:i+3) = &
                            relay%cplx(j+1:j+3, i+1:i+3) + transpose(Tpp_c)
                    end if
                else
                    relay%re(i+1:i+3, j+1:j+3) = &
                        relay%re(i+1:i+3, j+1:j+3) + Tpp
                    if (i_atom /= j_atom) then
                        relay%re(j+1:j+3, i+1:i+3) = &
                            relay%re(j+1:j+3, i+1:i+3) + transpose(Tpp)
                    end if
                end if
            end do ! j_atom
        end do ! i_atom
    end do ! i_G_vector
    if (have_plock('add_dipole_ewald')) then
        if (present(k_point)) then
            call sync_sum(relay%cplx)
        else
            call sync_sum(relay%re)
        end if
    end if
    do i_atom = 1, size(geom%coords, 1) ! self energy
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            if (present(k_point)) then
                relay%cplx(i, i) = relay%cplx(i, i) - 4*ewald_alpha**3/(3*sqrt(pi))
            else
                relay%re(i, i) = relay%re(i, i) - 4*ewald_alpha**3/(3*sqrt(pi))
            end if
        end do
    end do
    do_surface = .true.
    if (present(k_point)) then
        k_sq = sum(k_point**2)
        if (sqrt(k_sq) > 1.d-15) then
            do_surface = .false.
            do i_atom = 1, size(geom%coords, 1)
                do j_atom = 1, i_atom
                    do i_xyz = 1, 3
                        do j_xyz = 1, 3
                            i = 3*(i_atom-1)+i_xyz
                            j = 3*(j_atom-1)+j_xyz
                            elem = &
                                4*pi/volume*k_point(i_xyz)*k_point(j_xyz)/k_sq &
                                * exp(-k_sq/(4*ewald_alpha**2))
                            if (present(k_point)) then
                                relay%cplx(i, j) = relay%cplx(i, j) + elem
                                if (i_atom /= j_atom) then
                                    relay%cplx(j, i) = relay%cplx(j, i) + elem
                                end if
                            else
                                relay%re(i, j) = relay%re(i, j) + elem
                                if (i_atom /= j_atom) then
                                    relay%re(j, i) = relay%re(j, i) + elem
                                end if
                            end if ! present(k_point
                        end do ! j_xyz
                    end do ! i_xyz
                end do ! j_atom
            end do ! i_atom
        end if ! k_sq >
    end if ! k_point present
    if (do_surface) then ! surface energy
        do i_atom = 1, size(geom%coords, 1)
            do j_atom = 1, i_atom
                do i_xyz = 1, 3
                    i = 3*(i_atom-1)+i_xyz
                    j = 3*(j_atom-1)+i_xyz
                    if (present(k_point)) then
                        relay%cplx(i, j) = relay%cplx(i, j) + 4*pi/(3*volume)
                        relay%cplx(j, i) = relay%cplx(i, j)
                    else
                        relay%re(i, j) = relay%re(i, j) + 4*pi/(3*volume)
                        relay%re(j, i) = relay%re(i, j)
                    end if
                end do ! i_xyz
            end do ! j_atom
        end do ! i_atom
    end if
    call stop_clock('add_ewald_dipole_parts')
end subroutine

subroutine run_scs(ctx, geom, alpha, damping, alpha_full, lambda, k_point)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha(:)
    type(Damping_t), intent(in) :: damping
    type(Relay_t), intent(out) :: alpha_full
    real(8), intent(in), optional :: lambda
    real(8), intent(in), optional :: k_point(3)

    integer :: i_atom, i_xyz, i

    if (present(k_point)) then
        alpha_full%cplx(:, :) = cmplx(0.d0, 0.d0, 8)
    else
        alpha_full%re(:, :) = 0.d0
    end if
    call add_dipole_matrix(ctx, geom, damping, alpha_full, k_point)
    if (present(lambda)) then
        if (present(k_point)) then
            alpha_full%cplx = lambda*alpha_full%cplx
        else
            alpha_full%re = lambda*alpha_full%re
        end if
    end if
    do i_atom = 1, size(geom%coords, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            if (present(k_point)) then
                alpha_full%cplx(i, i) = alpha_full%cplx(i, i) + 1.d0/alpha(i_atom)
            else
                alpha_full%re(i, i) = alpha_full%re(i, i) + 1.d0/alpha(i_atom)
            end if
        end do
    end do
    call start_clock('invert')
    if (present(k_point)) then
        call invert(alpha_full%cplx)
    else
        call invert(alpha_full%re)
    end if
    call stop_clock('invert')
end subroutine

function get_screened_alpha(ctx, geom, alpha, damping) result(alpha_screened)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha(0:, :)
    type(Damping_t), intent(in) :: damping
    real(8) :: alpha_screened(0:ctx%n_freq, size(geom%coords, 1))

    type(Relay_t) :: relay
    real(8), target :: alpha_full(3*size(geom%coords, 1), 3*size(geom%coords, 1))
    integer :: i_freq

    call get_plock('screened_alpha')
    alpha_screened(:, :) = 0.d0
    relay%re => alpha_full
    do i_freq = 0, ctx%n_freq
        if (have_plock('screened_alpha')) then
            if (my_task /= modulo(i_freq, n_tasks)) cycle
        end if
        call run_scs(ctx, geom, alpha(i_freq+1, :), damping, relay)
        alpha_screened(i_freq, :) = contract_polarizability(alpha_full)
        if (i_freq == 0) call mute()
    end do
    call unmute()
    if (have_plock('screened_alpha')) then
        call sync_sum(alpha_screened)
    end if
    call release_plock('screened_alpha')
end function

subroutine run_mbd(ctx, geom, alpha_0, omega, damping, results)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha_0(:), omega(:)
    type(Damping_t), intent(in) :: damping
    type(MBDResults_t), intent(inout) :: results

    type(Relay_t) :: alpha
    integer :: i_atom

    if (ctx%do_rpa) then
        allocate (alpha%re(0:ctx%n_freq, size(alpha_0)))
        forall (i_atom = 1:size(alpha_0))
            alpha%re(:, i_atom) = alpha_dynamic_ts( &
                ctx%freq_grid, alpha_0(i_atom), omega(i_atom) &
            )
        end forall
    end if
    if (.not. geom%is_periodic) then
        if (ctx%do_rpa) then
            call run_single_rpa(ctx, geom, alpha%re, damping, results)
        else
            call run_single_mbd(ctx, geom, alpha_0, omega, damping, results)
        end if
    else if (.not. ctx%do_rpa) then
        if (ctx%do_reciprocal) then
            call run_reciprocal_mbd(ctx, geom, alpha_0, omega, damping, results)
        else
            call run_supercell_mbd(ctx, geom, alpha_0, omega, damping, results)
        end if
    end if
    if (ctx%do_rpa) then
        deallocate (alpha%re)
    end if
end subroutine

subroutine run_single_mbd(ctx, geom, alpha_0, omega, damping, results)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha_0(:), omega(:)
    type(Damping_t), intent(in) :: damping
    type(MBDResults_t), intent(inout) :: results

    integer :: i_atom, j_atom, i_xyz, i, j
    integer :: n_negative_eigs
    type(Relay_t) :: relay
    real(8), pointer :: eigs(:)
    character(len=1) :: mode

    if (associated(results%omegas)) then
        eigs => results%omegas
    else
        allocate (eigs(3*size(geom%coords, 1)))
    end if
    if (associated(results%modes)) then
        relay%re => results%modes
    else
        allocate (relay%re(3*size(geom%coords, 1), 3*size(geom%coords, 1)))
    end if
    call add_dipole_matrix(ctx, geom, damping, relay)
    do i_atom = 1, size(geom%coords, 1)
        do j_atom = 1, size(geom%coords, 1)
            i = 3*(i_atom-1)
            j = 3*(j_atom-1)
            relay%re(i+1:i+3, j+1:j+3) = &  ! relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                * sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                relay%re(i+1:i+3, j+1:j+3)
        end do
    end do
    do i_atom = 1, size(geom%coords, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            relay%re(i, i) = &  ! relay = w^2+sqrt(a*a)*w*w*T
                relay%re(i, i) + omega(i_atom)**2
        end do
    end do
    if (associated(relay%re, results%modes)) then
        mode = 'V'
    else
        mode = 'N'
    end if
    call sdiagonalize(mode, relay%re, eigs)
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call print_warning( &
            'CDM Hamiltonian has ' // trim(tostr(n_negative_eigs)) &
            // ' negative eigenvalues' &
        )
        results%has_error = .true.
        return
    end if
    eigs = sqrt(eigs)
    results%energy = 1.d0/2*sum(eigs)-3.d0/2*sum(omega)
    if (.not. associated(eigs, results%omegas)) then
        deallocate (eigs)
    end if
    if (.not. associated(relay%re, results%modes)) then
        deallocate (relay%re)
    end if
end subroutine

subroutine run_single_recip_mbd(ctx, geom, alpha_0, omega, k_point, damping, results)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha_0(:), omega(:)
    real(8), intent(in) :: k_point(3)
    type(Damping_t), intent(in) :: damping
    type(MBDResults_t), intent(inout) :: results

    integer :: i_atom, j_atom, i_xyz, i, j
    integer :: n_negative_eigs
    type(Relay_t) :: relay
    real(8), pointer :: eigs(:)
    character(len=1) :: mode

    if (associated(results%omegas)) then
        eigs => results%omegas
    else
        allocate (eigs(3*size(geom%coords, 1)))
    end if
    if (associated(results%modes_im)) then
        relay%cplx => results%modes_im
    else
        allocate (relay%cplx(3*size(geom%coords, 1), 3*size(geom%coords, 1)))
    end if
    relay%cplx(:, :) = cmplx(0.d0, 0.d0, 8)
    call add_dipole_matrix(ctx, geom, damping, relay, k_point)
    do i_atom = 1, size(geom%coords, 1)
        do j_atom = 1, size(geom%coords, 1)
            i = 3*(i_atom-1)
            j = 3*(j_atom-1)
            relay%cplx(i+1:i+3, j+1:j+3) = & ! relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                * sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                relay%cplx(i+1:i+3, j+1:j+3)
        end do
    end do
    do i_atom = 1, size(geom%coords, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            relay%cplx(i, i) = &  ! relay = w^2+sqrt(a*a)*w*w*T
                relay%cplx(i, i) + omega(i_atom)**2
        end do
    end do
    if (associated(relay%cplx, results%modes_im)) then
        mode = 'V'
    else
        mode = 'N'
    end if
    call sdiagonalize(mode, relay%cplx, eigs)
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call print_warning( &
            'CDM Hamiltonian has ' // trim(tostr(n_negative_eigs)) &
            // ' negative eigenvalues' &
        )
        results%has_error = .true.
    end if
    eigs = sqrt(eigs)
    results%energy = 1.d0/2*sum(eigs)-3.d0/2*sum(omega)
    if (.not. associated(eigs, results%omegas)) then
        deallocate (eigs)
    end if
    if (.not. associated(relay%cplx, results%modes_im)) then
        deallocate (relay%cplx)
    end if
end subroutine

subroutine run_reciprocal_mbd(ctx, geom, alpha_0, omega, damping, results)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha_0(:), omega(:)
    type(Damping_t), intent(in) :: damping
    type(MBDResults_t), intent(inout) :: results

    integer :: i_kpt
    real(8) :: energy

    call get_plock('mbd_recip')
    energy = 0.d0
    if (associated(results%bandomegas)) then
        allocate (results%omegas(3*size(geom%coords, 1)))
    end if
    if (associated(results%bandmodes)) then
        allocate (results%modes_im(3*size(geom%coords, 1), 3*size(geom%coords, 1)))
    end if
    do i_kpt = 1, size(geom%k_points, 1)
        if (have_plock('mbd_recip')) then
            if (my_task /= modulo(i_kpt, n_tasks)) cycle
        end if
        call run_single_recip_mbd( &
            ctx, geom, alpha_0, omega, geom%k_points(i_kpt, :), damping, results &
        )
        energy = energy + results%energy
        if (associated(results%omegas)) then
            results%bandomegas(i_kpt, :) = results%omegas(:)
        end if
        if (associated(results%modes_im)) then
            results%bandmodes(i_kpt, :, :) = results%modes_im(:, :)
        end if
        if (i_kpt == 1) call mute()
    end do ! k_point loop
    call unmute()
    if (have_plock('mbd_recip')) then
        call sync_sum(energy)
        if (associated(results%bandomegas)) then
            call sync_sum(results%bandomegas)
        end if
        if (associated(results%bandmodes)) then
            call sync_sum(results%bandmodes)
        end if
    end if
    results%energy = energy/size(geom%k_points, 1)
    call release_plock('mbd_recip')
end subroutine

subroutine run_supercell_mbd(ctx, geom, alpha_0, omega, damping, results)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in) :: alpha_0(:), omega(:)
    type(Damping_t), intent(in) :: damping
    type(MBDResults_t), intent(inout) :: results

    type(Geometry_t) :: geom_super
    type(Damping_t) :: damping_super
    real(8) :: &
        alpha_0_super(product(geom%supercell)), &
        omega_super(product(geom%supercell))
    real(8) :: R_cell(3)
    integer :: i_atom, i
    integer :: i_cell
    integer :: idx_cell(3), n_cells

    geom_super = geom
    damping_super = damping
    n_cells = product(geom%supercell)
    allocate (geom_super%coords(n_cells*size(geom%coords, 1), 3))
    if (associated(damping%alpha_0)) then
        allocate (damping_super%alpha_0(n_cells*size(damping%alpha_0)))
    end if
    if (associated(damping%R_vdw)) then
        allocate (damping_super%R_vdw(n_cells*size(damping%R_vdw)))
    end if
    if (associated(damping%C6)) then
        allocate (damping_super%C6(n_cells*size(damping%C6)))
    end if
    forall (i = 1:3)
        geom_super%lattice_vectors(i, :) = geom%lattice_vectors(i, :)*geom%supercell(i)
    end forall
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, n_cells
        call shift_cell(idx_cell, (/ 0, 0, 0 /), geom%supercell-1)
        R_cell = matmul(idx_cell, geom%lattice_vectors)
        do i_atom = 1, size(geom%coords, 1)
            i = (i_cell-1)*size(geom%coords, 1)+i_atom
            geom_super%coords(i, :) = geom%coords(i_atom, :)+R_cell
            alpha_0_super(i) = alpha_0(i_atom)
            omega_super(i) = omega(i_atom)
            if (associated(damping%alpha_0)) then
                damping_super%alpha_0(i) = damping%alpha_0(i_atom)
            end if
            if (associated(damping%R_vdw)) then
                damping_super%R_vdw(i) = damping%R_vdw(i_atom)
            end if
            if (associated(damping%C6)) then
                damping_super%C6(i) = damping%C6(i_atom)
            end if
        end do
    end do
    call run_single_mbd( &
        ctx, geom_super, alpha_0_super, omega_super, damping_super, results &
    )
    results%energy = results%energy/n_cells
    deallocate (geom_super%coords)
    if (associated(damping_super%alpha_0)) then
        deallocate (damping_super%alpha_0)
    end if
    if (associated(damping_super%R_vdw)) then
        deallocate (damping_super%R_vdw)
    end if
    if (associated(damping_super%C6)) then
        deallocate (damping_super%C6)
    end if
end subroutine

subroutine run_single_rpa(ctx, geom, alpha, damping, results)
    type(Context_t), intent(in) :: ctx
    type(Geometry_t), intent(in) :: geom
    real(8), intent(in), target :: alpha(0:, :)
    type(Damping_t), intent(in) :: damping
    type(MBDResults_t), intent(inout) :: results

    type(Relay_t) :: relay
    type(Damping_t) :: damping_freq
    real(8) :: AT(3*size(geom%coords, 1), 3*size(geom%coords, 1))
    complex(8) :: eigs(3*size(geom%coords, 1))
    integer :: i_atom, i_xyz, i_freq, i
    integer :: n_order, n_negative_eigs
    real(8) :: energy

    if (ctx%n_freq > n_tasks) then
        call get_plock('rpa')
    end if
    allocate (relay%re(3*size(geom%coords, 1), 3*size(geom%coords, 1)))
    damping_freq = damping
    energy = 0.d0
    if (associated(results%rpa_orders)) then
        results%rpa_orders(:) = 0.d0
    end if
    do i_freq = 0, ctx%n_freq
        if (have_plock('rpa')) then
            if (my_task /= modulo(i_freq, n_tasks)) cycle
        end if
        damping_freq%alpha_0 => alpha(i_freq, :)
        call add_dipole_matrix(ctx, geom, damping_freq, relay)
        do i_atom = 1, size(geom%coords, 1)
            do i_xyz = 1, 3
                i = (i_atom-1)*3+i_xyz
                relay%re(i, :) = &  ! relay = alpha*T
                    alpha(i_freq, i_atom)*relay%re(i, :)
            end do
        end do
        AT = relay%re
        do i = 1, 3*size(geom%coords, 1)
            relay%re(i, i) = 1.d0+relay%re(i, i) ! relay = 1+alpha*T
        end do
        call diagonalize('N', relay%re, eigs)
        n_negative_eigs = count(dble(eigs(:)) < 0)
        if (n_negative_eigs > 0) then
            call print_warning( &
                '1+AT matrix has ' // trim(tostr(n_negative_eigs)) &
                // ' negative eigenvalues' &
            )
        end if
        energy = &
            energy + 1.d0/(2*pi)*sum(log(dble(eigs)))*ctx%freq_grid_w(i_freq)
        if (associated(results%rpa_orders)) then
            call diagonalize('N', AT, eigs)
            do n_order = 2, size(results%rpa_orders)
                results%rpa_orders(n_order) = &
                    results%rpa_orders(n_order) &
                    + (-1.d0/(2*pi)*(-1)**n_order &
                    * sum(dble(eigs)**n_order)/n_order) &
                    * ctx%freq_grid_w(i_freq)
            end do
        end if
        if (i_freq == 0) call mute()
    end do
    call unmute()
    if (have_plock('rpa')) then
        call sync_sum(energy)
        if (associated(results%rpa_orders)) then
            call sync_sum(results%rpa_orders)
        end if
    end if
    results%energy = energy
    call release_plock('rpa')
end subroutine

function eval_mbd_nonint_density(pts, coords, charges, masses, omegas) result(rho)
    real(8), intent(in) :: &
        pts(:, :), &
        coords(:, :), &
        charges(size(coords, 1)), &
        masses(size(coords, 1)), &
        omegas(size(coords, 1))
    real(8) :: rho(size(pts, 1))

    integer :: i_pt, i_atom, n_atoms
    real(8), dimension(size(coords, 1)) :: pre, kernel, rsq

    pre = charges*(masses*omegas/pi)**(3.d0/2)
    kernel = masses*omegas
    n_atoms = size(coords, 1)
    rho(:) = 0.d0
    do i_pt = 1, size(pts, 1)
        if (my_task /= modulo(i_pt, n_tasks)) cycle
        forall (i_atom = 1:n_atoms)
            rsq(i_atom) = sum((pts(i_pt, :)-coords(i_atom, :))**2)
        end forall
        rho(i_pt) = sum(pre*exp(-kernel*rsq))
    end do
    call sync_sum(rho)
end function

function eval_mbd_int_density(pts, coords, charges, masses, omegas, modes) result(rho)
    real(8), intent(in) :: &
        pts(:, :), &
        coords(:, :), &
        charges(size(coords, 1)), &
        masses(size(coords, 1)), &
        omegas(3*size(coords, 1)), &
        modes(3*size(coords, 1), 3*size(coords, 1))
    real(8) :: rho(size(pts, 1))

    integer :: i_pt, i_atom, n_atoms, i, i_xyz, j_xyz
    integer :: self(3), other(3*(size(coords, 1)-1))
    real(8) :: &
        pre(size(coords, 1)), &
        factor(size(coords, 1)), &
        rdiffsq(3, 3), &
        omegas_p(3*size(coords, 1), 3*size(coords, 1)), &
        kernel(3, 3, size(coords, 1)), &
        rdiff(3)

    omegas_p = matmul(matmul(modes, diag(omegas)), transpose(modes))
    n_atoms = size(coords, 1)
    kernel(:, :, :) = 0.d0
    pre(:) = 0.d0
    do i_atom = 1, n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        self(:) = (/ (3*(i_atom-1)+i, i = 1, 3) /)
        other(:) = (/ (i, i = 1, 3*(i_atom-1)),  (i, i = 3*i_atom+1, 3*n_atoms) /)
        kernel(:, :, i_atom) = &
            masses(i_atom)*(omegas_p(self, self) &
                - matmul(matmul(omegas_p(self, other), inverted(omegas_p(other, other))), &
                    omegas_p(other, self)))
        pre(i_atom) = charges(i_atom)*(masses(i_atom)/pi)**(3.d0/2) &
            *sqrt(product(omegas)/product(sdiagonalized(omegas_p(other, other))))
    end do
    call sync_sum(kernel)
    call sync_sum(pre)
    rho(:) = 0.d0
    do i_pt = 1, size(pts, 1)
        if (my_task /= modulo(i_pt, n_tasks)) cycle
        do i_atom = 1, n_atoms
            rdiff(:) = pts(i_pt, :)-coords(i_atom, :)
            forall (i_xyz = 1:3, j_xyz = 1:3)
                rdiffsq(i_xyz, j_xyz) = rdiff(i_xyz)*rdiff(j_xyz)
            end forall
            factor(i_atom) = sum(kernel(:, :, i_atom)*rdiffsq(:, :))
        end do
        rho(i_pt) = sum(pre*exp(-factor))
    end do
    call sync_sum(rho)
end function

! function get_single_reciprocal_rpa_ene(mode, version, xyz, alpha, k_point, &
!         geom%lattice_vectors, R_vdw, beta, a, overlap, C6, custom, &
!         custom_potential, rpa_orders) result(ene)
!     character(len=*), intent(in) :: mode, version
!     real(8), intent(in) :: &
!         xyz(:, :), &
!         alpha(:, :), &
!         k_point(3), &
!         geom%lattice_vectors(3, 3)
!     real(8), intent(in), optional :: &
!         R_vdw(size(xyz, 1)), &
!         beta, a, &
!         overlap(size(xyz, 1), size(xyz, 1)), &
!         C6(size(xyz, 1)), &
!         custom(size(xyz, 1), size(xyz, 1)), &
!         custom_potential(size(xyz, 1), size(xyz, 1), 3, 3)
!     real(8), intent(out), optional :: rpa_orders(20)
!     real(8) :: ene
!
!     complex(8), dimension(3*size(xyz, 1), 3*size(xyz, 1)) :: relay, AT
!     complex(8) :: eigs(3*size(xyz, 1))
!     integer :: i_atom, i_xyz, i_freq, i
!     integer :: n_order, n_negative_eigs
!     logical :: ctx%is_parallel, get_orders
!     character(len=1) :: muted
!
!     ctx%is_parallel = is_in('P', mode)
!     get_orders = is_in('O', mode)
!     if (is_in('M', mode)) then
!         muted = 'M'
!     else
!         muted = ''
!     end if
!
!     do i_freq = 0, n_grid_omega
!         ! MPI code begin
!         if (ctx%is_parallel) then
!             if (my_task /= modulo(i_freq, n_tasks)) cycle:
!         end if
!         ! MPI code end
!         relay(:, :) = 0.d0
!         call add_dipole_matrix( & ! relay = T
!             blanked('P', mode)//muted, &
!             version, &
!             xyz, &
!             alpha=alpha(i_freq+1, :), &
!             R_vdw=R_vdw, &
!             beta=beta, &
!             a=a, &
!             overlap=overlap, &
!             C6=C6, &
!             custom=custom, &
!             custom_potential=custom_potential, &
!             k_point=k_point, &
!             geom%lattice_vectors=geom%lattice_vectors, &
!             relay%cplx=relay)
!         do i_atom = 1, size(xyz, 1)
!             do i_xyz = 1, 3
!                 i = (i_atom-1)*3+i_xyz
!                 relay(i, :) = alpha(i_freq+1, i_atom)*relay(i, :)
!                 ! relay = alpha*T
!             end do
!         end do
!         AT = relay
!         do i = 1, 3*size(xyz, 1)
!             relay(i, i) = 1.d0+relay(i, i) ! relay = 1+alpha*T
!         end do
!         relay = sqrt(relay*conjg(relay))
!         call diagonalize('N', relay, eigs)
!         n_negative_eigs = count(dble(eigs) < 0)
!         if (n_negative_eigs > 0) then
!             call print_warning('1+AT matrix has ' &
!                 //trim(tostr(n_negative_eigs))//' negative eigenvalues')
!         end if
!         ene = ene + 1.d0/(2*pi)*sum(log(dble(eigs)))*ctx%freq_grid_w(i_freq)
!         if (get_orders) then
!             AT = 2*dble(AT)+AT*conjg(AT)
!             call diagonalize('N', AT, eigs)
!             do n_order = 2, param_rpa_order_max
!                 rpa_orders(n_order) = &
!                     rpa_orders(n_order) +(-1.d0/(2*pi)*(-1)**n_order &
!                     *sum(dble(eigs)**n_order)/n_order) &
!                     *ctx%freq_grid_w(i_freq)
!             end do
!         end if
!         muted = 'M'
!     end do
!     if (ctx%is_parallel) then
!         call sync_sum(ene)
!         if (get_orders) then
!             call sync_sum(rpa_orders)
!         end if
!     end if
! end function get_single_reciprocal_rpa_ene

function make_g_grid(n1, n2, n3) result(g_grid)
    integer, intent(in) :: n1, n2, n3
    real(8) :: g_grid(n1*n2*n3, 3)

    integer :: g_kpt(3), i_kpt, kpt_range(3), g_kpt_shifted(3)

    g_kpt = (/ 0, 0, -1 /)
    kpt_range = (/ n1, n2, n3 /)
    do i_kpt = 1, n1*n2*n3
        call shift_cell (g_kpt, (/ 0, 0, 0 /), kpt_range-1)
        g_kpt_shifted = g_kpt
        where (2*g_kpt > kpt_range) g_kpt_shifted = g_kpt-kpt_range
        g_grid(i_kpt, :) = dble(g_kpt_shifted)/kpt_range
    end do
end function make_g_grid

function make_k_grid(g_grid, uc) result(k_grid)
    real(8), intent(in) :: g_grid(:, :), uc(3, 3)
    real(8) :: k_grid(size(g_grid, 1), 3)

    integer :: i_kpt
    real(8) :: ruc(3, 3)

    ruc = 2*pi*inverted(transpose(uc))
    do i_kpt = 1, size(g_grid, 1)
        k_grid(i_kpt, :) = matmul(g_grid(i_kpt, :), ruc)
    end do
end function make_k_grid

function supercell_circum(uc, radius) result(sc)
    real(8), intent(in) :: uc(3, 3), radius
    integer :: sc(3)

    real(8) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = 2*pi*inverted(transpose(uc))
    forall (i = 1:3) layer_sep(i) = sum(uc(i, :)*ruc(i, :)/sqrt(sum(ruc(i, :)**2)))
    sc = ceiling(radius/layer_sep+0.5d0)
end function

subroutine shift_cell(ijk, first_cell, last_cell)
    integer, intent(inout) :: ijk(3)
    integer, intent(in) :: first_cell(3), last_cell(3)

    integer :: i_dim, i

    do i_dim = 3, 1, -1
        i = ijk(i_dim)+1
        if (i <= last_cell(i_dim)) then
            ijk(i_dim) = i
            return
        else
            ijk(i_dim) = first_cell(i_dim)
        end if
    end do
end subroutine

function contract_polarizability(alpha_full) result(alpha)
    real(8), intent(in) :: alpha_full(:, :)
    real(8) :: alpha(size(alpha_full, 1)/3)

    integer :: i_atom, i_xyz, j_xyz, dim_3n
    real(8) :: alpha_aniso(3, 3), alpha_diag(3)

    dim_3n = size(alpha_full, 1)
    do i_atom = 1, size(alpha)
        forall (i_xyz = 1:3, j_xyz = 1:3)
            alpha_aniso(i_xyz, j_xyz) &
                = sum(alpha_full(i_xyz:dim_3n:3, 3*(i_atom-1)+j_xyz))
        end forall
        alpha_diag = sdiagonalized(alpha_aniso)
        alpha(i_atom) = sum(alpha_diag)/3
    end do
end function contract_polarizability

subroutine get_frequency_grid(n, L, x, w)
    integer, intent(in) :: n
    real(8), intent(in) :: L
    real(8), intent(out) :: x(n), w(n)

    call gauss_legendre(n, x, w)
    w = 2*L/(1-x)**2*w
    x = L*(1+x)/(1-x)
    w = w(n:1:-1)
    x = x(n:1:-1)
end subroutine get_frequency_grid

subroutine gauss_legendre(n, r, w)
    integer, intent(in) :: n
    real(8), intent(out) :: r(n), w(n)

    integer, parameter :: q = 16
    integer, parameter :: n_iter = 1000
    real(q) :: x, f, df, dx
    integer :: k, iter, i
    real(q) :: Pk(0:n), Pk1(0:n-1), Pk2(0:n-2)

    if (n > 60) then
        call print_error( &
            'Cannot construct accurate Gauss-Legendre quadrature grids for n > 60.' &
        )
    end if
    if (n == 1) then
        r(1) = 0.d0
        w(1) = 2.d0
        return
    end if
    Pk2(0) = 1._q  ! k = 0
    Pk1(0:1) = (/ 0._q, 1._q /)  ! k = 1
    do k = 2, n
        Pk(0:k) = ((2*k-1)*(/ 0.0_q, Pk1(0:k-1) /)-(k-1)*(/ Pk2(0:k-2), 0._q, 0._q /))/k
        if (k < n) then
            Pk2(0:k-1) = Pk1(0:k-1)
            Pk1(0:k) = Pk(0:k)
        end if
    end do
    ! now Pk contains k-th Legendre polynomial
    do i = 1, n
        x = cos(pi*(i-0.25_q)/(n+0.5_q))
        do iter = 1, n_iter
            df = 0._q
            f = Pk(n)
            do k = n-1, 0, -1
                df = f + x*df
                f = Pk(k) + x*f
            end do
            dx = f/df
            x = x-dx
            if (abs(dx) < 10*epsilon(dx)) exit
        end do
        r(i) = dble(x)
        w(i) = dble(2/((1-x**2)*df**2))
    end do
end subroutine

real(8) elemental function alpha_dynamic_ts(freq, alpha_0, omega, C6) result(alpha)
    real(8), intent(in) :: freq, alpha_0
    real(8), intent(in), optional :: omega, C6

    if (present(omega)) then
        alpha = alpha_0/(1+(freq/omega)**2)
    else if (present(C6)) then
        alpha = alpha_0/(1+(freq/omega_eff(C6, alpha_0))**2)
    end if
end function

real(8) elemental function combine_C6(C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6)
    real(8), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j

    C6 = 2*C6_i*C6_j/(alpha_0_j/alpha_0_i*C6_i+alpha_0_i/alpha_0_j*C6_j)
end function

real(8) elemental function V_to_R(V) result(R)
    real(8), intent(in) :: V

    R = (3.d0*V/(4.d0*pi))**(1.d0/3)
end function

real(8) elemental function omega_eff(C6, alpha)
    real(8), intent(in) :: C6, alpha

    omega_eff = 4.d0/3*C6/alpha**2
end function

real(8) elemental function get_sigma_selfint(alpha) result(sigma)
    real(8), intent(in) :: alpha

    sigma = (sqrt(2.d0/pi)*alpha/3.d0)**(1.d0/3)
end function

real(8) function get_C6_from_alpha(ctx, alpha) result(C6)
    type(Context_t), intent(in) :: ctx
    real(8), intent(in) :: alpha(0:)

    C6 = 3.d0/pi*sum((alpha**2)*ctx%freq_grid_w)
end function

subroutine set_damping_parameters(damping, method, xc)
    type(Damping_t), intent(inout) :: damping
    character(len=*), intent(in) :: method, xc

    select case (method)
    case ('ts')
        damping%d = 20.d0
    case ('mbd@rsscs', 'mbd@ts')
        damping%a = 6.d0
    case ('mbd@scs')
        damping%beta = 1.d0
    end select
    select case (trim(method) // '@' // trim(xc))
    case ('ts@pbe')
        damping%s_R = 0.94d0
    case ('ts@pbe0', 'ts@hse')
        damping%s_R = 0.96d0
    case ('ts@blyp')
        damping%s_R = 0.62d0
    case ('ts@b3lyp')
        damping%s_R = 0.84d0
    case ('ts@revpbe')
        damping%s_R = 0.60d0
    case ('ts@am05')
        damping%s_R = 0.84d0
    case ('mbd@scs@pbe')
        damping%a = 2.56d0
    case ('mbd@ts@pbe')
        damping%beta = 0.81d0
    case ('mbd@rsscs@pbe')
        damping%beta = 0.83d0
    case ('mbd@scs@pbe0', 'mbd@scs@hse')
        damping%a = 2.53d0
    case ('mbd@ts@pbe0', 'mbd@ts@hse')
        damping%beta = 0.83d0
    case ('mbd@rsscs@pbe0', 'mbd@rsscs@hse')
        damping%beta = 0.85d0
    end select
end subroutine

function T_bare(rxyz) result(T)
    real(8), intent(in) :: rxyz(3)
    real(8) :: T(3, 3)

    integer :: i, j
    real(8) :: r_sq, r_5

    r_sq = sum(rxyz(:)**2)
    r_5 = sqrt(r_sq)**5
    do i = 1, 3
        T(i, i) = (3.d0*rxyz(i)**2-r_sq)/r_5
        do j = i+1, 3
            T(i, j) = 3.d0*rxyz(i)*rxyz(j)/r_5
            T(j, i) = T(i, j)
        end do
    end do
    T = -T
end function

real(8) function B_erfc(r, a) result(B)
    real(8), intent(in) :: r, a

    B = (erfc(a*r)+(2*a*r/sqrt(pi))*exp(-(a*r)**2))/r**3
end function

real(8) elemental function C_erfc(r, a) result(C)
    real(8), intent(in) :: r, a

    C = (3*erfc(a*r)+(2*a*r/sqrt(pi))*(3.d0+2*(a*r)**2)*exp(-(a*r)**2))/r**5
end function

function T_erfc(rxyz, alpha) result(T)
    real(8), intent(in) :: rxyz(3), alpha
    real(8) :: T(3, 3)

    integer :: i, j
    real(8) :: r, B, C

    r = sqrt(sum(rxyz(:)**2))
    B = B_erfc(r, alpha)
    C = C_erfc(r, alpha)
    do i = 1, 3
        do j = i, 3
            T(i, j) = -C*rxyz(i)*rxyz(j)
            if (i /= j) T(j, i) = T(i, j)
        end do
        T(i, i) = T(i, i) + B
    end do
end function

function damping_fermi(r, sigma, a) result(f)
    real(8), intent(in) :: r, sigma, a
    real(8) :: f

    f = 1.d0/(1+exp(-a*(r/sigma-1)))
end function

function damping_erf(r, sigma, a) result(f)
    real(8), intent(in) :: r, sigma, a
    real(8) :: f

    f = erf((r/sigma)**a)
end function

function damping_1mexp(r, sigma, a) result(f)
    real(8), intent(in) :: r, sigma, a
    real(8) :: f

    f = 1-exp(-(r/sigma)**a)
end function

function damping_overlap(r, overlap, C6, beta, a) result(f)
    real(8), intent(in) :: r, overlap, C6, beta, a
    real(8) :: f

    f = 1.d0-terf(-overlap/(erf(r/6)**6*C6/r**6), beta, a)
end function

function T_overlap_coulomb(rxyz, overlap, C6, beta, a) result(T)
    real(8), intent(in) :: rxyz(3), overlap, C6, beta, a
    real(8) :: T(3, 3)

    real(8) :: zeta_1, zeta_2
    real(8) :: r, erff, exp36, qene, qenep, qenepp

    r = sqrt(sum(rxyz**2))
    erff = erf(r/6)
    exp36 = exp(r**2/36)
    qene = overlap*r**6/(C6*erff**6)
    qenep = 2.d0*overlap*r**5*(-(1.d0/exp36)*r/sqrt(pi)+3.d0*erff)/(C6*erff**7)
    qenepp = (1.d0/exp36**2)*overlap*r**4/(9*C6*pi*erff**8) &
        * (42*r**2+exp36*sqrt(pi)*r*(-216+r**2)*erff &
        + 270.d0*exp36**2*pi*erff**2)
    zeta_1 = 1.d0/2*(2.d0-erf(a*(beta-qene))+erf(a*(beta+qene)) &
        + 2*a*r*qenep/sqrt(pi)*(-exp(-a**2*(beta-qene)**2) &
        - exp(-a**2*(beta+qene)**2)))
    zeta_2 = 1.d0/sqrt(pi)*a*exp(-a**2*(beta+qene)**2)*r**2 &
        * (2*a**2*qenep**2*(beta*(-1.d0+exp(4*a**2*beta*qene)) &
        - qene*(1.d0+exp(4*a**2*beta*qene))) &
        + qenepp*(1.d0+exp(4*a**2*beta*qene)))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function

function T_fermi_coulomb(rxyz, sigma, a) result(T)
    real(8), intent(in) :: rxyz(3), sigma, a
    real(8) :: T(3, 3)

    real(8) :: r_sigma, d_r_sigma_m_1, zeta_1, zeta_2

    r_sigma = sqrt(sum(rxyz**2))/sigma
    d_r_sigma_m_1 = a*(r_sigma-1)
    zeta_1 = 1.d0/(1.d0+exp(-d_r_sigma_m_1)) &
        - a/2.d0*r_sigma/(1.d0+cosh(-d_r_sigma_m_1))
    zeta_2 = 2.d0*a**2*r_sigma**2/sinh(-d_r_sigma_m_1)**3 &
        * sinh(-d_r_sigma_m_1/2.d0)**4
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function

function T_erf_coulomb(rxyz, sigma, a) result(T)
    real(8), intent(in) :: rxyz(3), sigma, a
    real(8) :: T(3, 3)

    real(8) :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = erf(r_sigma)-2.d0/sqrt(pi)*a*r_sigma*exp(-r_sigma**2)
    zeta_2 = -2.d0/sqrt(pi)*a*r_sigma*exp(-r_sigma**2) &
        * (1.d0+a*(-1.d0+2.d0*r_sigma**2))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function

function T_1mexp_coulomb(rxyz, sigma, a) result(T)
    real(8), intent(in) :: rxyz(3), sigma, a
    real(8) :: T(3, 3)

    real(8) :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = 1.d0-exp(-r_sigma)-a*r_sigma*exp(-r_sigma)
    zeta_2 = -r_sigma*a*exp(-r_sigma)*(1+a*(-1+r_sigma))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function

elemental function terf(r, r0, a)
    real(8), intent(in) :: r, r0, a
    real(8) :: terf

    terf = 0.5d0*(erf(a*(r+r0))+erf(a*(r-r0)))
end function

end module mbd
