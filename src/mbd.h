struct Damping {
    char *label;
    double d, s_R, a, beta;
    double *alpha_0, *C6, *R_vdw, *overlap, *custom, *custom_potential;
    int n_atoms;
};

struct Context {
    _Bool do_rpa, do_reciprocal, do_ewald;
    double ts_energy_accuracy, ts_shell_thickness,
           dipole_low_dim_cutoff, mayer_scaling,
           ewald_real_cutoff_scaling, ewald_rec_cutoff_scaling;
    double *freq_grid, *freq_grid_w;
    int n_freq;
};

struct Geometry {
    double *coords;
    _Bool is_periodic, vacuum_axes[3];
    double lattice_vectors[9], k_points;
    int supercell[3];
    int n_atoms, n_kpoints;
};

struct MBDResults {
    double energy;
    _Bool has_error;
    double *omegas, *modes, *bandomegas, *rpa_orders;
    double *modes_im, *bandmodes;  // _Complex in fact
    int n_atoms, n_orders, n_kpoints;
};

struct Relay {
    double *re;
    double *cplx; // _Complex in fact
    int n_atoms;
};

extern int my_task, n_tasks;
extern double bohr;

void _init_grid(struct Context *ctx, int n);
void _destroy_grid(struct Context *ctx);
void _run_scs(struct Context ctx, struct Geometry geom, double * alpha,
             struct Damping damping, struct Relay alpha_full,
             double lambda, double *k_point);
double _get_C6_from_alpha(struct Context ctx, double *alpha);
void _omega_eff(double *C6, double *alpha, double *omega);
void _run_mbd(struct Context ctx, struct Geometry geom, double *alpha_0,
             double *omega, struct Damping damping,
             struct MBDResults results);
