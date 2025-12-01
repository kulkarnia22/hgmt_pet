#ifndef PORE_GEOMETRY_H
#define PORE_GEOMETRY_H

#include <stddef.h>  // size_t
#include "linear_algebra.h"

#define PORE_NUM_LAYERS 12

//layer inner and outer radii tuples in mm(layers are T = 1 inch thick)
#define PORE_LAYER_RADII_INIT { {450.0, 475.4}, {500.0, 525.4}, {555.0, 575.4},\
{600.0, 625.4}, {650.0, 675.4}, {700, 725.4}, {750, 775.4}, {800, 825.4},\
{850.0, 875.4}, {900.0, 925.4}, {950.0, 975.4}, {1000.0, 1025.4}}


#define STR(x)   #x
#define SHOW_DEFINE(x) printf("%s=%s\n", #x, STR(x))

#define ALPHA 5 //microns
#define BETA 50 //microns
#define GAMMA 1000 //microns
#define TAU 51 //microns

// tolerances for AXIS_SLAB / axis_slab
#define EPS_U 1e-12      // treat |u| < EPS_U as parallel (units: same as u, typically mm per t)
#define EPS_X 1e-9       // positional slack (units: mm if your (u,v) are mm)
#define EPS_T 1e-12      // param interval slack (unitless, t in [0,1])

// =============================
// --------- Types -------------
// =============================

typedef struct {
  double Rin, Rout; // mm
} PG_Radii;

extern const PG_Radii PG_LAYER_RADII[PORE_NUM_LAYERS];

typedef struct {
  // ---- inputs (configure per layer) ----
  double Rin, Rout;         // mm
  double p_phi_mm;          // tangential pitch (center-to-center) in mm
  double p_z_mm;            // axial pitch (center-to-center) in mm
  double phi0;              // anchor pore angle (rad) for (i,j)=(0,0)
  double z0;                // anchor pore z (mm)  for (i,j)=(0,0)
  double w_phi_mm;          // pore width along e_phi (mm) = α
  double w_z_mm;            // pore height along z (mm)    = β
  double z_min, z_max;      // active axial span (mm)

  // ---- derived (computed once) ----
  double Rstar;             // (Rin+Rout)/2
  double a_r;               // radial half-size = (Rout-Rin)/2
  double a_phi, a_z;        // half-sizes of opening (α/2, β/2)
  int    Nphi;              // pores around circumference
  double dphi;              // angular pitch = 2π / Nphi
  double p_phi_eff_mm;      // effective mm pitch = Rstar * dphi
} LatticeSpec;

extern LatticeSpec specs[PORE_NUM_LAYERS];  // declaration

// Search window (signed, per-axis)
typedef struct {
  int di_minus, di_plus;    // i ∈ [i0 - di_minus, i0 + di_plus]
  int dj_minus, dj_plus;    // j ∈ [j0 - dj_minus, j0 + dj_plus]
  double smax_eff;          // effective range used (may be radially capped)
} PG_SearchWindow;

// =============================
// ------- API functions -------
// =============================

// --- Small utilities (double-precision friendly) ---
double pg_wrap2pi(double a);                 // -> [0, 2π)
double pg_angdiff(double a, double b);       // signed smallest diff in (-π, π]
long   pg_nearest_int(double x);             // half-away-from-zero rounding
int    pg_imod(int a, int m);                // positive modulus for indices


void pg_local_dir_components_at_point(vec3d r0, vec3d uhat,
                                             double* ur, double* uphi, double* uzloc);
                                        
// --- Lattice setup ---
void pg_compute_pitches_from_dims(LatticeSpec* L,
                                  double alpha, double tau,
                                  double beta, double gamma);
// (sets L->p_phi_mm = alpha + tau, L->p_z_mm = beta + gamma)

void pg_lattice_init(LatticeSpec* L);        // computes derived fields

// Initialize one layer from raw inputs (α,β,τ,γ in microns)
void pg_init_layer_from_dims(LatticeSpec* L,
                             double Rin_mm, double Rout_mm,
                             double phi0_rad, double z0_mm,
                             double z_min_mm, double z_max_mm);

// --- Index/center conversions ---
void pg_pore_center_cartesian_from_ij(const LatticeSpec* L, long i, long j,
                                      vec3d* C, vec3d* er,
                                      vec3d* ephi, vec3d* ez);

// Nearest (i,j) to a world point (x,y,z)
void pg_nearest_indices_from_point(const LatticeSpec* L,
                                   vec3d scatter_point,
                                   long* i0, long* j0);

// --- Search window sizing (from dir + range) ---
void pg_compute_search_window(const LatticeSpec* L,
                              vec3d scatter_point,
                              vec3d scatter_dir, // unit
                              double smax_mm,                  // range in Kapton
                              int cap_by_radial_exit,          // 0/1
                              PG_SearchWindow* W);


// --- Intersection: does line segment enter pore (i,j)? ---
int pg_line_hits_pore_voxel(const LatticeSpec* L, long i, long j,
                            vec3d scatter_point,
                            vec3d scatter_dir, // unit
                            double smax_mm,
                            double* t_enter_mm);             // nullable

// --- Kapton electron range (user-supplied table) ---
void   pg_set_kapton_range_table(size_t n,
                                 const double* E_keV,   // strictly increasing
                                 const double* R_mm);   // mm
double pg_kapton_range_mm(double E_keV);                // linear interp; 0 if unset
double pg_kapton_energy_keV(double R_mm);

#endif