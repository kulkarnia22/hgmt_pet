// pore_geometry.c
#include "pore_geometry.h"
#include <math.h>
#include <float.h>
#include <string.h>

#define M_PI 3.14159265358979323846

//==============================
//  Layer radii (definition)
//==============================
const PG_Radii PG_LAYER_RADII[PORE_NUM_LAYERS] = PORE_LAYER_RADII_INIT;

LatticeSpec specs[PORE_NUM_LAYERS];  // the single definition

//==============================
//  Small helpers (internal)
//==============================
static inline double pg_hypot2(double x, double y){ return sqrt(x*x + y*y); }

//==============================
//  Public utils (match .h)
//==============================
double pg_wrap2pi(double a){
  double t = fmod(a, 2.0*M_PI);
  return (t < 0.0) ? t + 2.0*M_PI : t;
}

double pg_angdiff(double a, double b){
  double d = pg_wrap2pi(a) - pg_wrap2pi(b);
  if (d >  M_PI) d -= 2.0*M_PI;
  if (d <= -M_PI) d += 2.0*M_PI;
  return d;
}

long pg_nearest_int(double x){
  return (x >= 0.0) ? (long)floor(x + 0.5) : (long)(-floor(-x + 0.5));
}

int pg_imod(int a, int m){
  int r = a % m;
  return (r < 0) ? r + m : r;
}

static inline double um_to_mm(double u){ return u * 1e-3; }

//==============================
//  Lattice setup
//==============================
void pg_compute_pitches_from_dims(LatticeSpec* L,
                                  double alpha,  // mm pore opening along φ
                                  double tau,    // mm wall along φ
                                  double beta,   // mm pore opening along z
                                  double gamma)  // mm wall along z
{  
  L->w_phi_mm = alpha;
  L->w_z_mm   = beta;
  L->p_phi_mm = (alpha + tau);    // τ along φ
  L->p_z_mm   = (beta  + gamma);  // γ along z
}

void pg_lattice_init(LatticeSpec* L){
  //set pitches
  pg_compute_pitches_from_dims(L, um_to_mm(ALPHA), um_to_mm(BETA), um_to_mm(GAMMA), um_to_mm(TAU));  
  // Derived sizes
  L->Rstar = 0.5*(L->Rin + L->Rout);
  L->a_r   = 0.5*(L->Rout - L->Rin);
  L->a_phi = 0.5*L->w_phi_mm;
  L->a_z   = 0.5*L->w_z_mm;

  // Choose integer Nphi that tiles ring; derive exact dphi & effective pitch
  double Nphi_real = (2.0*M_PI * L->Rstar) /L->p_phi_mm;
  int Nphi = (int) llround(Nphi_real);
  if (Nphi < 1) Nphi = 1;
  L->Nphi = Nphi;
  //printf("Nphi = %i\n", L->Nphi);
  L->dphi = 2.0*M_PI / (double)Nphi;
  L->p_phi_eff_mm = L->Rstar * L->dphi;
}

void pg_init_layer_from_dims(LatticeSpec* L,
                             double Rin_mm, double Rout_mm,
                             double phi0_rad, double z0_mm,
                             double z_min_mm, double z_max_mm)
{
    L->Rin   = Rin_mm;
    L->Rout  = Rout_mm;
    L->phi0  = phi0_rad;
    L->z0    = z0_mm;
    L->z_min = z_min_mm;
    L->z_max = z_max_mm;

    pg_lattice_init(L);
}

//==============================
//  (i,j) -> center + local basis
//==============================
void pg_pore_center_cartesian_from_ij(const LatticeSpec* L, long i, long j,
                                      vec3d* C, vec3d* er, vec3d* ephi, vec3d* ez)
{
  //long iw = (L->Nphi > 0) ? pg_imod((int)i, L->Nphi) : 0;
  long iw = pg_imod((int)i, L->Nphi);
  double phic = pg_wrap2pi(L->phi0 + iw * L->dphi);
  double zc   = L->z0 + j * L->p_z_mm;

  double c = cos(phic), s = sin(phic);
  if (C)    { C->x = L->Rstar*c; C->y = L->Rstar*s; C->z = zc; }
  if (er)   { er->x = c; er->y = s; er->z = 0.0; }
  if (ephi) { ephi->x = -s; ephi->y =  c; ephi->z = 0.0; }
  if (ez)   { ez->x = 0.0; ez->y = 0.0; ez->z = 1.0; }
}

//==============================
//  cartesian -> nearest (i,j)
//==============================
void pg_nearest_indices_from_point(const LatticeSpec* L,
                                   vec3d scatter_point,
                                   long* i0, long* j0)
{
  double x = scatter_point.x, y = scatter_point.y, z = scatter_point.z;

  //double phi_pt = pg_wrap2pi(atan2(y, x));
  double phi_pt = pg_wrap2pi(atan2(y, x));
  //double pphi   = (L->p_phi_eff_mm > 0.0) ? L->p_phi_eff_mm : 1.0;
  double pphi = L->p_phi_eff_mm;
  double di  = (L->Rstar * pg_angdiff(phi_pt, L->phi0)) / pphi;
  //double dj     = (z - L->z0) / ((L->p_z_mm > 0.0) ? L->p_z_mm : 1.0);
  double dj = (z - L->z0) / L->p_z_mm;

  if (i0) *i0 = pg_nearest_int(di);
  if (j0) *j0 = pg_nearest_int(dj);
}

//==============================
//  Search window sizing
//==============================
void pg_local_dir_components_at_point(vec3d r0, vec3d uhat,
                                             double* ur, double* uphi, double* uzloc)
{
  double rxy = pg_hypot2(r0.x, r0.y);
  double rx = r0.x / rxy;
  double ry = r0.y / rxy;
  double phix = -ry, phiy = rx;

  if (ur)    *ur    = uhat.x*rx  + uhat.y*ry;     // u_r
  if (uphi)  *uphi  = uhat.x*phix+ uhat.y*phiy;   // u_phi
  if (uzloc) *uzloc = uhat.z;                     // u_z
}

static double pg_time_to_radial_exit_mm(const LatticeSpec* L, vec3d r0, double ur){
  double r0mag = pg_hypot2(r0.x, r0.y);
  if (fabs(ur) < 1e-16) return DBL_MAX;
  if (ur > 0.0){
    double dr = L->Rout - r0mag;
    return (dr / ur);
  } else {
    double dr = r0mag - L->Rin;
    return (dr / (-ur));
  }
}

void pg_compute_search_window(const LatticeSpec* L,
                              vec3d scatter_point,
                              vec3d scatter_dir,
                              double smax_mm,
                              int cap_by_radial_exit,
                              PG_SearchWindow* W)
{
  // Normalize direction (defensive)
  double un = vec_mag(scatter_dir);
  vec3d u = scatter_dir;
  if (un > 0.0){ u = vec_norm(u); }

  double ur, uphi, uzloc;
  pg_local_dir_components_at_point(scatter_point, u, &ur, &uphi, &uzloc);

  double smax_eff = smax_mm;
  if (cap_by_radial_exit){
    double t_rad = pg_time_to_radial_exit_mm(L, scatter_point, ur);
    if (t_rad < smax_eff) smax_eff = t_rad;
  }
  if (smax_eff < 0.0) smax_eff = 0.0;

  // per-axis budgets (signed split)
  double sphi_plus  = (uphi  > 0.0) ?  smax_eff* uphi  : 0.0;
  double sphi_minus = (uphi  < 0.0) ? -smax_eff* uphi  : 0.0;
  double sz_plus    = (uzloc > 0.0) ?  smax_eff* uzloc : 0.0;
  double sz_minus   = (uzloc < 0.0) ? -smax_eff* uzloc : 0.0;

  double aphi = L->a_phi, az = L->a_z;
  double pphi = (L->p_phi_eff_mm > 0.0) ? L->p_phi_eff_mm : 1.0;
  double pz   = (L->p_z_mm        > 0.0) ? L->p_z_mm        : 1.0;

  int di_plus  = (int)ceil((sphi_plus  + aphi)/pphi - 1e-12);
  int di_minus = (int)ceil((sphi_minus + aphi)/pphi - 1e-12);
  int dj_plus  = (int)ceil((sz_plus    + az  )/pz   - 1e-12);
  int dj_minus = (int)ceil((sz_minus   + az  )/pz   - 1e-12);

  if (di_plus  < 0) di_plus  = 0;
  if (di_minus < 0) di_minus = 0;
  if (dj_plus  < 0) dj_plus  = 0;
  if (dj_minus < 0) dj_minus = 0;

  W->di_plus  = di_plus;
  W->di_minus = di_minus;
  W->dj_plus  = dj_plus;
  W->dj_minus = dj_minus;
  W->smax_eff = smax_eff;
}


//==============================
//  Intersection (slab test)
//==============================
static inline int axis_slab(double p0, double u, double a,
                            double *tmin, double *tmax)
{
    if (fabs(u) < EPS_U) {
        // Parallel to the slab planes: reject if we're outside by more than EPS_X
        if (fabs(p0) > (a + EPS_X)) {
            // fprintf(stderr, "axis_slab: parallel reject |p0|=%g > %g\n", fabs(p0), a + EPS_X);
            return 0;
        }
        // parallel but inside: no tightening on [tmin,tmax]
    } else {
        double t1 = (-(a) - p0) / u;
        double t2 = ( +(a) - p0) / u;
        if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }

        if (t1 > *tmin) {*tmin = t1;}
        if (t2 < *tmax) {*tmax = t2;}
        //if (*tmin*u < *tmax){printf("promise"); return 0;}
        if (*tmin > *tmax) {
            //fprintf(stderr, "axis_slab: empty interval tmin=%g > tmax=%g\n", *tmin, *tmax);
            return 0; 
        }
    }
    return 0; // keep going
}
//This function contains the main logic of the detection scheme and may need to be rewritten
static int pg_line_hits_box_local(double p0r, double p0p, double p0z,
                                  double ur, double up, double uzl,
                                  double a_r, double a_p, double a_z,
                                  double smax,
                                  double* t_enter)
{

  double tmin = 0.0;
  double tmax = smax;
  
  if (!axis_slab(p0r, ur,  a_r,  &tmin, &tmax)) {/*printf("here!");*/ return 0;}
  if (!axis_slab(p0p, up,  a_p,  &tmin, &tmax)) {/*printf("here!");*/ return 0;}
  if (!axis_slab(p0z, uzl, a_z,  &tmin, &tmax)) {/*printf("here!");*/ return 0;}
  
  return 1;
}

int pg_line_hits_pore_voxel(const LatticeSpec* L, long i, long j,
                            vec3d scatter_point,
                            vec3d scatter_dir,
                            double smax_mm,
                            double* t_enter_mm)
{
  // Normalize direction (defensive)
  double un = vec_mag(scatter_dir);
  if (un == 0.0) return 0;
  vec3d u = vec_norm(scatter_dir);

  // Build pore center + basis
  vec3d C, er, ephi, ez;
  pg_pore_center_cartesian_from_ij(L, i, j, &C, &er, &ephi, &ez);

  // Local start coords
  vec3d d = vec_sub(scatter_point,C);
  /*printf("difference: %f, ", vec_mag(d));
  printf("max range: %f, ", smax_mm);*/
  //if(vec_mag(d) < smax_mm){printf("here!");}

  double p0r = vec_dot(d,er);
  double p0p = vec_dot(d,ephi);
  double p0z = vec_dot(d,ez);

  // Local direction
  double ur = vec_dot(u,er);
  double up = vec_dot(u,ephi);
  double uz = vec_dot(u,ez);

  return pg_line_hits_box_local(p0r, p0p, p0z, ur, up, uz,
                                L->a_r, L->a_phi, L->a_z,
                                smax_mm, t_enter_mm);
}


//==============================
//  Kapton range table (interp)
//==============================
static const double* g_range_E_keV = NULL;
static const double* g_range_R_mm  = NULL;
static size_t        g_range_N     = 0;

void pg_set_kapton_range_table(size_t n,
                               const double* E_keV,
                               const double* R_mm)
{
  g_range_E_keV = E_keV;
  g_range_R_mm  = R_mm;
  g_range_N     = n;
}

double pg_kapton_range_mm(double E_keV)
{
  if (!g_range_E_keV || !g_range_R_mm || g_range_N < 2) return 0.0;
  if (E_keV <= g_range_E_keV[0])              return g_range_R_mm[0];
  if (E_keV >= g_range_E_keV[g_range_N-1])    return g_range_R_mm[g_range_N-1];

  // binary search
  size_t lo = 0, hi = g_range_N - 1;
  while (hi - lo > 1){
    size_t mid = lo + (hi - lo)/2;
    if (E_keV < g_range_E_keV[mid]) hi = mid; else lo = mid;
  }
  double x0 = g_range_E_keV[lo], x1 = g_range_E_keV[hi];
  double y0 = g_range_R_mm[lo],  y1 = g_range_R_mm[hi];
  double t  = (E_keV - x0) / (x1 - x0);
  return y0 + t*(y1 - y0);
}



