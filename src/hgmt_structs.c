#include <math.h>
#include <stdio.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "linear_algebra.h"
#include "pore_geometry.h"

hit event_to_hit(event *single_event) {
  vec_mag(three_vec(single_event->position.x, single_event->position.y, 0.0));
  hit new_hit;
  new_hit.source = single_event;
  vec3d z_hat = three_vec(0.0, 0.0, 1.0);
  vec3d circ_hat = vec_norm(vec_cross(z_hat, single_event->position));
  new_hit.position =
      vec_add(single_event->position, vec_scale(z_hat, gaussian(LONG_UNC)));
  vec3d offset = vec_add(vec_scale(z_hat, gaussian(LONG_UNC)),
                         vec_scale(circ_hat, gaussian(CIRC_UNC)));
  if (DETECTOR_SEGMENTATION) {
    // we move the radial component to the midpoint of the detector which it hit
    double rad_dist = radial_dist(single_event->position);
    single_event->position = radial_scale(
        single_event->position, (detector_positions[single_event->detector_id] +
                                 DETECTOR_THICKNESS / 2) /
                                    rad_dist);
  } else {
    vec3d r_hat = vec_norm(
        three_vec(single_event->position.x, single_event->position.y, 0));
    offset = vec_add(offset, vec_scale(r_hat, gaussian(RAD_UNC)));
  }
  new_hit.tof = single_event->tof + gaussian(TIME_UNC);
  new_hit.position = vec_add(single_event->position, offset);
  return new_hit;
}
bool is_detected(event *single_event, double eff_by_energy[COLS]) {
  /*if (single_event->energy > 250){
    printf("Energy: %.3f\n", single_event->energy);
    printf("det_prob: %.3f\n", 
        linear_interpolation(eff_by_energy, E_MIN, E_MAX, single_event->energy));
  }*/
  //printf("%f, ", single_event->energy);  
  return (single_event->detector_id != -1 &&
          drand48() < linear_interpolation(eff_by_energy, E_MIN, E_MAX,
                                           single_event->energy));
}

static inline double um_to_mm(double u){ return u * 1e-3; }

static inline double wrap_pi(double a) {
    a = fmod(a + PI, 2.0*PI);
    if (a < 0) a += 2.0*PI;
    return a - PI;  // in (-π, π]
}

static inline double wrap_2pi(double phi)
{
    phi = fmod(phi, 2.0*PI);
    if (phi < 0)
        phi += 2.0*PI;
    return phi;
}

bool plane_crossingv2(event *single_event){
    vec3d scatter_pos = vec_scale(single_event->position, 10);
    vec3d u = vec_norm(single_event->direction);
    double R_scatter = radial_dist(scatter_pos);
    double phi_scatter = atan2(scatter_pos.y, scatter_pos.x);
    double tau = um_to_mm(TAU);
    double alpha = um_to_mm(ALPHA);
    double i = pg_nearest_int((R_scatter*phi_scatter)/(tau + alpha));//index of nearest pore

    //now I need to see if I am moving away from or towards this pore
    double phi_pore = wrap_pi(i*(tau + alpha)/R_scatter);
    //need to check if I am inside pore plane
    if(R_scatter*fabs(wrap_pi(phi_scatter - phi_pore)) < alpha/2){
        return 0;
    }

    double phi_diff = wrap_pi(phi_pore - phi_scatter);
    vec3d ephi = three_vec(-sin(phi_scatter), cos(phi_scatter), 0);
    double uphi = vec_dot(u, ephi);
    
    //adjusting pore to look at based on direction of scatter
    if(uphi * phi_diff < 0){
        if(phi_pore > phi_scatter){
            i -= 1;
            phi_pore = wrap_pi(i*(tau + alpha)/R_scatter);
        }
        else{
            i += 1;
            phi_pore = wrap_pi(i*(tau + alpha)/R_scatter);
        }
    }
    //now we check if the electron reaches a pore
    double scatter_energy = single_event->energy;
    double s_max_mm = pg_kapton_range_mm(scatter_energy);
    //first, we need to make adjustments to s_max_mm in case electron leaves first
    //I will deal with that once I get things working
    vec3d final_pos = vec_add(scatter_pos, vec_scale(u, s_max_mm));
    double r1 = R_scatter;
    double r2 = radial_dist(final_pos);
    double final_phi = atan2(final_pos.y, final_pos.x);
    phi_diff = wrap_pi(final_phi - phi_scatter);
    double pore_scatter_diff = wrap_pi(phi_pore - phi_scatter);
    double r = (r1 + r2)/2;
    double alpha_diff = wrap_pi(alpha/(2*r));
    double angle_to_pore = wrap_pi(pore_scatter_diff - alpha_diff);
    if (fabs(r*phi_diff) > fabs(r*angle_to_pore)){
        return 1;
    }
    return 0;

    //commented out bottom due to problematic edge case that I try to fix above
    /*double phi1_mm = r1*phi_scatter;
    double phi2_mm = r2*atan2(final_pos.y,final_pos.x); 
    double phi_travelled = phi2_mm - phi1_mm;
    if (phi_pore > phi_scatter){
        if (fabs(phi_travelled) > fabs(((r1+r2)/2)*(phi_pore - phi_scatter) - alpha/2)){
            return 1;
        }
        else{
            return 0;
        }
    }
    if(phi_pore < phi_scatter){
        if (fabs(phi_travelled) > fabs(((r1+r2)/2)*(phi_pore - phi_scatter) + alpha/2)){
            return 1;
        }
    }
    return 0;*/
}

int pores_crossed(event *first_hit){
    /*Here I will do the analysis on the path of our first hits
    This will tell me how many pores are crossed by scatters resulting in a hit.
    In another function that I will likely write after this, I will get the 
    energy of my electron at the entrance of each pore by playing with my
    kapton range table. 

    Importantly, this function will rerun much of the logic in plan_crosingv2.
    It is really silly how much I am just repeating. Ideally I just do everything in 
    plane_crossingv2 but to interpret that I would have to change the structure of the
    post processing flow too much.
    */

    vec3d hit_pos = vec_scale(first_hit->position, 10);
    vec3d u = vec_norm(first_hit->direction);
    double hit_R = sqrt(hit_pos.x*hit_pos.x + hit_pos.y*hit_pos.y);
    double hit_phi = wrap_2pi(atan2(hit_pos.y, hit_pos.x));
    double tau = um_to_mm(TAU);
    double alpha = um_to_mm(ALPHA);
    int i = pg_nearest_int((hit_R*hit_phi)/(tau + alpha));//index of nearest pore

    double pitch_phi = (tau + alpha)/hit_R;
    double phi_pore  = i * pitch_phi; 
    vec3d ephi = three_vec(-sin(hit_phi), cos(hit_phi), 0);
    double uphi = vec_dot(u, ephi);
    int dir = (uphi >= 0) ? 1 : -1;
    double dphi_fwd = (dir > 0) ? wrap_2pi(phi_pore - hit_phi): wrap_2pi(hit_phi - phi_pore);
    if (dphi_fwd > 0.5*pitch_phi){
        i += dir;
        phi_pore = i*pitch_phi;
    }
    double phi_diff = (dir > 0) ? wrap_2pi(phi_pore - hit_phi): wrap_2pi(hit_phi - phi_pore);

    double scatter_energy = first_hit->energy;
    double s_max_mm = pg_kapton_range_mm(scatter_energy);
    vec3d final_pos = vec_add(hit_pos, vec_scale(u, s_max_mm));
    double r1 = hit_R;
    double r2 = radial_dist(final_pos);
    double final_phi = wrap_2pi(atan2(final_pos.y, final_pos.x));

    //once I'm in mm I'm good
    //but because of angle trouble I need to stay in phi coords until the very end
    //technically electron won't lose energy in the pore
    double r = (r1 + r2)/2;
    double phi_alpha_diff = alpha/(2*r);
    double new_phi = hit_phi + dir * phi_diff - dir * phi_alpha_diff;
    double new_phi_diff = (dir > 0) ? wrap_2pi(final_phi - wrap_2pi(new_phi))
                                 : wrap_2pi(wrap_2pi(new_phi) - final_phi);
    int num_pores_crossed = 1 + floor(fabs(r*new_phi_diff/tau));
    return num_pores_crossed;

}

double min_energy(event *first_hit, int num_crosses){
    /*for this I just need to start with the energy of my scatter
    and then see how much energy it loses by the time it reaches the first pore.
    After that, for each new pore it crosses, I just need to decrease energy by
    tau mm equivalent.*/
    vec3d hit_pos = vec_scale(first_hit->position, 10);
    vec3d u = vec_norm(first_hit->direction);
    double hit_R = sqrt(hit_pos.x*hit_pos.x + hit_pos.y*hit_pos.y);
    double hit_phi = wrap_2pi(atan2(hit_pos.y, hit_pos.x));
    double tau = um_to_mm(TAU);
    double alpha = um_to_mm(ALPHA);
    double i = pg_nearest_int((hit_R*hit_phi)/(tau + alpha));//index of nearest pore

    double phi_pore = i*(tau + alpha)/hit_R;

    double phi_diff = wrap_2pi(phi_pore - hit_phi);
    vec3d ephi = three_vec(-sin(hit_phi), cos(hit_phi), 0);
    double uphi = vec_dot(u, ephi);
    
    //adjusting pore to look at based on direction of scatter
    if(uphi * phi_diff < 0){
        if(phi_pore > hit_phi){
            i -= 1;
            phi_pore = wrap_2pi(i*(tau + alpha)/hit_R);
        }
        else{
            i += 1;
            phi_pore = wrap_2pi(i*(tau + alpha)/hit_R);
        }
    }

    double scatter_energy = first_hit->energy;
    double s_max_mm = pg_kapton_range_mm(scatter_energy);
    vec3d final_pos = vec_add(hit_pos, vec_scale(u, s_max_mm));
    double r1 = hit_R;
    double r2 = radial_dist(final_pos);
    //double final_phi = atan2(final_pos.y, final_pos.x);

    //now I need to find distance to first pore
    double r = (r1 + r2)/2;
    phi_diff = wrap_2pi(phi_pore - hit_phi);
    double phi_alpha_diff = wrap_2pi(alpha/(2*r));
    double phi_prime = hit_phi + phi_diff - phi_alpha_diff;
    //first need to find change of phi to first pore
    double t = (hit_pos.x*sin(phi_prime) - hit_pos.y*cos(phi_prime))/(u.y*cos(phi_prime) - u.x*sin(phi_prime));
    double new_range = s_max_mm - t;
    //double min_energy = scatter_energy - pg_kapton_energy_keV(new_range);
    double min_energy = pg_kapton_energy_keV(new_range);
    int num_pores = num_crosses - 1;
    int dir = (uphi >= 0) ? +1 : -1;
    double phi_per_pore = dir * (tau/r);
    double t_old = t;
    if (t < 0){
        printf("event detected?: %d\n", first_hit->detected);
        printf("s_max_mm: %f\n", s_max_mm);
        printf("initial pos pi: (%f, %f, %f)\n", hit_pos.x, hit_pos.y, hit_pos.z);
        printf("unit direction u: (%f, %f, %f)\n", u.x, u.y, u.z);
        printf("hit_phi: %f\n", hit_phi);
        printf("phi_pore: %f\n", phi_pore);
        printf("t_first value: %f\n\n", t);
    }
    for (int i = 0; i < num_pores; i ++){
        phi_prime += phi_per_pore;
        //phi_prime = wrap_pi(phi_prime);
        t = (hit_pos.x*sin(phi_prime) - hit_pos.y*cos(phi_prime))/(u.y*cos(phi_prime) - u.x*sin(phi_prime)) - t_old;
        new_range -= t;
        t_old = (hit_pos.x*sin(phi_prime) - hit_pos.y*cos(phi_prime))/(u.y*cos(phi_prime) - u.x*sin(phi_prime));
        if (t < 0){
            printf("event detected?: %d\n", first_hit->detected);
            printf("s_max_mm: %f\n", s_max_mm);
            printf("initial pos pi: (%f, %f, %f)\n", hit_pos.x, hit_pos.y, hit_pos.z);
            printf("unit direction u: (%f, %f, %f)\n", u.x, u.y, u.z);
            printf("hit_phi: %f\n", hit_phi);
            printf("phi_pore: %f\n", phi_pore);
            printf("t value: %f\n\n", t);
        }
        //small cases where t is negative. will have to work through that. That doesn't make any sense tbh.
        //min_energy -= pg_kapton_energy_keV(new_range);
        min_energy = pg_kapton_energy_keV(new_range);
    }
    return min_energy;    
}

bool plane_crossing(event *single_event){
    /*first check: is scatter position in alpha beta space?
    to do this I need to find the i coord of the nearest pore, then
    check if the phi coord of my scatter is within the phi space of my voxel
    */

   vec3d scatter_pos = vec_scale(single_event->position, 10);
   double scatter_R = sqrt(scatter_pos.x*scatter_pos.x + scatter_pos.y*scatter_pos.y); 
   //printf("scatter pos = (%f,%f,%f)\n", scatter_pos.x, scatter_pos.y, scatter_pos.z);

   if (single_event->detector_id == -1){
    return 0;
   }
   LatticeSpec L = specs[single_event->detector_id];
   long i0;
   long j0;
   pg_nearest_indices_from_point(&L, scatter_pos, &i0, &j0);
   double phi_center = pg_wrap2pi(L.phi0 + i0*L.dphi);
   double phi_scatter = pg_wrap2pi(atan2(scatter_pos.y,scatter_pos.x));

   //check if scatter is in alpha/beta voxel space
   //printf("in pore check: %f\n", fabs(scatter_R * (pg_angdiff(phi_center, phi_scatter))));
   if(fabs(scatter_R * (pg_angdiff(phi_center, phi_scatter))) < L.a_phi){
    //count1 ++;
    //printf("edge case");
    //printf("in pore\n");
    return 0;
   }
   
   /*ok now I need to determine the plane crossing direction to check*/
   vec3d C;
   vec3d er;
   vec3d ephi;
   vec3d ez;
   pg_pore_center_cartesian_from_ij(&L, i0, j0, &C, &er, &ephi, &ez);

   vec3d u = vec_norm(single_event->direction);
   double uphi_comp = vec_dot(u, ephi);
   //check if scatter is moving away from the nearest pore
   if(uphi_comp < 0){
    if(pg_wrap2pi(atan2(C.y,C.x)) > pg_wrap2pi(atan2(scatter_pos.y,scatter_pos.x))){
        i0 -= 1;
        pg_pore_center_cartesian_from_ij(&L, i0, j0, &C, &er, &ephi, &ez);
    }
    else{
        i0 += 1;
        pg_pore_center_cartesian_from_ij(&L, i0, j0, &C, &er, &ephi, &ez);
    }
   }

   /*Ok, at this point I know which plane to look at. Now I need to determine
   if my scatter crosses that plane. I will do this with the range and direction 
   of my scatter. I also need to check for radial or axial exit.*/
   double scatter_energy = single_event->energy;
   double s_max_mm = pg_kapton_range_mm(scatter_energy);
   double ur;
   double uphi;
   double uzloc;
   pg_local_dir_components_at_point(scatter_pos, u, &ur, &uphi, &uzloc);
   double R, z;
   if(ur > 0){
    R = L.Rout;
   }
   else{
    R = L.Rin; 
   }

   if(uzloc < 0){
    z = L.z_min;
   }
   else{
    z = L.z_max;
   }
   //checking for early radial or axial exit
   double t_exit = fmin(fabs(R-scatter_R)/fabs(ur), fabs(z-scatter_pos.z)/fabs(uzloc));
   if(t_exit < s_max_mm){
    //count1 ++;
    //printf("edge case");
    s_max_mm = t_exit;
   }

   //now I can check if the scatter crosses its plane
   double Cphi = scatter_R*pg_wrap2pi(atan2(C.y,C.x));
   double scatter_phi = scatter_R*pg_wrap2pi(atan2(scatter_pos.y,scatter_pos.x));
   if(Cphi < scatter_phi){
    if((scatter_R*fabs(s_max_mm*uphi)) > fabs(Cphi + L.a_phi - scatter_phi)){
        return 1;
   }
    else if((scatter_R*fabs(s_max_mm*uphi)) > fabs(Cphi - L.a_phi - scatter_phi)){
        return 1;
    }
   }
   return 0;

}
//extern LatticeSpec specs[PORE_NUM_LAYERS];
// helper: reject scatters that *start inside* the nearest pore
static inline int starts_inside_nearest_pore(const LatticeSpec* L, vec3d P){
    long i0, j0;
    pg_nearest_indices_from_point(L, P, &i0, &j0);

    vec3d C, er, ephi, ez;
    pg_pore_center_cartesian_from_ij(L, i0, j0, &C, &er, &ephi, &ez);
    //printf("%f, ", vec_mag(vec_sub(C, P)));

    double dx = P.x - C.x, dy = P.y - C.y, dz = P.z - C.z;
    double pr = dx*er.x   + dy*er.y;      // er.z = 0
    double pp = dx*ephi.x + dy*ephi.y;    // ephi.z = 0
    double pz = dz;                       // ez = (0,0,1)

    const double eps = 1e-9 * fmax(L->a_r, fmax(L->a_phi, L->a_z));
    return (fabs(pr) <= L->a_r+eps) && (fabs(pp) <= L->a_phi+eps) && (fabs(pz) <= L->a_z+eps);
}

// minimal boolean detector (geometry only)
bool is_detected_geom(event* e){
    if (!e || e->detector_id < 0 || e->detector_id >= PORE_NUM_LAYERS){
        //printf("first, ");
        return false;
    }

    const LatticeSpec* L = &specs[e->detector_id];

    // energy -> Kapton range (mm)
    double smax = pg_kapton_range_mm(e->energy);
    //printf("%f, ", e->energy);
    if (smax <= 0.0){
        //printf("smax <= 0, ");    
        return false;
    }

    // unit direction
    vec3d u = e->direction;
    double un = vec_mag(u);
    if (un <= 0.0){
        //printf("un <= 0, ");
        return false;
    }
    u = vec_norm(u);

    // discard if starting inside a pore
    if (starts_inside_nearest_pore(L, e->position)){
        //printf("inside a pore, ");
        return false;
    }

    // seed and window
    long i0, j0;
    pg_nearest_indices_from_point(L, e->position, &i0, &j0);

    PG_SearchWindow W;
    pg_compute_search_window(L, e->position, u, smax, /*cap_by_radial_exit=*/1, &W);
    if (W.smax_eff <= 0.0){
        return false;
    }

    // scan candidate pores; first true hit wins
    //const double T_EPS = 1e-12;  // ignore t≈0 (start-inside noise)
    for (int dj = -W.dj_minus; dj <= W.dj_plus; dj++){
        for (int di = -W.di_minus; di <= W.di_plus; di++){
            long i = i0 + di, j = j0 + dj;
            double t_enter = 0.0;
            //printf("here");
            if (pg_line_hits_pore_voxel(L, i, j, e->position, u, W.smax_eff, &t_enter)){
                //printf("here!");
                return true;
            }
        }
    }
    //printf("last, ");
    return false;
}

int get_detector(vec3d position) {
  double rad_dist = radial_dist(position);
  for (int i = 0; i < sizeof(detector_positions) / sizeof(double); i++) {
    if (rad_dist > detector_positions[i] &&
        rad_dist < detector_positions[i] + DETECTOR_THICKNESS) {
      return i;
    }
  }
  return -1;
}
void free_annihilation(annihilation *annihil) {
  free(annihil->events);
  free(annihil->hits);
  free(annihil->photon1.events);
  free(annihil->photon1.hits);
  free(annihil->photon2.events);
  free(annihil->photon2.hits);
}
