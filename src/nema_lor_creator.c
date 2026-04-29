// including standard files
#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "linear_algebra.h"
#include "read_write.h"
#include "pore_geometry.h"

#define THREADS 20

nema_lor apply_hit_uncertainty(nema_lor *n) {
    nema_lor out = *n;

    vec3d z_hat = three_vec(0.0, 0.0, 1.0);
    // Endpoint 1 uncertainty
    {
        vec3d p = out.pos1;

        vec3d r_hat = vec_norm(three_vec(p.x, p.y, 0.0));
        vec3d rad_offset = vec_scale(r_hat, gaussian(RAD_UNC));
        p = vec_add(p, rad_offset);

        vec3d circ_hat = vec_norm(vec_cross(z_hat, p));
        vec3d tangential_offset = vec_scale(circ_hat, gaussian(CIRC_UNC));
        vec3d longitudinal_offset = vec_scale(z_hat, gaussian(LONG_UNC));

        p = vec_add(p, vec_add(longitudinal_offset, tangential_offset));

        out.pos1 = p;
        out.tof1 = out.tof1 + gaussian(TIME_UNC);
    }
    // Endpoint 2 uncertainty
    {
        vec3d p = out.pos2;

        vec3d r_hat = vec_norm(three_vec(p.x, p.y, 0.0));
        vec3d rad_offset = vec_scale(r_hat, gaussian(RAD_UNC));
        p = vec_add(p, rad_offset);

        vec3d circ_hat = vec_norm(vec_cross(z_hat, p));
        vec3d tangential_offset = vec_scale(circ_hat, gaussian(CIRC_UNC));
        vec3d longitudinal_offset = vec_scale(z_hat, gaussian(LONG_UNC));

        p = vec_add(p, vec_add(longitudinal_offset, tangential_offset));

        out.pos2 = p;
        out.tof2 = out.tof2 + gaussian(TIME_UNC);
    }

    return out;
}

sym_matrix hit_covariance(vec3d position) {
  sym_matrix covariance;
  double rad = vec_mag(position);
  double cos = position.x / rad;
  double sin = position.y / rad;
  covariance.xx = RAD_VAR * cos * cos + CIRC_VAR * sin * sin;
  covariance.xy = (RAD_VAR - CIRC_VAR) * sin * cos;
  covariance.xz = 0;
  covariance.yy = CIRC_VAR * cos * cos + RAD_VAR * sin * sin;
  covariance.yz = 0;
  covariance.zz = LONG_VAR;
  return covariance;
}

lor create_lor_from_nema(nema_lor *n) {

  vec3d a = n->pos1;
  vec3d b = n->pos2;

  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_scale(vec_add(a, b), 0.5);
  vec3d c_hat = vec_norm(c);

  // Same convention as old create_lor:
  // delta_t = -(tof1 - tof2) = tof2 - tof1
  double delta_t = n->tof2 - n->tof1;

  vec3d displacement_from_center =
      vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);

  vec3d annihilation_loc =
      vec_add(geometric_center, displacement_from_center);

  lor new_lor;
  new_lor.center = annihilation_loc;

  // Non-collinearity covariance, transverse to the LOR direction
  vec3d u = vec_norm(vec_sub(b, a));
  sym_matrix transverse_projector =
      sym_add(sym_id(1.0), sym_scale(sym_proj(u), -1.0));

  sym_matrix nc_cov =
      sym_scale(transverse_projector, NONCOLLINEARITY_VARIANCE);

  // Endpoint spatial covariance contribution
  new_lor.covariance =
      sym_scale(sym_add(hit_covariance(a),
                        hit_covariance(b)),
                0.25);

  // TOF covariance contribution along LOR direction
  new_lor.covariance =
      sym_add(new_lor.covariance,
              sym_scale(sym_proj(c_hat),
                        0.5 * SPD_LGHT * SPD_LGHT * TIME_VAR));

  if (NONCOLLINEARITY) {
    new_lor.covariance = sym_add(new_lor.covariance, nc_cov);
  }

  return new_lor;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [input_nema_binary] [output_lor]\n", argv[0]);
    return 1;
  }

  FILE *input_file = fopen(argv[1], "rb");
  if (!input_file) {
    perror("Could not open input file");
    return 1;
  }

  FILE *output_file = fopen(argv[2], "wb");
  if (!output_file) {
    perror("Could not open output file");
    fclose(input_file);
    return 1;
  }

  nema_lor raw_nema;
  unsigned long long count = 0;

  while (read_nema_lor(&raw_nema, input_file)) {
    nema_lor smeared = apply_hit_uncertainty(&raw_nema);
    lor out_lor = create_lor_from_nema(&smeared);

    print_lor(&out_lor, output_file);

    count++;
    if (count % 1000000000 == 0) {
        printf("Processed %llu LORs\n", count);
    }
  }
  fclose(input_file);
  fclose(output_file);

  printf("Read NEMA LORs: %llu\n", count);

  return 0;
}