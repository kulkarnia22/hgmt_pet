#include "read_write.h"
#include "hgmt_structs.h"
#include "linear_algebra.h"

void read_event(event *new_event, FILE *source, double eff_by_energy[COLS]) {
  fread(&new_event->tof, sizeof(double), 1, source);
  fread(&new_event->energy, sizeof(double), 1, source);
  read_vec(&new_event->position, source);
  read_vec(&new_event->direction, source);
  fread(&new_event->primary, sizeof(uint), 1, source);
  new_event->detector_id = get_detector(new_event->position);
  //new_event->detected = is_detected(new_event, eff_by_energy);
  //new_event->detected = plane_crossing(new_event);
  new_event->detected = plane_crossingv2(new_event);
  new_event->parent_gamma_energy = -1; //placeholder
  new_event->parent_gamma_dir = (vec3d){0,0,0};
  //printf("scatter direction = (%f, %f, %f)\n", new_event->direction.x, new_event->direction.y, new_event->direction.z);
  //printf("scatter energy = %f\n", new_event->energy);
  //printf("detector layer = %d\n", new_event->detector_id);
  //printf("gamma = %d\n", new_event->primary);
  //printf("detected = %d\n\n", new_event->detected);
}

bool read_annihilation(annihilation *annihil, FILE *source,
                       double eff_by_energy[COLS]) {
  *annihil = (annihilation){0};
  if (!fread(&annihil->time, sizeof(double), 1, source))
    return false;
  read_vec(&annihil->origin, source);
  read_vec(&annihil->center, source);
  fread(&annihil->num_events, sizeof(uint), 1, source);
  uint num_primary1 = 0;
  uint num_primary2 = 0;
  annihil->events = (event *)malloc(sizeof(event) * annihil->num_events);
  uint num_hits = 0;
  uint num_primary1_hits = 0;
  uint num_primary2_hits = 0;
  for (int i = 0; i < annihil->num_events; i++) {
    event *next_event = &annihil->events[i];
    read_event(next_event, source, eff_by_energy);

    if (next_event->primary == 1) {
      next_event->number = num_primary1;
      //printf("scatter number = %d\n\n", next_event->number);
      num_primary1++;
    } else if (next_event->primary == 2) {
      next_event->number = num_primary2;
      //printf("scatter number = %d\n\n", next_event->number);
      num_primary2++;
    }
    if (next_event->detected) {
      num_hits++;
      if (next_event->primary == 1)
        num_primary1_hits++;
      else if (next_event->primary == 2)
        num_primary2_hits++;
    }
  }
  annihil->hits = (hit *)malloc(sizeof(hit) * num_hits);
  annihil->photon1.events = (event **)malloc(sizeof(event *) * num_primary1);
  annihil->photon2.events = (event **)malloc(sizeof(event *) * num_primary2);
  annihil->photon1.hits = (hit **)malloc(sizeof(hit *) * num_primary1_hits);
  annihil->photon2.hits = (hit **)malloc(sizeof(hit *) * num_primary2_hits);
  //I'll try to reconstruct parent gamma energies and directions here
  for (int i = 0; i < annihil->num_events; i++) {
    event *current_event = &annihil->events[i];
    if (current_event->primary == 1) {
      annihil->photon1.events[annihil->photon1.num_events] = current_event;
      annihil->photon1.num_events++;
      //first scatter case of compton reconstruction
      if(annihil->photon1.num_events == 1){
        //setting parent_gamma_energy
        current_event->parent_gamma_energy = 511.0;
        //getting direction vector of incoming/parent gamma
        vec3d parent_gamma_dir = vec_norm(vec_sub(current_event->position, annihil->center));
        current_event->parent_gamma_dir = parent_gamma_dir;
      }
      else{
        //getting info from previous scatter
        event *prev_event = annihil->photon1.events[annihil->photon1.num_events - 2];
        //setting parent_gamma_energy
        current_event->parent_gamma_energy = prev_event->parent_gamma_energy - prev_event->energy;
        //getting direction vector of incoming/parent gamma
        vec3d prev_parent_gamma_dir = prev_event->parent_gamma_dir;
        double prev_parent_gamma_energy = prev_event->parent_gamma_energy;
        vec3d prev_electron_dir = prev_event->direction;
        double prev_electron_kin = prev_event->energy;
        //defined ME_C2_KEV as electron rest energy in macros. Treat speed of light as 1
        vec3d p_gamma_i = vec_scale(prev_parent_gamma_dir, prev_parent_gamma_energy);
        double p_electron_mag = sqrt(prev_electron_kin*(prev_electron_kin + 2*ME_C2_KEV));
        vec3d p_electron_f = vec_scale(prev_electron_dir, p_electron_mag);
        vec3d outgoing_gamma_dir = vec_norm(vec_sub(p_gamma_i, p_electron_f));
        current_event->parent_gamma_dir = outgoing_gamma_dir;
      }
    } else if (current_event->primary == 2) {
      annihil->photon2.events[annihil->photon2.num_events] = current_event;
      annihil->photon2.num_events++;
      if(annihil->photon2.num_events == 1){
        //setting parent_gamma_energy
        current_event->parent_gamma_energy = 511.0;
        //getting direction vector of incoming/parent gamma
        vec3d parent_gamma_dir = vec_norm(vec_sub(current_event->position, annihil->center));
        current_event->parent_gamma_dir = parent_gamma_dir;
      }
      else{
        //getting info from previous scatter
        event *prev_event = annihil->photon2.events[annihil->photon2.num_events - 2];
        //setting parent_gamma_energy
        current_event->parent_gamma_energy = prev_event->parent_gamma_energy - prev_event->energy;
        //getting direction vector of incoming/parent gamma
        vec3d prev_parent_gamma_dir = prev_event->parent_gamma_dir;
        double prev_parent_gamma_energy = prev_event->parent_gamma_energy;
        vec3d prev_electron_dir = prev_event->direction;
        double prev_electron_kin = prev_event->energy;
        //defined ME_C2_KEV as electron rest energy in macros. Treat speed of light as 1
        vec3d p_gamma_i = vec_scale(prev_parent_gamma_dir, prev_parent_gamma_energy);
        double p_electron_mag = sqrt(prev_electron_kin*(prev_electron_kin + 2*ME_C2_KEV));
        vec3d p_electron_f = vec_scale(prev_electron_dir, p_electron_mag);
        vec3d outgoing_gamma_dir = vec_norm(vec_sub(p_gamma_i, p_electron_f));
        current_event->parent_gamma_dir = outgoing_gamma_dir;
      }
    }
    if (current_event->detected) {
      annihil->hits[annihil->num_hits] = event_to_hit(current_event);
      hit *current_hit = &annihil->hits[annihil->num_hits];
      annihil->num_hits++;
      if (current_event->primary == 1) {
        annihil->photon1.hits[annihil->photon1.num_hits] = current_hit;
        annihil->photon1.num_hits++;
      } else if (current_event->primary == 2) {
        annihil->photon2.hits[annihil->photon2.num_hits] = current_hit;
        annihil->photon2.num_hits++;
      }
    }
  }
  return true;
}
void print_lor(lor *new_lor, FILE *output) {
  print_vec(new_lor->center, output);
  print_sym_matrix(&new_lor->covariance, output);
}
void read_eff(FILE *source, double eff_by_energy[COLS]) {
  // loops through all the entries in a row
  for (int i = 0; i < COLS; i++) {
    fscanf(source, "%lf,", &eff_by_energy[i]);
  }
}

void read_ranges_and_energies(FILE *source,
                              double ranges_mm[RANGE_COLS],
                              double energies_keV[RANGE_COLS])
{
    // Skip first 8 lines (header)
    for (int i = 0; i < 8; i++) {
        int ch;
        while ((ch = fgetc(source)) != EOF && ch != '\n') { /* skip to EOL */ }
        if (ch == EOF) return; // file shorter than header
    }

    // Read RANGE_COLS data rows: Energy(MeV)  Range(g/cm^2)
    for (int i = 0; i < RANGE_COLS; i++) {
        double E_MeV, R_gcm2;
        if (fscanf(source, " %lf %lf", &E_MeV, &R_gcm2) != 2) break;

        // Convert units inline
        energies_keV[i] = 1000.0 * E_MeV;                         // MeV → keV
        ranges_mm[i]     = (R_gcm2 / KAPTON_RHO_G_CM3) * 10.0;    // g/cm² → mm
    }
}

void print_double(double numb, FILE *output) {
  fwrite(&numb, sizeof(double), 1, output);
}
void print_int(int numb, FILE *output) {
  fwrite(&numb, sizeof(int), 1, output);
}
void print_vec(vec3d vec, FILE *output) {
  fwrite(&vec.x, sizeof(double), 1, output);
  fwrite(&vec.y, sizeof(double), 1, output);
  fwrite(&vec.z, sizeof(double), 1, output);
}
void print_sym_matrix(sym_matrix *mat, FILE *output) {
  fwrite(&mat->xx, sizeof(double), 1, output);
  fwrite(&mat->xy, sizeof(double), 1, output);
  fwrite(&mat->xz, sizeof(double), 1, output);
  fwrite(&mat->yy, sizeof(double), 1, output);
  fwrite(&mat->yz, sizeof(double), 1, output);
  fwrite(&mat->zz, sizeof(double), 1, output);
}
void print_lower_matrix(lower_matrix *mat, FILE *output) {
  fwrite(&mat->l11, sizeof(double), 1, output);
  fwrite(&mat->l21, sizeof(double), 1, output);
  fwrite(&mat->l22, sizeof(double), 1, output);
  fwrite(&mat->l31, sizeof(double), 1, output);
  fwrite(&mat->l32, sizeof(double), 1, output);
  fwrite(&mat->l33, sizeof(double), 1, output);
}
bool read_vec(vec3d *vec, FILE *source) {
  bool worked = fread(&vec->x, sizeof(double), 1, source);
  fread(&vec->y, sizeof(double), 1, source);
  fread(&vec->z, sizeof(double), 1, source);
  return worked;
}
bool read_sym(sym_matrix *mat, FILE *source) {
  bool worked = fread(&mat->xx, sizeof(double), 1, source);
  fread(&mat->xy, sizeof(double), 1, source);
  fread(&mat->xz, sizeof(double), 1, source);
  fread(&mat->yy, sizeof(double), 1, source);
  fread(&mat->yz, sizeof(double), 1, source);
  fread(&mat->zz, sizeof(double), 1, source);
  return worked;
}
bool read_lor(lor *new_lor, FILE *input) {
  read_vec(&new_lor->center, input);
  return read_sym(&new_lor->covariance, input);
}

void print_prim_lor(primitive_lor *prim_lor, FILE *output) {}
void print_history(debug_context context, FILE *output) {
  print_vec(context.annihil->center, output);
  fwrite(&context.annihil->num_events, sizeof(uint), 1, output);
  for (int i = 0; i < context.annihil->num_events; i++) {
    print_vec(context.annihil->events[i].position, output);
    fwrite(&context.annihil->events[i].energy, sizeof(double), 1, output);
    fwrite(&context.annihil->events[i].primary, sizeof(uint), 1, output);
    fwrite(&context.annihil->events[i].detected, sizeof(bool), 1, output);
  }
  bool made_lor = context.prim_lor != NULL;
  fwrite(&made_lor, sizeof(bool), 1, output);
  if (made_lor) {
    print_vec(context.lor->center, output);
    lower_matrix transform = cholesky(&context.lor->covariance);
    print_lower_matrix(&transform, output);
    print_vec(context.prim_lor->hit1.position, output);
    print_vec(context.prim_lor->hit2.position, output);
  }
}
