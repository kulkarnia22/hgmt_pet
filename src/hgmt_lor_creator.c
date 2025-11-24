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
double eff_by_energy[COLS];
double energies_keV[RANGE_COLS];
double ranges_mm[RANGE_COLS]; 
// params
bool writing_to_lor = true;
uint vis_events = 0;
uint counter = 0;
uint hit_first_counter = 0;
uint hit_second_counter = 0;
uint hit_third_counter = 0;
uint chain_count = 0;
uint lor_count = 0;
uint first_scatter_count = 0;
uint second_scatter_count = 0;
uint third_scatter_count = 0;
// global variables
#define NUM_CUTS 6
#define NUM_DEBUG_OPTIONS 5
// cuts are {occured, interacted with something, wasn't inpatient, detected,
// first or second detected}
char *cut_descriptions[] = {
    "occured",
    "interected with something",
    "something detected",
    "wasn't inpatient",
    "first scatter detected",
    "first scatter identified (requires lor to be made)"};
uint cuts[NUM_CUTS] = {0};
// dual cuts are the same but require both to happen in an annihilation
uint num_annihilations = 0;
uint max_annihilations = 0;
uint dual_cuts[NUM_CUTS] = {0};
uint num_scatters = 0;
uint num_hits = 0;
uint first_scatters = 0;
uint first_hits = 0;
event *first_event;
bool debug_options[NUM_DEBUG_OPTIONS];
FILE *debug[NUM_DEBUG_OPTIONS];
FILE *visualization;
FILE *lor_output;
FILE *eff_table_file;
FILE *energy_and_range_file;
FILE *input_file;
FILE *first_detected;
FILE *lor_config_output;
FILE *lor_layer_decomp;
FILE *compton_check;
FILE *incoming_photons;
FILE *incoming_photons_firstscattered;
FILE *incoming_photons_first_firstid;
FILE *incoming_photons_first_secondid;
FILE *incoming_photons_first_thirdid;
FILE *incoming_photons_secondscattered;
FILE *incoming_photons_thirdscattered;
FILE *incoming_photons_fourthscattered;
FILE *incoming_photons_fifthscattered;
FILE *incoming_photons_sixthscattered;
FILE *debug_test;
FILE *first_scatter_ranges;
FILE *second_scatter_ranges;
FILE *third_scatter_ranges;
FILE *pores_crossed_first_hit;
FILE *min_energy_pores;
FILE *scatter_layers;
//FILE *angular_output;
//FILE *det_angle;
//FILE *energy_dist;
//FILE *energy_det;
//FILE *det_check;
pthread_mutex_t read_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t write_lock = PTHREAD_MUTEX_INITIALIZER;
primitive_lor create_prim_lor(hit_split split) {
  primitive_lor prim_lor;
  // hit *hit1 = initial_by_neural_network(split.hits1, split.num_hits1);
  // hit *hit2 = initial_by_neural_network(split.hits2, split.num_hits2);
  // hit *hit1 = initial_by_truth(split.hits1, split.num_hits1);
  // hit *hit2 = initial_by_truth(split.hits2, split.num_hits2);
  // hit *hit1 = initial_by_least_time(split.hits1, split.num_hits1);
  // hit *hit2 = initial_by_least_time(split.hits2, split.num_hits2);
  //hit *hit1 = initial_by_truth(split.hits1, split.num_hits1);
  //hit *hit2 = initial_by_truth(split.hits2, split.num_hits2);
  hit *hit1 = initial_by_least_radial(split.hits1, split.num_hits1);
  hit *hit2 = initial_by_least_radial(split.hits2, split.num_hits2);
  prim_lor.hit1 = *hit1;
  prim_lor.hit2 = *hit2;
  return prim_lor;
}

double impact_parameter(vec3d loc1, vec3d loc2, double tof1, double tof2,
                        vec3d true_center) {
  vec3d c = vec_sub(loc1, loc2);
  vec3d geometric_center = vec_add(loc2, vec_scale(c, 0.5));
  vec3d c_hat = vec_norm(c);
  double delta_t = tof2 - tof1;
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d estimated_loc = vec_add(geometric_center, displacement_from_center);
  return vec_mag(vec_rejection(vec_sub(estimated_loc, true_center), c));
}

/*void angle_from_gamma(annihilation annihil){

  printing to angular output file
  need to figure out angle between outgoing electron and incoming gamma
  do I have information for the direction of the incoming gamma? Yes, I can 
  find this out from the ground truth lor. If I have position of event[0] from
  first photon and event[0] from second photon, then I can draw a direction vector 
  from one to the other. Photon1 direction vector will be give by photon1_pos - photon2_pos
  and vice versa for Photon2 direction vector. In that case this should be done in
  debug_annihilation.
  print_double(path->events[0]->direction, );

    photon_path path1 = annihil.photon1;
    photon_path path2 = annihil.photon2;

    vec3d p1 = path1.events[0]->position; 
    vec3d p2 = path2.events[0]->position;
    vec3d first_dir = vec_sub(p1,p2);
    vec3d second_dir = vec_sub(p2,p1);

    double p1_angle = vec_angle(first_dir, path1.events[0]->direction);
    double p2_angle = vec_angle(second_dir, path2.events[0]->direction);

    if (path1.events[0]->detected){
      print_double(p1_angle, det_angle);
    }

    if (path2.events[0]->detected){
      print_double(p2_angle, det_angle);
    }

    print_double(p1_angle, angular_output);
    print_double(p2_angle, angular_output);
}*/

// gets the detector an event happened in. return -1 if it didn't happen in a
// detector
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
lor create_lor(primitive_lor *prim_lor) {

  vec3d a = prim_lor->hit1.position;
  vec3d b = prim_lor->hit2.position;
  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_scale(vec_add(a, b), 0.5);
  vec3d c_hat = vec_norm(c);
  double delta_t = -(prim_lor->hit1.tof - prim_lor->hit2.tof);
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d annihilation_loc = vec_add(geometric_center, displacement_from_center);

  lor new_lor;
  new_lor.center = annihilation_loc;
  new_lor.covariance =
      sym_scale(sym_add(hit_covariance(prim_lor->hit1.position),
                        hit_covariance(prim_lor->hit2.position)),
                0.25);
  new_lor.covariance =
      sym_add(new_lor.covariance,
              sym_scale(sym_proj(c_hat), 0.5 * SPD_LGHT * SPD_LGHT * TIME_VAR));
  new_lor.covariance = sym_add(new_lor.covariance, sym_id(DIFFUSION_VARIANCE));
  return new_lor;
}

// provide debug statistics
int debug_path(photon_path *path, debug_context context) {
  if (path->num_events == 0)
    return 0;

  if (path->num_hits != 0){
    chain_count ++;
    //here I will do all the work for counting # of pore crossings
    hit *first_hit = path->hits[0];
    event *single_event = first_hit->source;
    int num_pores_crossed = pores_crossed(single_event);
    double energy = min_energy(single_event, num_pores_crossed);
    
    fwrite(&num_pores_crossed, sizeof(int), 1, pores_crossed_first_hit);
    fwrite(&num_pores_crossed, sizeof(int), 1, min_energy_pores);
    fwrite(&energy, sizeof(double), 1, min_energy_pores); 
  }

  // getting all the important statistics
  if (debug_options[1]){
    if(path->num_events >= 1){
        first_scatter_count ++;
        double range1 = pg_kapton_range_mm(path->events[0]->energy);
        fwrite(&range1, sizeof(double), 1, first_scatter_ranges);
    }
    if(path->num_events >= 2){
        second_scatter_count ++;
        double range2 = pg_kapton_range_mm(path->events[1]->energy);
        fwrite(&range2, sizeof(double), 1, second_scatter_ranges);
    }
    if(path->num_events >= 3){
        third_scatter_count ++;
        double range3 = pg_kapton_range_mm(path->events[2]->energy);
        fwrite(&range3, sizeof(double), 1, third_scatter_ranges);
    }
    print_double(path->events[0]->detector_id, debug[1]);
    first_scatters++;
    if (path->events[0]->detected && first_detected != NULL){
        first_hits++;
        print_double(path->events[0]->detector_id, first_detected);  // Detected only
    }
  }

  //printing to energy distribution file
  /*double energy_check = path->events[0]->energy;
  print_double(energy_check, energy_dist);
  if (path->events[0]->detected){
    print_double(energy_check, energy_det);
  }*/

 //figure out which cut the photon got to, format is: if (not cut n) cut=n-1
  int cut;
  if (path->num_hits == 0)
    cut = 1;
  else if (path->events[0]->detector_id == -1)
    cut = 2;
  else if (!path->events[0]->detected)
    cut = 3;
  else if (context.prim_lor == NULL ||
           (path->events[0] != context.prim_lor->hit1.source &&
            path->events[0] != context.prim_lor->hit2.source))
    cut = 4;
  else
    cut = 5;
  return cut;
}

int debug_annihilation(debug_context context) {
  annihilation *annihil = context.annihil;
  num_scatters += annihil->num_events;
  num_hits += annihil->num_hits;
  for(int i = 0; i < annihil->num_events; i ++){
    print_double(annihil->events[i].detector_id, scatter_layers);
  }
  if (debug_options[0])
    for (int j = 0; j < annihil->num_events; j++)
      print_double(annihil->events[j].detector_id, debug[0]);

  if (vis_events > 0) {
    print_history(context, visualization);
    vis_events--;
  }
  int cut1 = debug_path(&annihil->photon1, context);
  int cut2 = debug_path(&annihil->photon2, context);

  int cut = MIN(cut1, cut2);
  cuts[cut1]++;
  cuts[cut2]++;
  dual_cuts[cut]++;


  if (debug_options[4] && cut >= 1) {
    for (int i = 0; i < MIN(annihil->photon1.num_events, 4); i++)
      for (int j = 0; j < MIN(annihil->photon2.num_events, 4); j++) {
        vec3d true_center = annihil->center;
        vec3d loc1 = annihil->photon1.events[i]->position;
        vec3d loc2 = annihil->photon2.events[j]->position;
        double tof1 = annihil->photon1.events[i]->tof;
        double tof2 = annihil->photon2.events[j]->tof;
        if (i < j) {
          print_int(i + 1, debug[4]);
          print_int(j + 1, debug[4]);
        } else {
          print_int(j + 1, debug[4]);
          print_int(i + 1, debug[4]);
        }
        print_double(impact_parameter(loc1, loc2, tof1, tof2, true_center),
                     debug[4]);
      }
    int separator = -1;
    fwrite(&separator, sizeof(int), 1, debug[4]);
  }

  return cut;
}
void debug_lor(debug_context context) {
  if (debug_options[2] || debug_options[3]) {
    lor *new_lor = context.lor;
    vec3d truecenter = context.annihil->center;
    vec3d dir = vec_norm(vec_sub(context.prim_lor->hit1.position,
                                 context.prim_lor->hit2.position));
    if (debug_options[2])
      print_double(
          vec_mag(vec_rejection(vec_sub(new_lor->center, truecenter), dir)),
          debug[2]);
    if (debug_options[3])
      print_double(
          vec_mag(vec_projection(vec_sub(new_lor->center, truecenter), dir)),
          debug[3]);
  }

}
void debug_prim_lor(debug_context context) {}
void debug_all(debug_context context) {
  debug_annihilation(context);
  num_annihilations++;
  printm(num_annihilations, 1000000);
  if (context.prim_lor != NULL) {
    debug_prim_lor(context);
    debug_lor(context);
  }
}
void *worker(void *arg) {
  annihilation annihil;
  while (true) {
    pthread_mutex_lock(&read_lock);
    bool worked =
        read_annihilation(&annihil, input_file, eff_by_energy) &&
        (max_annihilations == 0 || num_annihilations < max_annihilations);
    pthread_mutex_unlock(&read_lock);

    if (!worked)
      return NULL;
    //hit_split split = create_hit_split(annihil.hits, annihil.num_hits);
    hit_split split;
    split.num_hits1 = annihil.photon1.num_hits;
    split.num_hits2 = annihil.photon2.num_hits;
    split.hits1 = annihil.photon1.hits;
    split.hits2 = annihil.photon2.hits;
    for (int i = 0; i < split.num_hits1; i ++){
        hit *current_hit = split.hits1[i];
        int scatter_num = current_hit->source->number;
        if (scatter_num == 0){
            hit_first_counter ++;
        } 
        else if(scatter_num == 1){
            hit_second_counter ++;
        }
        else if(scatter_num == 2){
            hit_third_counter ++;
        }
    }
    for (int i = 0; i < split.num_hits2; i ++){
        hit *current_hit = split.hits2[i];
        int scatter_num = current_hit->source->number;
        if (scatter_num == 0){
            hit_first_counter ++;
        } 
        else if(scatter_num == 1){
            hit_second_counter ++;
        }
        else if(scatter_num == 2){
            hit_third_counter ++;
        }
    }
    /*if (annihil.photon1.num_events > 0 &&
    annihil.photon1.events != NULL &&
    annihil.photon1.events[0] != NULL &&
    annihil.photon2.num_events > 0 &&
    annihil.photon2.events != NULL &&
    annihil.photon2.events[0] != NULL){
        angle_from_gamma(annihil);
    }*/
    

    debug_context context = {0};
    context.annihil = &annihil;
    context.split = &split;
    primitive_lor prim_lor;
    lor new_lor;
    if (split.num_hits1 >= 1 && split.num_hits2 >= 1) { 
      lor_count ++;
      prim_lor = create_prim_lor(split);
      new_lor = create_lor(&prim_lor);
      context.prim_lor = &prim_lor;
      context.lor = &new_lor;
    }
    pthread_mutex_lock(&write_lock);
    debug_all(context);
    if (context.prim_lor != NULL) {
      double hit1_photon = context.prim_lor->hit1.source->parent_gamma_energy;
      double hit2_photon = context.prim_lor->hit2.source->parent_gamma_energy;
      //printf("double: %f/n", hit1_photon);
      if (hit1_photon > 0 && hit2_photon > 0){
        //printf("here");
        fwrite(&hit1_photon, sizeof(double), 1, incoming_photons);
        fwrite(&hit2_photon, sizeof(double), 1, incoming_photons); 
        //code for hit1
        if (context.prim_lor->hit1.source->number == 0){
            fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_firstscattered);
            if (context.prim_lor->hit1.source->detector_id == 1){
                fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_first_firstid);
            }
            if (context.prim_lor->hit1.source->detector_id == 2){
                fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_first_secondid);
                
            }
            if (context.prim_lor->hit1.source->detector_id == 3){
                fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_first_thirdid);
                
            }
        }
        if (context.prim_lor->hit1.source->number == 1){
            fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_secondscattered);
        }
        if (context.prim_lor->hit1.source->number == 2){
            fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_thirdscattered);
        }
        if (context.prim_lor->hit1.source->number == 3){
            fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_fourthscattered);
        }
        if (context.prim_lor->hit1.source->number == 4){
            fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_fifthscattered);
        }
        if (context.prim_lor->hit1.source->number == 5){
            fwrite(&hit1_photon, sizeof(double), 1, incoming_photons_sixthscattered);
        }
        //code for hit2
        if (context.prim_lor->hit2.source->number == 0){
            fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_firstscattered);
            if (context.prim_lor->hit2.source->detector_id == 1){
                fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_first_firstid);
            }
            if (context.prim_lor->hit2.source->detector_id == 2){
                fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_first_secondid);
                
            }
            if (context.prim_lor->hit2.source->detector_id == 3){
                fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_first_thirdid);
                
            }
        }
        if (context.prim_lor->hit2.source->number == 1){
            fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_secondscattered);
        }
        if (context.prim_lor->hit2.source->number == 2){
            fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_thirdscattered);
        }
        if (context.prim_lor->hit2.source->number == 3){
            fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_fourthscattered);
        }
        if (context.prim_lor->hit2.source->number == 4){
            fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_fifthscattered);
        }
        if (context.prim_lor->hit2.source->number == 5){
            fwrite(&hit2_photon, sizeof(double), 1, incoming_photons_sixthscattered);
        }
      }
      if (writing_to_lor){
        //writing configs and lor structs
        uint first_num = context.prim_lor->hit1.source->number;
        uint second_num = context.prim_lor->hit2.source->number;
        uint first_id = context.prim_lor->hit1.source->detector_id;
        uint second_id = context.prim_lor->hit2.source->detector_id;
        //printf("check = %i\n", second_id);
        fwrite(&first_id, sizeof(int), 1, lor_layer_decomp);
        fwrite(&second_id, sizeof(int), 1, lor_layer_decomp);
        fwrite(&first_num, sizeof(int), 1, lor_config_output);
        fwrite(&second_num, sizeof(int), 1, lor_config_output);
        print_lor(&new_lor, lor_output);
      }
    }
    pthread_mutex_unlock(&write_lock);
    free_annihilation(&annihil);
    //free_hit_split(&split);
  }
}
int main(int argc, char **argv) { 
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // handling all flags and arguments
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_lor_creator [TOPAS_file_position (not the phsp)] "
             "[efficiency_table_position.csv] [output_directory]\n");
      printf("-d: run in debug mode, do not write to lor file\n");
      printf("-v#: visualize # events\n");
      printf("-m#: only read a maximum of # annihilations\n");
      printf("-e#: run with debug option #\n");
      printf("\t0: histogram of detector vs number of scatters\n");
      printf("\t1: histogram of detector vs number of first scatters\n");
      printf("\t2: lor reconstruction error to real center (transverse)\n");
      printf("\t3: lor reconstruction error to real center (longitudinal)\n");
      printf("\t4: scatter number vs reconstruction error truth study\n");
      exit(0);
    } else if (strcmp(flags[i], "-d") == 0) {
      printf("running in debug mode, won't write to a lor file\n");
      writing_to_lor = false;
    } else if (strncmp(flags[i], "-e", 2) == 0) {
      uint debug_option;
      sscanf(flags[i], "-e%u", &debug_option);
      debug_options[debug_option] = true;
    } else if (strncmp(flags[i], "-m", 2) == 0) {
      sscanf(flags[i], "-m%u", &max_annihilations);
    } else if (strncmp(flags[i], "-v", 2) == 0) {
      sscanf(flags[i], "-v%u", &vis_events);
      printf("outputting data to visualize %u events\n", vis_events);
    }
  }
  
  //macros check
  //SHOW_DEFINE(TAU);
  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, three arguments required.\n");
    printf("Use the -h command to get options.\n");
    exit(1);
  }
  // opens file for debgug test
  char *test_filename;
  asprintf(&test_filename, "%sdebug_test.data", args[2]);
  debug_test = fopen(test_filename, "wb");
  free(test_filename);

  // opens files for debug output
  printf("running with debug options:"); // Output: 42
  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i]) {
      printf(" %i", i);
      char *filename;
      asprintf(&filename, "%sdebug%d.data", args[2], i);
      debug[i] = fopen(filename, "wb");
      free(filename);
    }

  // opens scatter file for detected first scatters
  if (debug_options[1]) {
    char *fd_filename;
    asprintf(&fd_filename, "%sfirst_detected.data", args[2]);
    first_detected = fopen(fd_filename, "wb");
    free(fd_filename);
  }
  printf("\n");

  // reads in efficiency table into 2D array called eff_by_ang
  /* Above comment is probably old. I believe lines bellow reads in
  a 1D array called eff_by_energy */
  printf("HGMT LOR Creator\n\nLoading in '%s' as efficiencies table...\n",
         args[1]);
  eff_table_file = fopen(args[1], "r");
  read_eff(eff_table_file, eff_by_energy);
  fclose(eff_table_file);

  //reading in energy,range vals in kapton
  energy_and_range_file = fopen("simulation_materials/kapton_ranges.csv", "r");
  read_ranges_and_energies(energy_and_range_file, ranges_mm, energies_keV);
  fclose(energy_and_range_file);
  pg_set_kapton_range_table(RANGE_COLS, energies_keV, ranges_mm);


  //initializing detector layer specs
  for(int i = 0; i < PORE_NUM_LAYERS; i++){
    LatticeSpec L;
    pg_init_layer_from_dims(&L, PG_LAYER_RADII[i].Rin, PG_LAYER_RADII[i].Rout,\
    0.0, 0.0, -1000.0, 1000.0);
    specs[i] = L;
  }

  /*for (int i = 0; i < COLS; i++){
    printf("eff_by_energy[%d] = %.3f\n", i, eff_by_energy[i]);
  }*/

  //making angle info files
  /*char *angular_file_loc;
  asprintf(&angular_file_loc, "%sangular_ouput.data", args[2]);
  angular_output = fopen(angular_file_loc, "wb");
  free(angular_file_loc);

  char *det_angle_loc;
  asprintf(&det_angle_loc, "%sdet_angle_output.data", args[2]);
  det_angle = fopen(det_angle_loc, "wb");
  free(det_angle_loc);

  //making energy files
  char *energy_output_loc;
  asprintf(&energy_output_loc, "%senergy_dist.data", args[2]);
  energy_dist = fopen(energy_output_loc, "wb");
  free(energy_output_loc);

  char *energy_det_loc;
  asprintf(&energy_det_loc, "%senergy_det.data", args[2]);
  energy_det = fopen(energy_det_loc, "wb");
  free(energy_det_loc);*/

  // opens up a .lor file to output each LOR into
  // can change HGMTDerenzo.lor to HGMTPoint.lor and vice versa

  
  //opening up file for pores crossed
  char *pores_crossed_loc;
  asprintf(&pores_crossed_loc, "%snum_pores_crossed.data", args[2]);
  pores_crossed_first_hit = fopen(pores_crossed_loc, "wb");
  free(pores_crossed_loc);

  //open up file for min energy
  char *min_energy_loc;
  asprintf(&min_energy_loc, "%smin_energy.data", args[2]);
  min_energy_pores = fopen(min_energy_loc, "wb");
  free(min_energy_loc);

  //open file for scatter layer counts
  char *scatter_layer_loc;
  asprintf(&scatter_layer_loc, "%sscatter_layers.data", args[2]);
  scatter_layers = fopen(scatter_layer_loc, "wb");
  free(scatter_layer_loc);

  //opens scatter ranges files
  char *first_ranges_loc;
  asprintf(&first_ranges_loc, "%sfirst_scatter_ranges.data", args[2]);
  first_scatter_ranges = fopen(first_ranges_loc, "wb");
  free(first_ranges_loc);

  char *second_ranges_loc;
  asprintf(&second_ranges_loc, "%ssecond_scatter_ranges.data", args[2]);
  second_scatter_ranges = fopen(second_ranges_loc, "wb");
  free(second_ranges_loc);

  char *third_ranges_loc;
  asprintf(&third_ranges_loc, "%sthird_scatter_ranges.data", args[2]);
  third_scatter_ranges = fopen(third_ranges_loc, "wb");
  free(third_ranges_loc);

  //opens up an incoming photon .data output file
  char *incoming_photon_file;
  asprintf(&incoming_photon_file, "%sIncoming_Photons.data", args[2]);
  incoming_photons = fopen(incoming_photon_file, "wb");
  free(incoming_photon_file);

  char *incoming_photon_first_file;
  asprintf(&incoming_photon_first_file, "%sIncoming_Photons_First.data", args[2]);
  incoming_photons_firstscattered = fopen(incoming_photon_first_file, "wb");
  free(incoming_photon_first_file);

  char *incoming_photon_first_firstid_file;
  asprintf(&incoming_photon_first_firstid_file, "%sIncoming_Photons_Firstid.data", args[2]);
  incoming_photons_first_firstid = fopen(incoming_photon_first_firstid_file, "wb");
  free(incoming_photon_first_firstid_file);

  char *incoming_photon_first_secondid_file;
  asprintf(&incoming_photon_first_secondid_file, "%sIncoming_Photons_Secondid.data", args[2]);
  incoming_photons_first_secondid = fopen(incoming_photon_first_secondid_file, "wb");
  free(incoming_photon_first_secondid_file);

  char *incoming_photon_first_thirdid_file;
  asprintf(&incoming_photon_first_thirdid_file, "%sIncoming_Photons_Thirdid.data", args[2]);
  incoming_photons_first_thirdid = fopen(incoming_photon_first_thirdid_file, "wb");
  free(incoming_photon_first_thirdid_file);

  char *incoming_photon_second;
  asprintf(&incoming_photon_second, "%sIncoming_Photons_Second.data", args[2]);
  incoming_photons_secondscattered = fopen(incoming_photon_second, "wb");
  free(incoming_photon_second);

  char *incoming_photon_third;
  asprintf(&incoming_photon_third, "%sIncoming_Photons_Third.data", args[2]);
  incoming_photons_thirdscattered = fopen(incoming_photon_third, "wb");
  free(incoming_photon_third);

  char *incoming_photon_fourth;
  asprintf(&incoming_photon_fourth, "%sIncoming_Photons_Fourth.data", args[2]);
  incoming_photons_fourthscattered = fopen(incoming_photon_fourth, "wb");
  free(incoming_photon_fourth);

  char *incoming_photon_fifth;
  asprintf(&incoming_photon_fifth, "%sIncoming_Photons_Fifth.data", args[2]);
  incoming_photons_fifthscattered = fopen(incoming_photon_fifth, "wb");
  free(incoming_photon_fifth);

  char *incoming_photon_sixth;
  asprintf(&incoming_photon_sixth, "%sIncoming_Photons_Sixth.data", args[2]);
  incoming_photons_sixthscattered = fopen(incoming_photon_sixth, "wb");
  free(incoming_photon_sixth);

  if (writing_to_lor) {
    char *lor_file_loc;
    asprintf(&lor_file_loc, "%sHGMTPointVac.lor", args[2]);
    lor_output = fopen(lor_file_loc, "wb");
    free(lor_file_loc);
  }

  //opens lor configuration file
  char *lor_config_loc;
  asprintf(&lor_config_loc, "%slor_config.data", args[2]);
  lor_config_output = fopen(lor_config_loc, "wb");
  free(lor_config_loc);

  char *lor_layer_decomp_loc;
  asprintf(&lor_layer_decomp_loc, "%slor_layer.data", args[2]);
  lor_layer_decomp = fopen(lor_layer_decomp_loc, "wb");
  free(lor_layer_decomp_loc);

  //opens compton check file
  char *compton_check_loc;
  asprintf(&compton_check_loc, "%scompton_check.data", args[2]);
  compton_check = fopen(compton_check_loc, "wb");
  free(compton_check_loc);

  if (vis_events) {
    char *filename;
    asprintf(&filename, "%svisualization.data", args[2]);
    visualization = fopen(filename, "w");
    free(filename);
  }
  input_file = fopen(args[0], "rb");
  printf("Loading in '%s' as the annihilations...\n", args[0]);

  printf("Constructing the hits...\n\n");
  // everythings happens here, process and read the annihilations
  pthread_t threads[THREADS];
  for (int i = 0; i < THREADS; i++)
    pthread_create(&threads[i], NULL, worker, NULL);
  for (int i = 0; i < THREADS; i++)
    pthread_join(threads[i], NULL);
  printf("\n");
  // fixing cuts formating to be cumulative
  for (int i = NUM_CUTS - 2; i >= 0; i--) {
    cuts[i] += cuts[i + 1];
    dual_cuts[i] += dual_cuts[i + 1];
  }
  printf("total annihilations: %u\n", num_annihilations);
  printf("total scatters: %u\n", num_scatters);
  printf("total hits: %u\n", num_hits);
  printf("total first scatters: %u\n", first_scatters);
  printf("total first hits: %u\n", first_hits);
  printf("total chains: %u\n", chain_count);
  printf("total number of lors: %u\n", lor_count);
  printf("first scatter hits: %u\n", hit_first_counter);
  printf("first scatter count: %u\n", first_scatter_count);
  printf("hits/first scatter: %d\n", hit_first_counter/first_scatter_count);
  printf("hits/second scatter: %d\n", hit_second_counter/second_scatter_count);
  printf("hits/third scatters: %d\n", hit_third_counter/third_scatter_count);
  printf(
      "(DUAL)CUT 'N': 'num' 'percent passing' 'cumulative percent passing'\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("%u: %s\n", i, cut_descriptions[i]);
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("CUT %u: %u %lf %lf\n", i, cuts[i], (double)cuts[i] / cuts[i - 1],
           (double)cuts[i] / cuts[0]);
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("DUALCUT %u: %u %lf %lf\n", i, dual_cuts[i],
           (double)dual_cuts[i] / dual_cuts[i - 1],
           (double)dual_cuts[i] / dual_cuts[0]);
  // closing stuff out
  pthread_mutex_destroy(&read_lock);
  pthread_mutex_destroy(&write_lock);
  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i]){
      fclose(debug[i]);
    }
  if (first_detected != NULL){
    fclose(first_detected);
  }
  if (visualization != NULL){
    printf("Closing visualization file.\n");
    fclose(visualization);
  }
  if (lor_output != NULL){
    fclose(lor_output);
  }
  if (lor_config_output != NULL){
    fclose(lor_config_output);
  }
  if (lor_layer_decomp != NULL){
    fclose(lor_layer_decomp);
  }
  if(compton_check != NULL){
    fclose(compton_check);
  }
  if (incoming_photons != NULL){
    fclose(incoming_photons);
  } 
  if(incoming_photons_firstscattered != NULL){
    fclose(incoming_photons_firstscattered);
  } 
  if(incoming_photons_first_firstid != NULL){
    fclose(incoming_photons_first_firstid);
  }
  if(incoming_photons_first_secondid != NULL){
    fclose(incoming_photons_first_secondid);
  }
  if(incoming_photons_first_thirdid != NULL){
    fclose(incoming_photons_first_thirdid);
  }
  if(incoming_photons_secondscattered != NULL){
    fclose(incoming_photons_secondscattered);
  }
  if(incoming_photons_thirdscattered!=NULL){
    fclose(incoming_photons_thirdscattered);
  }
  if(incoming_photons_fourthscattered!=NULL){
    fclose(incoming_photons_fourthscattered);
  }
  if(incoming_photons_fifthscattered!=NULL){
    fclose(incoming_photons_fifthscattered);
  }
  if(incoming_photons_sixthscattered!=NULL){
    fclose(incoming_photons_sixthscattered);
  }
  if(first_scatter_ranges != NULL){
    fclose(first_scatter_ranges);
  }
  if(second_scatter_ranges != NULL){
    fclose(second_scatter_ranges);
  }
  if(third_scatter_ranges != NULL){
    fclose(third_scatter_ranges);
  }
  if(pores_crossed_first_hit != NULL){
    fclose(pores_crossed_first_hit);
  }
  if(min_energy_pores != NULL){
    fclose(min_energy_pores);
  }
  if(scatter_layers != NULL){
    fclose(scatter_layers);
  }
  /*if (angular_output != NULL){
    fclose(angular_output);
  }
  if (det_angle != NULL){
    fclose(det_angle);
  }
  if (energy_dist != NULL){
    fclose(energy_dist);
  }
  if (energy_det != NULL){
    fclose(energy_det);
  }*/
  return 0;
}
