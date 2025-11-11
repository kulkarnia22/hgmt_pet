#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: %s [lor_file] [config_file]\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        return 1;
    }

    const char *config_name = argv[2];
    FILE *config = fopen(config_name, "rb");
    if (!config){
        perror("Error opening config file");
        return 1;
    }

    // Get file size
    fseek(file, 0L, SEEK_END);
    long filesize = ftell(file);
    fclose(file);

    fseek(config, 0L, SEEK_END);
    long configsize = ftell(config);
    fclose(config);

    // Each LOR is 72 bytes (9 doubles)
    const int bytes_per_lor = 9 * sizeof(double);
    long num_lors = filesize / bytes_per_lor;

    // Each config is 8 bytes (2 ints)
    const int bytes_per_config = 2 *sizeof(int);
    long num_configs = configsize / bytes_per_config; 

    // Total annihilations simulated
    const long total_annihilations = 1000000;

    double sensitivity = (double)num_lors / total_annihilations;

    printf("Number of LORs: %ld\n", num_lors);
    printf("Number of configs: %ld\n", num_configs);
    printf("Total annihilations: %ld\n", total_annihilations);
    printf("Sensitivity: %.6f\n", sensitivity);

    return 0;
}
