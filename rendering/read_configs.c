#include <stdio.h>
#include <stdlib.h>
FILE *lor_input;
FILE *config_input;

FILE *lor_11;
FILE *lor_12;
FILE *lor_21;
FILE *lor_22;
FILE *missed;

int main() {
    //open input files for reading
    lor_input = fopen("data/HGMTPointVac.lor", "rb");
    config_input = fopen("data/lor_config.data", "rb");

    //open output files for writing
    lor_11 = fopen("data/HGMTPointVac11.lor", "wb");
    lor_12 = fopen("data/HGMTPointVac12.lor", "wb");
    lor_21 = fopen("data/HGMTPointVac21.lor", "wb");
    lor_22 = fopen("data/HGMTPointVac22.lor", "wb");
    missed = fopen("data/Missed.lor", "wb");

    int first, second;
    double lor_data[9];

    while(fread(&first, sizeof(int), 1, config_input) == 1 && 
        fread(&second, sizeof(int), 1, config_input) == 1 && 
        fread(lor_data, sizeof(double), 9, lor_input) == 9) {

            if (!(first == 0 && second == 0)) {
                fwrite(lor_data, sizeof(double), 9, missed);
            }
            
            if (first == 0 && second == 0) {
                fwrite(lor_data, sizeof(double), 9, lor_11);
            } else if (first == 0 && second == 1) {
                fwrite(lor_data, sizeof(double), 9, lor_12);
            } else if (first == 1 && second == 0) {
                fwrite(lor_data, sizeof(double), 9, lor_21);
            } else if (first == 1 && second == 1) {
                fwrite(lor_data, sizeof(double), 9, lor_22);
            }
        }
    fclose(lor_input);
    fclose(config_input);
    fclose(lor_11);
    fclose(lor_12);
    fclose(lor_21);
    fclose(lor_22);
    fclose(missed);

    return 0;
}