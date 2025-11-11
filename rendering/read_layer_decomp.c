#include <stdio.h>
#include <stdlib.h>
FILE *lor_input;
FILE *layer_input;

FILE *lor_layer1;
FILE *lor_layer2;
FILE *lor_layer3;
FILE *lor_layer4;
FILE *lor_layer5;
FILE *lor_layer6;

static inline int max(int a, int b) {
        return (a > b) ? a : b;
    }

int main(){
    //open input files for reading
    lor_input = fopen("data/HGMTPointVac", "rb");
    layer_input = fopen("data/lor_layer.data", "rb");

    lor_layer1 = fopen("data/layer1.lor", "wb");
    lor_layer2 = fopen("data/layer2.lor", "wb");
    lor_layer3 = fopen("data/layer3.lor", "wb");
    lor_layer4 = fopen("data/layer4.lor", "wb");
    lor_layer5 = fopen("data/layer5.lor", "wb");
    lor_layer6 = fopen("data/layer6.lor", "wb");

    int first, second;
    double lor_data[9];

    while(fread(&first, sizeof(int), 1, layer_input) == 1 && 
        fread(&second, sizeof(int), 1, layer_input) == 1 && 
        fread(lor_data, sizeof(double), 9, lor_input) == 9) {
            
        if (max(first, second) == 1) {
            //printf("here");
            fwrite(lor_data, sizeof(double), 9, lor_layer1);
        } else if (max(first, second) == 2) {
            //printf("here");
            fwrite(lor_data, sizeof(double), 9, lor_layer2);
        } else if (max(first, second) == 3) {
            //printf("here");
            fwrite(lor_data, sizeof(double), 9, lor_layer3);
        } else if (max(first, second) == 4) {
            //printf("here");
            fwrite(lor_data, sizeof(double), 9, lor_layer4);
        } else if (max(first, second) == 5) {
            //printf("here");
            fwrite(lor_data, sizeof(double), 9, lor_layer5);
        } else if (max(first, second) == 6){
            //printf("here");
            fwrite(lor_data, sizeof(double), 9, lor_layer6);
        }
    }
    fclose(lor_input);
    fclose(layer_input);
    fclose(lor_layer1);
    fclose(lor_layer2);
    fclose(lor_layer3);
    fclose(lor_layer4);
    fclose(lor_layer5);
    fclose(lor_layer6);

    return 0;
}