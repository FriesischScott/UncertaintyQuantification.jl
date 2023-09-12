#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{
    FILE *input = fopen(argv[1], "r");

    double x;

    fscanf(input, "%lf", &x);

    fclose(input);

    double y = pow(x, 2);

    FILE *output = fopen("out-squared.txt", "w");

    fprintf(output, "%.16f", y);

    fclose(output);

    return 0;
}