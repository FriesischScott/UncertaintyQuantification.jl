#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{
    FILE *input = fopen(argv[1], "r");

    double x[2];
    int i;

    for (i = 0; i < 2; i++)
    {
        fscanf(input, "%lf", &x[i]);
    }

    fclose(input);

    double y = sqrt(pow(x[0], 2) + pow(x[1], 2));

    FILE *output = fopen("out.txt", "w");

    fprintf(output, "%.16f", y);

    fclose(output);

    return 0;
}