// Dan Wolf
// main.cpp

#include <string>
#include <cstring>
#include "SimilarityMatrix.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        std::printf("Got %d arguments, expected a filename with no extension and true or false to set printing.\n",
                    argc - 1);
        return 1;
    }
    auto baseFilename = argv[1];
    auto printFlag = argv[2];
    char filenameOut[256] = {0};
    char filename1[256] = {0};
    char filename2[256] = {0};
    sprintf(filename1, "%s.1.txt", baseFilename);
    sprintf(filename2, "%s.2.txt", baseFilename);
    if (strcmp(printFlag, "true") == 0) {
        sprintf(filenameOut, "%s.out.txt", baseFilename);
    }
    SimilarityMatrix matrix(filename1, filename2, filenameOut);
    matrix.Load();
    double time = matrix.Fill();
    printf("time to fill matrix: %fns\n", time);

    return 0;
}
