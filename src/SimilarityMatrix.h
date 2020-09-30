// Dan Wolf
// SimilarityMatrix.h

#ifndef DNA_LOCALITY_SIMILARITYMATRIX_H
#define DNA_LOCALITY_SIMILARITYMATRIX_H

#include <fstream>
#include <vector>
#include "LocatedScore.h"

std::vector<char> ReadFile(std::ifstream inFile);



class SimilarityMatrix {
private:
    std::ifstream inFile1;
    std::ifstream inFile2;
    std::ofstream outFile;
    std::vector<char> strandT;
    std::vector<char> strandS;
    std::vector<std::vector<int>> matrix;
    std::vector<std::vector<Direction>> directions;
    LocatedScore best;
    bool printOutput;

    void score(int i, int j);
    void printOut();

public:
    explicit SimilarityMatrix(const char *filename1, const char *filename2, const char *filenameOut);

    ~SimilarityMatrix();

    void Load();

    double Fill();

    std::string MatchString();

    LocatedScore At(Coordinates location);
};


#endif //DNA_LOCALITY_SIMILARITYMATRIX_H
