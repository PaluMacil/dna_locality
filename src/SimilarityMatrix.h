// Dan Wolf
// SimilarityMatrix.h

#ifndef DNA_LOCALITY_SIMILARITYMATRIX_H
#define DNA_LOCALITY_SIMILARITYMATRIX_H

#include <fstream>
#include <vector>

std::vector<char> ReadFile(std::ifstream inFile);

enum Direction {
    None, North, West, Northwest
};

struct Score {
    Direction direction;
    int Value;
};

class SimilarityMatrix {
private:
    std::ifstream inFile1;
    std::ifstream inFile2;
    std::ofstream outFile;
    std::vector<char> strandT;
    std::vector<char> strandS;
    std::vector<std::vector<int>> matrix;
    std::vector<std::vector<Direction>> directions;
    bool printOutput;

    void score(int i, int j);
    void printOut();

public:
    explicit SimilarityMatrix(const char *filename1, const char *filename2, const char *filenameOut);

    ~SimilarityMatrix();

    void Load();

    double Fill();

    std::string Calculate();
};


#endif //DNA_LOCALITY_SIMILARITYMATRIX_H
