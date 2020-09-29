// Dan Wolf
// SimilarityMatrix.cpp

#include <cstring>
#include <chrono>
#include "SimilarityMatrix.h"

Score MaxScore(const std::vector<Score> &scores) {
    Score best{None, 0};
    for (auto score : scores) {
        if (score.Value > best.Value) {
            best = {score.direction, score.Value};
        }
    }

    return best;
}

SimilarityMatrix::SimilarityMatrix(const char *filename1, const char *filename2, const char *filenameOut)
        : inFile1(filename1, std::ifstream::binary),
          inFile2(filename2, std::ifstream::binary),
          outFile() {
    if (!inFile1) {
        throw "could not open inFile1";
    }
    if (!inFile2) {
        throw "could not open inFile2";
    }
    printOutput = strcmp(filenameOut, "") != 0;
    if (printOutput) {
        outFile.open(filenameOut, std::ios::trunc);
        if (!outFile) {
            throw "could not open outfile";
        }
    }
}

SimilarityMatrix::~SimilarityMatrix() = default;

void SimilarityMatrix::Load() {
    strandT = ReadFile(std::move(inFile1));
    strandS = ReadFile(std::move(inFile2));
}

double SimilarityMatrix::Fill() {
    // start time
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // i is the row index and j is the column index.
    // With strandT across the top and strandS going down the matrix,
    // the number of rows is the size of strand 2,
    // and the number of columns is the size of strandT.

    // expand i vector (rows) by size of strand 2
    matrix.resize(strandS.size());
    directions.resize(strandS.size());
    for (int i = 0; i < strandS.size(); i++) {
        // expand j vector (columns) by size of strand 1, initialize all scores to 0
        matrix[i].resize(strandT.size(), 0);
        directions[i].resize(strandT.size(), None);
        for (int j = 0; j < strandT.size(); j++) {
            // if in a cell initialized to 0 (first column and row), skip other calculations
            if (i == 0 || j == 0) {
                continue;
            }
            score(i, j);
        }
    }

    // end time
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    printOut();

    return duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
}

std::string SimilarityMatrix::Calculate() {
    return std::string();
}

void SimilarityMatrix::score(int i, int j) {
    std::vector<Score> scores;

    // F(s[1..i],t[1..j-1])-2
    Score westScore = {West, matrix[i][j - 1] - 2};
    scores.push_back(westScore);

    // F(s[1..i-1],t[1..j-1])+(if s[i]=t[j] then 1 else - 1)
    // TODO: compare i+1 and j+1 in strands to account for offset from corner of grid
    int predicate = strandS[i] == strandT[j] ? 1 : -1;
    Score northwestScore = {Northwest, matrix[i - 1][j - 1] + predicate};
    scores.push_back(northwestScore);

    // F(s[1..i-1],t[1..j])-2
    Score northScore = {North, matrix[i - 1][j] - 2};
    scores.push_back(northScore);

    auto max = MaxScore(scores);
    matrix[i][j] = max.Value;
    directions[i][j] = max.direction;
}

void SimilarityMatrix::printOut() {
    if (printOutput) {
        // print first row
        outFile << "  ";
        for (auto protein : strandT) {
            auto proteinChar = protein == 'I' ? ' ' : protein;
            outFile << proteinChar << ' ';
        }
        outFile << '\n';

        // print each row
        for (int i = 0; i < strandS.size(); i++) {
            auto proteinChar = strandS[i] == 'I' ? ' ' : strandS[i];
            outFile << proteinChar << " ";
            for (int j = 0; j < strandT.size(); j++) {
                outFile << matrix[i][j] << " ";
            }
            outFile << '\n';
        }
    }
}

std::vector<char> ReadFile(std::ifstream inFile) {
    std::vector<char> strand;
    // The initial (I) row and column are all initialized to a 0 score
    strand.push_back('I');

    for (char currentChar; !inFile.eof(); inFile.get(currentChar)) {
        // only grab valid characters and the wildcard (?)
        if (currentChar == 'A' || currentChar == 'G' || currentChar == 'T' ||
            currentChar == 'C' || currentChar == '?') {
            strand.push_back(currentChar);
        }
    }

    return strand;
}