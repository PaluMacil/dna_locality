// Dan Wolf
// SimilarityMatrix.cpp

#include <cstring>
#include <chrono>
#include <iterator>
#include <sstream>
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
    best = {None, 0, 0, 0};
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
#pragma omp parallel
            {
#pragma omp task shared(i, strandS, strandT, j) depend(in:this->matrix[i - 1][j], this->matrix[i][j - 1], this->matrix[i - 1][j - 1]) depend(out:this->matrix[i][j])
                {
                    // if in a cell initialized to 0 (first column and row), skip other calculations
                    if (!(i == 0 || j == 0)) {
                        score(i, j);
                    }
                }
            }
        }
    }

    // end time
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    printOut();

    return duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
}

std::string SimilarityMatrix::MatchString() {
    std::vector<char> finalT;
    std::vector<char> finalS;

    // sequence after match of t
    for (int i = strandS.size() - 1; i > best.Location.i; i--) {
        finalS.insert(finalS.begin(), strandS[i]);
    }

    // keep track of first matched place in sequence
    int firstMatchedIndex = best.Location.i;
    // match space
    for (LocatedScore s = best; s.Value > 0; s = At(s.Next())) {
        switch (s.direction) {
            case North:
                finalT.insert(finalT.begin(), '-');
                finalS.insert(finalS.begin(), strandS[s.Location.i]);
                firstMatchedIndex = s.Location.i;
                break;
            case Northwest:
                finalT.insert(finalT.begin(), strandT[s.Location.j]);
                finalS.insert(finalS.begin(), strandS[s.Location.i]);
                firstMatchedIndex = s.Location.i;
                break;
            case West:
                finalT.insert(finalT.begin(), strandT[s.Location.j]);
                finalS.insert(finalS.begin(), '-');
                break;
            case None:
                break;
        }
    }

    // sequence before match of t
    for (int i = firstMatchedIndex - 1; i > 0; i--) {
        finalS.insert(finalS.begin(), strandS[i]);
        finalT.insert(finalT.begin(), ' ');
    }

    std::stringstream ss;
    ss << "Sequence (s):\t";
    for (auto p: finalS) { ss << p << ' '; }
    ss << '\n';
    ss << "Unknown (t):\t";
    for (auto p: finalT) { ss << p << ' '; }
    ss << '\n';

    return ss.str();
}

void SimilarityMatrix::score(int i, int j) {
    std::vector<Score> scores;

    // F(s[1..i],t[1..j-1])-2
    Score westScore = {West, matrix[i][j - 1] - 2};
    scores.push_back(westScore);

    // F(s[1..i-1],t[1..j-1])+(if s[i]=t[j] then 1 else - 1)
    // '?' is a wildcard and always matches the other value
    int predicate = (strandS[i] == strandT[j] || strandS[i] == '?' || strandS[j] == '?')
                    ? 1 : -1;
    Score northwestScore = {Northwest, matrix[i - 1][j - 1] + predicate};
    scores.push_back(northwestScore);

    // F(s[1..i-1],t[1..j])-2
    Score northScore = {North, matrix[i - 1][j] - 2};
    scores.push_back(northScore);

    auto max = MaxScore(scores);
    matrix[i][j] = max.Value;
    directions[i][j] = max.direction;
    // save if best score in matrix (replace if equal, as furthest bottom right value wins)
#pragma omp critital
    {
        if (max.Value >= best.Value) {
            best.Location = {.i = i, .j = j};
            best.direction = max.direction;
            best.Value = max.Value;
        }
    };
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
        outFile << '\n' << "max = " << best.Value << ",\ti = " << best.Location.i << ",\tj = " << best.Location.j
                << '\n';
        outFile << MatchString();
    }
}

LocatedScore SimilarityMatrix::At(Coordinates location) {
    LocatedScore locatedScore{};
    locatedScore.Location = location;
    locatedScore.direction = directions[location.i][location.j];
    locatedScore.Value = matrix[location.i][location.j];

    return locatedScore;
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