// Dan Wolf
// LocatedScore.h

#ifndef DNA_LOCALITY_LOCATEDSCORE_H
#define DNA_LOCALITY_LOCATEDSCORE_H

struct Coordinates { int i, j; };

enum Direction {
    None, North, West, Northwest
};

class Score {
public:
    Direction direction;
    int Value;
};

class LocatedScore : public Score {
public:
    Coordinates Location;

    // Next returns the coordinates in the direction of the this Score
    Coordinates Next();
};


#endif //DNA_LOCALITY_LOCATEDSCORE_H
