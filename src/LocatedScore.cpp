//
// Created by dan on 2020-09-29.
//

#include "LocatedScore.h"

Coordinates LocatedScore::Next() {
    Coordinates next;
    switch (this->direction){
        case North:
            next.i = this->Location.i - 1;
            next.j = this->Location.j;
            break;
        case Northwest:
            next.i = this->Location.i - 1;
            next.j = this->Location.j - 1;
            break;
        case West:
            next.i = this->Location.i;
            next.j = this->Location.j - 1;
            break;
        case None:
            break;
    }
    return next;
}
