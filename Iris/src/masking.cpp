#include <algorithm>
#include <iostream>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

const double PI_ = 3.141592653589793;


int** getCircleCoordinates(int radius, int points){
    double slice = 2 * PI_ / points;
    //empty container for the coordinates
    int** coords = 0;
    coords = new int *[points];

    //calculate all coordinates
    for (int i = 0; i < points; i++){
        coords[i] = new int[2];
        double angle = slice * i;
        coords[i][0] = (int)(radius * cos(angle));
        coords[i][1] = (int)(radius * sin(angle));
    }
    return coords;
}

//standard function to delete a dynamic array
void deleteArray(int **coords, int radius){
    for (int i=0; i < radius*4; i++){
        delete [] coords[i];
    }
    delete [] coords;
}

NumericMatrix expand(NumericMatrix mask, int** coords, int points){
    bool found;
    int xAdj, yAdj;
    NumericMatrix filled(mask.nrow(),mask.ncol());
    for (int x=0;x<mask.nrow();x++){
        for (int y=0;y<mask.ncol();y++){
            found = false;
            for (int i=0;!found && i < points; i++){
                xAdj = x + coords[i][0];
                yAdj = y + coords[i][1];
                //only check the coordinate if its within the image boundary
                if (xAdj >= 0 && xAdj < mask.nrow() &&
                    yAdj >= 0 && yAdj < mask.ncol()){
                    if (mask(xAdj, yAdj) != 0){
                        filled(x, y) = 1;
                        found = true;
                    }
                }
            }
        }
    }
    return filled;
}

// [[Rcpp::export]]
NumericMatrix growMarginC(NumericMatrix mask, int distance){
    //number of points around the coordinate we need
    int points = distance * 4;

    //get all the circumference coordinates
    int** coords = getCircleCoordinates(distance, points);

    //expand the margin
    NumericMatrix filled = expand(mask, coords, points);

    deleteArray(coords, distance);
    return filled;
}


