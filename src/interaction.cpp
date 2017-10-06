#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <Rcpp.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Watershed and cleanup//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

//removes a given cell id from the map
void remove_id(NumericMatrix &padded_map, int cell_id){
    int nrow = padded_map.nrow();
    int ncol = padded_map.ncol();
    for (int i=0; i < nrow; i++){
        for (int j=0; j < ncol; j++){
            if (padded_map(i,j) == cell_id){
                padded_map(i,j) = 0;
            }
        }
    }
}

int fill_in_cell(NumericMatrix &padded_map, int x, int y, int cell_id, NumericMatrix &cell_coords){
    int exit = 1;
    int sid = 0;
    NumericVector x_stack(2000);
    NumericVector y_stack(2000);
    x_stack[sid] = x;
    y_stack[sid] = y;

    //fill in the id
    padded_map(x,y) = cell_id;

    //as long as the stack is not empty
    while (sid >= 0){
        if (sid >= 999){
            //this means that the centroid of the cell was outside of a cell boundary,
            //which makes it impossible to tell where it really belongs
            sid = -1;
            exit = 0;
            //remove old results
            remove_id(padded_map, cell_id);
            //remove cell from the list
            cell_coords(cell_id-1,1) = -1;
            cell_coords(cell_id-1,2) = -1;
        }else{
            //pop coordinate from stack
            x = x_stack[sid];
            y = y_stack[sid];
            sid -= 1;

            //then check all of the direct neighbors, if they are 0 push that coordinate onto the stack
            if (padded_map(x-1,y) == 0){
                padded_map(x-1,y) = cell_id;
                sid += 1;
                x_stack[sid] = x-1;
                y_stack[sid] = y;
            }
            if (padded_map(x,y-1) == 0){
                padded_map(x,y-1) = cell_id;
                sid += 1;
                x_stack[sid] = x;
                y_stack[sid] = y-1;
            }
            if (padded_map(x+1,y) == 0){
                padded_map(x+1,y) = cell_id;
                sid += 1;
                x_stack[sid] = x+1;
                y_stack[sid] = y;
            }
            if (padded_map(x,y+1) == 0){
                padded_map(x,y+1) = cell_id;
                sid += 1;
                x_stack[sid] = x;
                y_stack[sid] = y+1;
            }
        }
    }
    return(exit);
}

//sums all the neighborhood values
int check_neighbors(NumericMatrix &padded_map, int x, int y, int dis){
    int sum = padded_map(x-dis,y);
    sum += padded_map(x,y-dis);
    sum += padded_map(x+dis,y);
    sum += padded_map(x,y+dis);
    return(sum);
}

//only run this if exit is 0 - meaning the cell was not already filled in
//if that coordinate is 0 and it was already filled in it means that the
//filling is ambiguous and we just remove that coordinate
int check_and_fill(NumericMatrix &padded_map,
                   int x,
                   int y,
                   int cell_type,
                   NumericMatrix &cell_coords){
    int exit = fill_in_cell(padded_map, x, y, cell_type, cell_coords);
    cell_coords(cell_type-1,1) = x;
    cell_coords(cell_type-1,2) = y;
    return(exit);
}

//removes a particular coordinate that couldn't be resolved
void remove_coord(NumericMatrix &padded_map, int cell_type, NumericMatrix &cell_coords){
    cell_coords(cell_type-1,1) = -1;
    cell_coords(cell_type-1,2) = -1;
    remove_id(padded_map, cell_type);
}

//fills the open neighborhood with a given cell type
//if it encounters two possibilities it removes the coordinate since it's ambiguous
int fill_neighborhood(NumericMatrix &padded_map, int x, int y, int dis, int cell_type, NumericMatrix &cell_coords){
    int exit = 0;
    bool finished = false;
    // only run this if there is one unique solution
    if (padded_map(x-dis,y) == 0){
        exit = check_and_fill(padded_map, x-dis, y, cell_type, cell_coords);
        finished = true;
    }

    if(padded_map(x,y-dis)==0){
        if (!finished){
            exit = check_and_fill(padded_map, x, y-dis, cell_type, cell_coords);
            finished = true;
        }else{
            exit=0;
        }
    }

    if(padded_map(x+dis,y)==0){
        if (!finished){
            exit = check_and_fill(padded_map, x+dis, y, cell_type, cell_coords);
            finished = true;
        }else{
            exit=0;
        }
    }

    if(padded_map(x,y+dis)==0){
        if (!finished){
            exit = check_and_fill(padded_map, x, y+dis, cell_type, cell_coords);
        }else{
            exit=0;
        }
    }
    if (exit == 0){
        remove_coord(padded_map, cell_type, cell_coords);
    }
    return exit;
}

//Since inForm sometimes provides coordinates that are outside of the original cell we use
//this function to fix that issue as good as possible
void fix_coordinate_issues(NumericMatrix &padded_map, std::list<std::pair <int,int> > problems, NumericMatrix &cell_coords){
    int first_cell_type, first_x, first_y, second_cell_type, second_x, second_y;
    int issues = 0;
    int solved = 0;
    int exit;

    for (std::list<std::pair <int,int> >::iterator it = problems.begin(); it != problems.end(); it++){
        exit = 0;
        //get the current coordinate
        first_cell_type = (*it).first;
        first_x = cell_coords((first_cell_type-1),1);
        first_y = cell_coords((first_cell_type-1),2);

        //if the coordinate is directly on top of a membrane
        if ((*it).second == -1){
            issues++;
            // look only at neighborhood of 1
            exit = fill_neighborhood(padded_map, first_x, first_y, 1, first_cell_type, cell_coords);
            solved += exit;

        }else{
            issues+=2;
            second_cell_type = (*it).second;
            second_x = cell_coords((second_cell_type-1),1);
            second_y = cell_coords((second_cell_type-1),2);

            //start with checking the neighborhood then increase up to a distance of 3
            bool checking = true;
            for (int dis=1;dis<4 && checking;dis++){
                //check if there is a neighboring membrane within the current distance
                int prim = check_neighbors(padded_map,first_x,first_y,dis);
                int sec = check_neighbors(padded_map,second_x,second_y,dis);

                //if both hit membranes at the same time we just drop both coordinates
                if (prim<0 && sec <0){
                    checking = false;
                    cell_coords(first_cell_type-1,1) = -1;
                    cell_coords(first_cell_type-1,2) = -1;
                    cell_coords(second_cell_type-1,1) = -1;
                    cell_coords(second_cell_type-1,2) = -1;
                }

                //if the second coordinate hits a membrane first
                if (prim == 0 && sec < 0){
                    //fill in the first one, exit
                    exit = fill_in_cell(padded_map, first_x, first_y, first_cell_type, cell_coords);
                    // and then see if we can also fill in the second one
                    exit += fill_neighborhood(padded_map, second_x, second_y, dis+1, second_cell_type, cell_coords);
                    checking = false;
                }
                //if the first coordinate hits a membrane first
                if (prim < 0 && sec == 0){
                    //fill in the second one
                    exit = fill_in_cell(padded_map, second_x, second_y, second_cell_type, cell_coords);
                    //and then see if we can fill in the first one also
                    exit += fill_neighborhood(padded_map, first_x, first_y, dis+1, first_cell_type, cell_coords);
                    checking = false;
                }
                solved += exit;
            }
            // if both cells are really close in the middle we need to remove them
            if (checking){
                cell_coords(first_cell_type-1,1) = -1;
                cell_coords(first_cell_type-1,2) = -1;
                cell_coords(second_cell_type-1,1) = -1;
                cell_coords(second_cell_type-1,2) = -1;
            }
        }
    }
    Rprintf("%d solved of %d issues\n",solved,issues);
}


// [[Rcpp::export]]
List watershedC(NumericMatrix padded_map, NumericMatrix cell_coords){
    int num_cells = cell_coords.nrow();
    int x, y;

    //list for cell ids with issues
    std::list<std::pair <int,int> > problems;

    for (int idx = 0; idx < num_cells; idx++) {
        // make sure to account for the coordinate shift between R an C++
        x = cell_coords(idx,1);
        y = cell_coords(idx,2);
        int cell_id = cell_coords(idx,0);

        if (padded_map(x,y) == 0){
            fill_in_cell(padded_map, x, y, cell_id, cell_coords);
        }else{
            //This is the case where there is already a color filled in or the
            //coordinate is on the membrane
            if (padded_map(x,y)==(-1)){
                problems.push_back(std::make_pair(cell_id,(int)-1));
            }else{
                problems.push_back(std::make_pair(cell_id,(int)padded_map(x,y)));
                remove_id(padded_map, (int)padded_map(x,y));
            }
        }
    }
    fix_coordinate_issues(padded_map,problems,cell_coords);
    //return padded_map;
    return List::create(padded_map,cell_coords);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Marker masks //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void outline_cells(NumericMatrix &marker_map, NumericMatrix &padded_map, int x, int y){
    int sid = 0;
    NumericVector x_stack(2000);
    NumericVector y_stack(2000);
    x_stack[sid] = x;
    y_stack[sid] = y;

    //fill in the starting point
    padded_map(x,y) = 1;

    //as long as the stack is not empty
    while (sid >= 0){
        //pop coordinate from stack
        x = x_stack[sid];
        y = y_stack[sid];
        sid -= 1;

        if (padded_map(x-1,y)==-1 ||
            padded_map(x,y-1)==-1 ||
            padded_map(x+1,y)==-1 ||
            padded_map(x,y+1)==-1 ||
            padded_map(x-1,y-1)==-1 ||
            padded_map(x-1,y+1)==-1 ||
            padded_map(x+1,y-1)==-1 ||
            padded_map(x+1,y+1)==-1 ){
            marker_map(x,y) = 1;
        }

        //then check all of the direct neighbors, if they are 0 push that coordinate onto the stack
        if (padded_map(x-1,y) == 0){
            padded_map(x-1,y) = 1;
            sid += 1;
            x_stack[sid] = x-1;
            y_stack[sid] = y;
        }
        if (padded_map(x,y-1) == 0){
            padded_map(x,y-1) = 1;
            sid += 1;
            x_stack[sid] = x;
            y_stack[sid] = y-1;
        }
        if (padded_map(x+1,y) == 0){
            padded_map(x+1,y) = 1;
            sid += 1;
            x_stack[sid] = x+1;
            y_stack[sid] = y;
        }
        if (padded_map(x,y+1) == 0){
            padded_map(x,y+1) = 1;
            sid += 1;
            x_stack[sid] = x;
            y_stack[sid] = y+1;
        }
    }
}

// [[Rcpp::export]]
NumericMatrix generate_maskC(NumericMatrix marker_map, NumericMatrix padded_map, NumericMatrix cell_coords){
    int num_cells = cell_coords.nrow();
    int x, y;
    for (int idx = 0; idx < num_cells; idx++) {
        x = cell_coords(idx,0);
        y = cell_coords(idx,1);
        //since we removed some coordinates
        if (x>0){
            outline_cells(marker_map, padded_map, x, y);
        }
    }
    return marker_map;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Marker masks //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List getInteractionsC(NumericMatrix filled_map){

    std::map<int,std::vector<int> > interactions;
    int nrow = filled_map.nrow();
    int ncol = filled_map.ncol();

    for (int x=2;x<nrow-2;x++){
        for (int y=2;y<ncol-2;y++){
            if (filled_map(x-2,y)>0 && filled_map(x-2,y)!=filled_map(x,y)){
                interactions[filled_map(x,y)].push_back(filled_map(x-2,y));
            }
            if (filled_map(x,y-2)>0 && filled_map(x,y-2)!=filled_map(x,y)){
                interactions[filled_map(x,y)].push_back(filled_map(x,y-2));
            }
            if (filled_map(x+2,y)>0 && filled_map(x+2,y)!=filled_map(x,y)){
                interactions[filled_map(x,y)].push_back(filled_map(x+2,y));
            }
            if (filled_map(x,y+2)>0 && filled_map(x,y+2)!=filled_map(x,y)){
                interactions[filled_map(x,y)].push_back(filled_map(x,y+2));
            }
        }
    }

    for(std::map<int, std::vector<int> >::iterator it = interactions.begin(); it != interactions.end(); it++) {
        std::vector<int> vec = it->second;
        std::sort( vec.begin(), vec.end() );
        vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
        it->second = vec;
    }
    return List::create(interactions);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// fill in marker masks///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_mask(NumericMatrix &mask, NumericMatrix &padded, int xx, int yy){
    int sid = 0;
    int x,y;
    NumericVector x_stack(2000);
    NumericVector y_stack(2000);
    x_stack[sid] = xx;
    y_stack[sid] = yy;

    //fill in the starting point
    padded(xx,yy) = 1;

    //as long as the stack is not empty
    while (sid >= 0){
        if (sid > 999){
            Rprintf("STACK issue that shouldn't exist!\n");
            Rprintf("x: %d y: %d\n",xx,yy);
            break;
        }

        //pop coordinate from stack
        x = x_stack[sid];
        y = y_stack[sid];
        sid -= 1;

        //fill in the id
        mask(x,y) = 2;

        //then check all of the direct neighbors, if they are 0 push that coordinate onto the stack
        if (padded(x-1,y) == 0){
            padded(x-1,y) = 1;
            sid += 1;
            x_stack[sid] = x-1;
            y_stack[sid] = y;
        }
        if (padded(x,y-1) == 0){
            padded(x,y-1) = 1;
            sid += 1;
            x_stack[sid] = x;
            y_stack[sid] = y-1;
        }
        if (padded(x+1,y) == 0){
            padded(x+1,y) = 1;
            sid += 1;
            x_stack[sid] = x+1;
            y_stack[sid] = y;
        }
        if (padded(x,y+1) == 0){
            padded(x,y+1) = 1;
            sid += 1;
            x_stack[sid] = x;
            y_stack[sid] = y+1;
        }
    }
}

// [[Rcpp::export]]
NumericMatrix fillMaskC(NumericMatrix mask, NumericMatrix padded, NumericMatrix cell_coords){
    int num_cells = cell_coords.nrow();
    int x, y;

    for (int idx = 0; idx < num_cells; idx++) {
        // make sure to account for the coordinate shift between R an C++
        x = cell_coords(idx,0);
        y = cell_coords(idx,1);
        if (x>0){
            fill_mask(mask, padded, x, y);
        }
    }
    return mask;
}

