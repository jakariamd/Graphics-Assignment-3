#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stack>
#include <vector>
#include <iomanip>
#include<math.h>
#include <stdlib.h>

#include "bitmap_image.hpp"
#define PI 3.14159265
using namespace std;

ifstream infile1, infile2;
ofstream outfile;
int screen_width, screen_height;
double x_left_limit, y_bottom_limit;
double z_front_limit, z_rear_limit;

struct Point
{
    double x,y,z;
};

struct Triangle{
        Point points[3];
        int color[3];
        Triangle *next;
};

Triangle *head, *tail;

void get_config(){
    infile1 >> screen_width >> screen_height;
    infile1 >> x_left_limit;
    infile1 >> y_bottom_limit;
    infile1 >> z_front_limit >> z_rear_limit;
    //cout << screen_width << " " << screen_height << endl << x_left_limit << endl << y_bottom_limit << endl << z_front_limit << " " << z_rear_limit << endl;
}


void get_triangles(){
    head = tail = NULL;
    while(!infile2.eof()){
        Triangle *triangle = new Triangle;
        triangle->next = NULL;
        infile2 >> triangle->points[0].x >> triangle->points[0].y >> triangle->points[0].z ;
        infile2 >> triangle->points[1].x >> triangle->points[1].y >> triangle->points[1].z ;
        infile2 >> triangle->points[2].x >> triangle->points[2].y >> triangle->points[2].z ;

        triangle->color[0] = rand()%256;
        triangle->color[1] = rand()%256;
        triangle->color[2] = rand()%256;

        if(head == NULL){
            head = tail = triangle;
        }
        else{
            tail->next = triangle;
            tail = triangle;
        }
    }

    /*std::cout << '\n' << endl;
    Triangle *triangle = head;
    while (triangle != NULL){
        cout << triangle->points[0].x << " " << triangle->points[0].y << " " << triangle->points[0].z << endl ;
        cout << triangle->points[1].x << " " << triangle->points[1].y << " " << triangle->points[1].z << endl ;
        cout << triangle->points[2].x << " " << triangle->points[2].y << " " << triangle->points[2].z << endl;
        cout << triangle->color[0] << " " << triangle->color[1] << " " << triangle->color[2] << endl;
        std::cout << '\n';
        triangle = triangle->next;
    }*/

}


double top_y, left_x;
double dx, dy;
bitmap_image image;

vector<vector<double> > z_buffer;
void init_buffers(){
    dx = (-2*x_left_limit)/screen_width;
    dy = (-2*y_bottom_limit)/screen_height;
    top_y = -y_bottom_limit - dy/2 ;
    left_x = x_left_limit + dx/2 ;

    for(int i=0; i<screen_width; i++){
        vector<double> it;
        for(int j=0; j<screen_height; j++){
            it.push_back(z_rear_limit);
        }
        z_buffer.push_back(it);
    }

    image.setwidth_height(screen_width,screen_height);
    for(int i=0;i<screen_width;i++){
        for(int j=0;j<screen_height;j++){
            image.set_pixel(i,j,0,0,0);
        }
    }

    //image.save_image("1.bmp");
}

Point* find_top(Triangle *t){
    Point *top_point = &(t->points[0]);
    if(t->points[1].y > top_point->y)
        top_point = &(t->points[1]);
    if(t->points[2].y > top_point->y)
        top_point = &(t->points[2]);
    return top_point;
}

Point* find_bottom(Triangle *t){
    Point *top_point = &(t->points[0]);
    if(t->points[1].y < top_point->y)
        top_point = &(t->points[1]);
    if(t->points[2].y < top_point->y)
        top_point = &(t->points[2]);
    return top_point;
}

Point* find_mid(Triangle *t, Point *top, Point *bottom){
    Point *mid_point = &(t->points[0]);
    if(mid_point != top && mid_point != bottom)
        return mid_point;
    mid_point = &(t->points[1]);
    if(mid_point != top && mid_point != bottom)
        return mid_point;
    return &(t->points[2]);
}

double find_top_scan(double top_point_y){
    if(top_point_y >= top_y)
        return top_y;
    if(top_point_y <= -top_y)
        return -top_y;
    return top_point_y;
}

double find_left_scan(double left_point_x){
    if(left_point_x <= left_x)
        return left_x;
    if(left_point_x >= -left_x)
        return -left_x;
    return left_point_x;
}


void apply_procedure(){
    Triangle *triangle = head;
    while(triangle){
        Point *top_point= find_top(triangle);
        Point *bottom_point= find_bottom(triangle);
        Point *mid_point= find_mid(triangle, top_point,bottom_point);
        double top_scan = find_top_scan(top_point->y);
        double bottom_scan = find_top_scan(bottom_point->y);

        int top_pix = screen_height - ((top_scan + top_y) / dy);
        int bottom_pix = screen_height - ((bottom_scan + top_y) / dy);

        for(int j = top_pix ; j < bottom_pix ; j++){
            if((bottom_point->y - top_point->y)==0) continue;
            double row1 = -y_bottom_limit - ((j+0.5)*(-2*y_bottom_limit)/screen_height);
            double left_scan =  top_point->x + (row1  - top_point->y)*(bottom_point->x - top_point->x)/(bottom_point->y - top_point->y);
            double z1 = top_point->z + (row1  - top_point->y)*(bottom_point->z - top_point->z)/(bottom_point->y - top_point->y);
            double right_scan;
            double z2;
            if(row1  >= mid_point->y){
                if((mid_point->y - top_point->y)==0) continue;
                right_scan = top_point->x + (row1  - top_point->y)*(mid_point->x - top_point->x)/(mid_point->y - top_point->y);
                z2 = top_point->z + (row1  - top_point->y)*(mid_point->z - top_point->z)/(mid_point->y - top_point->y);
            }
            else{
                if((mid_point->y - bottom_point->y)==0) continue;
                right_scan = mid_point->x + (row1  - mid_point->y)*(mid_point->x - bottom_point->x)/(mid_point->y - bottom_point->y);
                z2 = mid_point->z + (row1  - mid_point->y)*(mid_point->z - bottom_point->z)/(mid_point->y - bottom_point->y);
            }
            if(left_scan > right_scan){
                double temp = right_scan;
                double temp2 = z2;
                right_scan = left_scan;
                z2 = z1;
                left_scan = temp;
                z1 = temp2;
            }

            left_scan = find_left_scan(left_scan);
            right_scan = find_left_scan(right_scan);

            int left_pix = floor((left_scan - left_x) / dx) ;
            int right_pix = ceil((right_scan - left_x) / dx) ;

            for(int i = left_pix; i < right_pix; i++){
                double col1 = (i*(-2*x_left_limit)/screen_width)+x_left_limit;
                double z = z1 + (z1-z2)*(col1 - left_scan)/(left_scan-right_scan);

                if(z < z_buffer[i][j] && z >= z_front_limit){
                    z_buffer[i][j] = z;
                    image.set_pixel(i,j,triangle->color[0],triangle->color[1],triangle->color[2]);
                }
            }
        }

        triangle = triangle->next;
    }

}

void save(){
    image.save_image("1.bmp");
    outfile.open("z_buffer.txt");
    for(int i=0; i<z_buffer.size(); i++){
        for(int j=0; j<z_buffer[i].size(); j++){
            if(z_buffer[j][i] < z_rear_limit){
                 outfile << std::fixed << std::setprecision(6) << z_buffer[j][i] << '\t';
                 //cout << z_buffer[j][i] << "\t";
            }
        }
        outfile << '\n';
    }
    z_buffer.clear();
    outfile.close();
}

void free_memory(){
    head = tail = NULL;
    for(int i=0; i<z_buffer.size(); i++){
        z_buffer[i].clear();
    }
    z_buffer.clear();
}

int main(){
    infile1.open("config.txt");
    infile2.open("stage3.txt");
    get_config();
    get_triangles();
    init_buffers();
    apply_procedure();
    save();
    free_memory();

    infile1.close();
    infile2.close();
}
