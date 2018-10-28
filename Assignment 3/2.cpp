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

struct Triangle
{
    Point points[3];
    int color[3];
    Triangle *next;
};

Triangle *head, *tail;

void get_config()
{
    infile1 >> screen_width >> screen_height;
    infile1 >> x_left_limit;
    infile1 >> y_bottom_limit;
    infile1 >> z_front_limit >> z_rear_limit;
    //cout << screen_width << " " << screen_height << endl << x_left_limit << endl << y_bottom_limit << endl << z_front_limit << " " << z_rear_limit << endl;
}


void get_triangles()
{
    head = tail = NULL;
    while(!infile2.eof())
    {
        Triangle *triangle = new Triangle;
        triangle->next = NULL;
        infile2 >> triangle->points[0].x >> triangle->points[0].y >> triangle->points[0].z ;
        infile2 >> triangle->points[1].x >> triangle->points[1].y >> triangle->points[1].z ;
        infile2 >> triangle->points[2].x >> triangle->points[2].y >> triangle->points[2].z ;

        triangle->color[0] = rand()%256;
        triangle->color[1] = rand()%256;
        triangle->color[2] = rand()%256;

        if(head == NULL)
        {
            head = tail = triangle;
        }
        else
        {
            tail->next = triangle;
            tail = triangle;
        }
    }
}


double top_y, left_x;
double dx, dy;


void init_buffers()
{
    dx = (-2*x_left_limit)/screen_width;
    dy = (-2*y_bottom_limit)/screen_height;
    top_y = -y_bottom_limit - dy/2 ;
    left_x = x_left_limit + dx/2 ;

}

struct Polygon
{
    int id;
    double plane[4];
    int color[3];
    bool in_out;
    Polygon *next;
};

struct Edge
{
    int x_ymin;
    int ymax;
    double del_x;
    Edge *next;
    Polygon *polygon;
};


Point* find_top(Triangle *t)
{
    Point *top_point = &(t->points[0]);
    if(t->points[1].y > top_point->y)
        top_point = &(t->points[1]);
    if(t->points[2].y > top_point->y)
        top_point = &(t->points[2]);
    return top_point;
}

Point* find_bottom(Triangle *t)
{
    Point *top_point = &(t->points[0]);
    if(t->points[1].y < top_point->y)
        top_point = &(t->points[1]);
    if(t->points[2].y < top_point->y)
        top_point = &(t->points[2]);
    return top_point;
}

Point* find_mid(Triangle *t, Point *top, Point *bottom)
{
    Point *mid_point = &(t->points[0]);
    if(mid_point != top && mid_point != bottom)
        return mid_point;
    mid_point = &(t->points[1]);
    if(mid_point != top && mid_point != bottom)
        return mid_point;
    return &(t->points[2]);
}

double find_top_scan(double top_point_y)
{
    if(top_point_y >= top_y)
        return top_y;
    if(top_point_y <= -top_y)
        return -top_y;
    return top_point_y;
}

double find_left_scan(double left_point_x)
{
    if(left_point_x <= left_x)
        return left_x;
    if(left_point_x >= -left_x)
        return -left_x;
    return left_point_x;
}

vector<vector<Edge> > edge_table_vect;


void insert_edge_table2(Point *p1, Point *p2, Polygon *polygon)   //scanline
{
    Edge edge ;

    double top_scan = find_top_scan(p1->y);
    double right_scan = (p1->x);
    double left_scan = (p2->x);
    double bottom_scan = find_top_scan(p2->y);

    int top_pix =  (top_scan + top_y) / dy;
    int bottom_pix = (bottom_scan + top_y) / dy;
    int left_pix = ceil((left_scan - left_x) / dx);
    int right_pix = ceil((right_scan - left_x) / dx);


    int top_pix1 = p1->y / dy;
    int bottom_pix1 = p2->y / dy;
    int left_pix1 = p2->x / dx;
    int right_pix1 = p1->x / dx;


    edge.x_ymin = left_pix;
    edge.ymax = top_pix;
    edge.del_x = (dy/dx)*(p1->x - p2->x)/(p1->y - p2->y);
    double del2 = (top_pix - bottom_pix)/( right_pix- left_pix);
    edge.polygon = polygon;


    edge.next = NULL;

    if(bottom_pix >=0 and bottom_pix < screen_height)
    {
        edge_table_vect[bottom_pix].push_back(edge);
    }

    for(int i=0; i< edge_table_vect[bottom_pix].size(); i++){
        for(int j=i+1; j< edge_table_vect[bottom_pix].size(); j++){
            if(edge_table_vect[bottom_pix][j].x_ymin < edge_table_vect[bottom_pix][i].x_ymin){
                Edge e = edge_table_vect[bottom_pix][j];
                edge_table_vect[bottom_pix][j] = edge_table_vect[bottom_pix][i];
                edge_table_vect[bottom_pix][i] = e;
            }
        }

    }
}

void calc_plane_eqn(Polygon *polygon,Point *top_point, Point *mid_point, Point *bottom_point){

    Point p1,p2;
    {
        p1.x = top_point->x - bottom_point->x;
        p1.y = top_point->y - bottom_point->y;
        p1.z = top_point->z - bottom_point->z;

        p2.x = top_point->x - mid_point->x;
        p2.y = top_point->y - mid_point->y;
        p2.z = top_point->z - mid_point->z;
    }
    double a = p1.y * p2.z - p1.z * p2.y;
    double b = p1.z * p2.x - p1.x * p2.z;
    double c = p1.x * p2.y - p1.y * p2.x;
    double d = a * top_point->x + b * top_point->y + c * top_point->z ;

    polygon->plane[0] = a;
    polygon->plane[1] = b;
    polygon->plane[2] = c;
    polygon->plane[3] = d;
}


void init_edge_table_and_polygon_table2()
{
    for(int i=0; i<screen_height; i++)
    {
        vector<Edge> it;
        edge_table_vect.push_back(it);
    }

    int id_stack = 1;

    Triangle *triangle = head;
    while(triangle){
        Point *top = find_top(triangle);
        Point *bottom = find_bottom(triangle);
        Point *mid = find_mid(triangle,top,bottom);

        Polygon *polygon = new Polygon;
        polygon->id = id_stack; id_stack++;
        polygon->color[0] = triangle->color[0];
        polygon->color[1] = triangle->color[1];
        polygon->color[2] = triangle->color[2];
        calc_plane_eqn(polygon,top,mid,bottom);
        polygon->in_out = false;
        polygon->next = NULL;

        //cout << top->y << " " << mid->y << " " << bottom->y << endl;
        if(top->y != mid->y){
            insert_edge_table2(top,mid,polygon);
        }
        if(top->y != bottom->y){
            insert_edge_table2(top,bottom,polygon);
        }
        if(bottom->y != mid->y){
            insert_edge_table2(mid,bottom,polygon);
        }

        triangle = triangle->next;
    }

}


struct Active_Edge
{
    int ymax;
    double x_A;
    double del_x;
    Polygon *polygon;
    Active_Edge *next;
};

Active_Edge *active_edge_table;
Polygon *active_polygon_table;
bitmap_image image;

vector<Active_Edge>active_edge_table_vect;
vector<Polygon*>active_polygon_table_vect;

void update_active_edges3(int scan_y){  //scanline
    for(int i =0; i < edge_table_vect[scan_y].size(); i++){
        Active_Edge edge ;
        Edge temp1 = edge_table_vect[scan_y][i];
        edge.del_x = temp1.del_x;
        edge.ymax = temp1.ymax;
        edge.x_A = temp1.x_ymin;
        edge.polygon = temp1.polygon;
        edge.next = NULL;

        active_edge_table_vect.push_back(edge);
    }
    return;

}



void remove_active_edges3(int y_scan){ //scanline
    vector<Active_Edge> temp;
    for(int i=0; i<active_edge_table_vect.size();i++){
        if(active_edge_table_vect[i].ymax > y_scan){
            temp.push_back(active_edge_table_vect[i]);
        }
    }
    active_edge_table_vect.clear();
    for(int i=0; i<temp.size();i++){
        active_edge_table_vect.push_back(temp[i]);
    }
}




void sort_active_edges(int y_scan){  //scanline
    if(active_edge_table_vect.size()==0) return;

    for(int i=0; i<active_edge_table_vect.size(); i++){
        for(int j=i+1; j<active_edge_table_vect.size(); j++){
            if(active_edge_table_vect[j].x_A < active_edge_table_vect[i].x_A){
                Active_Edge e = active_edge_table_vect[j];
                active_edge_table_vect[j] = active_edge_table_vect[i];
                active_edge_table_vect[i] = e;
            }
        }
    }
}


void update_active_polygon3(Polygon *polygon){ //scanline
    active_polygon_table_vect.push_back(polygon);
}


void remove_active_poligon3(Polygon *polygon){ //scanline
    if(active_polygon_table_vect.size() == 0) return;
    for(int i = 0; i< active_polygon_table_vect.size(); i++){
        if(active_polygon_table_vect[i] == polygon){
            active_polygon_table_vect.erase(active_polygon_table_vect.begin()+i);
            break;
        }
    }
}




void apply_procedure3()
{
    image.setwidth_height(screen_width,screen_height);
    for(int i=0; i<screen_width; i++){
        for(int j=0; j<screen_height; j++)
        {
            image.set_pixel(i,j,0,0,0);
        }
    }

    int y_scan;
    for(y_scan=0; y_scan<screen_height; y_scan++){
        if(edge_table_vect[y_scan].size() != 0)break;
    }

    for(int j = y_scan ; j < screen_height ; j++)
    {
        update_active_edges3(j);
        sort_active_edges(j);

        if(active_edge_table_vect.size()==0)continue;
        for(int p=0 ; p < active_polygon_table_vect.size(); p++){
            active_polygon_table_vect[p]->in_out = false;
        }
        active_polygon_table_vect.clear();

        for(int e = 0; e < active_edge_table_vect.size()-1 ; e++){
            Active_Edge edge = active_edge_table_vect[e];
            edge.polygon->in_out = !edge.polygon->in_out;
            if(edge.polygon->in_out) update_active_polygon3(edge.polygon);
            else remove_active_poligon3(edge.polygon);

            Polygon *min_p = NULL ;
            double z_min_val = z_rear_limit;

            double left_scan = (edge.x_A*(-2*x_left_limit)/screen_width)+x_left_limit;
            double row = y_bottom_limit + ((j+0.5)*(-2*y_bottom_limit)/screen_height);

            for(int k = 0; k < active_polygon_table_vect.size(); k++){
                Polygon *p = active_polygon_table_vect[k];

                if(p->plane[2] == 0) {cout << "c=0\n"; continue;}
                double z = (p->plane[3] - p->plane[0] * left_scan - p->plane[1] * row ) / p->plane[2]  ;

                if(z < z_min_val && z >= z_front_limit){
                    min_p = p;
                    min_p->id = p->id;
                    min_p->color[0] = p->color[0];
                    min_p->color[1] = p->color[1];
                    min_p->color[2] = p->color[2];
                    z_min_val = z;
                }
            }

            if(min_p){
                double exa = edge.x_A, nexa = active_edge_table_vect[e+1].x_A ;
                if(exa < 0) exa = 0; if(nexa >= screen_width) nexa = screen_width-1;
                for(int i=exa; i<nexa;i++){
                    image.set_pixel(i,screen_height- j,min_p->color[0],min_p->color[1],min_p->color[2]);
                }
            }
        }

        remove_active_edges3(j);
        for(int k=0; k<active_edge_table_vect.size();k++){
            active_edge_table_vect[k].x_A += active_edge_table_vect[k].del_x;
        }
        sort_active_edges(j);
    }

}




void save()
{
    image.save_image("2.bmp");
}

void free_memory()
{
    head = tail = NULL;
    active_edge_table_vect.clear();
    active_polygon_table_vect.clear();
}

int main()
{
    infile1.open("config.txt");
    infile2.open("stage3.txt");
    get_config();
    get_triangles();
    init_buffers();

    init_edge_table_and_polygon_table2();
    apply_procedure3();
    save();
    free_memory();

    infile1.close();
    infile2.close();
}

