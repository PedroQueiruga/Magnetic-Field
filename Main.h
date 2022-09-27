// System Headers

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "Mathematics.h"
#include "Physics.h"
//#include "Atmosphere.h"

//Global variables

//char ch;
//char string[8],line[100];
//double altitude[1481],temperature[1481],pressure[1481],density[1481];

const double q=1.602*pow(10,-19); //carga de um próton
double m=1.672*pow(10.0, -27); //massa de um próton
const double UA=1.495978707*pow(10,11); //unidades astronômicas
const int c=299792458; //Light velocity
const double mi0=12.5663706143592*pow(10,-7);
const double mi=8.00*pow(10,22);
const double rs=2.469*pow(10,20); //radius of sun's orbit
const double kp=3.24078*pow(10,-20); //meters to kiloparsec

int i,j,k,l,counter=0;
double t,w,T;
double modr,modv,modB,modF,modp;

double O[3]={0,0,0};
double F[3],B[3],B0[3],BD[3],BG[3],p[3];
double r[3],r0[3]={1E16,1E17,0},v[3],v0[3],a[3],angle,modvper,modvpar;
double m_0,mr;
double E,gama;
//double r0[3]={1E16,1E17,0};

double vpar[3],vper[3];
double x_min,y_min,z_min,x_max,y_max,z_max;
double prod[3],rL,rL0;
double vol;















