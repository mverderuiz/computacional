#include <iostream>
#include <fstream>
#include <cmath>
#define PI 3.1415

using namespace std;

float AlgoritmoVerlet2D(double p[][6], double m[], int Np, double h, float t);
float newtonsol2D(double p2[][6], double m[], int Np);

float AlgoritmoVerlet2D(double p[][6], double m[], int Np, float h, float t)
{
    int i;
    double aux[Np][2];
    //Le damos valor al vector w_i para cada planeta    
    for (i=0; i<Np; i++)
    {
    aux[i][0]=p[i][2]+h*p[i][4]/2;
    aux[i][1]=p[i][3]+h*p[i][5]/2; 
    }
    // Evaluamos el nuevo valor de la posición 
    for (i=0; i<Np; i++)
    {
    p[i][0]=p[i][0]+h*aux[i][0];
    p[i][1]=p[i][1]+h*aux[i][1];
    }
    //Evaluamos las nuevas acelaraciones aplicando la ley de Newton
    newtonsol2D(p,m,Np);
    //Evaluamos las nuevas velocidades
    for (i=0; i<Np; i++)
    {
    p[i][2]=aux[i][0]+h*p[i][4]/2;
    p[i][3]=aux[i][1]+h*p[i][5]/2;
    }

    return t+h;
}

float newtonsol2D(double p2[][6], double m[], int Np)
{
    int i,k;
    float modulo, posicion[2];

    //Reiniciamos todos las aceleraciones
    for (i=0;i<Np;i++)
    {
        p2[i][4]=0.;
        p2[i][5]=0.;
    }
    //Calculamos las nuevas aceleraciones
    for (i=0;i<Np;i++)
    {
        if (i!=Np-1)
        {
        for(k=i+1;k<Np;k++)
        {
            modulo=pow(sqrt( (p2[i][0]-p2[k][0])*(p2[i][0]-p2[k][0])+(p2[i][1]-p2[k][1])*(p2[i][1]-p2[k][1]) ), 3);
            posicion[0]=p2[i][0]-p2[k][0];
            posicion[1]=p2[i][1]-p2[k][1];
            p2[i][4]=p2[i][4]-m[k]*posicion[0]/modulo;
            p2[i][5]=p2[i][5]-m[k]*posicion[1]/modulo;
            p2[k][4]=p2[k][4]+m[i]*posicion[0]/modulo;
            p2[k][5]=p2[k][5]+m[i]*posicion[1]/modulo;
        }
        }
        
    }
    return 0;    
}

int main (void)
{
    //Enteros para ciclos y numero de planetas
    int i,k,z;
    int Np=9;
    //Matriz con los ocho planetas y sus posicones, velocidads y aceleración
    double planetas[9][6];
    //Variable del tiempo y paso
    float t=0.;
    float h=0.01;
    //Archivo donde se van a escribir las posiciones
    ofstream planpos("descripcion_geocentrica.txt");

    //Vector que contiene las masas de los planetas y del sol en masas solares
    double masas[9];
    masas[0]=1; //Masa del sol
    masas[1]=1.66E-7; //Masa de Mercurio
    masas[2]=2.44E-6; //Masa de Venus
    masas[3]=3.02E-6; //Masa de La Tierra
    masas[4]=3.23E-7; //Masa de Marte
    masas[5]=9.56E-4; //Masa de Júpiter
    masas[6]=2.96E-4; //Masa de Saturno
    masas[7]=4.38E-5; //Masa de Urano
    masas[8]=5.30E-5; //Masa de Neptuno

    //Definimos las posiciones (en UA) y velocidades iniciales (en m/s) de los planetas
    for (k=0;k<9;k++)
    {
        planetas[k][1]=0.;
        planetas[k][2]=0.;
    }
        planetas[0][0]=0.; //Posicion del Sol
        planetas[0][3]=0.; //Velocidad del Sol
        planetas[1][0]=0.38; //Posicion incial de Mercurio
        planetas[1][3]=48.92E3; //Velocidad inical de Mercurio
        planetas[2][0]=0.72; //Posicion incial de Venus
        planetas[2][3]=35.02E3; //Velocidad inical de Venus
        planetas[3][0]=1.; //Posicion incial de La Tierra
        planetas[3][3]=29.78E3; //Velocidad inical de La Tierra
        planetas[4][0]=1.52; //Posicion incial de Marte
        planetas[4][3]=24.07E3; //Velocidad inical de Marte
        planetas[5][0]=5.21; //Posicion incial de Jupiter
        planetas[5][3]=13.05E3; //Velocidad inical de Jupiter
        planetas[6][0]=9.54; //Posicion incial de Saturno
        planetas[6][3]=9.64E3; //Velocidad inical de Saturno
        planetas[7][0]=19.18; //Posicion incial de Urano
        planetas[7][3]=6.81E3; //Velocidad inical de Urano
        planetas[8][0]=30.11; //Posicion incial de Neptuno
        planetas[8][3]=5.43E3; //Velocidad inical de Neptuno
    //Reescalamos las velocidades y posiciones iniciales (Para que coincidan con las usadas en la aceleración)
    for (k=1;k<9;k++)
    {
        planetas[k][3]=planetas[k][3]*3.36E-5;
    }
    //A partir de las posiciones calculamos sus acelraciones inciales
        newtonsol2D(planetas, masas, Np);
    //Escribimos las posiciones iniciales en el fichero
    for(k=0;k<9;k++) planpos << planetas[k][0]-planetas[3][0] <<  ", " << planetas[k][1]-planetas[3][1] << endl; planpos << endl;
        
    //Aplicamos el algoritmo de Verlet
    for(i=1;i<10000;i++)
    {
        //Calculamos las nuevas posiciones con Verlet.
        t=AlgoritmoVerlet2D(planetas, masas, Np, h, t);
        if (i%5==0)
        {
            for(k=0;k<8;k++)  planpos << planetas[k][0]-planetas[3][0] <<  ", " << planetas[k][1]-planetas[3][1] << endl;
            planpos << endl;
        }
    }    

    planpos.close();
    return 0;
}