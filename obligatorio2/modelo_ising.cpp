#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#define PI 3.1415

using namespace std;
float energiasist(int red[][100], int N);

int main (void)
{
    //Inicializamos todas las variables que vamos a utilizar
    int N, i ,k, j, red[100][100], escrib;
    float T, aux, energia, deltaE;
    bool aleatorio;
    ofstream sistema("ising_data.txt");
    int xmas, ymas, xmenos, ymenos, xal, yal;
    int npasos, iteraciones;
    float p, expo, epsilon;

    srand (time(NULL));

    //Definimos el nº de nudos de la red bi-dimensionar
    N=64;
    //Introducimos la temperatura del sistema
    T=0.0001;
    //Numero de pasos del programa
    npasos = 1000;
    //Elementos de la red aleatorios o constantes
    aleatorio=true;

    //Inicilizamos los elementos de la red   
    for (i=0; i<100; i++)
    {
        for (k=0; k<100; k++) red[i][k] = 0;
    }    
    if (aleatorio==true)
    {
        for (i=0; i<N; i++)
        {
            for (k=0; k<N; k++) 
            {
                aux=(rand() % 100);
                if (aux<50) red[i][k]=1;
                else red[i][k]=-1;
            }
        }
    }
    else
    {
        for (i=0; i<N; i++)
        {
            for (k=0; k<N; k++) red[i][k]=1;
        }
    }
    //Escribimos en un fichero los espines.
    for (i=0; i<N; i++)
    {
        for (k=0; k<N-1; k++) 
        {
            sistema << red[i][k] << ", " ;
        }
        sistema << red[i][N-1] << endl;
    }
    sistema << endl;   

    //Inicializamos el algoritmo
    iteraciones = npasos*N*N;
    for (i=0; i<iteraciones; i++)
    {
        xal = rand() % N;
        yal = rand() % N;

            xmas=xal+1;
            xmenos=xal-1;
            ymas=yal+1;
            ymenos=yal-1;
            if (xmenos==-1) xmenos=N-1;
            if (xmas==N) xmas=0;
            if (ymenos==-1) ymenos=N-1;
            if (ymas==N) ymas=0;
        deltaE = (float) 2*red[xal][yal]*(red[xmas][yal]+red[xmenos][yal]+red[xal][ymas]+red[xal][ymenos]);
        expo= exp(-deltaE/T);
        if (expo<1) p=expo;
        else p=1;
        epsilon = (rand() % 1001)/1000.;
        if (epsilon<p) red[xal][yal]=-red[xal][yal];
        //Escribimos la nueva matriz
        escrib=N*N;
        if (i%escrib==0)
        {
            for (j=0; j<N; j++)
            {
                for (k=0; k<N-1; k++) 
                {
                    sistema << red[j][k] << ", " ;
                }
                sistema << red[j][N-1] << endl;
            }
            sistema << endl;
        }
    }

    sistema.close();
    return 0;
}

float energiasist(int red[][100], int N)
{
    int imas, kmas, imenos, kmenos;
    float energia;
    int i,k;

    for(i=1;i<N-1;i++)
    {
        for(k=1;k<N-1;k++)
        {
            kmas=k+1;
            kmenos=k-1;
            imas=i+1;
            imenos=i-1;
            if (imenos==-1) imenos=N-1;
            if (imas==N) imas=0;
            if (kmenos==-1) kmenos=N-1;
            if (kmas==N) kmas=0;
            energia=energia+red[i][k]*(red[i][kmas]+red[i][kmenos]+red[imas][k]+red[imenos][k]);
        }
    }
    energia=energia/2.;

    return energia;
}