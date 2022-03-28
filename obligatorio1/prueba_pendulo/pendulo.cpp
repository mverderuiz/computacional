#include <iostream>
#include <fstream>
#include <cmath>
#define PI 3.1415

using namespace std;

int main (void)
{

    double h, teta, l, g, v, a, t, aux;
    int i;
    ofstream posicion("posicion.txt");
    
    // Se inicializan constantes y valores predeterminados
    t=0.;
    g=9.8;
    aux=0.;

    //Se le otorgan valores iniciales a las variables del sistema
    h=0.01;
    v=0.;
    teta=30.*PI/180;
    l=1.;
    a=-g*sin(teta)/l;
    
    if (!posicion)
    {
    cout << "Error al abrir ejemplo.dat\n";
    exit(EXIT_FAILURE);
    }

    posicion << l*sin(teta) << ", " << -l*cos(teta) << endl << endl;
    for (i=1; i<=1000; i++)
    {
        aux=v+h*a/2;
        teta=teta+h*aux;
        a=-g*sin(teta)/l;
        v=aux+h*a/2;
        posicion << l*sin(teta) << ", " << -l*cos(teta) << endl << endl;
    }
    posicion.close();
        
    return 0;
}