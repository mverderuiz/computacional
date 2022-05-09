//Comprobar que la norma se conserva
//Modificar el valor de lambda

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <complex>
#define PI 3.14159265

using namespace std;

int main (void)
{
    //Variables modificables
    int N=1000;
    int nciclos=50;
    int duracion=500;
    float lambda=0.3;
    float s=0.01; //Discretización del tiempo
    float h=0.01; //Discretización del espacio

    ofstream datos("datos.txt");
    ofstream modulotxt("modulo.txt");
    //Otras variables
    complex<float> phi[N+1], xi[N+1], aux3;
    complex<float> A0[N], b[N], beta[N], alpha[N], gamma[N];
    float modulo;
    float Aplus, Aminus;
    float k0virg, V[N+1], svirg;
    float aux1, aux2;
    int j, n;    

    k0virg=2.*PI*nciclos/N;
    svirg=1/(4*k0virg*k0virg);

    //Calculo del potencial
    for (j=0; j<N+1; j++)
    {
        aux1=(float) j;
        aux2=(float)N/5.;
        if ((aux1>=2*aux2)&&(aux1<=3*aux2)) V[j]=lambda*k0virg*k0virg;
        else V[j]=0.;
    }

    //Cálculo de la función de onda inical
    phi[0]=complex<float>(0.,0.);
    phi[N]=complex<float>(0.,0.);
    for (j=1; j<N; j++)
    {
        aux1=exp(-8.*(4*j-N)*(4*j-N)/(N*N));
        phi[j]=complex<float>(cos(k0virg*j),sin(k0virg*j))*aux1;
    }
    //Normalizacion de la funcion de onda
    modulo=0.;
    for (j=1; j<N; j++) modulo= modulo + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
    modulo=sqrt(modulo);
    for (j=1; j<N; j++) phi[j]=phi[j]/modulo;

    //Calculamos los coeficientes gamma y alpha
    Aplus=1.;
    Aminus=1.;
    for (j=1; j<N; j++) A0[j]=complex<float>(-2.-V[j],2*Aminus/svirg);
    beta[N-1]=0.;
    alpha[N-1]=0.;
    gamma[N-1]=Aminus/A0[N-1];
    for (j=N-2; j>=0; j--) 
    {
        alpha[j]=-Aminus*gamma[j+1];
        gamma[j]=Aplus/(A0[j]+alpha[j]);
    }

    xi[0]=0;
    for (j=0; j<N+1; j++) datos << j*h << ", " << norm(phi[j]) << "\n"; datos << endl;

    //Iniciamos el bucle de pasos
    for (n=0; n<duracion; n++)
    {
        //Escribimos el modulo
        modulo=0.;
        for (j=1; j<N; j++) modulo= modulo + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
        modulo=sqrt(modulo);
        modulotxt << n << "\t" << modulo << endl;

        for (j=1; j<N; j++) b[j]=complex<float>(0, 4)*phi[j]/svirg;
        for (j=N-2; j>=0; j--) beta[j]=gamma[j+1]*(b[j+1]-beta[j+1]);    

        for (j=1; j<N; j++)
        {
            xi[j]=alpha[j-1]*xi[j-1]+beta[j-1];
            phi[j]=xi[j]-phi[j];
        } 

        for (j=0; j<N+1; j++) datos << j*h << ", " << norm(phi[j]) << "\n";
        datos << endl;
    }

    modulotxt.close();
    datos.close();
    return 0;
}