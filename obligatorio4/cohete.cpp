#include <iostream>
#include <cmath>
#include <fstream>
#define PI 3.14159265

using namespace std;
double rpunto (double y1, double y2, double y3, double y4);
double phipunto (double y1, double y2, double y3, double y4);
double prpunto (double y1, double y2, double y3, double y4 , double delta, double mu, double omega, double t);
double pphipunto (double y1, double y2, double y3, double y4 , double delta, double mu, double omega, double t);


int main (void)
{

    double nave[4], theta0, phi0, v0, r0;
    double t=0., h;
    double m=4.6768E4;
    double k[4][4];
    ofstream posiciones("datos.txt");
    int i, j, npasos;

    //Constantes
    double omega=2.6617E-6;
    double dtl=3.844E8;
    double Rt=6.378160E6;
    double Rl=1.7374E6;
    double Mt=5.9736E24;
    double Ml=0.07349E24;
    double G=6.67E-11;
    double delta, mu;
    delta = G*Mt/(dtl*dtl*dtl);
    mu = Ml/Mt;

    //Caracteristicas de la situacion inicial de la nave y modificables
    //theta0 =0.4548;
    //phi0 =0.4548;
    //v0 = 1.119E4/dtl;
    theta0 =0.0;
    phi0 =PI/2;
    v0 = 1.07E4/dtl;
    r0 = Rt/dtl;
    h=1.;
    npasos = 300000;

    //Situacion inicial de la nave
    nave[0]=r0;
    nave[1]=phi0;
    nave[2]=v0*cos(theta0-phi0);
    nave[3]=v0*r0*sin(theta0-phi0);

    //Escribimos en un documento la distribución inicial
    posiciones << 0 << ", " << 0 << endl << cos(omega*t) << ", " << sin(omega*t) << endl;
    posiciones << nave[0]*cos(nave[1]) << ", " << nave[0]*sin(nave[1]) << endl << endl;

    //Inicializimos el algoritmo de rungekuta
    for (i=0; i<=npasos; i++)
    {
        k[0][0]=h*rpunto(nave[0], nave[1], nave[2], nave[3]);
        k[0][1]=h*phipunto(nave[0], nave[1], nave[2], nave[3]);
        k[0][2]=h*prpunto(nave[0], nave[1], nave[2], nave[3], delta, mu, omega, t);
        k[0][3]=h*pphipunto(nave[0], nave[1], nave[2], nave[3], delta, mu, omega, t);

        k[1][0]=h*rpunto(nave[0]+k[0][0]/2, nave[1]+k[0][1]/2, nave[2]+k[0][2]/2, nave[3]+k[0][3]/2);
        k[1][1]=h*phipunto(nave[0]+k[0][0]/2, nave[1]+k[0][1]/2, nave[2]+k[0][2]/2, nave[3]+k[0][3]/2);
        k[1][2]=h*prpunto(nave[0]+k[0][0]/2, nave[1]+k[0][1]/2, nave[2]+k[0][2]/2, nave[3]+k[0][3]/2, delta, mu, omega, t+h/2);
        k[1][3]=h*pphipunto(nave[0]+k[0][0]/2, nave[1]+k[0][1]/2, nave[2]+k[0][2]/2, nave[3]+k[0][3]/2, delta, mu, omega, t+h/2);     
        
        k[2][0]=h*rpunto(nave[0]+k[1][0]/2, nave[1]+k[1][1]/2, nave[2]+k[1][2]/2, nave[3]+k[1][3]/2);
        k[2][1]=h*phipunto(nave[0]+k[1][0]/2, nave[1]+k[1][1]/2, nave[2]+k[1][2]/2, nave[3]+k[1][3]/2);
        k[2][2]=h*prpunto(nave[0]+k[1][0]/2, nave[1]+k[1][1]/2, nave[2]+k[1][2]/2, nave[3]+k[1][3]/2, delta, mu, omega, t+h/2);
        k[2][3]=h*pphipunto(nave[0]+k[1][0]/2, nave[1]+k[1][1]/2, nave[2]+k[1][2]/2, nave[3]+k[1][3]/2, delta, mu, omega, t+h/2);

        k[3][0]=h*rpunto(nave[0]+k[2][0], nave[1]+k[2][1], nave[2]+k[2][2], nave[3]+k[2][3]);
        k[3][1]=h*phipunto(nave[0]+k[2][0], nave[1]+k[2][1], nave[2]+k[2][2], nave[3]+k[2][3]);
        k[3][2]=h*prpunto(nave[0]+k[2][0], nave[1]+k[2][1], nave[2]+k[2][2], nave[3]+k[2][3], delta, mu, omega, t+h);
        k[3][3]=h*pphipunto(nave[0]+k[2][0], nave[1]+k[2][1], nave[2]+k[2][2], nave[3]+k[2][3], delta, mu, omega, t+h);    
        
        for (j=0; j<4; j++) nave[j] = nave[j] + 1./6*(k[0][j]+2*k[1][j]+2*k[2][j]+k[3][j]);

        t=t+h;
        if (i%750==0)
        {
            posiciones << 0 << ", " << 0 << endl << cos(omega*t) << ", " << sin(omega*t) << endl;
            posiciones << nave[0]*cos(nave[1]) << ", " << nave[0]*sin(nave[1]) << endl << endl;
        }
    }
    
    posiciones.close();
    return 0;
}

double rpunto (double y1, double y2, double y3, double y4)
{
    double k;

    k=y3;

    return k;
}

double phipunto (double y1, double y2, double y3, double y4)
{
    double k;

    k=y4/(y1*y1);
    
    return k;
}

double prpunto (double y1, double y2, double y3, double y4 , double delta, double mu, double omega, double t)
{
    double k;
    double rprim;

    rprim=sqrt(y1*y1+1-2*y1*cos(y2-omega*t));
    k=y4*y4/(y1*y1*y1)-delta*(1/(y1*y1)+mu/(rprim*rprim*rprim)*(y1-cos(y2-omega*t)));
    
    return k;
}

double pphipunto (double y1, double y2, double y3, double y4 , double delta, double mu, double omega, double t)
{
    double k;
    double rprim;

    rprim=sqrt(y1*y1+1-2*y1*cos(y2-omega*t));
    k=-delta*mu*y1/(rprim*rprim*rprim)*sin(y2-omega*t);

    return k;
}