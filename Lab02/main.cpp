#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

ofstream out1;
ofstream out2;

void Trapez_Picard(int N,double beta, double gamma, int tmax,double deltaT,double u0,double TOL, string fname1,string fname2) 
{
    out1.open(fname1);
    out2.open(fname2);
    double alfa=beta*N-gamma;
    double un=u0;
    double un1;
    int i=0;
    
    for(double t=0;t<tmax;t+=deltaT)
    {
        do{
            un=un1;
            un1=u0+(deltaT/2)*((alfa*u0-beta*u0*u0)+(alfa*un-beta*un*un));
            i++;
        }while(abs(un1-un)<TOL && i<=20);

        u0=un1;
        out1<<t<<" "<<u0<<endl;
        out2<<t<<" "<<N-u0<<endl;
    }
    out1.close();
    out2.close(); 
}

void Trapez_Newton(int N,double beta, double gamma, int tmax,double deltaT,double u0,double TOL, string fname1,string fname2)
{
    out1.open(fname1);
    out2.open(fname2);
    double alfa=beta*N-gamma;
    int i=0;
    double un=u0;
    double un1;

    for(double t=0;t<tmax;t+=deltaT)
    {
        do{
            un=un1;
            double p1=alfa*u0-beta*u0*u0;
            double p2=alfa*un-beta*un*un;
            double p3=alfa-2*beta*un;
            un1=un-((un-u0-(deltaT/2)*(p1+p2))/(1-(deltaT/2)*p3));
            i++;
        }while(abs(un1-un)<TOL && i<=20);
         u0=un1;
         out1<<t<<" "<<u0<<endl;
         out2<<t<<" "<<N-u0<<endl;
    }
    out1.close();
    out2.close(); 
}


void niewjawnaRK2(int N,double beta, double gamma, int tmax,double deltaT,double u0,double TOL, string fname1,string fname2)
{
    out1.open(fname1);
    out2.open(fname2);

    double alfa=beta*N-gamma;

    double A[2][2];
    A[1][1]=1/4;
    A[1][2]=1/4-sqrt(3)/6;
    A[2][1]=1/4+sqrt(3)/6;
    A[2][2]=1/4;

    double b[2]={0.5,0.5};
    double M[2][2];
    for(double t=0;t<tmax;t+=deltaT)
    {
        double U1=u0;
        double U2=u0;
        for(int i=0;i<20;i++)
        {
            double F1=U1-u0-deltaT*(A[1][1]*(alfa*U1-beta*U1*U1)+A[1][2]*(alfa*U2-beta*U2*U2));
            double F2=U2-u0-deltaT*(A[2][1]*(alfa*U1-beta*U1*U1)+A[2][2]*(alfa*U2-beta*U2*U2));

            M[1][1]=1-deltaT*A[1][1]*(alfa-2*beta*U1);
            M[1][2]=-deltaT*A[1][2]*(alfa-2*beta*U2);
            M[2][1]=-deltaT*A[2][1]*(alfa-2*beta*U1);
            M[2][2]=1-deltaT*A[2][2]*(alfa-2*beta*U2);

            double deltaU1=(F2*M[1][2]-F1*M[2][2])/(M[1][1]*M[2][2]-M[1][2]*M[2][1]);
            double deltaU2=(F1*M[2][1]-F2*M[1][1])/(M[1][1]*M[2][2]-M[1][2]*M[2][1]);

            U1=U1+deltaU1;
            U2=U2+deltaU2;
            if(abs(deltaU1)<TOL && abs(deltaU2)<TOL){
                break;
            }
        }

        u0=u0+deltaT*((b[1]*(alfa*U1-beta*U1*U1))+(b[2]*(alfa*U2-beta*U2*U2)));
        out1<<t<<" "<<u0<<endl;
        out2<<t<<" "<<N-u0<<endl;
    }
    out1.close();
    out2.close();
}

int main()
{
    int N=500;
    double beta=0.001;
    double gamma=0.1;
    int tmax=100;
    double deltaT=0.1;
    double u0=1;
    double TOL=pow(10,-6);
    //Metoda Trapezow z metoda picarda 
    Trapez_Picard(N,beta,gamma,tmax,deltaT,u0,TOL,"uP.txt","zP.txt");
    //Metoda Trapezow z iteracja Newtona 
    Trapez_Newton(N,beta,gamma,tmax,deltaT,u0,TOL,"uN.txt","zN.txt");
    //Metoda niejawna RK2
    niewjawnaRK2(N,beta,gamma,tmax,0.1,u0,TOL,"uRK2.txt","zRK2.txt");

    return 0;
}