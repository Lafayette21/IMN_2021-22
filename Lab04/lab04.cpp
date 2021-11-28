#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

ofstream out1, out2,out3;



double gestosc(int x, double xmax,double sigmaX,int y,double ymax,double sigmaY,double delta)
{
    double ro_1=exp(-(pow(x*delta-0.35*xmax,2))/pow(sigmaX,2)-(pow(y*delta-0.5*ymax,2))/pow(sigmaY,2));
    double ro_2=-exp(-(pow(x*delta-0.65*xmax,2))/pow(sigmaX,2)-(pow(y*delta-0.5*ymax,2))/pow(sigmaY,2));
    return ro_1+ro_2; 
}

void globalRelax(double epsilon,double delta, int nx, int ny,double V1,double V2,double TOL,double omegaG,string fmane1,string fname2,string fname3)
{
    out1.open(fmane1);
    out2.open(fname2);
    out3.open(fname3);
    //Nowa tablica potencjalu
    long double Vn[nx+1][ny+1]={0};
    //Stara tablica potencjalu
    long double Vs[nx+1][ny+1]={0};
    
    //ilosc iteracji
    int it=0;
    double xmax=nx*delta;
    double ymax=ny*delta;
    double sigmaX=0.1*xmax;
    double sigmaY=0.1*ymax;
    //Ustawienie poczatku i konca
    for(int i=0;i<=nx;i++){
        Vn[i][0]=V1;
        Vs[i][0]=V1;
        Vn[i][ny]=V2;
        Vs[i][ny]=V2;
    }

    //Poczatkowe S
    double S=0;
    double S1=0;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
            S+=pow(delta,2)*(0.5*pow((Vs[i+1][j]-Vs[i][j])/delta,2)+0.5*pow((Vs[i][j+1]-Vs[i][j])/delta,2)-ro*Vs[i][j]);
        }
    }
    //Glowna petla
    do{
        for(int i=1;i<nx;i++){
            for(int j=1;j<ny;j++){
                double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
                Vn[i][j]=0.25*(Vs[i+1][j]+Vs[i-1][j]+Vs[i][j+1]+Vs[i][j-1]+(pow(delta,2)/epsilon)*ro);
            }
        }

        for(int j=1;j<ny;j++){
            Vn[0][j]=Vn[1][j];
            Vn[nx][j]=Vn[nx-1][j];
        }
        for(int i = 0; i <= nx; i++){
            for(int j = 1; j < ny; j++) {
                Vs[i][j] = (1.0 - omegaG) * Vs[i][j] + omegaG * Vn[i][j];
            }
        }

        S1=S;
        S=0;
        for(int i=0;i<nx;i++){
            for(int j=0;j<ny;j++){
                double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
                S+=pow(delta,2)*(0.5*pow((Vs[i+1][j]-Vs[i][j])/delta,2)+0.5*pow((Vs[i][j+1]-Vs[i][j])/delta,2)-ro*Vs[i][j]);
            }
        }
        it++;
        out1<<it<<" "<<S1<<endl;

    }
    while(fabs((S-S1)/S1)>=TOL);
    //Blad
    long double blad[nx+1][ny+1]={0}; 
    for(int i=1;i<nx;i++){
        for(int j=1;j<ny;j++){
            double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
            blad[i][j]=((Vn[i+1][j]-2*Vn[i][j]+Vn[i-1][j])/pow(delta,2))+(Vn[i][j+1] - 2*Vn[i][j] + Vn[i][j-1])/(pow(delta,2)) + (ro/epsilon);
            out2<<blad[i][j]<<"\t";
        }
        out2<<endl;
    }
    //Tablica potencjalu
    for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            out3<<Vn[i][j]<<"\t";
        }
        out3<<endl;
    }
    out1.close();
    out2.close();
    out3.close();
}

void localRelax(double epsilon,double delta, int nx, int ny,double V1,double V2,double TOL,double omegaL,string fmane1,string fname2)
{
    out1.open(fmane1);
    out2.open(fname2);
    //Tablica potencjalu
    long double V[nx+1][ny+1]={};

    int it=0;
    double xmax=nx*delta;
    double ymax=ny*delta;
    double sigmaX=0.1*xmax;
    double sigmaY=0.1*ymax;
    //Ustawienie poczatku i konca
    for(int i=0;i<=nx;i++){
        V[i][0]=V1;
        V[i][ny]=V2;
    }
    //Poczatkowe S
    double S=0;
    double S1=0;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
            S+=pow(delta,2)*(0.5*pow((V[i+1][j]-V[i][j])/delta,2)+0.5*pow((V[i][j+1]-V[i][j])/delta,2)-ro*V[i][j]);
        }
    }
    cout<<S<<endl;
    //Glowna petla
    do{
        for(int i=1;i<nx;i++){
            for(int j=1;j<ny;j++){
                double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
                V[i][j]=(1-omegaL)*V[i][j]+(omegaL/4)*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]+(pow(delta,2)/epsilon)*ro);
            }
        }

        for(int j=1;j<ny;j++){
            V[0][j]=V[1][j];
            V[nx][j]=V[nx-1][j];
        }
        S1=S;
        S=0;
        for(int i=0;i<nx;i++){
            for(int j=0;j<ny;j++){
                double ro=gestosc(i,xmax,sigmaX,j,ymax,sigmaY,delta);
                S+=pow(delta,2)*(0.5*pow((V[i+1][j]-V[i][j])/delta,2)+0.5*pow((V[i][j+1]-V[i][j])/delta,2)-ro*V[i][j]);
            }
        }
        out1<<it<<" "<<S1<<endl;
        it++;


    }while(fabs((S-S1)/S1)>=TOL);
    //Tablica potencjalow
    for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            out2<<V[i][j]<<"\t";
        }
        out2<<endl;
    }
    out1.close();
    out2.close();
}
int main()
{  
    double epsilon=1;
    double delta=0.1;
    int nx=150;
    int ny=100;
    double V1=10;
    double V2=0;
    double TOL=pow(10,-8);
    //global relax omegaG=0.6
    //globalRelax(epsilon,delta,nx,ny,V1,V2,TOL,0.6,"glob_relax_0.6.txt","blad_0.6.txt","Gpotencjal_0.6.txt");
    //global relax omegaG=1.0
    //globalRelax(epsilon,delta,nx,ny,V1,V2,TOL,1.0,"glob_relax_1.0.txt","blad_1.0.txt","Gpotencjal_1.0.txt");
    
    //Local relax omegaL=1.0
    localRelax(epsilon,delta,nx,ny,V1,V2,TOL,1.0,"loc_relax_1.0.txt","Lpotencjal_1.0.txt");
    //Local relax omegaL=1.4
    localRelax(epsilon,delta,nx,ny,V1,V2,TOL,1.4,"loc_relax_1.4.txt","Lpotencjal_1.4.txt");
    //Local relax omegaL=1.8
    localRelax(epsilon,delta,nx,ny,V1,V2,TOL,1.8,"loc_relax_1.8.txt","Lpotencjal_1.8.txt");
    //Local relax omegaL=1.9
    localRelax(epsilon,delta,nx,ny,V1,V2,TOL,1.9,"loc_relax_1.9.txt","Lpotencjal_1.9.txt");
    
    return 0;

}