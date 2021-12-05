#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
//S(it)
ofstream outs16,outs8,outs4,outs2,outs1;
//V
ofstream outV16,outV8,outV4,outV2,outV1;


void poisson(double delta,int nx,int ny,double TOL)
{
    //Pliki
    outs16.open("s16.txt");
    outs8.open("s8.txt");
    outs4.open("s4.txt");
    outs2.open("s2.txt");
    outs1.open("s1.txt");
    outV16.open("V16.txt");
    outV8.open("V8.txt");
    outV4.open("V4.txt");
    outV2.open("V2.txt");
    outV1.open("V1.txt");

    double xmax=delta*nx;
    double ymax=delta*ny;
    int it=0; 

    //Tablica potencjalu
    long double V[nx+1][ny+1]={0};
    //Warunki brzegowe
    for(int i=0;i<=ny;i++){
        //VB1
        V[0][i]=sin(M_PI*(i*delta/ymax));
        //Vb3
        V[nx][i]=sin(M_PI*(i*delta/ymax));
    }
    for(int i=0;i<=nx;i++){
        //VB2
        V[i][ny]=-sin(2*M_PI*(i*delta/xmax));
        //VB4
        V[i][0]=sin(2*M_PI*(i*delta/xmax));
    }
    //Poczatkowa wartosc S
    double S=0;
    double S1=0;
    int k=16;
    for(int i=0;i<=nx-k;i+=k){
        for(int j=0;j<=ny-k;j+=k){
            S+=(pow(k*delta,2)/2)*(pow(((V[i+k][j]-V[i][j])/(2*k*delta))+((V[i+k][j+k]-V[i][j+k])/(2*k*delta)),2)+pow(((V[i][j+k]-V[i][j])/(2*k*delta))+((V[i+k][j+k]-V[i+k][j])/2*k*delta),2));
        }        
    }

    //Glowna petla po k
    for(k=16;k>=1;k/=2)
    {   
        it=0;
        do
        {
            for(int i=k;i<=nx-k;i+=k){
                for(int j=k;j<=ny-k;j+=k){
                    V[i][j]=0.25*(V[i+k][j]+V[i-k][j]+V[i][j+k]+V[i][j-k]);
                }
            }
            S1=S;
            S=0;

            for(int i=0;i<=nx-k;i+=k){
                for(int j=0;j<=ny-k;j+=k){
                    S+=(pow(k*delta,2)/2)*(pow(((V[i+k][j]-V[i][j])/(2*k*delta))+((V[i+k][j+k]-V[i][j+k])/(2*k*delta)),2)+pow(((V[i][j+k]-V[i][j])/(2*k*delta))+((V[i+k][j+k]-V[i+k][j])/2*k*delta),2));
                }
            }
            if(k==16) outs16<<it<<" "<<S1<<endl;
            if(k==8)  outs8<<it<<" "<<S1<<endl;
            if(k==4)  outs4<<it<<" "<<S1<<endl;
            if(k==2)  outs2<<it<<" "<<S1<<endl;
            if(k==1)  outs1<<it<<" "<<S1<<endl;
            it++;
        }
        while(fabs((S - S1)/S1) >= TOL);
        //Zageszczanie sieci
        if(k!=1)
        {
            for(int i=0;i<=nx-k;i+=k){
                for(int j=0;j<=ny-k;j+=k){
                    V[i+k/2][j+k/2]=0.25*(V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k]);
                    if(i+k!=nx){
                        V[i+k][j+k/2]=0.5*(V[i+k][j]+V[i+k][j+k]);
                    }
                    if(j+k!=ny){
                        V[i+k/2][j+k]=0.5*(V[i][j+k]+V[i+k][j+k]);
                    }
                    if(j!=0){
                        V[i+k/2][j]=0.5*(V[i][j]+V[i+k][j]);
                    }
                    if(i!=0){
                        V[i][j+k/2]=0.5*(V[i][j]+V[i][j+k]);

                    }
                }
            }
        }
        
        //Zapisywanie do pliku
        for(int i=0;i<=nx;i+=k){
            for(int j=0;j<=ny;j+=k){
                
                if(k==16){
                    outV16<<V[i][j]<<"\t";
                }
                if(k==8)
                    outV8<<V[i][j]<<"\t";
                if(k==4)
                    outV4<<V[i][j]<<"\t";
                if(k==2)
                    outV2<<V[i][j]<<"\t";
                if(k==1)
                    outV1<<V[i][j]<<"\t";
            }
            if(k==16) outV16<<endl;
            if(k==8)  outV8<<endl;
            if(k==4)  outV4<<endl;
            if(k==2)  outV2<<endl;
            if(k==1)  outV1<<endl;
        }
    }
    outs16.close();
    outs8.close();
    outs4.close();
    outs2.close();
    outs1.close();
    outV16.close();
    outV8.close();
    outV4.close();
    outV2.close();
    outV1.close();

}


int main()
{
    double delta=0.2;
    int nx=128;
    int ny=128;
    double TOL=pow(10,-8);
    poisson(delta,nx,ny,TOL);

}