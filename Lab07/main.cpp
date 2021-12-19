#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define nx 200
#define ny 90

ofstream out1,out2,out3,out4;

void WBpsi(double delta,double psi[nx+1][ny+1], int i1,int j1,double mi,double Qwe,double Qwy){
    double Yny=ny*delta;
    double Yj1=j1*delta;
    //Brzeg A (wejscie)
    for(int j=j1;j<=ny;j++){
        double y=j*delta;
        psi[0][j]=(Qwe/(2*mi))*(pow(y,3)/3-(pow(y,2)/2)*(Yj1+Yny)+y*Yj1*Yny);
    }
    //Brzeg C (wyjscie)
    for(int j=0;j<=ny;j++){
        double y=j*delta;
        psi[nx][j]=(Qwy/(2*mi))*(pow(y,3)/3-(pow(y,2)/2)*Yny)+(Qwe*pow(Yj1,2)*(-Yj1+3*Yny))/(12*mi);
    }
    //Brzeg B
    for(int i=1;i<=nx-1;i++){
        psi[i][ny]=psi[0][ny];
    }
    //brzeg D
    for(int i=i1;i<=nx-1;i++){
        psi[i][0]=psi[0][j1];
    }
    //brzeg E
    for(int j=1;j<=j1;j++){
        psi[i1][j]=psi[0][j1];
    }
    //brzeg F
    for(int i=1;i<=i1;i++){
        psi[i][j1]=psi[0][j1];
    }
}

void WBdzeta(double delta,double dzeta[nx+1][ny+1],double psi[nx+1][ny+1], int i1,int j1,double mi,double Qwe,double Qwy)
{
    double Yny=ny*delta;
    double Yj1=j1*delta;
    //brzeg A (wejscie)
    for(int j=j1;j<=ny;j++){
        double y=j*delta;
        dzeta[0][j]=(Qwe/(2*mi))*(2*y-Yj1-Yny);
    }
    //brzeg C (wyjscie)
    for(int j=0;j<=ny;j++){
        double y=j*delta;
        dzeta[nx][j]=(Qwy/(2*mi))*(2*y-Yny);
    }
    //brzeg B
    for(int i=1;i<=nx-1;i++){
        dzeta[i][ny]=(2/pow(delta,2))*(psi[i][ny-1]-psi[i][ny]);
    }
    //brzeg D
    for(int i=i1+1;i<=nx-1;i++){
        dzeta[i][0]=(2/pow(delta,2))*(psi[i][1]-psi[i][0]);
    }
    //brzeg E
    for(int j=1;j<=j1-1;j++){
        dzeta[i1][j]=(2/pow(delta,2))*(psi[i1+1][j]-psi[i1][j]);
    }
    //brzeg F
    for(int i=1;i<=i1;i++){
        dzeta[i][j1]=(2/pow(delta,2))*(psi[i][j1+1]-psi[i][j1]);
    }
    //wierzcholek E/F
    dzeta[i1][j1]=0.5*(psi[i1-1][j1]+psi[i1][j1-1]);

}

void rozw(double Qwe,string fname1,string fname2,string fname3,string fname4)
{
    double delta=0.01;
    double ro=1;
    double mi=1;
    int i1=50;
    int j1=55;
    int IT_MAX=20000;

    double psi[nx+1][ny+1]={0};
    double dzeta[nx+1][ny+1]={0};
    double u[nx+1][ny+1]={0};
    double v[nx+1][ny+1]={0};

    /*for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            cout<<psi[i][j]<<" ";
        }
        cout<<endl;
    }*/

    double Yny=ny*delta;
    double Yj1=j1*delta;
    double Qwy=Qwe*((pow(Yny,3)-pow(Yj1,3)-3*Yj1*pow(Yny,2)+3*pow(Yj1,2)*Yny)/pow(Yny,3));
    double omega=0;

   WBpsi(delta,psi,i1,j1,mi,Qwe,Qwy);
   for(int it=1;it<=IT_MAX;it++)
   {
        if(it<2000){
           omega=0;
        }
        else{
           omega=1;
        }
        for(int i=1;i<=nx-1;i++){
            for(int j=1;j<=ny-1;j++){
                if(i>i1 || j>j1 ){
                    psi[i][j]=0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-pow(delta,2)*dzeta[i][j]);
                    dzeta[i][j]=0.25*(dzeta[i+1][j]+dzeta[i-1][j]+dzeta[i][j+1]+dzeta[i][j-1])-omega*(ro/(16*mi))*((psi[i][j+1]-psi[i][j-1])*(dzeta[i+1][j]-dzeta[i-1][j])-(psi[i+1][j]-psi[i-1][j])*(dzeta[i][j+1]-dzeta[i][j-1]));
                    u[i][j] =(psi[i][j+1] - psi[i][j-1])/(2*delta);
                    v[i][j] = -(psi[i+1][j] - psi[i-1][j])/(2*delta);
                }
                else{
                    u[i][j] = 0;
                    v[i][j] = 0;
                }
            }
        }    
        WBdzeta(delta,dzeta,psi,i1,j1,mi,Qwe,Qwy);
        //blad
        double gamma=0;
        int j2=j1+2;
        for(int i =1; i<=nx-1; i++){
            gamma += (psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4.0*psi[i][j2] - pow(delta,2)*dzeta[i][j2]);
        }
    }
    //Zapisywanie do plikow
    out1.open(fname1);
    out2.open(fname2);
    out3.open(fname3);
    out4.open(fname4);
    for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            out1<<psi[i][j]<<"\t";
            out2<<dzeta[i][j]<<"\t";
            out3<<u[i][j]<<"\t";
            out4<<v[i][j]<<"\t";
        }
        out1<<endl;
        out2<<endl;
        out3<<endl;
        out4<<endl;
    }
    out1.close();
    out2.close();
    out3.close();
    out4.close();
}


int main()
{
    //Qwe=-1000
    rozw(-1000,"psiQ-1000.txt","dzetaQ-1000.txt","uQ-1000.txt","vQ-1000.txt");
    //Qwe=-4000
    rozw(-4000,"psiQ-4000.txt","dzetaQ-4000.txt","uQ-4000.txt","vQ-4000.txt");
    //Qwe=4000
    rozw(4000,"psiQ4000.txt","dzetaQ4000.txt","uQ4000.txt","vQ4000.txt");
    
    return 0;
}