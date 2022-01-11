#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::string;

std::ofstream out1, out2;

int deltaKroneckera(double x, double xF){
    if(fabs(x-xF) < pow(10, -8)){
        return 1;
    }
    else{
        return 0;
    }
}


void rozwiazanie(double alfa,double beta,string fname1,string fname2){
    //Pliki
    out1.open(fname1);
    out2.open(fname2);
    //Parametry
    int nx=150;
    int nt=1000;
    double delta=0.1;
    double deltaT=0.05;
    double xA=7.5;
    double sigma=0.5;
    double tmax=nt*deltaT;
    double xF=2.5;
    //Tablice (od razu są tu warunki brzegowe)
    double u0[nx+1]={0};
    double u[nx+1]={0};
    double v[nx+1]={0};
    double vp[nx+1]={0};
    double a[nx+1]={0};
    //Warunki początkowe
    for(int i=1;i<=nx-1;i++){
        u[i]=exp(-pow(i*delta-xA,2)/(2*pow(sigma,2)));
        v[i]=0;
    }
    //Zachowanie wyniku
    for(int i=0;i<=nx;i++){
        u0[i]=u[i];
    }
    //Wyznaczenie a
    for(int i=1;i<=nx-1;i++){
        a[i]=(u[i+1]-2*u[i]+u[i-1])/pow(delta,2)-beta*((u[i]-u0[i])/deltaT)+alfa*cos(50*0.0/tmax)*deltaKroneckera(delta*i,xF);
    }
    //Algorytm
    for(int n=1;n<=nt;n++){
        for(int i=0;i<=nx;i++){
            vp[i]=v[i]+(deltaT/2)*a[i];
        }

        for(int i=0;i<=nx;i++){
            u0[i]=u[i];
        }

        for(int i=0;i<=nx;i++){
            u[i]=u[i]+deltaT*vp[i];
        }
        for(int i=1;i<=nx-1;i++){
            a[i]=a[i]=(u[i+1]-2*u[i]+u[i-1])/pow(delta,2)-beta*((u[i]-u0[i])/deltaT)+alfa*cos(50*(deltaT*n)/tmax)*deltaKroneckera(delta*i,xF);
        }
        //Liczenie E
        double E=0;
        E = delta/4.0*(pow((u[1]-u[0])/delta,2)+pow((u[nx]-u[nx-1])/delta,2));
        for(int i=1;i<=nx-1;i++){
            E+=delta/2.0*(pow(v[i],2)+pow(((u[i+1]-u[i-1])/(2*delta)),2));
        }
        //Wynik
        out1<<deltaT*n<<"\t"<<E<<endl;
        for(int i=0;i<=nx;i++){
            out2<<u[i]<<"\t";
        }
        out2<<endl;

    }
    out1.close();
    out2.close();
}



int main(){
    //alfa=0 beta=0
    rozwiazanie(0,0,"E1.txt","u1.txt");
    //alfa=0 beta=0.1
    rozwiazanie(0,0.1,"E2.txt","u2.txt");
    //alfa=0 beta=1
    rozwiazanie(0,1,"E3.txt","u3.txt");
    //alfa=1 beta=1
    rozwiazanie(1,1,"E4.txt","u4.txt");
    return 0;
}