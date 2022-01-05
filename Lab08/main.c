#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nx 400
#define ny 90
#define i1 200
#define i2 210
#define j1 50
#define delta 0.01
#define sigma 10*delta
#define xA 0.45
#define yA 0.45
#define IT_MAX 10000

void rozw(double D, FILE *fileCX, FILE *filemap)
{
    //Zmienne
    
    //Tablice
    double U0[nx+1][ny+1];
    double U1[nx+1][ny+1];
    double Vx[nx+1][ny+1];
    double Vy[nx+1][ny+1];
    double psi[nx+1][ny+1];

    //Wczytanie funkcji strumienia
    FILE* psifile=fopen("psi.dat","r");
    int i,j;
    double psit;
    while(fscanf(psifile,"%d%d%lf",&i,&j,&psit)==3)
    {
        psi[i][j]=psit;
    }
    fclose(psifile);
    //pole predkosci
    for(int i=1;i<=nx-1;i++){
        for(int j=1;j<=ny-1;j++){
            Vx[i][j]=(psi[i][j+1]-psi[i][j-1])/(2.0*delta);
            Vy[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2.0*delta);
        }
    }
    //Na zastawce
    for(int i=i1;i<=i2;i++){
        for(int j=0;j<=j1;j++){
            Vx[i][j]=0.0;
            Vy[i][j]=0.0;
        }
    }
    //Na dolnym i lewym brzegu
    for(int i=1;i<=nx-1;i++){
        Vx[i][0]=0.0;
        Vy[i][ny]=0.0;
    }
    //Na lewym i prawym brzegu
    for(int j=0;j<=ny;j++){
        Vx[0][j]=Vx[1][j];
        Vx[nx][j]=Vx[nx-1][j];
    }
    //Zapis pol predkosci
    FILE* fileVx=fopen("Vx.txt","w");
    FILE* fileVy=fopen("Vy.txt","w");
    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            fprintf(fileVx,"%lf\t",Vx[i][j]);
            fprintf(fileVy,"%lf\t",Vy[i][j]);
        }
        fprintf(fileVx,"\n");
        fprintf(fileVy,"\n");
    }
    fclose(fileVx);
    fclose(fileVy);
    //Wyznaczenie vmax
    double Vmax=sqrt(Vx[0][0]*Vx[0][0] + Vy[0][0]*Vy[0][0]);
    for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            double Vp=sqrt(Vx[i][j]*Vx[i][j] + Vy[i][j]*Vy[i][j]);
            if(Vp>Vmax){
                Vmax=Vp;
            }
        }
    }
    //Wyznaczenie deltaT
    double deltaT=delta/(4.0*Vmax);
    printf("%g %g \n\n",Vmax,deltaT);
    //inicjalizacja gestosci
    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++)
            U0[i][j] = 1.0/(2.0*M_PI*pow(sigma, 2)) * exp( -(pow((delta*i) - xA, 2) + pow((delta*j) - yA, 2))/(2.0*pow(sigma, 2)));
    }

    for(int it=0;it<=IT_MAX;it++){
        printf("%d \n",it);
        //ustawienie U1=U0
        for(int i=0;i<=nx;i++){
            for(int j=0;j<=ny;j++){
                U1[i][j] = U0[i][j];
            }
        }
        //Start iteracji picarda
        for(int k=1;k<=20;k++){
            for(int i=0;i<nx;i++){
                for(int j=1;j<=ny-1;j++){
                    if(i<i1 || i>i2 || j>j1){
                        if(i==0){
                            U1[i][j] = (1.0/(1.0+((2.0*D*deltaT) / pow(delta, 2)))) * (U0[i][j] - (deltaT/2.0) * Vx[i][j] *
                            (((U0[i+1][j] - U0[nx][j])/(2.0*delta)) + (U1[i+1][j] - U1[nx][j])/(2.0*delta)) - (deltaT / 2.0) * Vy[i][j] *
                             ((U0[i][j+1] - U0[i][j-1])/(2.0*delta) + (U1[i][j+1] - U1[i][j-1])/(2.0*delta)) + (deltaT/2.0) * D *
                             ((U0[i+1][j] + U0[nx][j] + U0[i][j+1] + U0[i][j-1] - 4.0*U0[i][j])/pow(delta,2) + (U1[i+1][j] + U1[nx][j] + U1[i][j+1] + U1[i][j-1] )/pow(delta,2)));
                        }
                        else if(i==nx){
                            U1[i][j] = (1.0/(1.0+( (2.0*D*deltaT) / pow(delta, 2)))) * (U0[i][j] - (deltaT/2.0) * Vx[i][j] *
                            (((U0[0][j] - U0[i-1][j])/(2.0*delta)) + (U1[0][j] - U1[i-1][j])/(2.0*delta)) - (deltaT / 2.0) * Vy[i][j] *
                             ((U0[i][j+1] - U0[i][j-1])/(2.0*delta) + (U1[i][j+1] - U1[i][j-1])/(2.0*delta)) + (deltaT/2.0) * D *
                             ((U0[0][j] + U0[i-1][j] + U0[i][j+1] + U0[i][j-1] - 4.0*U0[i][j])/pow(delta,2) + (U1[0][j] + U1[i-1][j] + U1[i][j+1] + U1[i][j-1])/pow(delta,2)));
                        }
                        else{
                            U1[i][j] = (1.0/(1.0+((2.0*D*deltaT )/ pow(delta, 2)))) * (U0[i][j] - (deltaT/2.0) * Vx[i][j] *
                            (((U0[i+1][j] - U0[i-1][j])/(2.0*delta)) + (U1[i+1][j] - U1[i-1][j])/(2.0*delta)) - (deltaT / 2.0) * Vy[i][j] *
                            ((U0[i][j+1] - U0[i][j-1] )/(2.0*delta) + (U1[i][j+1] - U1[i][j-1])/(2.0*delta)) + (deltaT/2.0) * D *
                            ((U0[i+1][j] + U0[i-1][j] + U0[i][j+1] + U0[i][j-1] - 4.0*U0[i][j])/pow(delta,2) + (U1[i+1][j] + U1[i-1][j] + U1[i][j+1] + U1[i][j-1])/pow(delta,2)));
                        }
                    }
                }
            }
        }
        //Zachowanie rozwiazania
        for(int i=0;i<=nx;i++){
            for(int j=0;j<=ny;j++){
                U0[i][j]=U1[i][j];
            }
        }
        //Wyznaczenie c i x_sr
        double c=0.0;
        double x_sr=0.0;
        for(int i=0;i<=nx;i++){
            for(int j=0;j<=ny;j++){
                c+=U0[i][j];
                x_sr += (i*delta)*U0[i][j];
            }
        }
        c*=pow(delta,2);
        x_sr*=pow(delta,2);
        fprintf(fileCX,"%lf\t%lf\t%lf\n",deltaT*it, c, x_sr);
        printf("%lf\t%lf\t%lf\n",deltaT*it, c, x_sr);

        if(it%700==0 && it<=3500){
             for (int i = 0; i <= nx; i++) {
                for (int j = 0;j <= ny; j++) {
                    fprintf(filemap,"%lf\t",U1[i][j]);
                }
                fprintf(filemap,"\n");
            }
            fprintf(filemap,"\n\n\n");
        }

    }
}

int main()
{
    //D=0.0
    FILE* CX=fopen("CXsrD0.txt","w");
    FILE* MAP=fopen("mapa.txt","w");
    rozw(0.0,CX,MAP);
    fclose(CX);
    fclose(MAP);
    //D=0.1
    FILE* CX1=fopen("CXsrD1.txt","w");
    FILE* MAP1=fopen("mapa1.txt","w");
    //rozw(0.1,CX1,MAP1);
    fclose(CX1);
    fclose(MAP1);
    return 0;
}

