# include <stddef.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include "mgmres.c"

void wypiszDOUBLE(double* tab,int N){
    for(int i=0;i<N;i++){
        printf("%lf ",*(tab+i));
    }
    printf("\n");
}
void wypiszINT(int* tab,int N){
    for(int i=0;i<N;i++){
        printf("%d ",*(tab+i));
    }
    printf("\n");
}


void rozw(double delta,int nx,int ny,double V1,double V2, double V3,double V4,double eps1,double eps2,double ro1, double ro2, FILE* file1,FILE* file2,FILE* file3)
{
    int i=0,j=0;
    int N=(nx+1)*(ny+1);

    int nz_num = 0;
    double xmax=delta*nx;
    double ymax=delta*ny;
    double sigma=xmax/10;
    //Parametry do pmgmres_ilu_cr
    int itr_max=500;
    int mr=500;
    double tol_abs = pow(10,-8);
    double tol_rel = pow(10,-8);
    //Wektory definiujace macierz ukladu
    //a
    double a[5*N];
    for(int i=0;i<5*N;i++){
        a[i]=0.0;
    }
    //ja
    int ja[5*N];
    //ia
    int ia[N+1];
    for(int w=0;w<=N;w++){
        ia[w]=-1;
    }
    //b
    double b[N];
    //V
    double V[N];
    for(int w=0;w<N;w++){
        V[w]=0.0;
        b[w]=0.0;
    }

    double epsL[N];

    for(int l=0;l<N;l++){
        j = l/(nx+1);
        i=l-j*(nx+1);
        if (i <= nx / 2)
            epsL[l]=eps1;
        else
            epsL[l]=eps2;
     }
    
    //Algorytm wypelniania macierzy rzadkiej 
    int k=-1;
    for(int l=0;l<N;l++){
        j= (int)(l/(nx+1));
        i= l-j*(nx+1);

        int brzeg=0;
        double vb=0.0;
        if(i==0){ //lewy brzeg
            brzeg=1;
            vb=V1;
        }
        if(j==ny){ //gorny brzeg
            brzeg=1;
            vb=V2;
        }
        if(i==nx){ //prawy brzeg
            brzeg=1;
            vb=V3;
        }
        if(j==0){ //dolny brzeg
            brzeg=1;
            vb=V4;
        }
        //wypelnienie wektora wyrazow wolnych
        double ro=0.0;
        if( ro1 == -1.0 && ro2 ==-1.0){
            ro = exp(-(pow(delta*i - 0.25*xmax, 2)/pow(sigma,2))- (pow(delta*j - 0.5*ymax, 2)/pow(sigma,2)))+
            (-exp(-(pow(delta*i - 0.75*xmax, 2)/pow(sigma,2))- (pow(delta*j - 0.5*ymax, 2)/pow(sigma,2)))) ;
        }
        else{
            ro=ro1-ro2;
        }

        b[l]= -ro;
        //wartosc potencjalu na brzegu
        if(brzeg==1){ 
            b[l]=vb;
        }
        //Wypelnianie elementow macierzy A
        ia[l]=-1;
        //lewa skrajna przekatna
        if(l-nx-1>=0 && brzeg==0){
            k++;
            if(ia[l]<0){ 
                ia[l]=k;
            }
            a[k]=epsL[l]/pow(delta,2);
            ja[k]=l-nx-1;
        }
        //poddiagonala
        if(l-1>=0 && brzeg==0){
            k++;
            if(ia[l]<0){ 
                ia[l]=k;
            }
            a[k]= epsL[l]/pow(delta,2);
            ja[k] = l - 1;
        }
        //diagonala
        k++;
        if(ia[l]<0){ 
            ia[l]=k;
        }
        if(brzeg==0){
            a[k]=-(2*epsL[l]+epsL[l+1]+epsL[l+nx+1])/pow(delta,2);
        }
        else{
            a[k]=1;
        }
        ja[k]=l;
        //naddiagonala
        if(l<N && brzeg==0){
            k++;
            a[k] = epsL[l+1]/pow(delta,2);
            ja[k]=l+1;
        }
        //prawa skrajna przekatna
        if(l<N-nx-1 && brzeg==0){
            k++;
            a[k] = epsL[l+nx+1]/pow(delta,2);
            ja[k] = l + nx +1;
        }
 
    }
    if(nx==4){
        for(int l=0;l<5*N;l++){
            j= (int)(l/(nx+1));
            i= l-j*(nx+1);
            fprintf(file1,"%d\t%d\t%d\t%g\n",l,i,j,a[l]);
        }
        for(int l=0;l<N;l++){
            j= (int)(l/(nx+1));
            i= l-j*(nx+1);
            fprintf(file2,"%d\t%d\t%d\t%g\n",l,i,j,b[l]);
            
        }
    }
    


    fprintf(file2, "\n\n\n");
    fprintf(file1, "\n\n\n");

    nz_num = k+1;
    printf("%d ",nz_num);
    ia[N] = nz_num;
    
    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);
    for(int l = 0; l < N; l++){
        fprintf(file3, "%g\t", V[l]);
        if((l+1)%(nx+1) == 0){
            fprintf(file3, "\n");
        }
    }
}


int main()
{
    //void rozw(double delta,int nx,int ny,double V1,double V2, double V3,double V4,double eps1,double eps2,double ro1, double ro2, FILE* file1,FILE* file2,FILE* file3)
    FILE* file_A;
    file_A = fopen("macierzA4.txt","w");
    FILE* file_b;
    file_b = fopen("vektorb4.txt", "w");
    //nx=4 ny=4 eps1=1 eps2=1 V1=V3=10 V2=V4=-10 ro1=ro1=0
    FILE* file_V4;
    file_V4 = fopen("V4.txt", "w");
    rozw(0.1,4,4,10,-10,10,-10,1.0,1.0,0.0,0.0,file_A,file_b,file_V4);
    fclose(file_V4);
    //nx=50 ny=50 eps1=1 eps2=1 V1=V3=10 V2=V4=-10 ro1=ro2=0 
    FILE* file_V50;
    file_V50=fopen("V50.txt","w");
    rozw(0.1,50,50,10,-10,10,-10,1.0,1.0,0.0,0.0,file_A,file_b,file_V50);
    fclose(file_V50);
    //nx=100 ny=100 eps1=1 eps2=1 V1=V3=10 V2=V4=-10 ro1=ro2=0
    FILE* file_V100;
    file_V100=fopen("V100.txt","w");
    rozw(0.1,100,100,10,-10,10,-10,1.0,1.0,0.0,0.0,file_A,file_b,file_V100);
    fclose(file_V100);
    //nx=200 ny=200 eps1=1 eps2=1 V1=V3=10 V2=V4=-10 ro1=ro2=0
    FILE* file_V200;
    file_V200=fopen("V200.txt","w");
    rozw(0.1,200,200,10,-10,10,-10,1.0,1.0,0.0,0.0,file_A,file_b,file_V200);
    fclose(file_V200);
    //nx=100 ny=100 V1=V2=V3=V4=0 eps1=1.0 eps2=1.0 ro1=ro2=-1.0
    FILE* file_Veps11;
    file_Veps11=fopen("V_eps11.txt","w");
    rozw(0.1,100,100,0.0,0.0,0.0,0.0,1.0,1.0,-1.0,-1.0,file_A,file_b,file_Veps11);
    fclose(file_Veps11);
    //nx=100 ny=100 V1=V2=V3=V4=0 eps1=1.0 eps2=2.0 ro1=ro2=-1.0
    FILE* file_Veps12;
    file_Veps12=fopen("V_eps12.txt","w");
    rozw(0.1,100,100,0.0,0.0,0.0,0.0,1.0,2.0,-1.0,-1.0,file_A,file_b,file_Veps12);
    fclose(file_Veps12);
    //nx=100 ny=100 V1=V2=V3=V4=0 eps1=1.0 eps2=10.0 ro1=ro2=-1.0
    FILE* file_Veps110;
    file_Veps110=fopen("V_eps110.txt","w");
    rozw(0.1,100,100,0.0,0.0,0.0,0.0,1.0,10.0,-1.0,-1.0,file_A,file_b,file_Veps110);
    fclose(file_Veps110);

    fclose(file_A);
    fclose(file_b);
    return 0;
}