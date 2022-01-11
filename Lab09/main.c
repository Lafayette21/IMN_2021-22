#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

void dyfuzjaMacierzowa(){
    //parametry
    int nx=40;
    int ny=40;
    int N=(nx+1)*(ny+1);
    double delta=1.0;
    double deltaT=1.0;
    int TA=40;
    int TB=0;
    int TC=30;
    int TD=0;
    double kB=0.1;
    double kD=0.6;
    int IT_MAX=2000;
    //macierze
    gsl_matrix* A=gsl_matrix_calloc(N,N);
    gsl_matrix* B=gsl_matrix_calloc(N,N);
    gsl_vector* d=gsl_vector_calloc(N);
    gsl_vector* c=gsl_vector_calloc(N);

    //Do dekompozycji LU
    gsl_permutation *p=gsl_permutation_calloc(N);
    int signum=0;

    int l;
    
    for(int i=1;i<=nx-1;i++){
        for(int j=1;j<=ny-1;j++){
            l=i+j*(nx+1);
            //Wypelnienie A
            gsl_matrix_set(A, l, l-nx-1, deltaT/(2*pow(delta,2)));
            gsl_matrix_set(A, l, l-1,deltaT/(2*pow(delta,2)));
            gsl_matrix_set(A, l, l+1, deltaT/(2*pow(delta,2)));
            gsl_matrix_set(A, l, l+nx+1, deltaT/(2*pow(delta,2)));
            gsl_matrix_set(A, l, l, -((2*deltaT)/pow(delta,2))-1);
            //Wypelnienie B
            gsl_matrix_set(B, l, l-nx-1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l-1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l+1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l+nx+1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l, (2*deltaT)/pow(delta,2) - 1);
        }
    }
    //Warunki brzegowe Dirichleta
    //Lewy brzeg i=0
    for(int j=0;j<=ny;j++){
        l=0+j*(nx+1);
        //Macierz A
        gsl_matrix_set(A, l, l,1);
        //Macierz B
        gsl_matrix_set(B, l, l, 1);
        //Wektor c
        gsl_vector_set(c, l, 0);
    }
    //Prawy brzeg i=nx
    for(int j=0;j<=ny;j++){
        l=nx+j*(nx+1);
        //MacierzA
        gsl_matrix_set(A, l, l, 1);
        //MacierzB
        gsl_matrix_set(B, l, l, 1);
        //Wektor c
        gsl_vector_set(c, l, 0);
    }
    //WB von Neumanna na górnym brzegu
    for(int i=1;i<=nx-1;i++){
        l = i + ny*(nx+1);
        gsl_matrix_set(A, l, l-nx-1, - 1.0/(kB*delta)); 
        gsl_matrix_set(A, l, l, 1 + 1.0/(kB*delta)); 
        gsl_vector_set(c, l, TB);
        for(int k=0;k<N;k++){
            gsl_matrix_set(B, l, k, 0);
        }
    }
    //WB von Neumnanna na dolnym brzegu
    for(int i=1;i<=nx-1;i++){
        l = i + 0*(nx+1);
        gsl_matrix_set(A, l, l, 1 + 1.0/(kD*delta)); 
        gsl_matrix_set(A, l, l+nx+1, - 1.0/(kD*delta)); 
        gsl_vector_set(c, l, TD);
        for(int k=0;k<N;k++){
            gsl_matrix_set(B, l, k, 0);
        }
    }
    //Wektor startowy z warunkami początkowymi
    gsl_vector* T=gsl_vector_calloc(N);
    for(int i=0;i<=nx;i++){
        for(int j=0;j<=ny;j++){
            l=nx+j*(nx+1);
            if(i==0){
                gsl_vector_set(T,l,TA);
            }
            else if(i==nx){
                gsl_vector_set(T,l,TC);
            }
            else{
                gsl_vector_set(T,l,0);
            }
        }
    }
    //Pliki
    FILE* fileT;
    fileT=fopen("T.txt","w");
    FILE* filenT;
    filenT=fopen("nT.txt","w");


    //Rozklad LU
    gsl_linalg_LU_decomp(A,p,&signum);
    //Algorytm CN
    for(int it=0;it<=IT_MAX;it++){
        gsl_blas_dgemv( CblasNoTrans, 1, B, T, 0, d );
        gsl_blas_daxpy( 1, c, d );
        gsl_linalg_LU_solve( A, p, d, T );
        if( it==100 || it==200 || it==500 || it==1000 || it==2000 ){
            for(int i=0;i<N;i++ ){
                fprintf(fileT,"%g\t",gsl_vector_get(T, i));
    
                if((i+1)%(nx+1)==0){
                    fprintf(fileT, "\n");
                }   
            }

            for(int i=1;i<=nx-1;i++){
                for(int j=1;j<=ny-1;j++){
                    l=i+j*(nx+1);
                    fprintf(filenT,"%g\t",((gsl_vector_get(T,l+1) - 2*gsl_vector_get(T,l) + gsl_vector_get(T,l-1))/pow(delta,2))
                                      + ((gsl_vector_get(T,l+nx+1) - 2*gsl_vector_get(T,l) + gsl_vector_get(T,l-nx-1))/pow(delta,2)));
                }
                fprintf(filenT,"\n");
            }
            fprintf(filenT,"\n");
        }
    }
    fclose(filenT);
    fclose(fileT);

    

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(T);
    gsl_vector_free(c);
    gsl_permutation_free(p);
}

int main(){

    dyfuzjaMacierzowa();
    return 0;
}