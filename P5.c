#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define frand (rand()/((double) RAND_MAX+1))

#define PI 3.1415926
#define N_Data_MAX 500000000
#define it_TERM 3141592

//#define PRUEBAS
#define ENERGIAS
//#define TRAYECTORIAS
//#define DISTRIBUCIONES

//#define DEBUG


void BoxMuller(double *g1, double *g2);

double f(double p, double m);

double g(double x, double p, double m, double k, double eta);

double Z(double eta, double kT, double h);

void Euler(double *x_old, double *x_new, double *p_old, double *p_new, double m, double k, double eta, double kT, double h);

void Histograma(double *data,double *Hist, int N_data,int N_Intervalos, double *d,double *m, double *M);

double E_pot(double x, double k);

double E_cin(double p, double m);

void guardaTrayectoria(double x[], double p[], double h, double N_it, char nombre[]);

void guardaDistribucion(double x[], double p[], int N_Trayectoria, char nombre[]);

double gVer(double x, double k);

void Verlet (double *x_old, double *x_new, double *p_old, double *p_new, double m, double k, double eta, double kT, double h);

void RK(double *x_old, double *x_new, double *p_old, double *p_new, double m, double k, double eta, double kT, double h);

#ifdef PRUEBAS
int main()
{

 // PRUEBA DEL BOX MULLER
     double Data_example[N_Data_MAX],Hist[N_Data_MAX],aux1,aux2,min,max,dist;
     int i;

     FILE *fBoxMuller;

     fBoxMuller=fopen("Data_BoxMuller.txt","w");

     for(i=0;i<N_Data_MAX-1;i+=2)
     {
         BoxMuller(&aux1,&aux2);
         Data_example[i]=aux1;
         Data_example[i+1]=aux2;
     }

     for(i=0;i<N_Data_MAX;i++)
     {
         //printf("%d: %lf \n",i,Data_example[i]);
     }

     Histograma(Data_example,Hist,N_Data_MAX,20,&dist,&min,&max);

     printf("Min: %lf Max:%lf \n",min,max);

     for(i=0;i<N_Data_MAX;i++)
     {
         fprintf(fBoxMuller,"%lf \n",Data_example[i]);
     }
 // END PRUEBA DEL BOX MULLER

}
#endif // PRUEBAS

#ifdef ENERGIAS

int main()
{
    double x_old, x_new, p_old, p_new, m=1, k=1, eta=1, kT=1, h=0.001, Epot, Ecin;
    long long int it, nIt;

    FILE *fEnergias;
    //fEnergias=fopen("Energias_Euler_eta0_h0.0001.txt","w");

    nIt=2*PI*900/h;

    //CONDICIONES INICIALES
    //x_old=0;
    //p_old=1;
    //END CONDICIONES INICIALES
    /*
    for (it=0;it<it_MAX;it++)
    {
        Euler(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
        //Verlet();
        //RK4();

        // OBSERVABLES
        Epot=E_pot(x_old,k);
        Ecin=E_cin(p_old,m);

        fprintf(fEnergias,"%d %lf %lf %lf \n",it,Epot,Ecin,Epot+Ecin);
        // END OBSERVABLES

        x_old=x_new;
        p_old=p_new;

        //Control
        if(it % (it_MAX / 10) == 0)printf("Progreso: %d%%\n", (it * 100) / it_MAX);
    }
    */
    //guardaEnergia();
    equiparticion(1.0,0.0,k,m,eta,kT,h,nIt);
    //fclose(fEnergias);
}

#endif // ENERGIAS

#ifdef TRAYECTORIAS
int main()
{
    double x_old, x_new, p_old, p_new, m=1, k=1, eta=1, kT=1, h=0.1;
    int it_h, it_eta,n_data=10000000;
    long long int it, it_MAX;
    char filename[50];

    for (it_h=1; it_h<5; it_h++)
    {
        h=pow(0.1,it_h);
        it_MAX=n_data*pow(10,(it_h-1));

        for(it_eta=1; it_eta<4;it_eta++)
        {
            eta=pow(10,it_eta-2);

            FILE *fTrayectoria;
            sprintf(filename,"Trayectoria_RK_eta%.1lf_h%.4lf_data%d.txt", eta, h,n_data);

            fTrayectoria=fopen(filename,"w");
            printf("Generando %s\n",filename);

            //CONDICIONES INICIALES
            x_old=0;
            p_old=1;
            //END CONDICIONES INICIALES

            //TERMALIZACION
            for (it=0;it<it_TERM;it++)
            {
                //Euler(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
                //Verlet(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
                RK(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);



                x_old=x_new;
                p_old=p_new;

                        //Control
                if(it % (it_TERM / 10) == 0)printf("Progreso Termalizacion: %d%%\n", (it * 100) / it_TERM);
            }
            //END TERMALIZACION

            for (it=0;it<it_MAX;it++)
            {
                //Euler(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
                //Verlet(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
                RK(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);

                //OBSERVABLES
                if (fmod((double)it,pow(10,(it_h-1)))==0)
                {
                    fprintf(fTrayectoria,"%lf %lf %lf \n",it*h,x_old,p_old);
                }
                // END OBSERVABLES

                x_old=x_new;
                p_old=p_new;

                        //Control
                if(it % (it_MAX / 10) == 0)printf("Progreso: %d%%\n", (it * 100) / it_MAX);
            }

            fclose(fTrayectoria);
        }
    }


}
#endif // TRAYECTORIAS

#ifdef DISTRIBUCIONES
int main()
{
    double x_old, x_new, p_old, p_new, m=1, k=1, eta=0.1, kT=1, h=0.01;
    int it,n_Trayectoria;

    FILE *fDistribucion;

    fDistribucion=fopen("Distribucion_Euler_eta0.1_h0.01.txt","w");

    //CONDICIONES INICIALES
    x_old=0;
    p_old=1;
    //END CONDICIONES INICIALES

    for (it=0;it<it_MAX;it++)
    {
        Euler(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
        //Verlet(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
        //RK(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);

        // OBSERVABLES
        // END OBSERVABLES

        x_old=x_new;
        p_old=p_new;
    }
    fprintf(fDistribucion,"%d %lf %lf \n",n_Trayectoria,x_old,p_old);

    //Control
    if(n_Trayectoria % (N_Trayectoria_MAX /100) == 0)printf("Progreso: %d%%\n", (n_Trayectoria * 100) / N_Trayectoria_MAX);

//guardaDistribucion(x,p,N_Trayectoria_MAX,"Distribuciones_Euler_eta1_h0.001_it1000.txt");


}
#endif // DISTRIBUCIONES

void BoxMuller(double *g1, double *g2)
{
    double r1 = frand;
    double r2 = frand;

    while(r1==0 || r2==0)
    {
        r1 = frand;
        r2 = frand;
    }

    *g1 = -sqrt(-2*log(r1))*cos(2*PI*r2);
    *g2 = -sqrt(-2*log(r1))*sin(2*PI*r2);

}

double f(double p, double m)
{
    return p/m;
}

double g(double x, double p, double m, double k, double eta)
{
    return -eta*p/m-k*x;
}

double Z(double eta, double kT, double h)
{
    double g1,g2;
    BoxMuller(&g1,&g2);
    return sqrt(2*eta*kT*h)*g1;
}

void Euler(double *x_old, double *x_new, double *p_old, double *p_new, double m, double k, double eta, double kT, double h)
{
    *x_new=*x_old + f(*p_old,m)*h;
    *p_new=*p_old + g(*x_old,*p_old,m,k,eta)*h+Z(eta,kT,h);

}

void Histograma(double *data,double *Hist, int N_data,int N_Intervalos, double *d,double *m, double *M)
{
    /* Genera un histograma. Calcula, en funcion de los datos el minimo y maximo
    de los mismos para ajustar mejor los intervalos.
    *data-> input Datos sobre los que se genera el histograma
    *Hist-> output Histograma calculado
    N_data-> input Numero de datos
    N_Intervalos-> input Numero de intervalos del histograma
    *d -> output Medida de cada intervalo del histograma
    *m -> output Valor minimo de los datos
    *M -> output Valor maximo de los datos
    */
    int i,Indice;
    double del,min,max,Norm;
    for(i=0;i<N_Intervalos;i++)//Inicializo
    Hist[i]=0;
    min=10000000;
    max=-10000000;
    for(i=0;i<N_data;i++) //Calculo el minimo y el maximo
    {
    if(data[i]<min)min=data[i];
    if(data[i]>max)max=data[i];
    }
    del=(max-min)/N_Intervalos;
    if(del==0)
    {
    printf("Error: No se pueden calcular los intervalos; Max=%lf,Min=%lf\n",max,min);
    exit(1);
    }
    for(i=0;i<N_data;i++) //*********** Calculo el histograma: Nucleo del Programa**********
    {
    Indice=(data[i]-min)/del;
    if(Indice==N_Intervalos)
    Indice=N_Intervalos-1;
    Hist[Indice]++;
    #ifdef DEBUG
    printf("x=%f, H[%d]=%f\n",Indice*del+min,Indice,Hist[Indice]);
    #endif
    }
    // ************** Fin del nucleo del programa ***************
    *d=del;
    *m=min;
    *M=max;
    /* Ahora normalizo */
    // Recordar: 1=A*sum(h_i * delta_i)=A*delta*sum(h_i)=A*delta*N => A=1/(delta*N)
    Norm=1.0/(N_data*del);
    for(i=0;i<N_Intervalos;i++)
    Hist[i]*=Norm;
}

double E_pot(double x, double k)
{
    return 0.5*k*x*x;
}

double E_cin(double p, double m)
{
    return p*p/(2*m);
}

void guardaTrayectoria(double x[], double p[], double h, double N_it, char nombre[])
{
    int i;
    double t;
    FILE *fTrayectoria;

    fTrayectoria=fopen(nombre,"w");

    for(i=0;i<N_it;i++)
    {
        t=i*h;
        fprintf(fTrayectoria,"%lf %lf %lf \n",t,x[i],p[i]);
    }

    fclose(fTrayectoria);

}

void guardaDistribucion(double x[], double p[], int N_Trayectoria, char nombre[])
{
    int i;
    FILE *fDistribucion;

    fDistribucion=fopen(nombre,"w");

    for(i=0;i<N_Trayectoria;i++)
    {
        fprintf(fDistribucion,"%d %lf %lf \n",i,x[i],p[i]);
    }

    fclose(fDistribucion);

}


void RK(double *x_old, double *x_new, double *p_old, double *p_new, double m, double k, double eta, double kT, double h)
{
    double z = Z(eta,kT,h);
    double aux1,aux2,aux3;
    aux1= *p_old + z;
    double f1 = f(aux1,m);
    double g1 = g(*x_old,aux1,m,k,eta);

    aux2 = *p_old + h*g1;
    aux3 = *x_old+h*f1;

    double f2 = f(aux2,m);
    double g2 = g(aux3,aux2,m,k,eta);

    *x_new = *x_old + 0.5*h *(f1+f2);
    *p_new = *p_old + 0.5*h *(g1+g2)+z;
}

void equiparticion (double p0, double x0, double k, double m, double eta, double kT, double h, int nIt){
    double p_old, x_old, p_new, x_new, sumaKin, sumaV, T, V, promKin, promV, promTV, r;
    long long int i = 0;
    FILE *f;
    f = fopen("EnergiasEuler_eta0.1_h0.1.txt", "w");
    if (f == NULL) {
        printf("Error al abrir el archivo\n");
        return;
    }
    else{
        p_old = p0;
        x_old = x0;
        sumaKin = 0;
        sumaV = 0;
        //fprintf(f, "Datos introducidos: p0=%lf, x0=%lf, k=%lf, m=%lf, eta=%lf, kT=%lf, h=%lf, nIt=%d\n",p0,x0,k,m,eta,kT,h,nIt);
        //fprintf(f, "t\t p_old\t x_old\t EKin\t EV\t Etot\n");
        for(i = 0; i < nIt; i++){
            T = E_cin(p_old, m);
            V = E_pot(x_old, k);
            //fprintf(f, "%lf\t %lf \t%lf \t%lf \t%lf \t%lf \n", h*i, p_old, x_old, T, V, T+V);
            //Euler(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
            //RK(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
            Verlet(&x_old, &x_new, &p_old, &p_new, m, k, eta, kT, h);
            sumaKin += T;
            sumaV += V;
            x_old = x_new;
            p_old = p_new;
            promKin = sumaKin/(i+1);
            promV = sumaV/(i+1);
            promTV = promV + promKin;
            fprintf(f,"%lf\t %lf\t %lf\t %lf\n", h*i, promKin, promV, promTV);

            if(i % (nIt / 10) == 0)printf("Progreso: %d%%\n", (i * 100) / nIt);
            }
        }
        //fprintf(f, "El promedio de T es: %lf y el promedio de V es: %lf\n", sumaKin / nIt, sumaV / nIt);
        fclose(f);
    }


void Verlet (double *x_old, double *x_new, double *p_old, double *p_new, double m, double k, double eta, double kT, double h)
{
    double a,b, beta;

    a=(1-eta*h/(2*m))/(1+eta*h/(2*m));
    b=1/(1+eta*h/(2*m));
    beta=Z(eta,kT,h);

    *x_new=*x_old+*p_old*b*h/m+b*h*h*gVer(*x_old,k)/(2*m)+b*h/(2*m)*beta;
    *p_new=m*(* p_old/m*a+h/(2*m)*(a*gVer(*x_old,k)+gVer(*x_new,k))+b/m*beta);
}

double gVer(double x, double k)
{
    return -k*x;
}

