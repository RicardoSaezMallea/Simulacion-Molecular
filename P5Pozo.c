#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define frand (rand()/((double) RAND_MAX+1))

#define PI 3.1415926
#define it_MAX 10000000
#define it_TERM 0

#define ENERGIAS
//#define TRAYECTORIAS
//#define TIEMPOS
//#define DOS_ESTADOS


void BoxMuller(double *g1, double *g2);
double f(double p, double m);
double gVer(double x, double A, double F);
double Z(double eta, double kT, double h);
double E_pot(double x, double A,  double F);
double E_cin(double p, double m);
void Verlet (double *x_old, double *x_new, double *p_old, double *p_new, double m, double A, double F, double eta, double kT, double h);
void comparar (double x_new, double x_pozo, int estado_old, int *estado_new);

#ifdef ENERGIAS
int main()
{
    double p_old, x_old, p_new, x_new, suma_T, suma_V, T, V, prom_T, prom_V, prom_Etot, h=0.01, eta=0.01, A=0.5, F=0, m=1, kT=0.2;
    int it=0;
    char filename[50];
    FILE *f;

    for(A=0.5; A<3; A=A*2){
    printf("Empiezo A=%lf\n", A);
    eta = 0.01;
    while(eta<6){
    printf("Empiezo eta=%lf\n", eta);
    sprintf(filename,"EnergiasDP_A%.1f_eta%.2f_h%.4f_it%d.txt",A,eta,h,it_MAX);
    f = fopen(filename, "w");
    if (f == NULL)
    {
        printf("Error al abrir el archivo\n");
        return;
    }
    else
    {
        p_old = 0;
        x_old = 1;
        suma_T = 0;
        suma_V = 0;
        //fprintf(f, "Datos introducidos: p0=%lf, x0=%lf, k=%lf, m=%lf, eta=%lf, kT=%lf, h=%lf, nIt=%d\n",p0,x0,k,m,eta,kT,h,nIt);
        //fprintf(f, "t\t p_old\t x_old\t EKin\t EV\t Etot\n");
        for(it = 0; it < it_MAX; it++)
        {

            Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);

            //OBSERVABLES
            T = E_cin(p_old, m);
            V = E_pot(x_old, A,F);
            suma_T += T;
            suma_V += V;

            prom_T = suma_T/(it+1);
            prom_V = suma_V/(it+1);
            prom_Etot = prom_V + prom_T;
            fprintf(f,"%lf\t %lf\t %lf\t %lf\n", h*it, prom_T, prom_V, prom_Etot);
            //END OBSERVABLES

            x_old = x_new;
            p_old = p_new;

            //Control
            if(it % (it_MAX / 10) == 0)printf("Progreso: %d%%\n", (it * 100) / it_MAX);

        }
    }
        //fprintf(f, "El promedio de T es: %lf y el promedio de V es: %lf\n", sumaKin / nIt, sumaV / nIt);
        fclose(f);

        if(eta == 0.01) eta = 1;
        else eta = eta+1;
}
}
}

#endif // ENERGIAS

#ifdef TRAYECTORIAS
int main()
{
    double x_old, x_new, p_old, p_new, h=0.01, eta=1, A=1, F=0, m=1, kT=0.2;
    int it,it_h,it_eta,n_h=5,n_eta=4;
    char filename[50];

    for(it_h=1;it_h<n_h;it_h++)
    {
        h=pow(0.1,it_h);
        for(it_eta=1;it_eta<n_eta;it_eta++)
        {
            eta=pow(10,it_eta-2);
            sprintf(filename,"TrayectoriaDP_eta%.1lf_h%.4lf_it%d.txt",eta,h,it_MAX);
            FILE *fTrayectoria;
            fTrayectoria=fopen(filename,"w");
            printf("Generando %s\n",filename);

            //CONDICIONES INICIALES
            x_old=0;
            p_old=0.1;
            //END CONDICIONES INICIALES

            //TERMALIZACION
            for (it=0;it<it_TERM;it++)
            {
                Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);

                x_old=x_new;
                p_old=p_new;

                //Control
                if(it % (it_TERM / 10) == 0)printf("Progreso Termalizacion: %d%%\n", (it * 100) / it_TERM);
            }
            //END TERMALIZACION

            for (it=0;it<it_MAX;it++)
            {

                Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);


                // OBSERVABLES
                fprintf(fTrayectoria,"%lf %lf %lf \n",it*h,x_old,p_old);
                // END OBSERVABLES

                x_old=x_new;
                p_old=p_new;

                        //Control
                if(it % (it_MAX / 10) == 0)printf("Progreso: %d%%\n", (it * 100) / it_MAX);
            }
        }
    }
}
#endif // TRAYECTORIAS

#ifdef TIEMPOS
int main()
{
    double x_old, x_new, p_old, p_new, h=0.01, eta=1,eta_MIN=0.05,eta_MAX=1, A=1, F=0, m=1, kT=0.2, x_pozo=0.2, aux;
    int estado_old, estado_new;
    long long int it, it_A, it_eta;
    int contador = 0, n_eta=20;
    char filename[50];

    for(it_A=0;it_A<1;it_A++)
    {
        A=pow(0.5,(1-it_A));

        for(it_eta=0;it_eta<=n_eta;it_eta++)
        {
            eta=eta_MIN+(double)it_eta/n_eta*(eta_MAX-eta_MIN);

            FILE *fTiempos;
            sprintf(filename,"TiemposDP_A%.1f_eta%.4f_h%.4f_it%d.txt",A,eta,h,it_MAX);
            fTiempos=fopen(filename,"w");
            printf("Generando %s\n",filename);

            //CONDICIONES INICIALES
            contador=0;
            x_old=1;
            p_old=0.1;
            estado_old=1;
            //END CONDICIONES INICIALES

            //TERMALIZACION
            for (it=0;it<it_TERM;it++)
            {
                Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);

                x_old=x_new;
                p_old=p_new;

                //Control
                if(it % (it_TERM / 10) == 0)printf("Progreso Termalizacion: %d%%\n", (it * 100) / it_TERM);
            }
            //END TERMALIZACION

            for (it=0;it<it_MAX;it++)
            {

                Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);


                // OBSERVABLES
                comparar(x_new,x_pozo,estado_old,&estado_new);

                if(estado_new == estado_old)contador++;
                else
                {
                    aux = contador*h;
                    fprintf(fTiempos,"%lf\n", aux);
                    contador = 0;
                    estado_old = estado_new;
                }
                //fprintf(fTrayectoria,"%lf %lf %lf \n",it*h,x_old,p_old);
                // END OBSERVABLES

                x_old=x_new;
                p_old=p_new;

                        //Control
                if(it % (it_MAX / 100) == 0)printf("Progreso: %.1f%%\n",((double)it * 100) /it_MAX);
            }

            fclose(fTiempos);
        }
    }


}
#endif // TIEMPOS

#ifdef DOS_ESTADOS
int main()
{
    double x_old, x_new, p_old, p_new, h=0.01, eta=1, A=1, F_MAX=1,F, m=1, kT=0.2, x_pozo=0.2, aux;
    int estado_old, estado_new, n_F=20;
    long long int it,it_F, contador = 0;
    char filename[50];

    FILE *fDosEstados;
    sprintf(filename,"DosEstadosDP_A%.1f_F%.1f_eta%.2f_h%.4f_it%d.txt",A,F_MAX,eta,h,it_MAX);
    fDosEstados=fopen(filename,"w");
    for(it_F=0;it_F<=n_F;it_F++)
    {
        contador=0;
        F=-F_MAX+(double)it_F/n_F*2*F_MAX;
        printf("F: %lf\n",F);

        //CONDICIONES INICIALES
        x_old=1;
        p_old=0.1;
        estado_old=1;
        //END CONDICIONES INICIALES

        //TERMALIZACION
        for (it=0;it<it_TERM;it++)
        {
            Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);

            x_old=x_new;
            p_old=p_new;

            //Control
            if(it % (it_TERM / 10) == 0)printf("Progreso Termalizacion: %d%%\n", (it * 100) / it_TERM);
        }
        //END TERMALIZACION

        for (it=0;it<it_MAX;it++)
        {

            Verlet(&x_old, &x_new, &p_old, &p_new, m, A, F, eta, kT, h);


            // OBSERVABLES
            if(x_old>=0)contador++;
            // END OBSERVABLES

            x_old=x_new;
            p_old=p_new;

                    //Control
            if(it % (it_MAX / 10) == 0)printf("Progreso: %.1f%%\n",((double)it * 100) /it_MAX);
        }
        fprintf(fDosEstados,"%lf %lf %lf \n",F,(double)contador/it_MAX,1-(double)contador/it_MAX);

    }

    fclose(fDosEstados);
}
#endif // DOS_ESTADOS


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

double gVer(double x, double A, double F)
{
    return -4*A*x*(x*x-1)+F;
}

double Z(double eta, double kT, double h)
{
    double g1,g2;
    BoxMuller(&g1,&g2);
    return sqrt(2*eta*kT*h)*g1;
}

double E_pot(double x, double A, double F)
{
    return A*(x*x-1)*(x*x-1)-F*x ;
}

double E_cin(double p, double m)
{
    return p*p/(2*m);
}

void Verlet (double *x_old, double *x_new, double *p_old, double *p_new, double m, double A, double F, double eta, double kT, double h)
{
    double a,b, beta;

    a=(1-eta*h/(2*m))/(1+eta*h/(2*m));
    b=1/(1+eta*h/(2*m));
    beta=Z(eta,kT,h);


    *x_new=*x_old+*p_old*b*h/m+b*h*h*gVer(*x_old,A,F)/(2*m)+b*h/(2*m)*beta;
    *p_new=m*(* p_old/m*a+h/(2*m)*(a*gVer(*x_old,A,F)+gVer(*x_new,A,F))+b/m*beta);
}

void comparar (double x_new, double x_pozo, int estado_old, int *estado_new)
{
    if(estado_old == 1 && x_new < -x_pozo) *estado_new = 0;
    else if(estado_old == 0 && x_new > x_pozo) *estado_new = 1;
    else *estado_new = estado_old;
}
