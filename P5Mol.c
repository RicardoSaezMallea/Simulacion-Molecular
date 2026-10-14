#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define frand (rand()/((double) RAND_MAX+1))

#define PI 3.1415926
#define it_MAX 1000000
#define len 10
#define it_prom 1
#define it_TERM 1000

#define N 64

//#define TRAYECTORIA
//#define ENERGIA
//#define ENERGIA_ALT
#define ENERGIA_DISJUNTA
//#define FUERZA


double mod_diff(double X[3][N], int part1, int part2);

void V_muelle(double X[3][N], double *V, double k, double b, double F_const[3]);

double T_total(double P[3][N],double m);

double V_total(double X[3][N], double k, double b, double F_const[3]);

void F_part(double X[3][N], double F[3][N], double k, double b, double F_const[3]);

void BoxMuller(double *g1, double *g2);

double Z(double eta, double kT, double h);

void Verlet(double X_old[3][N], double X_new[3][N], double P_old[3][N], double P_new[3][N], double m, double k, double b_pot, double kT, double h, double eta, double F_const[3]);

void mover_vector(double T[], double var_nueva, int L);

void R_CM(double X[3][N], double RCM[]);

double R_G(double X[3][N]);

void R_ee(double X[3][N], double Ree[]);


#ifdef TRAYECTORIA
int main()
{
    //PRUEBA POTENCIAL

    /*double X[3][N],V[N], k=1, b=1;
    int i;

    for(i=0;i<3;i++)
    {
        X[0][i]=i;
        X[1][i]=i;
        X[2][i]=i;
    }

    V_muelle(X,V,k,b);

    printf("V[1] = %lf \n",V[1]);*/

    //PRUEBA FUERZAS

    /*double X[3][N],F[3][N], k=1, b=1;
    int i;

    for(i=0;i<3;i++)
    {
        X[0][i]=i;
        X[1][i]=i;
        X[2][i]=i;
    }

    F_part(X,F,k,b);
    printf("Fx= %lf Fy= %lf Fz= %lf \n",F[0][0],F[1][0],F[2][0]);*/


    double X_old[3][N], X_new[3][N], P_old[3][N], P_new[3][N], k=1, b=1, kT=1, h=0.01,eta=1, m=1, F_const[3];
    int i;
    long long int it;
    char filename [50];

    sprintf(filename,"Moleculas_Trayectoria_N%d_eta%.1f_h%.4f_k%.0f_b%.0f_it%d.txt",N,eta,h,k,b,it_MAX);
    printf("generando archivo %s\n",filename);

    FILE *fTrayectoria;
    fTrayectoria=fopen(filename,"w");

    //CONDICIONES INICIALES

    F_const[0] = 0;
    F_const[1] = 0;
    F_const[2] = 0;

    for(i=0; i<N; i++)
    {
        X_old[0][i]=i*2;
        X_old[1][i]=0;
        X_old[2][i]=0;
    }

    for(i=0; i<N; i++)
    {
        P_old[0][i]=0;
        P_old[1][i]=0;
        P_old[2][i]=0;
    }
    //END CONDICIONES INICIALES

    for(it=0; it<it_MAX; it++)
    {
        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta,F_const);

        //OBSERVABLES
        for(i=0; i<N; i++)
        {
            fprintf(fTrayectoria,"%lf %lf %lf %lf %lf %lf",X_old[0][i],X_old[1][i],X_old[2][i],P_old[0][i],P_old[1][i],P_old[2][i]);
        }
        fprintf(fTrayectoria,"\n");

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }
    }


}
#endif // TRAYECTORIA

#ifdef ENERGIA
int main()
{
    double X_old[3][N], X_new[3][N], P_old[3][N], P_new[3][N], T[len], V[len], k=1, b=1, kT=1, h=0.01,eta=1, m=1, suma_T=0, suma_V=0,  F_const[3];
    int i, n_prom, cont;
    long long int it;
    char filename [50];

    n_prom=it_MAX/it_prom;

    sprintf(filename,"Moleculas_Energia_N%d_eta%.1f_h%.4f_k%.0f_b%.0f_it%d.txt",N,eta,h,k,b,it_MAX);
    printf("generando archivo %s\n",filename);

    FILE *fEnergia;
    fEnergia=fopen(filename,"w");

    //CONDICIONES INICIALES

     F_const[0]=0;
     F_const[1]=0;
     F_const[2]=0;

    for(i=0; i<N; i++)
    {
        X_old[0][i]=i;
        X_old[1][i]=0;
        X_old[2][i]=0;
    }

    for(i=0; i<N; i++)
    {
        P_old[0][i]=0;
        P_old[1][i]=0;
        P_old[2][i]=0;
    }
    //END CONDICIONES INICIALES

    //TERMALIZACIÓN

    for(it=0; it<it_TERM; it++)
    {
        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta, F_const);

        //OBSERVABLES

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }
    }

    //INICIALIZACIÓN DE T V

    for(cont=0; cont<len; cont++)
    {
        T[cont]=T_total(P_old,m);
        V[cont]=V_total(X_old,k,b, F_const);
        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta,F_const);

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }
    }

    //PROMEDIO DE T V

    for (i=0;i<len;i++)
    {
        suma_T+=T[i];
        suma_V+=V[i];
    }

    suma_T=suma_T/len;
    suma_V=suma_V/len;

    fprintf(fEnergia,"%lf %lf\n",suma_T,suma_V);

    //MOVIMIENTO DE VENTANA
    for(int j=0; j<n_prom-1; j++)
    {
        for(it=0; it<it_prom; it++)
        {

            mover_vector(T, T_total(P_old, m), len);
            mover_vector(V, V_total(X_old, k,b, F_const), len);

            Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta,F_const);

            for(i=0; i<N; i++)
            {
                X_old[0][i]=X_new[0][i];
                X_old[1][i]=X_new[1][i];
                X_old[2][i]=X_new[2][i];

                P_old[0][i]=P_new[0][i];
                P_old[1][i]=P_new[1][i];
                P_old[2][i]=P_new[2][i];
            }
        }
        suma_T=0;
        suma_V=0;
        for (i=0;i<len;i++)
        {
            suma_T+=T[i];
            suma_V+=V[i];
        }

        suma_T=suma_T/len;
        suma_V=suma_V/len;

        fprintf(fEnergia,"%lf %lf \n",suma_T, suma_V);

        //Control
        if(j % ((n_prom-1) / 100) == 0)printf("Progreso: %d%%\n", (j * 100) / (n_prom-1));

    }


}

#endif // ENERGIA

#ifdef ENERGIA_ALT
int main()
{
    double X_old[3][N], X_new[3][N], P_old[3][N], P_new[3][N], T[len], V[len], k=1, b=1, kT=1, h=0.01,eta=1, m=1, suma_T=0, suma_V=0,  F_const[3];
    int i;
    long long int it;
    char filename [50];

    sprintf(filename,"Moleculas_Energia_N%d_eta%.1f_h%.4f_k%.0f_b%.0f_it%d.txt",N,eta,h,k,b,it_MAX);
    printf("generando archivo %s\n",filename);

    FILE *fEnergia;
    fEnergia=fopen(filename,"w");

    //CONDICIONES INICIALES

     F_const[0]=0;
     F_const[1]=0;
     F_const[2]=0;

    for(i=0; i<N; i++)
    {
        X_old[0][i]=i;
        X_old[1][i]=0;
        X_old[2][i]=0;
    }

    for(i=0; i<N; i++)
    {
        P_old[0][i]=0;
        P_old[1][i]=0;
        P_old[2][i]=0;
    }
    //END CONDICIONES INICIALES

    //TERMALIZACIÓN

    for(it=0; it<it_TERM; it++)
    {
        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta, F_const);

        //OBSERVABLES

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }
    }

    //INICIALIZACIÓN DE T V

    for(it=0; it<len; it++)
    {
        T[it]=T_total(P_old,m);
        V[it]=V_total(X_old,k,b, F_const);
        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta,F_const);

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }
    }

    //SIMULACIÓN

    for(it=0;it<it_MAX;it++)
    {
        if(it % it_prom == 0 && it< it_MAX-len)
        {
            //PROMEDIO DE T V
            suma_T=0;
            suma_V=0;

            for (i=0;i<len;i++)
            {
                suma_T+=T[i];
                suma_V+=V[i];
            }

            suma_T=suma_T/len;
            suma_V=suma_V/len;

            fprintf(fEnergia,"%d %lf %lf\n",it, suma_T,suma_V);
        }

        mover_vector(T, T_total(P_old, m), len);
        mover_vector(V, V_total(X_old, k,b, F_const), len);

        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta,F_const);

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }

        //Control
        if(it % ((it_MAX-len) / 100) == 0)printf("Progreso: %d%%\n", (it * 100) / (it_MAX-len));
    }


}
#endif // ENERGIA


#ifdef FUERZA
int main()
{
    double X_old[3][N], X_new[3][N], P_old[3][N], P_new[3][N], k=100, b=1, kT=1, h=0.01,eta=1, m=1, F_const[3], F_MIN=0, F_MAX=20, Ree[3], R_ee_x=0;
    int i, it_F, n_F=20;
    long long int it;
    char filename [50];

    sprintf(filename,"Moleculas_Fuerzas_Fmax%.0f_N%d_eta%.1f_h%.4f_k%.0f_b%.0f_it%d.txt",F_MAX, N,eta,h,k,b,it_MAX);
    printf("generando archivo %s\n",filename);

    FILE *fFuerza;
    fFuerza=fopen(filename,"w");

    //CONDICIONES INICIALES
    for(i=0; i<N; i++)
    {
        X_old[0][i]=i; // mas un random entre 0.05
        X_old[1][i]=0;
        X_old[2][i]=0;
    }

    for(i=0; i<N; i++)
    {
        P_old[0][i]=0;
        P_old[1][i]=0;
        P_old[2][i]=0;
    }
    //END CONDICIONES INICIALES

    for(it_F=0;it_F<n_F;it_F++) // F hasta 15 20
    {
        printf("bucle %d \n", it_F);
        F_const[0] = F_MIN+(double)it_F/n_F*(F_MAX-F_MIN);
        F_const[1] = 0;
        F_const[2] = 0;

        printf("Calculando F= %lf \n", F_const[0]);

        //TERMALIZACIÓN

        for(it=0; it<it_TERM; it++)
        {
            Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta, F_const);
            //OBSERVABLES

            for(i=0; i<N; i++)
            {
                X_old[0][i]=X_new[0][i];
                X_old[1][i]=X_new[1][i];
                X_old[2][i]=X_new[2][i];

                P_old[0][i]=P_new[0][i];
                P_old[1][i]=P_new[1][i];
                P_old[2][i]=P_new[2][i];
            }
        }

        //SIMULACIÓN DE R_ee

        for(it=0; it<it_MAX; it++)
        {
            Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta, F_const);

            //OBSERVABLES
            R_ee(X_old,Ree);
            R_ee_x += Ree[0];
            //END OBSERVABLES

            for(i=0; i<N; i++)
            {
                X_old[0][i]=X_new[0][i];
                X_old[1][i]=X_new[1][i];
                X_old[2][i]=X_new[2][i];

                P_old[0][i]=P_new[0][i];
                P_old[1][i]=P_new[1][i];
                P_old[2][i]=P_new[2][i];
            }
        }
        R_ee_x = R_ee_x/it_MAX;

        printf("Ree %lf\n" , R_ee_x);
        fprintf(fFuerza,"%lf %lf \n",F_const[0],R_ee_x);

        R_ee_x=0;

        /*//Control
        if(it_F % ( n_F / 100) == 0)printf("Progreso: %d%%\n", (it_F * 100) / n_F );
        */

    }

}
#endif // FUERZA

#ifdef ENERGIA_DISJUNTA
int main()
{
    double X_old[3][N], X_new[3][N], P_old[3][N], P_new[3][N], k=100, b=1, kT=1, h=0.01,eta=1, m=1, suma_T=0, suma_V=0,  F_const[3];
    int i, N_ventana=10, ind_ventana=0;
    long long int it;
    char filename [50];

    double T_ventana[N_ventana], V_ventana[N_ventana];
    {
    sprintf(filename,"Moleculas_Energia_N%d_eta%.1f_h%.4f_k%.0f_b%.0f_it%d.txt",N,eta,h,k,b,it_MAX);
    printf("generando archivo %s\n",filename);

    FILE *fEnergia;
    fEnergia=fopen(filename,"w");

    //CONDICIONES INICIALES

     F_const[0]=0;
     F_const[1]=0;
     F_const[2]=0;

    for(i=0; i<N; i++)
    {
        X_old[0][i]=i;
        X_old[1][i]=0;
        X_old[2][i]=0;
    }

    for(i=0; i<N; i++)
    {
        P_old[0][i]=0;
        P_old[1][i]=0;
        P_old[2][i]=0;
    }
    //END CONDICIONES INICIALES

    for (it=0; it < it_MAX; it ++)
    {
        Verlet(X_old,X_new,P_old,P_new,m,k,b,kT,h,eta, F_const);

        // OBSERVABLES
        fprintf(fEnergia,"%lf %lf \n",T_total(P_old,m),V_total(X_old,k,b,F_const));
        //END OBSERVABLES

        for(i=0; i<N; i++)
        {
            X_old[0][i]=X_new[0][i];
            X_old[1][i]=X_new[1][i];
            X_old[2][i]=X_new[2][i];

            P_old[0][i]=P_new[0][i];
            P_old[1][i]=P_new[1][i];
            P_old[2][i]=P_new[2][i];
        }

        //Control
        if(it % ((it_MAX) / 100) == 0)printf("Progreso: %d%%\n", (it * 100) / (it_MAX));
    }
    }
}
#endif // ENERGIA_DISJUNTA


double mod_diff(double X[3][N], int part1, int part2)
{
    return sqrt( (X[0][part1]-X[0][part2])*(X[0][part1]-X[0][part2]) + (X[1][part1]-X[1][part2])*(X[1][part1]-X[1][part2]) + (X[2][part1]-X[2][part2])*(X[2][part1]-X[2][part2]) );
}

void V_muelle(double X[3][N], double *V, double k, double b, double F_const[3])
{
    //FUNCIÓN ACTUALIZADA
    /*
        Devuelve un vector V de tamaño N-1, en lugar de N, porque realmente son los N-1 muelles los que almacenan el potencial elástico, es decir,
        son los delta_x_i=x_i-x_i+1 los grados de libertad.

        Distinguimos los extremos para añadir la fuerza constante F_const[3]

    */
    double mod_diff_dch=0;
    int muelle;

    for (muelle=0; muelle<N-1; muelle++)
    {
        mod_diff_dch=mod_diff(X,muelle,muelle+1);

        V[muelle]=0.5*k*(mod_diff_dch-b)*(mod_diff_dch-b);
        //printf("Potencial del muelle %d : %lf \n",muelle, V[muelle]);
    }

}

void F_part(double X[3][N], double F[3][N], double k, double b,  double F_const[3])
{
    double mod_diff_izq=0, mod_diff_dch=0;
    int part;
    //Comprobamos si estamos en uno de los extremos de la cadena. Calculamos la contribución al potencial de la partícula i de sus vecinas.

    //CALCULO LAS DE LOS EXTREMOS

    //EXTREMO IZQ
    mod_diff_dch = mod_diff(X,0,1);

    F[0][0] = -k*(mod_diff_dch-b)*(X[0][0]-X[0][1])/mod_diff_dch - (F_const[0]);
    F[1][0] = -k*(mod_diff_dch-b)*(X[1][0]-X[1][1])/mod_diff_dch - (F_const[1]);
    F[2][0] = -k*(mod_diff_dch-b)*(X[2][0]-X[2][1])/mod_diff_dch - (F_const[2]);

    //EXTREMO DCHA
    mod_diff_izq = mod_diff(X,N-1,N-2);

    F[0][N-1] = k*(mod_diff_izq-b)*(X[0][N-2]-X[0][N-1])/mod_diff_izq + (F_const[0]);
    F[1][N-1] = k*(mod_diff_izq-b)*(X[1][N-2]-X[1][N-1])/mod_diff_izq + (F_const[1]);
    F[2][N-1] = k*(mod_diff_izq-b)*(X[2][N-2]-X[2][N-1])/mod_diff_izq + (F_const[2]);


    for(part=1; part<N-1; part++)
    {
        mod_diff_dch = mod_diff(X,part,part+1);
        mod_diff_izq = mod_diff(X,part,part-1);

        F[0][part] = -k*(mod_diff_dch-b)*(X[0][part]-X[0][part+1])/mod_diff_dch + k*(mod_diff_izq-b)*(X[0][part-1]-X[0][part])/mod_diff_izq;
        F[1][part] = -k*(mod_diff_dch-b)*(X[1][part]-X[1][part+1])/mod_diff_dch + k*(mod_diff_izq-b)*(X[1][part-1]-X[1][part])/mod_diff_izq;
        F[2][part] = -k*(mod_diff_dch-b)*(X[2][part]-X[2][part+1])/mod_diff_dch + k*(mod_diff_izq-b)*(X[2][part-1]-X[2][part])/mod_diff_izq;
    }

}

void Verlet(double X_old[3][N], double X_new[3][N], double P_old[3][N], double P_new[3][N], double m, double k, double b_pot, double kT, double h,double eta, double F_const[3])
{
    double a,b, beta[3][N], F_old[3][N], F_new[3][N];

    int i,j;

    a=(1-eta*h/(2*m))/(1+eta*h/(2*m));
    b=1/(1+eta*h/(2*m));

    F_part(X_old,F_old,k,b_pot, F_const); //ineficiente



    for(i=0; i<3; i++)
    {
        for(j=0; j<N; j++)
        {
            beta[i][j]=Z(eta,kT,h);
        }
    }

    for(i=0; i<N; i++)
    {
        X_new[0][i] = X_old[0][i] + P_old[0][i]*b*h/m+b*h*h*F_old[0][i]/(2*m)+b*h/(2*m)*beta[0][i];
        X_new[1][i] = X_old[1][i] + P_old[1][i]*b*h/m+b*h*h*F_old[1][i]/(2*m)+b*h/(2*m)*beta[1][i];
        X_new[2][i] = X_old[2][i] + P_old[2][i]*b*h/m+b*h*h*F_old[2][i]/(2*m)+b*h/(2*m)*beta[2][i];
    }

    F_part(X_new,F_new,k,b_pot, F_const);

    for(i=0; i<N; i++)
    {
        P_new[0][i] = m * (P_old[0][i]*a/m + h/(2*m)*(a*F_old[0][i]+F_new[0][i]) + b/m*beta[0][i]);
        P_new[1][i] = m * (P_old[1][i]*a/m + h/(2*m)*(a*F_old[1][i]+F_new[1][i]) + b/m*beta[1][i]);
        P_new[2][i] = m * (P_old[2][i]*a/m + h/(2*m)*(a*F_old[2][i]+F_new[2][i]) + b/m*beta[2][i]);
    }


}

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

double Z(double eta, double kT, double h)
{
    double g1,g2;
    BoxMuller(&g1,&g2);
    return sqrt(2*eta*kT*h)*g1;
}

double T_total(double P[3][N],double m)
{
    double T=0;
    int i;

    for(i=0;i<N;i++)
    {
        T += P[0][i]*P[0][i]/(2*m) + P[1][i]*P[1][i]/(2*m) + P[1][i]*P[1][i]/(2*m);
    }

    return T;
}

double V_total(double X[3][N], double k, double b, double F_const[3])
{
    double V[N-1], V_tot=0;
    int i;

    V_muelle(X,V,k,b,F_const);

    for(i=0;i<N-1;i++)
    {
        V_tot += V[i];
    }
    //printf("Potencial total calculado para la cadena : %lf \n",V_tot);
    return V_tot;
}

void mover_vector(double T[], double var_nueva, int L)
{
    int i;

    for (i=0;i<L-1;i++)
    {
        T[i]=T[i+1];
    }

    T[L-1]=var_nueva;
}

void R_CM(double X[3][N], double RCM[])
{
    //INICIALIZACIÓN A 0

    for(int j=0;j<3;j++)
    {
        RCM[j]=0;
    }

    //CÁLCULO DEL CDM

    for(int i=0;i<N;i++)
    {
        RCM[0]+=X[0][i];
        RCM[1]+=X[1][i];
        RCM[2]+=X[2][i];
    }

    for(int j=0;j<3;j++)
    {
        RCM[j]=RCM[j]/N;
    }

}

double R_G(double X[3][N])
{
    double RCM[3], aux=0;

    R_CM(X, RCM);

    for(int i=0; i<N; i++)
    {
        aux+=(X[0][i]-RCM[0])*(X[0][i]-RCM[0])+(X[1][i]-RCM[1])*(X[1][i]-RCM[1])+(X[2][i]-RCM[2])*(X[2][i]-RCM[2]);
    }

    aux=sqrt(aux/N);

    return aux;

}

void R_ee(double X[3][N], double Ree[])
{
    Ree[0]=X[0][N-1]-X[0][0];
    Ree[1]=X[1][N-1]-X[1][0];
    Ree[2]=X[2][N-1]-X[2][0];
}

