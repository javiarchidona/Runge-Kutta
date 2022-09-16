#include<iostream>
#include<fstream>
#include<cmath>


#define G 6.67e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define dTL 3.844e8
#define w 2.6617e-6
#define RT 6.378160e6
#define RL 1.7374e6

#define theta 0.52   //Angulo en radianes (0º)
#define h 1.0
#define v 11150.0   //Velocidad inicial cohete m/s (VELOCIDAD orbital)

#define N 100000     //ITERACIONES
//comprobar que lanzando perp desde polo norte con vel orbital de radio terrestre se queda girando


double f(double r, double phi, double pr, double p_phi, double t, int i);


using namespace std;


int main(void){


    int i, j, k, cont;
    double cohete[4], luna[2], phi, p_phi, p_r, t, r, r_luna, H;
    double k1[4], k2[4], k3[4], k4[4];
    ofstream co, lu, ham;

    co.open("datos_cohete.txt");
    lu.open("datos_luna.txt");
    ham.open("hamiltoniano.txt");

    //Condiciones iniciales
    cohete[0]=RT;  //El cohete inicialmente se encuentra en la superficie terrestre (Posicion)
    cohete[1]=0.5;  //Phi
    cohete[2]=v*cos(theta-cohete[1]);  //Momento radial (p_r)
    cohete[3]=RT*v*sin(theta-cohete[1]);  //Momento angular (p_phi)


    //Reescalamos 
    cohete[0]=cohete[0]/dTL;
    cohete[2]=cohete[2]/dTL;
    cohete[3]=cohete[3]/(dTL*dTL);

    //Escribimos las condiciones iniciales de la luna y la reescalamaos
    luna[0]=(RT+RL+dTL)/dTL;
    luna[1]=0.0;

    t=0.0;
    cont=0;

    for(k=0;k<N;k++){

        //Escribimos los datos cada 100 iteraciones
        if(cont==100){
            co<<cohete[0]*cos(cohete[1])<<"\t"<<cohete[0]*sin(cohete[1])<<"\t"<<RT/dTL<<"\t"<<RT/dTL<<endl;
            lu<<luna[0]*cos(w*t)<<"\t"<<luna[0]*sin(w*t)<<endl;

            cont=0;
        }

        //Aplicamos Runge-Kutta

        for(j=0;j<4;j++){
            k1[j]=h*f(cohete[0], cohete[1], cohete[2], cohete[3], t, j);

        }

        for(j=0;j<4;j++){
            k2[j]=h*f(cohete[0]+(k1[0]/2.0), cohete[1]+(k1[1]/2.0), cohete[2]+(k1[2]/2.0), cohete[3]+(k1[3]/2.0), t+(h/2.0), j);
            
        }

        for(j=0;j<4;j++){
            k3[j]=h*f(cohete[0]+(k2[0]/2.0), cohete[1]+(k2[1]/2.0), cohete[2]+(k2[2]/2.0), cohete[3]+(k2[3]/2.0), t+(h/2.0), j);
            
        }

        for(j=0;j<4;j++){
            k4[j]=h*f(cohete[0]+k3[0], cohete[1]+k3[1], cohete[2]+k3[2], cohete[3]+k3[3], t+h, j);
            
        }

        //Calculamos las nuevas coordenadas

        for(j=0;j<4;j++){
            cohete[j]=cohete[j]+(1.0/6.0)*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j]);
        }


        //Reescalamos a las coordenadas originales (para calcular el hamiltoniano)
        r=cohete[0]*dTL;
        phi=cohete[1];
        p_r=cohete[2]*dTL;
        p_phi=cohete[3]*dTL*dTL;


        //Calculamos el Hamiltoniano del sistema para comrpobar que se trata de una cte del movimiento
        r_luna=sqrt((r*r)+(dTL*dTL)-(2.0*r*dTL*cos(phi-(w*t))));

        H=(p_r*p_r)/(2.0)+((p_phi*p_phi)/(2.0*r*r))-(G*MT/r)-(G*ML/r_luna);

        //Lo imprimimos por pantalla
        ham<<k<<"\t"<<H-(w*p_phi)<<endl;






        t=t+h;
        cont++;

    }

    co.close();
    lu.close();
    ham.close();


    return 0;
}


double f(double r, double phi, double pr, double pfi, double t, int i){

    double delta, mu, rp, derivada;

    delta=(G*MT)/(dTL*dTL*dTL);
    mu=ML/MT;
    rp=sqrt(1.0+r*r-2*r*cos(phi-w*t));

    
    if(i==0){
        derivada=pr;
    }

    if(i==1){
        derivada=pfi/(r*r);
    }

    if(i==2){
        derivada=((pfi*pfi)/(r*r*r))-delta*((1.0/(r*r))+(mu/(rp*rp*rp))*(r-cos(phi-w*t)));
    }

    if(i==3){
        derivada=-((delta*mu*r)/(rp*rp*rp))*sin(phi-w*t);
    }


    return derivada;
}