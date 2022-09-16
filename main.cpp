#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <list>
#include <fstream> // flujo de archivo


using namespace std;

double S = 1000.0; //tension puntos extremos
double E = 3.0e7;//modulo elasticidad
double I = 625.0; //momento inercia
double q = 100.0; //intensidad carga
double l = 120.0; //longitud


double f_y(double x, double w_i, double t){
    return -18*w_i*pow(x,-5);
}

double f_yp(double x, double w_i, double t){
    return 2*pow(x,-3)*t;
}

double f(double x, double w_i, double t){
    return pow(t,2)*pow(x,-3) - 9*pow(w_i,2)* pow(x,-5) + 4*x;
}

//double f_y(double x, double w_i, double t){
//    return( pow((1. + pow(t,2)), 3./2.) * S/(E*I) );
//}
//
//double f_yp(double x, double w_i, double t){
//    return 3*( (S/(E*I))*w_i + (x-l)*(q*x)/(2*E*I) )*pow((1. + pow(t,2)), 1./2.)*t;
//}
//
//double f(double x, double w_i, double t){
//    return pow((1. + pow(t,2)), 3./2.)*( (S/(E*I))*w_i + (x-l)*(q*x)/(2*E*I) );
//}

//Ejemplo del libro
//double f_y(double x, double w_i, double t){
//    return( -(1./8.)*t);
//}
//
//double f_yp(double x, double w_i, double t){
//    return -(1./8.)*(w_i);
//}
//
//double f(double x, double w_i, double t){
//    return (1./8.)*( 32 + pow(x,3) - (w_i)*t );
//}

//ecuaciones del péndulo simple no lineales.
//const double g=9.8;
//const double l=1.0;
//double f_y(double x, double w_i, double t){
//    return( -(g/l)*cos(x) );
//}
//
//double f_yp(double x, double w_i, double t){
//    return 0;
//}
//
//double f(double x, double w_i, double t){
//    return (-(g/l)*sin(x) );
//}

vector<vector<double>> dif_fin_nl(
//        [a,b] es el intervlo
        double a,
        double b,
//        alfa y beta son las condiciones de frontera en el intervalo
        double alpha,
        double beta,
//        N es el numero de intervalos en los que se dividira [a,b]
        int N,
//        M es el numero de iteraciones maxima
        int M,
//        tolerancia deseada
        double tol){


    double x;
//    t es la derivada
    double t;

    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;
    std::vector<double> D;
    std::vector<double> L;
    std::vector<double> U;
    std::vector<double> Z;
    std::vector<double> V(N+1);
    std::vector<double> X;
    vector<vector<double>> respuesta;

//    Paso 1
//    las condiciones iniciales se guardan en el vector W
    double h = (b - a) / (N + 1);
    std::vector<double> W;
    W.push_back(alpha); //posicion 0

//  Paso 2
    for (int i = 1; i <= N; ++i) {
        W.push_back(alpha + (i * h * (beta - alpha) / (b - a))); //posicion 1 hasta la N-1
    }

    W.push_back(beta); // posicion N+1


//    for ( double w: W){
//        cout << w << " ";
//    }
//    cout << endl;


//    Paso 3
    int k = 1;

//    Paso 4
    while (k <= M){
//        Paso 5
        x = a + h;
        t = (W[2] - alpha) / (2*h);
        A.push_back(2. + pow(h,2) * f_y(x, W[1],t)); //posicion 0
        B.push_back(-1 + h * 0.5 * f_yp(x, W[1],t)); //posicion 0
        D.push_back(-(2*W[1] - W[2] - alpha + pow(h,2)* f(x, W[1],t))); //posicion 0

//        Paso 6
        for(int i=1; i<N; i++){
            x = a + i*h;
            t = (W[i+1] - W[i-1]) / (2*h);
            A.push_back(2. + pow(h,2)* f_y(x, W[i],t)); //posiciones desde la 1 hasta la N-1
            B.push_back(-1 + h * 0.5 * f_yp(x, W[i],t)); //posiciones desde la 1 hasta la N-1
            C.push_back(-1 - h * 0.5 * f_yp(x, W[i],t)); //posiciones desde la 0 hasta la N-2
            D.push_back(-(2*W[i] - W[i+1] - W[i-1] + pow(h,2)* f(x, W[i],t))); //posiciones desde la 1 hasta la N-1

        }

//        Paso 7
        x = b - h;
        t = (beta - W[N-1]) / (2*h);
        A.push_back(2. + pow(h,2)* f_y(x, W[N],t)); // posicion N
        C.push_back(-1 - h * 0.5 * f_yp(x, W[N],t)); // posicion N-1
        D.push_back(-(2*W[N] - W[N-1] - beta + pow(h,2)* f(x, W[N],t))); //posicion N

//        Paso 8
        L.push_back(A[0]); //posicion 0
        U.push_back(B[0]/A[0]); //posicion 0
        Z.push_back(D[0]/L[0]); //posicion 0

//        Paso 9
        for(int i=1; i<N; i++){
            L.push_back(A[i] - C[i]*U[i-1]); //posiciones desde la 1 hasta la N-1
            U.push_back(B[i]/L[i]); //posiciones desde la 1 hasta la N-1
            Z.push_back((D[i] - C[i]*Z[i-1])/L[i]); //posiciones desde la 1 hasta la N-1
        }

//        Paso 10
        L.push_back(A[N] - C[N]*U[N-1]); // posicion N
        Z.push_back((D[N] - C[N]*Z[N-1])/L[N]); // posicion N

//        Paso 11
        V[N] = Z[N]; // posisicon N
        W[N] = W[N] + V[N]; // posisicon N, La posicion N+1 siempre sera beta

//        Paso 12
        for(int i=N-1; i>=0; i--){
            V[i] = Z[i] - U[i]*V[i+1]; // posisicon desde la N-1 hasta la 0
            W[i] = W[i] + V[i]; // posisicon desde la N-1 hasta la 0
        }

//        Paso 13
        double norma = 0;
        for(double i : V){
            norma += pow(i,2);
        }
        norma = sqrt(norma);
//        cout << norma << endl;

/*        for ( double v: V){
            cout << v << " ";
        }*/

        if (norma <= tol){
//            Paso 14
            for(int i=0; i<=N+1; i++){
                X.push_back(a + i*h);
            }
//            Paso 15
            respuesta.push_back(X);
            respuesta.push_back(W);
            return respuesta;
        }

        A.clear();
        B.clear();
        C.clear();
        D.clear();
        L.clear();
        U.clear();
        Z.clear();

        k++;
    }

//    Paso 16 y 17
    throw "Numero maximo de iteraciones alcanzado";

}


int main() {
    const double PI = 3.1416;

    ofstream archivo;
	
   archivo.open("programa.txt",ios::out); //le indicamos que queremos añadir texto
   if(archivo.fail()){
			cout<<"No se pudo abrir el archivo"<<endl;
			exit(1);
	}
    try{
        vector<vector<double>> respuestas;
//        double a = 0;
//        double b = PI * sqrt(l/g);
//    //        alfa y beta son las condiciones de frontera en el intervalo
//        double alpha = -PI/4;
//        double beta = PI/4;
//    //        N es el numero de intervalos en los que se dividira [a,b]
//        int N = 10;
//    //        M es el numero de iteraciones maxima
//        int M = 1000000;
//    //        tolerancia deseada
//        double tol = 1e-8;
//        respuestas = dif_fin_nl(a, b, alpha, beta, N, M, tol);

//    Ejemplo del libro
//        respuestas = dif_fin_nl(1,3,17,43./3.,50,100000, 1e-12); //a, b, alpha, beta, N, M, tol.

//        respuestas = dif_fin_nl(0, l, 0,0,50,100000, 1e-12); //a, b, alpha, beta, N, M, tol.

        respuestas = dif_fin_nl(1, 2, 0, log(256),50,100000, 1e-12); //a, b, alpha, beta, N, M, tol.

        for ( int i=0; i<respuestas[0].size(); i++ ){
            archivo << respuestas[0][i] << " " << respuestas[1][i] << endl;
        }
//        cout << endl;

//        for ( double w: respuestas[1]){
//            cout << w << ", ";
//        }
        cout << endl;
    }catch (const char* msg){
        std::cout<< msg << std::endl;
    }
    
    archivo.close();

    return 0;
}

