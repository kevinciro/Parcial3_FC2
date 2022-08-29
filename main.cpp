#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

double f_y(double x, double w_i, double t){
    return( -(1./8.)*t);
}

double f_yp(double x, double w_i, double t){
    return -(1./8.)*(w_i);
}

double f(double x, double w_i, double t){
    return (1./8.)*( 32 + pow(x,3) - (w_i)*t );
}

vector<vector<double>> dif_fin_nl(
        double a,
        double b,
        double alpha,
        double beta,
        int N,
        int M,
        double tol){
    double x;
    double t;
    std::vector<double> A(N+1);
    std::vector<double> B(N+1);
    std::vector<double> C(N+1);
    std::vector<double> D(N+1);
    std::vector<double> L(N+1);
    std::vector<double> U(N+1);
    std::vector<double> Z(N+1);
    std::vector<double> V(N+1);
    std::vector<double> X(N+1);
    vector<vector<double>> respuesta;

//    Paso 1
    double h = (b - a) / (N + 1);
    std::vector<double> W(N+1);
    W[0] = alpha;
    W[N+1] = beta;

//  Paso 2
    for (int i = 0; i < N; ++i) {
        W[i]=alpha + (i * h * (beta - alpha) / (b - a));
    }

//    Paso 3
    int k = 1;

//    Paso 4
    while (k <= M){
//        Paso 5
        x = a + h;
        t = (W[2] - alpha) / (2*h);
        A[0] = 2. + pow(h,2)* f_y(x, W[1],t);
        B[0] = -1 + h * 0.5 * f_yp(x, W[1],t);
        D[0] = -(2*W[1] - W[2] - alpha + pow(h,2)* f(x, W[1],t));

//        Paso 6
        for(int i=1; i<N; i++){
            x = a + i*h;
            t = (W[i+1] - W[i-1]) / (2*h);
            A[i] = 2. + pow(h,2)* f_y(x, W[i],t);
            B[i] = -1 + h * 0.5 * f_yp(x, W[i],t);
            C[i] = -1 - h * 0.5 * f_yp(x, W[i],t);
            D[i] = -(2*W[i] - W[i+1] - W[i-1] + pow(h,2)* f(x, W[i],t));

        }

//        Paso 7
        x = b - h;
        t = (beta - W[N-1]) / (2*h);
        A[N] = 2. + pow(h,2)* f_y(x, W[N],t);
        C[N] = -1 - h * 0.5 * f_yp(x, W[N],t);
        D[N] = -(2*W[N] - W[N-1] - beta + pow(h,2)* f(x, W[N],t));

//        Paso 8
        L[0] = A[0];
        U[0] = B[0]/A[0];
        Z[0] = D[0]/L[0];

//        Paso 9
        for(int i=1; i<N; i++){
            L[i] = A[i] - C[i]*U[i-1];
            U[i] = B[i]/L[i];
            Z[i] = (D[i] - C[i]*Z[i-1])/L[i];
        }

//        Paso 10
        L[N] = A[N] - C[N]*U[N-1];
        Z[N] = (D[N] - C[N]*Z[N-1])/L[N];

//        Paso 11
        V[N] = Z[N];
        W[N] = W[N] + V[N];

//        Paso 12
        for(int i=N; i>0; i--){
            V[i] = Z[i] - U[i]*V[i+1];
            W[i] = W[i] + V[i];
        }

//        Paso 13
        double norma = 0;
        for(double i: V){
            norma = norma + pow(V[i],2);
        }
        norma = sqrt(norma);

        if (norma <= tol){
//            Paso 14
            for(int i=0; i<=N+1; i++){
                X[i] = a + i*h;
            }
//            Paso 15
            respuesta.push_back(X);
            respuesta.push_back(W);
            return respuesta;
        }
        k++;
    }

//    Paso 16 y 17
    throw "Numero maximo de iteraciones alcanzado";
}


int main() {
    vector<vector<double>> respuestas;
    respuestas = dif_fin_nl(1,3,17,43./3.,19,100000, 1e-8);

    cout << respuestas[0].size() << endl;
    cout << respuestas[1].size() << endl;

    for ( double x: respuestas[0]){
        cout << x << " ";
    }
    cout << endl;

    for ( double w: respuestas[1]){
        cout << w << " ";
    }
    cout << endl;


    return 0;
}
