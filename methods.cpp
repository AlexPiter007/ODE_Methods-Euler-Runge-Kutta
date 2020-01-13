#include <iostream>
#include<cmath>
#include<iomanip>
#include<vector>
using namespace std;

    class Methods{
    public:

        double func(double x, double y){
            //double eta=5, k=1, m=0.6,n=1;
            //return eta*pow(y,m)-pow(y,n)*k;
            //return x*y;
            //return pow((y+x), 2);
            return x+y;
        }
        void Euler(double a,double b,double h,double yO ) {
            double n=(b-a)/h;
            double X[(int)n];
            double y[(int)n];
            X[0]=a;
            y[0]=yO;
            for(int i=0; i<n; i++){
                X[i+1]=X[i]+h;
                y[i+1] = (y[i]+h*func(X[i],y[i]));
            }
            for (int j = 0; j <=n; ++j) {
                printf("X = %.2f", X[j]);
                cout<<"\tY = "<<y[j]<<endl;
            }
        }
        void Runge_kutte(double a,double b, double h, double yO) {
            double n=(b-a)/h;
            double x[(int)n];
            double y[(int)n];
            double k1,k2,k3,k4;
            x[0]=a;
            y[0]=yO;
            for(int i=0; i<n; i++){
                x[i+1]=x[i]+h;
                k1=h*func(x[i],y[i]);
                k2=h*func(x[i]+h/2.0,y[i]+k1/2.0);
                k3=h*func(x[i]+h/2.0, y[i]+k2/2.0);
                k4=h*func(x[i]+h,y[i]+k3);
                y[i+1]=y[i]+(k1+2*k2+2*k3+k4)/6;

            }
            for (int j = 0; j <=n; ++j) {
                printf("X = %.2f", x[j]);
                cout<<"\tY = "<<y[j]<<endl;
            }

        }
        void Rkutm(double a,double b, double h, double yO,double e) {
            const double n=100000;
            double x[(int)n],y[(int)n];
            double k1,k2,k3,k4,k5,R;
            x[0]=a; y[0]=yO;
            for(int k=0; k<n; k++){
                x[k+1]=x[k]+h;
                k1=func(x[k],y[k]);
                k2=func(x[k]+h/3.0,y[k]+k1*h/3.0);
                k3=func(x[k]+h/3,y[k]+(k1*h/6)+(k2*h/3));
                k4=func(x[k]+h/2,y[k]+(k1*h/3)+(k3*(3*h)/8));
                k5=func(x[k]+h/3,y[k]+(k1*h/2)-(k3*(3*h)/2)+2*k3*h);
                y[k+1]=y[k]+(h*(k1+4*k4+k5))/6;
                R=(-2*k1+9*k3-8*k4+k1)/32;
                //cout<<"\tR = "<<fabs(R)<<endl;
                if((fabs(R)<=e)&&(fabs(R)>=(e/32))){
                    h=h;
                }
                if(fabs(R)>e){
                    cout<<"\tR = "<<fabs(R)<<" > "<<e<<endl;
                    h=h/2;
                    k=0;
                }
                if(fabs(R)<(e/32)){
                    cout<<"\tR = "<<fabs(R)<<" < "<<e/32<<endl;
                    h=h*2;
                    k=0;
                }
                if(x[k]>b)
                    break;
            }
            int k = 0;
            while (true) {
                printf("X = %.9f",x[k]);
                cout<<"\tY = "<<y[k]<<endl;
                if (x[k]>b)
                    break;
                k++;
            }


        }
        void inglenda(double a, double b, double h, double yO){
            double n=(b-a)/h;
            double x[(int) n];
            double y[(int) n];
            double k1,k2,k3,k4,k5,k6;
            x[0]=a;
            y[0]=yO;
            for(int i=0; i<=n; i++) {
                x[i + 1] = x[i] + h;
                k1 = func(x[0], y[0]);
                k2 = func(x[i] + h / 2, y[i] + k1 * (h / 2));
                k3 = func(x[i] + h / 2, y[i] + (h / 4) * (k1 + k2));
                k4 = func(x[i] + h, y[i] + h * (2 * k3 - k2));
                k5 = func(x[i] + (2 * h) / 3, y[i] + (h / 27) * (7 * k1 + 10 * k2 + k4));
                k6 = func(x[i] + h / 5, y[i] + (h / 625) * (28 * k1 - 125 * k2 + 546 * k3 + k4 * k4 - 378 * k5));
                y[i + 1] = y[i] + (1.0/ 336.0) * h * (14 * k1 + 35 * k4 + 162 * k5 + 125 * k6);

            }
            for (int j = 0; j <=n; ++j) {
                printf("X = %.2f", x[j]);
                cout<<"\tY = "<<y[j]<<endl;
            }


        }
        void felberg(double a, double b, double h, double yO){
            double n=(b-a)/h;
            double x[(int) n];
            double y[(int) n];
            double k1,k2,k3,k4,k5,k6;
            x[0]=a; y[0]=yO;
            for(int i=0; i<=n; i++){
                x[i+1]=x[i]+h;
                k1 = func(x[0],y[0]);
                k2 = func(x[i]+(h/4),y[i]+k1*(h/4));
                k3 = func(x[i]+((3*h)/8),y[i]+k1*((3*h)/32)+k2*((9*h)/32));
                k4 = func(x[i]+((12*h)/13), y[i]+k1*((1932*h)/2197)-k2*((7200*h)/2197)+k3*((7296*h)/2197));
                k5 = func(x[i]+h, y[i]+k1*((439*h)/2162)-8*h*k2+k3*((3680*h)/513)-k4*((845*h)/4104));
                k6 = func(x[i]+h/2,y[i]+k1*((8*h)/27)+2*h*k2-k3*((3544*h)/2565)+k4*((1859*h)/4104)-k5*((11*h)/40));
                y[i+1] = y[i]+h*((16.0/135)*k1+(6665.0/12825)*k3+(28561.0/56430)*k4+(9.0/50)*k5+(2.0/55)*k6);
            }
            for (int j = 0; j <=n; ++j) {
                printf("X = %.2f", x[j]);
                cout<<"\tY = "<<y[j]<<endl;
            }
        }
    };

int main() {
    cout<<setprecision(17);
    freopen("methods.txt", "w", stdout);
    double a=0;
    double b=1.0;
    double h=0.1;
    double yO=1;
    double eps=1e-3;
  Methods Methods;
  cout<<"\t\t\tNamudi Muodila\t***X+Y***"<<endl;
  cout<<"\ta = "<<a<<endl;
  cout << "\tb = " << b<<endl;
  cout << "\th = " << h << endl;
  cout << "\tY_0 = " << yO << endl;
  cout << "\tepsilon = " << eps << endl;


  cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
  cout << "\t\tEuler Methods\n";
  Methods.Euler(a, b, h, yO);
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
  cout << "\t\tRunge kutte Methods\n";
  Methods.Runge_kutte(a, b, h, yO);
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
  cout << "\t\tInglenda Methods\n";
  Methods.inglenda(a, b, h, yO);
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
  cout << "\t\tRungge kutte Merson Methods\n";
  Methods.Rkutm(a, b, h, yO, eps);
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
  cout << "\t\tFelberg Methods\n";
  Methods.felberg(a, b, h, yO);
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
 
  return 0;
}
