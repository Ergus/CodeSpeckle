#include "auxiliary.h"

void polint(double* xa, double* ya,int n, double x, double &y, double &dy){
    const int NMAX=10;
    int ns=0;
    double den, dif, dift, ho, hp, w, c[NMAX], d[NMAX];

    dif=fabs(x-xa[0]);
    for(int i=0; i<n; i++){
        dift=fabs(x-xa[i]);
        if (dift<dif){
            ns=i;
            dif=dift;
            }
        c[i]=ya[i];
        d[i]=ya[i];
        }
    y=ya[ns--];
    for(int m=0;m<n-1;m++){
        int bound=n-m-1;
        for (int i=0;i<bound;i++){
            ho=xa[i]-x;
            hp=xa[i+m+1]-x;
            w=c[i+1]-d[i];
            den=ho-hp;
            if(den==0){
                fprintf(stderr,"Failure in polint\n");
                exit (EXIT_FAILURE);
                }
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
            }
        dy=(2*(ns+1)<bound?c[ns+1]:d[ns--]);
        y+=dy;
        }
    }
//(C) Copr. 1986-92 Numerical Recipes Software #$!5,5.){2p491&&k"15.

//modified by S.P. from polin2
void polin3(double *x1a,double *x2a,double *x3a,double *ya,
            int l,int m,int n,
            double x1,double x2,double x3,
            double &y, double &dy){
    
    double yltmp[l], ymtmp[m], yntmp[n];

    /*The next loops are compleetly transversal over the memory
      but the order was not changed because that was the order in
      the fortran code*/
    
    for(int i=0;i<l;i++){
        for(int j=0;j<m;j++){
            for(int k=0;k<n;k++){            
                yntmp[k]=ya[k*m*l+j*l+i];
                }
            polint(x3a,yntmp,n,x3,ymtmp[j],dy);
            }
        polint(x2a,ymtmp,m,x2,yltmp[i],dy);
        }
    polint(x1a,yltmp,l,x1,y,dy);
    
    }

double bessj1(double x){
    double ax, z, xx,y,ans,ans1,ans2;
    
    if((ax=fabs(x)) <8.0){
        y=x*x;
        ans1=x*(72362614232.0 +
                y*(-7895059235.0 +
                   y*(242396853.1 +
                      y*(-2972611.439 +
                         y*(15704.48260 +
                            y*(-30.16036606))))));
        ans2=144725228442.0+
            y*(2300535178.0+
               y*(18583304.74+
                  y*(99447.43394+
                     y*(376.9991397+
                        y*1.0))));
        ans=ans1/ans2;
        }
    else{
        z=8.0/ax;
        y=z*z;
        xx=ax-2.356194491;
        ans1=1.0+
            y*(0.183105e-2+
               y*(-0.3516396496e-4+
                  y*(0.2457520174e-5+
                     y*(-0.240337019e-6))));

        ans2=0.04687499995+
            y*(-0.2002690873e-3+
               y*(0.8449199096e-5+
                  y*(-0.88228987e-6+
                     y*0.105787412e-6)));
        if(x<0) ans=-sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
        else ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
        }
    return ans;
    }

void dbg_printf(const char * format, ... ){
    #ifdef DEBUG
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    #endif
    }


void log_printf(const char * format, ... ){
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime(buffer,80,"%a_%F_%X",timeinfo);
    
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    }

