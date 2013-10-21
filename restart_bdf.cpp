#include <cmath>
#include <stdio.h>

static const double kd=0.001;


const static int N=600;								//spatial discretization
const static int Nt=70000;								//number of time steps
const static int Ni=1700;								//number of iterations
const static double h=2*M_PI/N;								//spacial step
const static double tau=0.0001;								//time step
const static int d=3;									//dimesionality
double amplit=1,prec=1E-22,koef=1;

typedef double matrix[N][N];
typedef double tensor[3][d][d][N][N];
matrix *F=new matrix[d];
matrix *F1=new matrix[d];
tensor *A=new tensor[3];

void eqn(bool first){
	double tstep;
	tstep=2*tau/3;
     	for(int j=0;j<N;++j){
		for(int k=0;k<N;++k){
			A[0][0][0][0][j][k]=(1-cos(-2*M_PI+(j+k)*h))/2-(1-cos(-2*M_PI+(j-k)*h))/2;							
			A[0][1][0][0][j][k]=-(1-cos(-2*M_PI+(j+k)*h))-(1-cos(-2*M_PI+(j-k)*h))-kd;
			A[0][2][0][0][j][k]=-(1-cos(-2*M_PI+(j+k)*h))/2+(1-cos(-2*M_PI+(j-k)*h))/2;							
			A[1][0][0][0][j][k]=-(1-cos(-2*M_PI+(j+k)*h))-(1-cos(-2*M_PI+(j-k)*h))-kd;
			A[1][1][0][0][j][k]=4*(1-cos(-2*M_PI+(j+k)*h))+4*(1-cos(-2*M_PI+(j-k)*h))+4*kd+(h*h)/tstep;
			A[1][2][0][0][j][k]=-(1-cos(-2*M_PI+(j+k)*h))-(1-cos(-2*M_PI+(j-k)*h))-kd;
			A[2][0][0][0][j][k]=-(1-cos(-2*M_PI+(j+k)*h))/2+(1-cos(-2*M_PI+(j-k)*h))/2;						
			A[2][1][0][0][j][k]=-(1-cos(-2*M_PI+(j+k)*h))-(1-cos(-2*M_PI+(j-k)*h))-kd;
			A[2][2][0][0][j][k]=(1-cos(-2*M_PI+(j+k)*h))/2-(1-cos(-2*M_PI+(j-k)*h))/2;	
				
                	for(int n=0;n<3;++n){
                            	for(int m=0;m<3;++m){
                                        A[n][m][1][1][j][k]=A[n][m][0][0][j][k];
                                        A[n][m][2][2][j][k]=A[n][m][0][0][j][k];
      		                }
                  	}
			

			A[0][0][0][1][j][k]=0;
			A[0][1][0][1][j][k]=-h*sin(-2*M_PI+(j+k)*h);
			A[0][2][0][1][j][k]=0;
			A[1][0][0][1][j][k]=h*sin(-2*M_PI+(j+k)*h);
			A[1][1][0][1][j][k]=0;
			A[1][2][0][1][j][k]=-h*sin(-2*M_PI+(j+k)*h);
			A[2][0][0][1][j][k]=0;
			A[2][1][0][1][j][k]=h*sin(-2*M_PI+(j+k)*h);
			A[2][2][0][1][j][k]=0;


			A[0][0][0][2][j][k]=0;
			A[0][1][0][2][j][k]=-h*sin(-2*M_PI+(j-k)*h);
			A[0][2][0][2][j][k]=0;
			A[1][0][0][2][j][k]=-h*sin(-2*M_PI+(j-k)*h);
			A[1][1][0][2][j][k]=0;
			A[1][2][0][2][j][k]=h*sin(-2*M_PI+(j-k)*h);
			A[2][0][0][2][j][k]=0;
			A[2][1][0][2][j][k]=h*sin(-2*M_PI+(j-k)*h);
			A[2][2][0][2][j][k]=0;

			A[0][0][1][0][j][k]=0;								
			A[0][1][1][0][j][k]=-2*h*sin(-2*M_PI+(j-k)*h);								
			A[0][2][1][0][j][k]=0;								
			A[1][0][1][0][j][k]=-2*h*sin(-2*M_PI+(j-k)*h);								
			A[1][1][1][0][j][k]=0;								
			A[1][2][1][0][j][k]=2*h*sin(-2*M_PI+(j-k)*h);								
			A[2][0][1][0][j][k]=0;								
			A[2][1][1][0][j][k]=2*h*sin(-2*M_PI+(j-k)*h);								
			A[2][2][1][0][j][k]=0;										

					

			A[0][0][1][2][j][k]=0;
			A[0][1][1][2][j][k]=0;
			A[0][2][1][2][j][k]=0;
			A[1][0][1][2][j][k]=0;
			A[1][1][1][2][j][k]=-4*cos(-2*M_PI+(j-k)*h)*h*h;
			A[1][2][1][2][j][k]=0;
			A[2][0][1][2][j][k]=0;
			A[2][1][1][2][j][k]=0;
			A[2][2][1][2][j][k]=0;
				
			A[0][0][2][0][j][k]=0;												
			A[0][1][2][0][j][k]=-2*h*sin(-2*M_PI+(j+k)*h);									
			A[0][2][2][0][j][k]=0;												
			A[1][0][2][0][j][k]=2*h*sin(-2*M_PI+(j+k)*h);									
			A[1][1][2][0][j][k]=0;									
			A[1][2][2][0][j][k]=-2*h*sin(-2*M_PI+(j+k)*h);									
			A[2][0][2][0][j][k]=0;												
			A[2][1][2][0][j][k]=2*h*sin(-2*M_PI+(j+k)*h);									
			A[2][2][2][0][j][k]=0;												
				
			A[0][0][2][1][j][k]=0;
			A[0][1][2][1][j][k]=0;
			A[0][2][2][1][j][k]=0;
			A[1][0][2][1][j][k]=0;
			A[1][1][2][1][j][k]=-4*cos(-2*M_PI+(j+k)*h)*h*h;
			A[1][2][2][1][j][k]=0;
			A[2][0][2][1][j][k]=0;
			A[2][1][2][1][j][k]=0;
			A[2][1][2][1][j][k]=0;

			for(int i=0;i<d;++i){
				for(int l=0;l<d;++l){
					for(int n=0;n<3;++n){
						for(int m=0;m<3;++m){
							A[n][m][i][l][j][k]*=tstep/(h*h);
						}
					}
				}
			}
		}
	}
	  
  
}
void bdf (int t){
	double a=0,dt=0,nr=0,nap=0,nr1=0,b=0,dt1=0,tmp=0,fnr=0,koef=1;
	matrix *r=new matrix[d];
	matrix *r1=new matrix[d];
	matrix *p=new matrix[d];
	matrix *u=new matrix[d];
	matrix *u1=new matrix[d];
	matrix *f=new matrix[d];
	FILE *log;

	int rx,lx,ry,ly;

	for(int i=0;i<d;++i){
		for(int j=0;j<N;++j){
			for(int k=0;k<N;++k){
				r[i][j][k]=0;
				r1[i][j][k]=0;
				p[i][j][k]=0;
				u[i][j][k]=2*F1[i][j][k]-F[i][j][k];
				f[i][j][k]=4*F1[i][j][k]/3-F[i][j][k]/3;
			}
		}
	}

	for(int i=0;i<d;++i){
		for(int j=0;j<N;++j){
			for(int k=0;k<N;++k){
					rx=j+1;
					lx=j-1;
					ry=k+1;
					ly=k-1;
					if (j==0) lx=N-1;
					if (j==N-1) rx=0;
					if (k==0) ly=N-1;
					if (k==N-1) ry=0;
					for(int n=0;n<d;++n){
						r[i][j][k]+=A[0][0][i][n][j][k]*u[n][lx][ly]+A[0][1][i][n][j][k]*u[n][lx][k]+A[0][2][i][n][j][k]*u[n][lx][ry]+
								A[1][0][i][n][j][k]*u[n][j][ly]+A[1][1][i][n][j][k]*u[n][j][k]+A[1][2][i][n][j][k]*u[n][j][ry]+
								A[2][0][i][n][j][k]*u[n][rx][ly]+A[2][1][i][n][j][k]*u[n][rx][k]+A[2][2][i][n][j][k]*u[n][rx][ry];
					}
					r[i][j][k]-=f[i][j][k];
					//nr+=r[i][j][k]*r[i][j][k];
					
			}
		}
	}

	log=fopen("param.log","a");
	for(int count=0;count<Ni;count++){
		
		for(int i=0;i<d;++i){
			for(int j=0;j<N;++j){
				for(int k=0;k<N;++k){
					p[i][j][k]=0;
					rx=j+1;
					lx=j-1;
					ry=k+1;
					ly=k-1;
					if (j==0) lx=N-1;
					if (j==N-1) rx=0;
					if (k==0) ly=N-1;
					if (k==N-1) ry=0;
					for(int n=0;n<d;++n){
						p[i][j][k]+=A[0][0][i][n][j][k]*r[n][lx][ly]+A[0][1][i][n][j][k]*r[n][lx][k]+A[0][2][i][n][j][k]*r[n][lx][ry]+
								A[1][0][i][n][j][k]*r[n][j][ly]+A[1][1][i][n][j][k]*r[n][j][k]+A[1][2][i][n][j][k]*r[n][j][ry]+
								A[2][0][i][n][j][k]*r[n][rx][ly]+A[2][1][i][n][j][k]*r[n][rx][k]+A[2][2][i][n][j][k]*r[n][rx][ry];
					}
					nap+=p[i][j][k]*p[i][j][k];
					nr1+=r[i][j][k]*p[i][j][k];
					nr+=r[i][j][k]*r[i][j][k];
			//		if (count==0&&t==0&&j<4&&k<4) printf("r(%d,%d,%d)=%e\n",i,j,k,p[i][j][k]);
				}
			}
		}
			//	if (t==0&&(count<100||count==Ni-1)) printf("count=%d nr=%e\n",count,nr);
		if ((t%10==0)&&(count ==1||count==Ni-1)) fprintf(log,"t=%d\t nr=%e\t",t,nr);
		a=nr1/nap;
		for(int i=0;i<d;++i){
			for(int j=0;j<N;++j){
				for(int k=0;k<N;++k){
					u[i][j][k]-=a*r[i][j][k];
					r[i][j][k]-=a*p[i][j][k];
					//nr1+=r[i][j][k]*r[i][j][k];
				}
			}
		}
		
		/*b=nr1/nr;
		
		for(int i=0;i<d;++i){
			for(int j=0;j<N;++j){
				for(int k=0;k<N;++k){
					p[i][j][k]=r[i][j][k]+b*p[i][j][k];
				}
			}
		}*/
		
		
		nr1=0;
		nap=0;
		if (nr<amplit*amplit*prec*koef){
			printf("t=%d\tcount=%d\n",t,count);
			if (t%10==0) fprintf(log,"count=%d  \tkoef=%f\tamplit=%f\n",count,koef,amplit);
			
			if (t%50==0){
				for(int i=0;i<d;++i){
					for(int j=0;j<N;++j){
						for(int k=0;k<N;++k){
								p[i][j][k]=0;
								rx=j+1;
								lx=j-1;
								ry=k+1;
								ly=k-1;
								if (j==0) lx=N-1;
								if (j==N-1) rx=0;
								if (k==0) ly=N-1;
								if (k==N-1) ry=0;
								for(int n=0;n<d;++n){
									p[i][j][k]+=A[0][0][i][n][j][k]*u[n][lx][ly]+A[0][1][i][n][j][k]*u[n][lx][k]+A[0][2][i][n][j][k]*u[n][lx][ry]+
											A[1][0][i][n][j][k]*u[n][j][ly]+A[1][1][i][n][j][k]*u[n][j][k]+A[1][2][i][n][j][k]*u[n][j][ry]+
											A[2][0][i][n][j][k]*u[n][rx][ly]+A[2][1][i][n][j][k]*u[n][rx][k]+A[2][2][i][n][j][k]*u[n][rx][ry];
								}
								p[i][j][k]-=f[i][j][k];
								fnr+=p[i][j][k]*p[i][j][k];
								
						}
					}
				}
				if (fnr<amplit*amplit*prec*koef*10) {
					fprintf(log,"fnr=%le\n",fnr);
					printf("fnr=%le\n",fnr);
					count=Ni;
				}
				else {
					fnr=0;
					koef*=.2;
					printf("\n\nKOEF multiplied\n\n");
				}
			}
			else {
				count=Ni;
			}
		}
		nr=0;
		
	}
	for(int i=0;i<d;++i){
		for(int j=0;j<N;++j){
			for(int k=0;k<N;++k){
				F[i][j][k]=F1[i][j][k];
				F1[i][j][k]=u[i][j][k];
				//printf("r(%d,%d,%d)=%f\n",i,j,k,p[i][j][k]);
			}
		}
	}
	
	if ((t%100==0)&&(1==0)){
		nr=0;
		for(int i=0;i<d;++i){
			for(int j=0;j<N;++j){
				for(int k=0;k<N;++k){
						r[i][j][k]=0;
						rx=j+1;
						lx=j-1;
						ry=k+1;
						ly=k-1;
						if (j==0) lx=N-1;
						if (j==N-1) rx=0;
						if (k==0) ly=N-1;
						if (k==N-1) ry=0;
						for(int n=0;n<d;++n){
							r[i][j][k]+=A[0][0][i][n][j][k]*u[n][lx][ly]+A[0][1][i][n][j][k]*u[n][lx][k]+A[0][2][i][n][j][k]*u[n][lx][ry]+
									A[1][0][i][n][j][k]*u[n][j][ly]+A[1][1][i][n][j][k]*u[n][j][k]+A[1][2][i][n][j][k]*u[n][j][ry]+
									A[2][0][i][n][j][k]*u[n][rx][ly]+A[2][1][i][n][j][k]*u[n][rx][k]+A[2][2][i][n][j][k]*u[n][rx][ry];
						}
						r[i][j][k]-=f[i][j][k];
						nr+=r[i][j][k]*r[i][j][k];
						
				}
			}
		}
	
	printf("final nr=%e\n",nr);
	fprintf(log,"final nr=%e\n",nr);
	
	}
	fclose(log);
	delete [] r; delete [] u; delete [] u1; delete [] r1;  delete [] f; delete [] p;
}



int main(int argc,char **argv){
	FILE *out[d];
	FILE *outF[d];
	FILE *in[d];
	FILE *log;
	char name[d][20];
	int trash,x,y,startt;
	double f,sum=0;
	
    if (argc==2) {
        sscanf(argv[1],"%d",&startt);
        sprintf(name[0],"kutta_a_%d.txt",startt);
        sprintf(name[1],"kutta_b_%d.txt",startt);
        sprintf(name[2],"kutta_c_%d.txt",startt);
        for(int i=0;i<d;++i){

            in[i]=fopen(name[i],"r");
            while (fscanf(in[i],"%d %d %le",&x,&y,&f)!=EOF){
                F1[i][x][y]=f;
                }
            fclose(in[i]);

        }
        sprintf(name[0],"kutta_a.txt");
        sprintf(name[1],"kutta_b.txt");
        sprintf(name[2],"kutta_c.txt");
        for(int i=0;i<d;++i){

            in[i]=fopen(name[i],"r");
            while (fscanf(in[i],"%d %d %le",&x,&y,&f)!=EOF){
                F[i][x][y]=f;
                }
            fclose(in[i]);

        }

    }
    else return -1;
	


	double tik=0;
	int nfile=startt;
	amplit=0;
	for(int i=0;i<d;++i){
		double min=F1[i][0][0],max=F1[i][0][0];
		for(int j=0;j<N;++j){
			for(int k=0;k<N;++k){				
						if (F1[i][j][k]>max) max=F1[i][j][k];
						if (F1[i][j][k]<min) min=F1[i][j][k];
						if (0<F1[i][j][k]) amplit+=F1[i][j][k]/(N*N*d);
						if (0>F1[i][j][k]) amplit+=-F1[i][j][k]/(N*N*d);			
			}
		}
			printf ("max=%f\t min=%f\n",max,min);
	}
	eqn(true);
	printf("amplit=%le\n",amplit);
	for(int t=(int)(nfile/(10*tau));t<(int)(nfile/(10*tau))+Nt;++t){
		//rk4(t);
		bdf(t);
		tik+=tau;
		if (tik>=0.1)
		{
			tik-=0.1;
			nfile++;
			trash=sprintf(name[0],"kutta_a_%d.txt",nfile);
			trash=sprintf(name[1],"kutta_b_%d.txt",nfile);
			trash=sprintf(name[2],"kutta_c_%d.txt",nfile);
			outF[0]=fopen("kutta_a.txt","w");
			outF[1]=fopen("kutta_b.txt","w");
			outF[2]=fopen("kutta_c.txt","w");
			for(int i=0;i<d;++i){
				out[i]=fopen(name[i],"w");
				for(int j=0;j<N;++j){
					for(int k=0;k<N;++k){
						fprintf(out[i],"%d %d %e\n",j, k, F1[i][j][k]);
						fprintf(outF[i],"%d %d %e\n",j, k, F[i][j][k]);
					}
					fprintf(out[i],"\n");
				}
				fclose(out[i]);
				fclose(outF[i]);
			}

		}
		if (t%10==0){
			amplit=0;
			double min,max;
			for(int i=0;i<d;++i){
				min=100;
				max=-100;
				for(int j=0;j<N;++j){
					for(int k=0;k<N;++k){
						if (F1[i][j][k]>max) max=F1[i][j][k];
						if (F1[i][j][k]<min) min=F1[i][j][k];
						if (0<F1[i][j][k]) amplit+=F1[i][j][k]/(N*N*d);
						if (0>F1[i][j][k]) amplit-=F1[i][j][k]/(N*N*d);
					}
				}
				printf ("max=%f\t min=%f\n",max,min);	
			}
			log=fopen("param.log","a");
			fprintf(log,"maxc=%f\n",max);
			fclose(log);
		}
		//printf("amplit=%le\n",amplit);		
	}

	delete [] F;
	delete [] F1;
	delete [] A;
	return 0;
}
