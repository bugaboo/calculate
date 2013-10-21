#include <cmath>
#include <stdio.h>

static const double kd=0.0;


const static int N=302; 									//lattice size
const static double T=30;
const static int Ni=1700;								//number of iterations
//const static double h=2*M_PI/N;								//spacial step
const static double tau=0.0001;								//time step
const static int Nt=(int)(T/tau);								//number of time steps
const static int d=3;									//dimesionality
double amplit=1,prec=1E-22,koef=1;

typedef double matrix[N+1][N+1];
typedef double tensor[3][d][d][N][N];
matrix *F=new matrix[d];
matrix *F1=new matrix[d];
tensor *A=new tensor[3];

void eqn(bool first){
	int qx,qy;
	double tstep;
	if (first) tstep=tau;
	else tstep=2*tau/3;
   	for(int j=0;j<N;++j){
		for(int k=0;k<N;++k){
			int x,y;
			x=-(int)(N/2)+j;
			y=-(int)(N/2)+k;
			for(int i=0;i<d;++i){
				for(int l=0;l<d;++l){
					for(int n=0;n<3;++n){
						for(int m=0;m<3;++m){
							A[n][m][i][l][j][k]=0;
						}
					}
				}
			}
			if ((j!=0)&&(k!=0)){
				for(qx=-1;qx<2;qx+=2){
					for(qy=-1;qy<2;qy+=2){
						int Dxy;								//velocity correlator
						Dxy=-qx*qy;
						A[1+qx][1+qy][0][0][j][k]=1-2*qx*(y*Dxy+x)+y*y+x*x-2*qx*qy*x*y;							
						A[1+qx][1+qy][0][1][j][k]=2*qx*qy-2*qy*(y*Dxy+x);
						A[1+qx][1+qy][0][2][j][k]=1;
					
						A[1+qx][1+qy][1][0][j][k]=Dxy-qx*(y+x*Dxy);
						A[1+qx][1+qy][1][1][j][k]=-2-qx*(y*Dxy+x)-qy*(y+x*Dxy)+y*y+x*x-2*qx*qy*x*y;
						A[1+qx][1+qy][1][2][j][k]=Dxy-qy*(y*Dxy+x);
					
						A[1+qx][1+qy][2][0][j][k]=1;
						A[1+qx][1+qy][2][1][j][k]=2*qx*qy-2*qx*(y+x*Dxy);
						A[1+qx][1+qy][2][2][j][k]=1-2*qy*(y+x*Dxy)+y*y+x*x-2*qx*qy*x*y;
					}
				}	
			
				A[1][1][0][0][j][k]=-2*(2+kd)*(y*y+x*x)-2/tstep;
				A[1][1][1][1][j][k]=-2*(2+kd)*(y*y+x*x)-2/tstep;
				A[1][1][2][2][j][k]=-2*(2+kd)*(y*y+x*x)-2/tstep;
			}
			for(int i=0;i<d;++i){
				for(int l=0;l<d;++l){
					for(int n=0;n<3;++n){
						for(int m=0;m<3;++m){
							A[n][m][i][l][j][k]*=-tstep/2;
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
		for(int j=0;j<N+1;++j){
			for(int k=0;k<N+1;++k){
				r[i][j][k]=0;
				r1[i][j][k]=0;
				p[i][j][k]=0;
				u[i][j][k]=2*F1[i][j][k]-F[i][j][k];
				if (t==0) {f[i][j][k]=F1[i][j][k]; u[i][j][k]=F1[i][j][k];}
				else f[i][j][k]=4*F1[i][j][k]/3-F[i][j][k]/3;
			}
		}
	}
	log=fopen("param.log","a");
	for(int count=0;count<Ni;count++){
		if (count%10==0) {
			for(int i=0;i<d;++i){
				for(int j=1;j<N;++j){
					for(int k=1;k<N;++k){
						r[i][j][k]=0;
						for(int n=0;n<3;++n){
							for(int m=0;m<3;++m){
								for(int l=0;l<d;++l){
									r[i][j][k]+=A[n][m][i][l][j][k]*u[l][j-1+n][k-1+m];
								}
							}
						}
						r[i][j][k]-=f[i][j][k];
					}
				}
			}
		}
		for(int i=0;i<d;++i){
			for(int j=1;j<N;++j){
				for(int k=1;k<N;++k){
					p[i][j][k]=0;
					for(int n=0;n<3;++n){
						for(int m=0;m<3;++m){
							for(int l=0;l<d;++l){
								p[i][j][k]+=A[n][m][i][l][j][k]*r[l][j-1+n][k-1+m];
							}
						}
					}
					nap+=p[i][j][k]*p[i][j][k];
					nr1+=r[i][j][k]*p[i][j][k];
				}
			}
		}
		//		if (count==0&&t==0&&j<4&&k<4) printf("r(%d,%d,%d)=%e\n",i,j,k,p[i][j][k]);			
		//	if (t==0&&(count<100||count==Ni-1)) printf("count=%d nr=%e\n",count,nr);
		if ((t%10==0)&&(count ==1||count==Ni-1)) fprintf(log,"t=%d\t nr=%e\t",t,nr);
		a=nr1/nap;
		nr=0;
		for(int i=0;i<d;++i){
			for(int j=0;j<N;++j){
				for(int k=0;k<N;++k){
					u[i][j][k]-=a*r[i][j][k];
					r[i][j][k]-=a*p[i][j][k];
					nr+=r[i][j][k]*r[i][j][k];
				}
			}
		}
		
		
		nr1=0;
		nap=0;
		if (nr<amplit*amplit*prec*koef){
			printf("t=%d\tcount=%d\n",t,count);
			if (t%10==0) fprintf(log,"count=%d  \tkoef=%f\tamplit=%f\n",count,koef,amplit);
			
			if (t%50==0){
				for(int i=0;i<d;++i){
					for(int j=1;j<N;++j){
						for(int k=1;k<N;++k){
							p[i][j][k]=0;
							for(int n=0;n<3;++n){
								for(int m=0;m<3;++m){
									for(int l=0;l<d;++l){
										p[i][j][k]+=A[n][m][i][l][j][k]*u[l][j-1+n][k-1+m];
									}
								}
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
					koef*=.2;
					printf("\n\nfnr=%le\tKOEF multiplied\n\n",fnr);
					fnr=0;
				}
			}
			else {
				count=Ni;
			}
		}
		
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
			for(int j=1;j<N;++j){
				for(int k=1;k<N;++k){
					r[i][j][k]=0;
					for(int n=0;n<3;++n){
						for(int m=0;m<3;++m){
							for(int l=0;l<d;++l){
								r[i][j][k]+=A[n][m][i][l][j][k]*u[l][j-1+n][k-1+m];
							}
						}
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
	FILE *in;
	FILE *log;
	char name[d][20];
	int trash,x,y;
	double f=0;
	
	if(out == NULL) /* ALWAYS CHECK RETURN VALUE !!! */
	{
		printf("Failed to open file for writing\n");
		return -1;
	}
	
	for(int i=0;i<d;++i){
		for(int j=0;j<N+1;++j){
			for(int k=0;k<N+1;++k){
				F[i][j][k]=0;
				F1[i][j][k]=0;
			}
		}
	}
	if(argc==2){
		sscanf(argv[1],"%d",&trash);
		if (trash==1){
			F1[0][(int)(N/2)][(int)(N/2)]=2;			
			F1[2][(int)(N/2)][(int)(N/2)]=2;			
		}
		for(int i=0;i<d;++i){
			for(int j=0;j<N;++j){
				for(int k=0;k<N;++k){
					if ((trash==1)&&(i!=1)){
						x=-(int)(N/2)+j;
						y=-(int)(N/2)+k;
						if (x*x+y*y==2)
							F1[i][j][k]=1;
					}
					if ((trash==2)&&(i!=1)){
						x=-(int)(N/2)+j;
						y=-(int)(N/2)+k;
						if (x*x+y*y==1)
							F1[i][j][k]=1;
					}
					
				}
			}
		}
	}
	double tik=0;
	int nfile=0;
	amplit=1;
	for(int i=0;i<d;++i){
		double min=F1[i][0][0],max=F1[i][0][0];
		for(int j=0;j<N;++j){
			for(int k=0;k<N;++k){
				if (j==k) 
					{
						if (F1[i][j][k]>max) max=F1[i][j][k];
						if (F1[i][j][k]<min) min=F1[i][j][k];
						if (amplit<F1[i][j][k]) amplit=F1[i][j][k];
						if (-amplit>F1[i][j][k]) amplit=-F1[i][j][k];
				}
			}
		}
			printf ("max=%f\t min=%f\n",max,min);
	}
	printf("amplit=%le\n",amplit);
	for(int t=0;t<Nt;++t){
		if (t==0) eqn(true);
		if (t==1) eqn(false);
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
				min=F1[i][0][0];
				max=F1[i][0][0];
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
