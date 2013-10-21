#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char **argv){
	char rname [3][30],lname [3][30], wname[3][30];
	FILE *rt[3],*in[3], *out[3];
	int N,n,m,rx,ry,lx,ly,count;
	double lambda,l,r,mean[3]={0},disp[3]={0};
	if (argc>=2){
		sscanf(argv[1],"%d",&n);
		sscanf(argv[2],"%d",&m);
	}
	else exit(EXIT_FAILURE);

	sprintf(wname[0],"lambda1.txt");
	sprintf(wname[1],"lambda2.txt");
	sprintf(wname[2],"lambda3.txt");
	sprintf(lname[0],"kutta_a_%d.txt",n);	
	sprintf(lname[1],"kutta_b_%d.txt",n);	
	sprintf(lname[2],"kutta_c_%d.txt",n);	
	sprintf(rname[0],"kutta_a_%d.txt",m);	
	sprintf(rname[1],"kutta_b_%d.txt",m);	
	sprintf(rname[2],"kutta_c_%d.txt",m);
	for (int i=0;i<3;i++){
		in[i]=fopen(lname[i],"r");
		rt[i]=fopen(rname[i],"r");
		out[i]=fopen(wname[i],"w");
		count=0;
		N=0;
		while((!feof(in[i]))||(!feof(rt[i])))
		{
			fscanf(in[i],"%d %d %le",&lx,&ly,&l);
			fscanf(rt[i],"%d %d %le",&rx,&ry,&r);
			if ((rx!=lx)||(ry!=ly)){
				printf("\nFUCK\n");
				break;
			}
			if ((r!=0)&&(l!=0)&&(r*l>1e-15)) {
				lambda=10*(log(r/l))/(m-n);
				if (argc>3){
					if (strcmp(argv[3],"dump")==0) fprintf(out[i],"%d %d %le\n",rx,ry,lambda);
					if ((strcmp(argv[3],"dump-small")==0)&&(rx%2==0)&&(ry%2==0)) fprintf(out[i],"%d %d %le\n",rx,ry,lambda);
				}
				count++;	
				mean[i]+=lambda;
				disp[i]+=lambda*lambda;
			}
			if(i==0) N++;
		}
		mean[i]/=count;
		disp[i]=disp[i]/count-mean[i]*mean[i];
		disp[i]=sqrt(disp[i]);
		printf("meanie#%d=%le\t\tdisp#%d=%le\t\tcount=%d",i,mean[i],i,disp[i],count);
		if (i==0) printf("\tN=%d\n",N);
		else printf("\n");
		fclose(in[i]);
		fclose(rt[i]);
		fclose(out[i]);
	}

  exit(EXIT_SUCCESS);
}
