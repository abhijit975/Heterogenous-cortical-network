// For homogeneous cluster.....................................
#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>

#include "ran2.cpp"
static int long seed = time(NULL);

#define ne 4000
#define ni 1000

#define K 800
#define pee 0.2
#define pei 0.5
#define pie 0.5
#define pii 0.5

#define taui 10.0
#define taue 15.0
#define	tauerise 1.0
#define	tauedecay 3.0
#define	tauirise 1.0
#define	tauidecay 2.0


#define	vre 0.

#define	threshe 1. 
#define	threshi 1.

#define	refrac 5.0
#define maxrate 100.

#define	muemin 1.1
#define	muemax 1.2
#define	muimin 1.0
#define	muimax 1.05

void initialize(int N, double array[])
{
	int i;
	for(i=0;i<N;i++){
		array[i] = 0.0;
	}	
}

int main()
{
FILE *fp,*fptr;
fptr = fopen("checkI.dat","w");
fp = fopen("spikeraster.dat","w");
printf("Creating Parameters..\n");

double jee, jei, jie, jii;
int Ncells = ne + ni;
double k = sqrt(K);


//defining connection strengths.............
jee = taui/(k*taue) ;
jie = (pee/pie)*(1./k);
jei = -1.2*(ne/ni)*(pee/pie)*(taui/taue)*(1./k) ;
jii = - (ne/ni)*(pee/pie)*(1./k) ;

printf("jee = %lf; jei = %lf; jie = %lf; jii = %lf;\n",jee,jei,jie,jii);

double dt = 0.1;
double T = 2000;  // runtime
int Nsteps = (int)round(T/dt);

//making arrays of bias voltages, tau and threshold values.................
double mu[Ncells], threshold[Ncells], tau[Ncells];
int i,j,m;

for(i=0;i<Ncells;i++)
	{
	if(i<ne)
		{
		mu[i] = (muemax-muemin)*ran2(&seed) + muemin ;	
		}	
	else
		{
		mu[i] = (muimax-muimin)*ran2(&seed) + muimin ;	
		}	
	}

for(i=0;i<Ncells;i++)
	{
	if(i<ne)
		{
		tau[i] = taue ;	
		}	
	else
		{
		tau[i] = taui ;	
		}	
	}

for(i=0;i<Ncells;i++){
	if(i<ne){
		threshold[i] = threshe ;	
		}	
	else{
		threshold[i] = threshi ;	
		}	
	}

//..................................................

//weights matrix definition.........................
double **weights = (double **)malloc(Ncells * sizeof(double *));
for(i=0;i<Ncells;i++)
	weights[i] = (double *)malloc(Ncells* sizeof(double));

for(i=0;i<Ncells;i++)
{
	for(j=0;j<Ncells;j++){
		weights[i][j] = 0;}}

for(i=0;i<Ncells;i++)
	{for(j=0;j<Ncells;j++)
		{
		if (i<ne && j<ne && ran2(&seed)<pee){
			weights[i][j] = jee ;	
			}	
		else if (i<ne && j>ne-1 && ran2(&seed)<pei){
			weights[i][j] = jei ; 
			}
		else if (i>ne-1 && j<ne && ran2(&seed)<pie){
			weights[i][j] = jie ;
			}
		else if (i>ne-1 && j>ne-1 && ran2(&seed)<pii){
			weights[i][j] = jii ;	
			}				
		}
	}	


for (int i = 0; i < Ncells; i++)
	weights[i][i] = 0.0;

int counts = 0;

for(i=0;i<ne;i++)
	for(j=0;j<ne;j++)
		if (weights[i][j]!=0)
			counts += 1; 

printf("%d\n",counts );		
//..................................................
double forwardInputsE[Ncells], forwardInputsI[Ncells], forwardInputsEPrev[Ncells], forwardInputsIPrev[Ncells];
double xerise[Ncells], xedecay[Ncells], xirise[Ncells], xidecay[Ncells];
int ns[Ncells];

initialize(Ncells, forwardInputsE);
initialize(Ncells, forwardInputsI);
initialize(Ncells, forwardInputsEPrev);
initialize(Ncells, forwardInputsIPrev);
initialize(Ncells, xerise);
initialize(Ncells, xedecay);
initialize(Ncells, xirise);
initialize(Ncells, xidecay);
for(i=0;i<Ncells;i++)
	ns[i]=0;

double v[Ncells], lastSpike[Ncells];
for(i=0;i<Ncells;i++)
	v[i] = ran2(&seed);

for(i=0;i<Ncells;i++)
	lastSpike[i] = -100;

double t=0.0;
double synInput=0;
int maxTimes = (int)round(maxrate*T/1000.0);

double **times = (double **)malloc(Ncells * sizeof(double *));
for(i=0;i<Ncells;i++)
	times[i] = (double *)malloc(maxTimes* sizeof(double));

for (i = 0; i < Ncells; i++){
	for (j = 0; j < maxTimes; j++){
		times[i][j]=0;
	}
}


printf("Starting Simulation\n");

for(m=0;m<Nsteps;m++){
	t += dt;
	initialize(Ncells,forwardInputsE);
	initialize(Ncells,forwardInputsI);
	for(i=0;i<Ncells;i++){
		xerise[i] += -dt*xerise[i]/tauerise + forwardInputsEPrev[i] ;
		xedecay[i] += -dt*xedecay[i]/tauedecay + forwardInputsEPrev[i] ;
		xirise[i] += -dt*xirise[i]/tauirise + forwardInputsIPrev[i] ;
		xidecay[i] += -dt*xidecay[i]/tauidecay + forwardInputsIPrev[i] ;

		synInput = (xedecay[i] - xerise[i])/(tauedecay - tauerise) + (xidecay[i] - xirise[i])/(tauidecay - tauirise);

		if(i==1)
		{
		fprintf(fptr,"%lf\t%lf\t%lf\t%lf\n",t,(xedecay[i] - xerise[i])/(tauedecay - tauerise),(xidecay[i] - xirise[i])/(tauidecay - tauirise),synInput);
		}
		

		if (t > lastSpike[i] + refrac){  //not in refractory period
			v[i] += dt*((1/tau[i])*(mu[i]-v[i]) + synInput);
		
		
			if (v[i] > threshold[i]) { //spike occurred
				v[i] = vre;
				lastSpike[i] = t;
				ns[i] += 1;
				if (ns[i] <= maxTimes){
					times[i][ns[i]-1] = t;
				}
				for(j=0;j<Ncells;j++){
					if (weights[j][i] > 0){ //E Synapse
						forwardInputsE[j] += weights[j][i];
					}
					else if (weights[j][i] < 0){  //I synapse
						forwardInputsI[j] += weights[j][i];	
					}
				}			
			}	
		}	
	}
	for(i=0;i<Ncells;i++)
		forwardInputsIPrev[i] = forwardInputsI[i];
	for(i=0;i<Ncells;i++)
		forwardInputsEPrev[i] = forwardInputsE[i];
}

double meanE=0.0;
for(i=0;i<ne;i++)
	meanE += ns[i];

meanE = meanE*1000/(ne*T);
printf("Average Excitatory Firing Rate = %lf Hz\n",meanE);
double meanI=0.0;
for(i=ne;i<Ncells;i++)
	meanI += ns[i];

meanI = meanI*1000/(ni*T);
printf("Average Inhibitory Firing Rate = %lf Hz\n",meanI);

for(i=0;i<ne;i++){
	for(j=0;j<maxTimes;j++){
		fprintf(fp, "%d\t%lf\n", i,times[i][j]);
	}
	//fprintf(fp,"\n");
}

/*double meanmue = 0.0;
double meanmui = 0.0;
for(i=0;i<ne;i++)
	meanmue += mu[i];

for(i=ne;i<Ncells;i++)
	meanmui += mu[i];

printf("MeanMuE = %lf\t MeanMuI = %lf\n",meanmue/ne,meanmui/ni);


printf("Average Ie = %lf Hz\n",jee*ne*pee*meanE + jei*ni*pei*meanI);
printf("Average Ii = %lf Hz\n",jie*ne*pie*meanE + jii*ni*pii*meanI);*/

fclose(fp);
fclose(fptr);
printf("Plotting data..\n");

system("gnuplot -p 'plot.plt'");

return(0);
}
