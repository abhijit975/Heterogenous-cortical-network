#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <vector>

#include "ran2.cpp"
static int long seed = time(NULL);  //seed for random number generator

//Define network size.....................................................................................................................

#define ne 4000 	//Excitatory neurons
#define ni 1000     //Inhibitory neurons
//.........................................................................................................................................

//Define connection probabilities..........................................................................................................

#define K 800      //Number of excitatory connections for each E-neuron
#define pee 0.2
#define pei 0.5
#define pie 0.5
#define pii 0.5
//.........................................................................................................................................

//Define Timescales........................................................................................................................

#define taui 10.0
#define taue 15.0
#define	tauerise 1.0
#define	tauedecay 3.0
#define	tauirise 1.0
#define	tauidecay 2.0
#define	refrac 5.0
//.........................................................................................................................................

//Define clustering parameters.............................................................................................................

#define ratiojee 1.9
#define ratiopee 4.0
//.........................................................................................................................................

//Define Voltages..........................................................................................................................

#define	vre 0.			// refractory potential
#define	threshe 1.      // Threshold potential for E-neuron
#define	threshi 1.		// Threshold potential for I-neuron
//.........................................................................................................................................

#define maxrate 100.

//Define voltage biases....................................................................................................................
#define	muemin 1.1
#define	muemax 1.2
#define	muimin 1
#define	muimax 1.05
//.........................................................................................................................................

//Function to initialize 1D arrays.........................................................................................................

void initialize(int N, double array[])
{
	int i;
	for(i=0;i<N;i++){
		array[i] = 0.0;
	}	
}
//.........................................................................................................................................

//Starting Main Program....................................................................................................................

int main()
{

//system("/Applications/Mathematica.app/Contents/MacOS/MathKernel -script pop.m");  //runs a mathematica script pop.m to generate a power law distributed population
//comment the top line to generate a population manually each time.

//Opening Files to write spike timings, E & I synaptic input current, and read the powerlaw distributed population.........................

FILE *fp,*fptr,*fptr1;
fp = fopen("./spike_eps.dat","w");
fptr1 = fopen("./current_power.dat","w");  //writes E & I synaptic input currents for a neuron in population 1,2, ncl-2 and ni
fptr = fopen("./popexp.dat","r"); 
//.........................................................................................................................................

//defining connection strengths............................................................................................................

printf("Creating Parameters..\n");
double jee_in, jee_out, jei, jie, jii, peeout, peein;
int Ncells = ne + ni;
double k = pow(1.0*K,0.5);

jee_out = 100.0*taui/(k*taue);
//jee_out = 2.357023;
jee_in = jee_out/1.8;
//jee_in = jee_out*ratiojee ;  			//uncomment this line and comment the above line to get an unbalanced network with jeein greater than jeeout
jie = 100.0*(pee/pie)*(1./(k));
//jei = -0.045255;
jei = -1.2*(ne/ni)*(pee/pie)*(taui/taue)*(1./k) ;
jii = - (ne/ni)*(pee/pie)*(1./k) ;

printf("jee = %lf\t jei = %lf\t jie = %lf\t jii = %lf\n",jee_out,jei,jie,jii);
//.........................................................................................................................................

// Runtime and timesteps...................................................................................................................
double dt = 0.1;
double T = 5000;  				// Runtime for the simulation
int Nsteps = (int)round(T/dt);

int Nstim = 100 ;
double stimstr = 0.0/taue ;  // use 0.1/taue for stimulation strength
double	stimstart = 500 ;
double	stimend = 2500 ;            //uncomment this part to provide external simulation
//.........................................................................................................................................

//making arrays of bias voltages, timescales and threshold voltage values..................................................................

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
//.........................................................................................................................................

//Defining the connection matrix...........................................................................................................

double **weights = (double **)malloc(Ncells * sizeof(double *));
for(i=0;i<Ncells;i++)
	weights[i] = (double *)malloc(Ncells* sizeof(double));

for(i=0;i<Ncells;i++)
{
	for(j=0;j<Ncells;j++){
		weights[i][j] = 0;}}

// Reading and storing populations distributed in clusters.................................................................................

int val,len = 0;
while(1>0)
{
fscanf(fptr,"%d\n",&val);
len = len+1;
if( feof(fptr) ) 
    break ;
}
fclose(fptr);
fptr = fopen("./popexp.dat","r");
int ncl = len;
int nec[ncl];
for(i=0;i<ncl;i++)
	fscanf(fptr,"%d\n",&nec[i]);

double **W = (double **)malloc((ncl+1) * sizeof(double *));
for(i=0;i<ncl+1;i++)
	W[i] = (double *)malloc((ncl+1)* sizeof(double));
fclose(fptr);

//fptr = fopen("./pop.dat","w"); //rewriting the population size file with pop. size and corresponding firing rates. 
//.........................................................................................................................................

int popstart, popend, temp = 0;
int kpop, counter = 0;

//E-E Connection matrix elements (Ne X Ne) matrix..........................................................................................

for(int p = 0; p < ncl; p++)
{
	popstart = temp ;
	popend = temp + nec[p];  	//popstart and popend indices of neuron for beginning and ending of a cluster
	temp += nec[p] ;
	
	peeout = K/(nec[p]*(ratiopee-1) + ne);
	//peeout = pee;				// uncomment this line and comment above line will give fixed p_ee value for each cluster but degree of each cluster will be different. 
	peein = ratiopee*peeout;

//E-E within cluster connection................................................................................................................

	for(i = popstart; i<popend; i++)
	{
		for(j = popstart; j<popend; j++)
		{
			if (ran2(&seed)<peein)
			{
				if(nec[p]<25)
				{
				weights[i][j]= 	jee_out*ratiojee/100.0;
				}
				else{
				weights[i][j] = jee_in*pee/(nec[p]*peeout);}
				//weights[i][j] = jee_in/nec[p];				//uncomment this line and comment above line if peeout = pee (fixed) or don't want symmetric matrix
			}
		}
	}
//.........................................................................................................................................

//Defining Weight matrix for eigenvalue calculation........................................................................................
/*
	W[p][p] = jee_in*pee*peein/peeout;
	//W[p][p] = (jee_in*peein);   		//uncomment this line and comment above line if using weights[i][j] = jee_in/nec[p]
	W[ncl][p] = jie*pie;
	W[p][ncl] = jei*pei*ni;
	W[ncl][ncl] = jii*pii*ni;

	for(int l = 0; l<ncl; l++)
	{
		if(l!=p){
			W[l][p] = jee_out*pee;
			//W[l][p] = jee_out*peeout;		//uncomment this line and comment above line if using weights[i][j] = jee_out/nec[p]
			}
	}*/
//.........................................................................................................................................

//E-E outside cluster connection...........................................................................................................

	counter = popend;
	for(i = popstart; i<popend; i++)
	{
		counter = popend;
		for(kpop = p+1; kpop < ncl; kpop++)
		{
			for(j=counter; j< counter + nec[kpop]; j++)
			{
				if (ran2(&seed)<peeout)
				{
					if(nec[kpop]<25)
					{
						weights[i][j] = jee_out/100.0;
					}
					else{
					weights[i][j] = jee_out*pee/(peeout*nec[kpop]);}
					//weights[i][j] = jee_out/nec[kpop];		//uncomment this line and comment above line if peeout = pee (fixed) or don't want symmetric matrix
				}	
			}
			counter += nec[kpop];		
		}
	}
	counter = 0;
	if(p!=0)
	{	
		for(i = popstart; i<popend; i++)
		{
			counter = 0;
			for(kpop = 0; kpop< p; kpop++)
			{
				for(j = counter; j < counter + nec[kpop]; j++)
				{
					if (ran2(&seed)<peeout)
					{
					if(nec[kpop]<25)
					{
						weights[i][j] = jee_out/100.0;
					}
					else{
					weights[i][j] = jee_out*pee/(peeout*nec[kpop]);}
					}	
				}
				counter += nec[kpop];	
			}		
		}
	}
}
//.........................................................................................................................................	
//End of E-E connection matrix elements definition.........................................................................................

//Defininig E-I and I-I connection elements................................................................................................

for(i=0;i<Ncells;i++)
{	
	for(j=0;j<Ncells;j++)
	{
		if (i<ne && j>ne-1 && ran2(&seed)<pei){
			weights[i][j] = jei ; 
		}
		else if (i>ne-1 && j>ne-1 && ran2(&seed)<pii){
			weights[i][j] = jii ;	
		}				
	}
}	
//.........................................................................................................................................

//Defining I-E connection elements.........................................................................................................

for(i=ne;i<Ncells;i++)
{
	counter = 0 ;
	for(kpop=0;kpop<ncl;kpop++)
	{
		for(j=counter; j< counter + nec[kpop]; j++)
		{
			if (ran2(&seed)<pie)
			{
				if(nec[kpop]<25)
				{
					weights[i][j] = jie/100.0;
				}
				else{
				weights[i][j] = jie/nec[kpop];}
			}
		}
	counter += nec[kpop] ; 	
	}
}
//.........................................................................................................................................

for (int i = 0; i < Ncells; i++)
	weights[i][i] = 0.0;			//Making diagonals zero.

int counts = 0;

for(i=0;i<ne;i++)
	for(j=0;j<ne;j++)
		if (weights[i][j]!=0)
			counts += 1; 		

printf("%d\n",counts );				//counting and printing the total number of excitatory connections (should be Ne*K if fixed K is used)
//.........................................................................................................................................

//Defining and initializing arrays required for the simulation.............................................................................

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
int maxTimes = (int)round(maxrate*T/1000.0);      		//Maximum number of spikes allowed

double **times = (double **)malloc(Ncells * sizeof(double *));
for(i=0;i<Ncells;i++)
	times[i] = (double *)malloc(maxTimes* sizeof(double));

for (i = 0; i < Ncells; i++){
	for (j = 0; j < maxTimes; j++){
		times[i][j]=0;								//This array writes the spike time of the neurons
	}
}
//.........................................................................................................................................

//Starting Integration of the differential equation........................................................................................

printf("Starting Simulation\n");

for(m=0;m<Nsteps;m++)
{
	t += dt;
	initialize(Ncells,forwardInputsE);
	initialize(Ncells,forwardInputsI);
	for(i=0;i<Ncells;i++)
	{
		xerise[i] += -dt*xerise[i]/tauerise + forwardInputsEPrev[i] ;
		xedecay[i] += -dt*xedecay[i]/tauedecay + forwardInputsEPrev[i] ;
		xirise[i] += -dt*xirise[i]/tauirise + forwardInputsIPrev[i] ;
		xidecay[i] += -dt*xidecay[i]/tauidecay + forwardInputsIPrev[i] ;

		synInput = (xedecay[i] - xerise[i])/(tauedecay - tauerise) + (xidecay[i] - xirise[i])/(tauidecay - tauirise);		//synaptic input for the neurons

// Writing E & I contribution of synaptic currents for different neurons in a data file....................................................

/*		if(i==1){
			fprintf(fptr1, "%lf\t%lf\t%lf\t%lf\t",t,(xedecay[i] - xerise[i])/(tauedecay - tauerise),(xidecay[i] - xirise[i])/(tauidecay - tauirise),synInput);
		}
		else if(i==nec[0]+50){
			fprintf(fptr1, "%lf\t%lf\t%lf\t",(xedecay[i] - xerise[i])/(tauedecay - tauerise),(xidecay[i] - xirise[i])/(tauidecay - tauirise),synInput);
		}
		else if(i==ne-nec[ncl-2]){
			fprintf(fptr1, "%lf\t%lf\t%lf\t",(xedecay[i] - xerise[i])/(tauedecay - tauerise),(xidecay[i] - xirise[i])/(tauidecay - tauirise),synInput);
		}
		else if(i==ne+1){
			fprintf(fptr1, "%lf\t%lf\t%lf\n",(xedecay[i] - xerise[i])/(tauedecay - tauerise),(xidecay[i] - xirise[i])/(tauidecay - tauirise),synInput);
		}*/   
//.........................................................................................................................................		

		if (i<250 && t > stimstart && t < stimend){ 
			synInput += stimstr; 
		}

		if(i>509 && i<634 && t>2500 && t<4000){
     			synInput += stimstr; 
		}			//uncomment this part if there is external simulation (***********************)
		int newcount = 0;		
		for(int q = 0;q<ncl;q++){
			newcount += nec[q];
			if(i==newcount){
				fprintf(fptr1, "%lf\t",synInput);
			}}
		if(i==ne+1){
			fprintf(fptr1, "%lf\t",synInput);}
		
		if (t > lastSpike[i] + refrac){  //not in refractory period
			v[i] += dt*((1/tau[i])*(mu[i]-v[i]) + synInput);
		
			if (v[i] > threshold[i]) { //spike occurred
				v[i] = vre;
				lastSpike[i] = t;
				ns[i] += 1;
				if (ns[i] <= maxTimes){
					times[i][ns[i]] = t;
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
	fprintf(fptr1,"\n");
}
//End of numerical integration.............................................................................................................

//Calculating the mean firing rates.........................................................................................................

double mean=0.0;
for(i=0;i<ne;i++)
	mean += ns[i];

mean = mean*1000/(ne*T);
printf("Average Excitatory Firing Rate = %lf Hz\n",mean);
mean=0.0;
for(i=ne;i<Ncells;i++)
	mean += ns[i];

mean = mean*1000/(ni*T);
printf("Average Inhibitory Firing Rate = %lf Hz\n",mean);
//.........................................................................................................................................

//Writing spike timings in data file and firing rates for individual cluster in another data files.........................................
fptr = fopen("ratepartial.dat","w");
counter = 0; double rate;
for(kpop = 0; kpop < ncl; kpop++)
{
	rate = 0.0;
	for(i = counter; i< counter + nec[kpop]; i++)
	{
		for(j = 0; j<maxTimes; j++)
		{
			if(times[i][j]!=0){
				fprintf(fp, "%d\t%lf\n", i,times[i][j]);
				rate += 1; 
			}
		}
	}
	counter += nec[kpop];
	fprintf(fp,"\n\n");
	rate = rate*1000/(nec[kpop]*T);
	fprintf(fptr,"%.3f\n",rate);
}
fprintf(fptr,"%.3f\n",mean);
for(i=ne;i<Ncells;i++)
{
	for(j=0;j<maxTimes;j++)
	{
		fprintf(fp, "%d\t%lf\n", i,times[i][j]);
	}
}
//.........................................................................................................................................

fclose(fp);
fclose(fptr);
fclose(fptr1);

//Plotting data and writing W matrix for eigenvalue calculation............................................................................
printf("Plotting data..\n");
system("gnuplot -p 'plot.plt'");           // To plot with alternate color use 'plot.plt'
/*printf("Writing balance matrix W to file... \n");
FILE *fp1;
fp1 = fopen("./balance_W.dat","w");   //Import this file to mathematica as 'Table' and then find eigenvalues
for(i=0;i<ncl+1;i++){
	for(j=0;j<ncl+1;j++){
		fprintf(fp1,"%.6f\t",W[i][j]);
	}
	fprintf(fp1,"\n");
}
fclose(fp1);*/
//.........................................................................................................................................
return(0);
}
//End of Program...........................................................................................................................
