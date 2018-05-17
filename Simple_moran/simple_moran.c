/*
A simple Moran process in continuous time

Author: Juvid Aryaman

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>	


#define RND drand48() //to generate random numbers

#define DUMMY


void Print_array(double *a, int num_elements){
	int i;
	for (i = 0; i < num_elements; i++)
	{
		printf("%d: %f\n", i, a[i]);
	}
	printf("\n");
}


int main(int argc, char *argv[]) {

clock_t begin = clock();

int initial_seed;


if(argc == 2){
	initial_seed = atoi(argv[1]); 	
	srand48(initial_seed);
}
else{
	printf("argc = %d\n", argc);
	printf("Usage: %s  initial_seed\n", argv[0]);	
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Sim Parameters //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double TFinal = 100; // total amount of time to run simulation (days)

double DELTAT_STORE = 10.0; // Time lattice on which to store events
int NREPEATS = 10000; // Number of repeats of the stochastic simulation

int n = 10000; // number of particles in the Moran model
double mu = 0.023; // death rate
double h0 = 0.4; // initial heteroplasmy

double totalrate = n*mu;
double nsq = (double)(n*n);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Define some variables ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

FILE *fpStat; 

int i, k, nStore, index, number, ref; 
double t, dt, t_store, r;

int nTimeSteps = (int)(TFinal/DELTAT_STORE)+1;



// Containers for data on discrete time lattice
int* m_Store;
m_Store = (int*)malloc(sizeof(int)*nTimeSteps*NREPEATS);


#ifdef TEST
	clock_t loop_time = clock();
#endif


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// PRINT OUTPUTS ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

char str[500];
//sprintf(str, "./output%07d.txt",seed_inc);

sprintf(str, "./output.txt");
fpStat = fopen(str, "w");

fprintf(fpStat,"NREPEATS = %d\n", NREPEATS);
fprintf(fpStat, "n = %d\n",n);
fprintf(fpStat, "mu = %f\n",mu);
fprintf(fpStat, "h0 = %f\n",h0);
fprintf(fpStat,"TFinal = %f\n", TFinal);
fprintf(fpStat, "t,rep,m,h\n");
fclose(fpStat); 

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// START SIM ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double h, thresh_pm_m;
int m;

for (i=0; i< NREPEATS; i++){
	if(i%100 == 0){
		printf("i = %d\n", i);
	}

	
	//INITIALIZE
	m = (int)round(n*h0);
	
	t = 0.0; dt = 0.0; 
	index = 0; //this counter increases by 1 every time we store a value.
	t_store = 0.0;//the next time at which we store m values
	nStore = 1; //keeps track of #times we stored values, start with 1.
	
	//START WHILE LOOP
	while (t < TFinal){

		//update time 
		r = RND; //random number between 0 and 1
		dt = -log(r)/totalrate; //time increment for Gillespie

		t += dt;


		//store state on time lattice
		if (t>(t_store - 1E-12)){
			number = floor((t - t_store)/(double)DELTAT_STORE);
			for (k=0;k<(number+1);k++){	
				if (index >= (int)round(((double)(TFinal+DELTAT_STORE)/DELTAT_STORE))) printf("wtf2, %d,%.4f, %d\n", index,t_store, k); //to debug.

				//store things. First time we record here is at t=0 (well, t is then larger but the state is still from t=0)
				ref = i*nTimeSteps + index; //where to print. Print for all runs.

				m_Store[ref] = m;
								
				h = m/((double)(n));

				// Append to file
				fpStat = fopen(str, "a");
				fprintf(fpStat, "%f,%d,%d,%f\n", index*DELTAT_STORE, i, m_Store[ref], h);
				fclose(fpStat); 

				#ifdef TEST
					loop_time = clock();
					printf("t=%f, simtime = %f\n", index*DELTAT_STORE, (double)(loop_time - begin) / CLOCKS_PER_SEC);
				#endif

				//update storing time and printing number
				index ++; //After this line, index gives the number of times we have stored.
				t_store = (nStore*DELTAT_STORE);
				nStore += 1;
				
				//this if statement below here is to make sure that the last time we store is really
				// at TFinal, so if we reach TFinal + DELTAT_STORE then we should break.
				if (t_store < TFinal + DELTAT_STORE + (double)DELTAT_STORE/10.0 && t_store > TFinal+DELTAT_STORE - (double)DELTAT_STORE/10.0) break;
			}
		}	

		
		//update state 
		thresh_pm_m = m*(n-m)/nsq;
		r = RND;
		if ( r < thresh_pm_m){
			m += 1;
		}
		else if (r < 2.0*thresh_pm_m){
			m -= 1;
		}
		// otherwise leave m as is
		
		
		//check for absorption
		if (m == 0){
			//first check whether we're already at the end. This can happen if above here we stored multiple things (i.e. number > 1) and we went extinct as well.
			if (t_store < TFinal + DELTAT_STORE + (double)DELTAT_STORE/10.0 && t_store > TFinal+DELTAT_STORE - (double)DELTAT_STORE/10.0) { 
				if (index != nTimeSteps) {printf("wtf6\n"); return -1;}
				break; //break out of while loop
			}				
			//fill up everything with zeros and go to next iteration
			t_store = TFinal + (double)DELTAT_STORE/10.0; //so that we print up until TFinal, make t_store slightly above TFinal
			number = floor((t_store - t)/(double)DELTAT_STORE) + 5; //just do some big number of times, have break statement later
			for (k=0;k<(number+1);k++){	
				if (index >= (int)round(((double)(TFinal+DELTAT_STORE)/DELTAT_STORE))){printf("wtf3, %d,%.7f, %d, %.10f, %d, %d, %.7f\n", index,t_store, k, t, (int)round(((double)(TFinal+DELTAT_STORE)/DELTAT_STORE)), number, (double)TFinal/DELTAT_STORE + 1);} //to debug...
				
				//store things
				ref = i*nTimeSteps + index; //where to print. Print for all runs.

				m_Store[ref] = 0;

				// Append to file
				fpStat = fopen(str, "a");
				fprintf(fpStat, "%f,%d,%d,%f\n", index*DELTAT_STORE, i, m_Store[ref], 0.0);				
				fclose(fpStat); 
				
				index ++;

				if (index >= (int)round((double)TFinal/DELTAT_STORE + 1)) {break;} 
			}
			break;	//break out of while loop, go to next iteration
		}//close check extinct loop


	} //close while loop
} //close NREPEATS for loop


clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

printf("Total time taken = %.3f\n", time_spent);

return 0; //return from main

} //close main