/*
Stochastic simulation of a two-species population process where each population has a physical state (singleton or fused), resulting in a 4D
population process. Replication is controlled via closed-loop linear-feedback control

Based on work by Charlotte Bowles and Iain Johnston.

Fused have reduced susceptibility to degradation
Quality control in the form of selective degradation of mutants

Author: Juvid Aryaman
Date: 02/10/17

Conventions:

// Fusion
// ws + ws --> wf + wf  (0)
// ms + ms --> mf + mf (1)
// wf + ws --> wf + wf (2)
// mf + ms --> mf + mf (3)

// Fission
// wf --> ws (4)
// mf --> ms (5)

// Replication
// wf --> 2wf (6)
// mf --> 2mf (7)
// ws --> 2wf (8) (NEW)
// ms --> 2mf (9) (NEW)

// Degradation
// wf --> 0 (10)
// mf --> 0 (11)
// ws --> 0 (12)
// ms --> 0 (13)

// Physical cross-processes (fusion)
// wf + ms --> wf + mf (14)
// mf + ws --> mf + wf (15)
// ws + ms --> wf + mf (16)

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define RND drand48() //to generate random numbers

//#define PRINTTIME

// Generate a single N(0,1) variable using Box-Muller transformation
double Rand_n(){ 
	double U1, U2;
	static double Z0, Z1; // static retains value between calls
	static unsigned int call = 0; // Variable to flip-flop upon each call
	//printf("call: %d\n", call);
	if(call == 1){ // spit out stored normal
		call = !call; // switch the state of call, so next time we calculate
		return Z1;
	}
	U1 = RND;
	U2 = RND;
	Z0 = sqrt( -2 * log(U1)) * cos(2*M_PI*U2);
	Z1 = sqrt( -2 * log(U1)) * sin(2*M_PI*U2);
	//printf("Generation:\n");
	//printf("%f, %f\n",Z0, Z1 );
	call = !call;
	return Z0;
}

double KahanSum_fly(double* Sum, double number, double* c){
	//I got this function from the Wikipedia page for the Kahan Summation algorithm. This function returns the sum of Sum and number
	//and it keeps track of erros bla through c and stuff. See the wikipedia page. 
	double y, t;
	y = number - *c;
	t = *Sum + y;
	*c = (t - *Sum) - y; 
	*Sum = t; 
	return *Sum;
}

double replicationRate(int w, int m, double kappa, double b, double delta, double mu){
	double reprate;
	reprate = mu + b*(kappa - w - delta*m);

	if (reprate < 0) {reprate = 0;}
	return reprate;
}

void Print_array(double *a, int num_elements){
	int i;
	for (i = 0; i < num_elements; i++)
	{
		printf("%d: %f\n", i, a[i]);
	}
	printf("\n");
}

int update_state(int j, int* w_f, int* m_f, int* w_s, int* m_s){

	// Physical processes
	if(j == 0){*w_s -= 2; *w_f += 2;} // ws + ws --> wf + wf 
	else if(j ==1){*m_s -= 2; *m_f += 2;} // ms + ms --> mf + mf
	else if(j == 2){*w_s -= 1; *w_f += 1;} // wf + ws --> wf + wf
	else if(j == 3){*m_s -= 1; *m_f += 1;} // mf + ms --> mf + mf
	else if(j == 4){*w_f -= 1; *w_s += 1;} // wf --> ws
	else if(j == 5){*m_f -= 1; *m_s += 1;} // mf --> ms

	// Genetic processes
	else if(j == 6){*w_f += 1;} // wf --> 2wf 		
	else if(j == 7){*m_f += 1;} // mf --> 2mf 		
	else if(j == 8){*w_s -= 1; *w_f += 2;} // ws --> 2wf 		
	else if(j == 9){*m_s -= 1; *m_f += 2;} // ms --> 2mf 
	//else if(j == 8){*w_s += 1;} // old
	//else if(j == 9){*m_s += 1;} // old
	else if(j == 10){*w_f -= 1;	} // wf --> 0		
	else if(j == 11){*m_f -= 1;} // mf --> 0		
	else if(j == 12){*w_s -= 1; } // ws --> 0		
	else if(j == 13){*m_s -= 1;} // ms --> 0

	// Physical cross-processes
	else if(j == 14){*m_s -= 1; *m_f += 1;} // wf + ms --> wf + mf
	else if(j == 15){*w_s -= 1; *w_f += 1;} // mf + ws --> mf + wf
	else if(j == 16){*w_s -= 1; *m_s -= 1; *w_f += 1; *m_f += 1;} // ws + ms --> wf + mf
	else{printf("j out of bounds\n"); return -1;}

	return 1; // exit ok
}

int main(int argc, char *argv[]) {

clock_t begin = clock();

int initial_seed;
int param_block;
int n_iters_per_block;
int seed_inc;

// Declare model parameters to be taken from command line
double gamma;
double beta;
double kappa;
double b;
double mu;
double xi;
double delta;
double Q;

int n_target;

double h_init_mean;

if(argc == 3){
	printf("TESTING MODE\n");
	initial_seed = atoi(argv[1]); 
	seed_inc = atoi(argv[2]); 

	printf("Base seed = %d\n", initial_seed);
	printf("Seed inc = %d\n", seed_inc);

	// See /Det_sim/Closed_loop_control/Singleton_birth_fus/sing_birth_fus.ipynb
	
	h_init_mean = 0.5;
	n_target = 1000;

	gamma = 0.03785142857142858;
	beta = 33.12; 
	kappa = 11.662903457629223;
	b = 1.2416523075924095e-05;
	mu = 0.023;
	xi = 0.0; // fused susceptibility to degradation
	delta = 1.0;
	Q = 1.0;

	param_block = 1;


	srand48(initial_seed+seed_inc);

}
else if(argc == 15){

	initial_seed = atoi(argv[1]); 
	param_block = atoi(argv[2]); 
	n_iters_per_block = atoi(argv[3]);
	seed_inc = atoi(argv[4]);

	h_init_mean = atof(argv[5]);
	n_target = atoi(argv[6]);
	
	xi = atof(argv[7]);
	beta = atof(argv[8]);
	gamma = atof(argv[9]);
	kappa = atof(argv[10]);
	b = atof(argv[11]);
	mu = atof(argv[12]);
	delta = atof(argv[13]);
	Q = atof(argv[14]);

	srand48(initial_seed+param_block*n_iters_per_block + seed_inc);

}
else{
	printf("argc = %d\n", argc);
	printf("Usage: %s  initial_seed param_block\n", argv[0]);
	printf("       %s  initial_seed param_block n_iters_per_block seed_inc h_init_mean n_target xi beta gamma kappa b mu delta  Q\n", argv[0]);
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Sim Parameters //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double TFinal = 5000; // total amount of time to run simulation (days)

double DELTAT_STORE = 50.0; // Time lattice on which to store events
int NREPEATS = 10; // Number of repeats of the stochastic simulation


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Define some variables ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//define some variables.

FILE *fpStat; //to write stats to file

int i, j, k, w_s, w_f, m_s, m_f, nStore; 
double t, dt, t_store;
double correct;

int nTimeSteps = (int)(TFinal/DELTAT_STORE)+1;

// ICs
int ws_init = (int)round(n_target*(1.0 - h_init_mean)/2.0);
int wf_init = (int)round(n_target*(1.0 - h_init_mean)/2.0);
int ms_init = (int)round(n_target*h_init_mean/2.0);
int mf_init = (int)round(n_target*h_init_mean/2.0);

int nreactions = 17;
double hazards[nreactions]; 

double totalrate, rate_sum, r1, r2;

int number, index, ref; 



printf("h_init_mean = %f\n",h_init_mean);
printf("w0_s = %d,m0_s = %d,w0_f = %d,m0_f = %d\n", ws_init, ms_init, wf_init, mf_init);
printf("xi = %f\n",xi);
printf("beta = %f\n",beta);
printf("gamma = %f\n",gamma);
printf("kappa = %f\n",kappa);
printf("b = %f\n",b);
printf("mu = %f\n",mu);
printf("delta = %f\n",delta);
printf("Q = %f\n",Q);


// Containers for data on discrete time lattice
int* wsStore;
int* msStore;
int* wfStore;
int* mfStore;

wsStore = (int*)malloc(sizeof(int)*nTimeSteps*NREPEATS);
msStore = (int*)malloc(sizeof(int)*nTimeSteps*NREPEATS);
wfStore = (int*)malloc(sizeof(int)*nTimeSteps*NREPEATS);
mfStore = (int*)malloc(sizeof(int)*nTimeSteps*NREPEATS);

int return_status;

//int debug_counter = 0;

int large_n_flag = 0;

double rep_rate;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// PRINT OUTPUTS ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//this stores all statistics we look at, except for the distributions in the bins and the times at which first hit certain h.
char str[500];
sprintf(str, "./output_%d_%d.txt",param_block,seed_inc);
fpStat = fopen(str, "w");

fprintf(fpStat,"NREPEATS = %d\n", NREPEATS);
fprintf(fpStat,"w0_s = %d,m0_s = %d,w0_f = %d,m0_f = %d\n", ws_init, ms_init, wf_init, mf_init);
fprintf(fpStat,"xi = %f\n",xi);
fprintf(fpStat,"beta = %f\n",beta);
fprintf(fpStat,"gamma = %f\n",gamma);
fprintf(fpStat,"kappa = %f\n",kappa);
fprintf(fpStat,"b = %f\n",b);
fprintf(fpStat,"mu = %f\n",mu);
fprintf(fpStat,"delta = %f\n",delta);
fprintf(fpStat,"Q = %f\n",Q);
fprintf(fpStat,"TFinal = %f\n", TFinal);
fprintf(fpStat, "t,rep,ws,ms,wf,mf,h,n\n");
fclose(fpStat); 

#ifdef PRINTTIME
	clock_t loop_time = clock();
#endif

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// START GILLESPIE /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double h;
int n;

for (i=0; i< NREPEATS; i++){
	large_n_flag = 0;

	printf("i = %d\n", i);
	//INITIALIZE

	w_s = ws_init;
	w_f = wf_init;
	m_s = ms_init;
	m_f = mf_init;

	t = 0.0; dt = 0.0; correct = 0; //correct is the correction term for adding times (using Kahan summation)
	index=0; //this is the counter of ws and ms. It increases by 1 every time we store a value.
	t_store = 0.0;//the next time at which we store w, m values
	nStore = 1; //keeps track of #times we stored values, start with 1.
	
	//START WHILE LOOP
	while (t < TFinal){

		//Update rates.

		totalrate = 0.0; // reset total rate
		
		// Fusion
		hazards[0] = w_s * (w_s - 1) * gamma/2.0; // ws + ws --> wf + wf 
		hazards[1] = m_s * (m_s - 1) * gamma/2.0; // ms + ms --> mf + mf
		hazards[2] = w_f * w_s * gamma; // wf + ws --> wf + wf
		hazards[3] = m_f * m_s * gamma; // mf + ms --> mf + mf

		// Fission
		hazards[4] = w_f * beta; // wf --> ws
		hazards[5] = m_f * beta; // mf --> ms

		// Replication
		rep_rate =replicationRate(w_f + w_s, m_f + m_s, kappa, b, delta, mu);
		hazards[6] = w_f * rep_rate; // wf --> 2wf
		hazards[7] = m_f * rep_rate; // mf --> 2mf
		hazards[8] = w_s * rep_rate; // ws --> 2wf
		hazards[9] = m_s * rep_rate; // ms --> 2mf

		// Degradation
		hazards[10] = w_f * xi * mu; // wf --> 0
		hazards[11] = m_f * xi * mu; // mf --> 0
		hazards[12] = w_s * mu; // ws --> 0
		hazards[13] = m_s * mu * Q; // ms --> 0

		// Physical cross-processes (fusion)
		hazards[14] = w_f * m_s * gamma; // wf + ms --> wf + mf
		hazards[15] = m_f * w_s * gamma; // mf + ws --> mf + wf
		hazards[16] = w_s * m_s * gamma; // ws + ms --> wf + mf

		for (j = 0; j < nreactions; j++){totalrate += hazards[j];} // Find total hazard
		 
		//update time 
		r1 = RND; //random number between 0 and 1
		dt = -log(r1)/totalrate; //time increment for Gillespie

		//printf("Mean time to next event (min) = %f\n", 60*24./totalrate);

		if (dt<0) {printf("Error, negative time intervals.., i = %d\ntotalrate is %.5f\n", i,totalrate); return -1;}
		t = KahanSum_fly(&t, dt, &correct);

		//store state on time lattice
		if (t>(t_store - 1E-12)){
			number = floor((t - t_store)/(double)DELTAT_STORE);
			for (k=0;k<(number+1);k++){	
				if (index >= (int)round(((double)(TFinal+DELTAT_STORE)/DELTAT_STORE))) printf("wtf2, %d,%.4f, %d\n", index,t_store, k); //to debug.

				//store things. First time we record here is at t=0 (well, t is then larger but the state is still from t=0)
				ref = i*nTimeSteps + index; //where to print. Print for all runs.

				wsStore[ref] = w_s;
				msStore[ref] = m_s;
				wfStore[ref] = w_f;
				mfStore[ref] = m_f;
				
				n = w_s + w_f + m_s + m_f;
				h = (m_s + m_f)/((double)(n));

				// Append to file
				fpStat = fopen(str, "a");
				fprintf(fpStat, "%f,%d,%d,%d,%d,%d,%f,%d\n", index*DELTAT_STORE, i, wsStore[ref], msStore[ref],wfStore[ref],mfStore[ref], h, n);
				fclose(fpStat); 

				#ifdef PRINTTIME
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

		//check times at which first hit a certain h. Use 't - dt' as time, since t also includes the next event, whereas hnow is from t - dt
		//Do this before we update the state, otherwise w+m might be zero and get undefined h. Now if w+m reaches zero we break before how is calculated here.
		if (w_f + m_f + w_s + m_s == 0) {printf("wtf5, we went extinct!\n"); return -1;}

		if ((w_f + m_f + w_s + m_s > 2000) && (large_n_flag == 0)) {printf("warning, n > 2000\n");large_n_flag =1;}

		//update state 
		r2 = RND;
		rate_sum = 0.0;
		for (j = 0; j < nreactions; j++){
			rate_sum += hazards[j];
			if(r2 <= rate_sum / totalrate){					

				return_status = update_state(j, &w_f, &m_f, &w_s, &m_s);

				/*printf("j = %d\n", j);
				debug_counter += 1;
				if(debug_counter > 50){
					printf("dbg = %d\n", debug_counter);
					return 0;
				}*/

				if(w_f < 0 || m_f < 0 || w_s < 0 || m_s < 0){
					printf("Negative species!\n");
					printf("j = %d\n", j);
					printf("w_f = %d, m_f = %d, w_s = %d, m_s = %d\n", w_f, m_f, w_s, m_s);
					printf("hazards:\n");
					Print_array(hazards, nreactions);
					return -1;
				}
				break;
			}
			if(j == nreactions - 1){printf("No reactions selected!\n");return -1;}
		}
		
		

		if (return_status < 0){printf("update_state returned error!\n"); return -1;}

		//check for extinction. If extinct, print zeros etc.. to all store arrays and break this while loop.
		if (w_f + m_f + w_s + m_s == 0){
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
				wsStore[ref] = 0;
				msStore[ref] = 0;
				wfStore[ref] = 0;
				mfStore[ref] = 0;

				// Append to file
				fpStat = fopen(str, "a");
				fprintf(fpStat, "%f,%d,%d,%d,%d,%d,%s,%d\n", index*DELTAT_STORE, i, wsStore[ref], msStore[ref],wfStore[ref],mfStore[ref], "nan", 0);
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