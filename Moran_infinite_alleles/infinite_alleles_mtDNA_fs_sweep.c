/*
Simulating the infinite alleles model with Moran steps

5.6 × 10− 7 mutations/bp/doubling
https://www.sciencedirect.com/science/article/pii/S0005272815001103#bb0440

Author: Juvid Aryaman

Date: 18/05/18
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#define RND gsl_rng_uniform(r) // to generate random numbers

// Number of individuals in the population
#define N_POP 1000
// Since we're not using dynamic arrays, set the maximum number of mutations allowed
#define MAX_MUTNS 1000
// LEN_GENOME sets the probability distribution of the number of mutations per individual per copying event
#define LEN_GENOME 16569

#define DUMMY

void Print_int_array(int *a, int num_elements){
	/* Print a 1D int array
	
	Inputs
	--------
	a, 2D array
	num_elements, number of elements
	*/
	for (int i = 0; i < num_elements; i++) printf("%d\n", a[i]);	
	printf("\n");
}

void Print_2d_int_array(int **a, int dm_1, int dm_2){	
	/* Print a 2D int array
	
	Inputs
	--------
	a, 2D array
	dm_1, number of rows (dimension 1)
	dm_2, number of columns (dimension 2)
	*/
	for (int i = 0; i < dm_1; i++){
		for (int j = 0; j < dm_2; j++){
			printf("%d ", a[i][j]);
			if (j == dm_2 - 1) printf("\n");			
		}
	}
	printf("\n");
}

void Copy_state(int **state, int index_born, int index_die){
	/* Copy the state of an individual born to an individual which will die

	Inputs
	--------
	state, the state of the system
	index_born, the index of the individual to be born
	index_die, the index of the individual to die
	*/
	if(index_born > N_POP-1 || index_die > N_POP-1){
		printf("index born/die exceeds bounds!\n");
		exit(99);
	}
	for (int i = 0; i < MAX_MUTNS; i++){
			if((state[i][index_die] == 0) && (state[i][index_born] == 0)){return;} // finished copying
			else{state[i][index_die] = state[i][index_born];} // overwrite the indiv to die
		}
}

void Moran_update_state(int **state, int index_born, int index_die){
	/* Update state according to a Moran step where an individual is chosen to be born and to die

	Inputs
	--------
	state, the state of the system
	index_born, the index of the individual to be born
	index_die, the index of the individual to die
	*/
	if(index_born == index_die) return;	
	else Copy_state(state, index_born, index_die);	
}

void Moran_update_state_mutate(int **state, int index_born, int index_die, int *mut_counter, int n_mutations){
	/* Update state according to a Moran step, where the individual overwritten (index_die) incurrs a number of mutations
	
	Inputs
	--------
	state, the state of the system (2D array)
	index_born, the index of the individual to be born
	index_die, the index of the individual to die
	mut_counter, points to an int which provides the identity of the last mutation which occurred
	n_mutations, number of mutations which will occur after copying the state	
	*/
	
	int num_mut_to_copy = state[0][index_born]; // first row is the number of mutants
	if(num_mut_to_copy >= MAX_MUTNS - n_mutations){
		printf("Indiv %d has %d to copy but MAX_MUTNS - n_mutations = %d. Increase MAX_MUTNS. Exiting.\n", 
			index_born, num_mut_to_copy, MAX_MUTNS - n_mutations);
		exit(99);
	}
	
	if(*mut_counter > INT_MAX - n_mutations){
		printf("mut_counter overflow! Exiting.\n");
		exit(99);; // return error code
	}

	if(index_born == index_die){
		for (int i = 1; i < n_mutations+1; i++){
			*mut_counter += 1; // identity of the mutation
			state[num_mut_to_copy+i][index_die] = *mut_counter; // store identity of mutation			
		}
		state[0][index_die] += n_mutations; // increment the total number of mutations on the newborn
		return; // return ok
	}
	else{
		Copy_state(state, index_born, index_die);
		for (int i = 1; i < n_mutations+1; i++){
			*mut_counter += 1; // identity of the mutation
			state[num_mut_to_copy+i][index_die] = *mut_counter; // store identity of mutation			
		}
		state[0][index_die] += n_mutations; // increment the total number of mutations on the newborn
		return; // return ok
	}	
}

void Add_mutant_to_lists(int* mutant_identities, int* mutant_counts, int mutant_idty){
	/* Add a mutant to an ordered list of mutants and a corresponding list of counts

	Inputs
	----------
	mutant_identities, a zero-padded ordered array of unique integers corresponding to the identities of each mutation to have occurred to the state. 
	mutant_counts, an array of integers corresponding to the counts of every entry of mutant_identities
	mutant_idty, the identity of the mutant to be added to mutant_identities and mutant_counts
	*/

	int temp_id_old = -1;
	int temp_id_new = -1;
	int temp_count_old = -1;
	int temp_count_new = -1;
	int discovery_flag = 0;
	
	int i;
	for(i = 0; i < MAX_MUTNS; i++){
		if(mutant_identities[i] == 0 && discovery_flag == 0){
			mutant_identities[i] = mutant_idty; // got to end of unique mutants, this index was largest
			mutant_counts[i] += 1;
			break;
		} 
		else if(mutant_identities[i] == mutant_idty){
			mutant_counts[i] += 1;
			break; // already in the list
		}
		else if(mutant_identities[i] > mutant_idty && discovery_flag == 0){ // Not in list, need to shift lists to right by 1
			temp_id_old = mutant_identities[i];
			temp_count_old = mutant_counts[i];

			mutant_identities[i] = mutant_idty;
			mutant_counts[i] = 1; // A new mutant is discovered, set its count to 1
			discovery_flag = 1;			
			continue; // go to next mutant
		} 		
		
		// Move previously-discovered mutants to the right by 1 to maintain an ordered list
		if(mutant_identities[i] == 0 && discovery_flag == 1){
		 	mutant_identities[i] = temp_id_old; // got to the end, assign temp
		 	mutant_counts[i] = temp_count_old;
		 	break;
		}
		else if (discovery_flag == 1){	
		 	temp_id_new = mutant_identities[i];
		 	temp_count_new = mutant_counts[i];

		 	mutant_identities[i] = temp_id_old;
		 	mutant_counts[i] = temp_count_old;

		 	temp_id_old = temp_id_new;
		 	temp_count_old = temp_count_new;
		}

	}
}

void Get_idtys_counts_mutants(int **state, int *mutant_identities, int *mutant_counts){
	/* Populate empty arrays of ints with the identity and the total counts of every mutation in the state

	Inputs
	---------
	state, the state of the system (2D array)
	mutant_identities, a list of 0's of length MAX_MUTNS
	mutant_counts, a list of 0's of length MAX_MUTNS

	*/
	int n_muts_j, mutant_idty; 
	
	for(int j = 0; j < N_POP; j++){ // for every individual
		n_muts_j = state[0][j]; // get number of mutations 
		if(n_muts_j>0){ 					
			for(int k = 0; k < n_muts_j; k++){
				mutant_idty = state[k+1][j];
				Add_mutant_to_lists(mutant_identities, mutant_counts, mutant_idty);
			}
		}		 
	}

}

void Print_summary_string(int **state, FILE *fpStat, int *mutant_identities, int *mutant_counts){
	/* Write a comma-separated string to file which has the identity and the total counts of every mutation in the state

	Inputs
	-----------
	state, the state of the system (2D array)
	fpStat, a pointer to the file where output will be written which is free for appending
	mutant_identities, a zero-padded ordered array of unique integers corresponding to the identities of each mutation to have occurred to the state. 
	mutant_counts, an array of integers corresponding to the counts of every entry of mutant_identities
	*/
	
	for(int i = 0; i < MAX_MUTNS; i++){
		if(mutant_identities[i]==0){
			if(i==0) fprintf(fpStat, ",");
			break;
		}
		else{
			fprintf(fpStat, ",%d,%d",mutant_identities[i],mutant_counts[i]);
		}
	}


}

void Print_n_mut_n_het(int **state, FILE *fpStat, int *mutant_counts){
	/* Write a comma-separated string to file which has the total number of mutations, and the total number of unfixated mutations

	Inputs
	-----------
	state, the state of the system (2D array)
	fpStat, a pointer to the file where output will be written which is free for appending
	mutant_counts, a zero-padded array of integers corresponding to the counts of every entry of mutant_identities
	*/

	int n_mutations_total = 0;
	int n_mutations_het = 0;

	for (int i = 0; i < MAX_MUTNS; i++){
		if(mutant_counts[i]==0) break;
		if(mutant_counts[i]<N_POP){
			n_mutations_total += 1;	
			n_mutations_het += 1;
		}
		else if(mutant_counts[i]==N_POP){
			n_mutations_total += 1;
		} 
		else{
			printf("mutant_counts > N_POP!\n");
			exit(99);
		}
	}
		
	fprintf(fpStat, ",%d,%d",n_mutations_total,n_mutations_het);
	

}

void mean_std_mutants_per_indiv(int **state, double *mean_mutations_indiv, double *std_mutations_indiv){
	/* Compute the mean and standard deviation of the the number of distinct mutations across all individuals for a particular state
	
	Inputs
	-------
	state, the state of the system (2D array)
	mean_mutations_indiv, a pointer to a double which will be updated with the mean 
	std_mutations_indiv, a pointer to a double which will be updated with the sample standard deviation
	*/
	double n_muts[N_POP];
	for(int i = 0; i < N_POP; i++){
		n_muts[i] = (double)state[0][i];
	}
	*mean_mutations_indiv = gsl_stats_mean(n_muts,1,N_POP);
	*std_mutations_indiv = gsl_stats_sd_m(n_muts,1,N_POP, *mean_mutations_indiv);
}

int main(int argc, char *argv[]) {

/* set up GSL RNG */
gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);	
/* end of GSL setup */

clock_t begin = clock();


int initial_seed, fs_index;
double mutation_rate;
double replication_rate;	
if(argc == 5){
	initial_seed = atoi(argv[1]); 
	replication_rate = atof(argv[2]);
	mutation_rate = atof(argv[3]);
	fs_index = atoi(argv[4]);
	gsl_rng_set(r,initial_seed);
}
else{
	printf("argc = %d\n", argc);
	printf("Usage: %s  initial_seed mutation_rate fs_index\n", argv[0]);
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Sim Parameters //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double TFinal = 100*365.0; // total amount of time to run simulation (days)

int print_interval_years = 5;
printf("Printing every %d years\n", print_interval_years);
double DELTAT_STORE = print_interval_years*365.0; // Time lattice on which to store events
int NREPEATS = 1000; // Number of repeats of the stochastic simulation

double fs_arr[4] = {0.2,0.4,0.6,0.8};
if (fs_index < 0 || fs_index > 3){
	printf("0 <= fs <= 3, but fs = %d. Exiting.\n", fs_index); 
	return 0;
}

double fs_val = fs_arr[fs_index];
replication_rate = fs_val*replication_rate; // rescale time by fraction of singletons

printf("fs = %.1f\n", fs_val);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// PRINT OUTPUTS ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//this stores all statistics we look at, except for the distributions in the bins and the times at which first hit certain h.
char str[500];

sprintf(str, "./output_%d.txt", fs_index);
FILE *fpStat;
fpStat = fopen(str, "w");

fprintf(fpStat,"NREPEATS = %d\n", NREPEATS);
fprintf(fpStat,"replication_rate = %f\n", replication_rate);
fprintf(fpStat,"mutation_rate = %f\n", mutation_rate);
fprintf(fpStat,"TFinal = %f\n", TFinal);
fprintf(fpStat,"fs = %f\n", fs_val);
//fprintf(fpStat, "t,rep,mut_str\n");
fprintf(fpStat, "t,rep,n_mutants,n_het_mutants,mean_mutn_indiv,std_mutn_indiv\n");
fclose(fpStat); 


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Setup ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int** state;
/*
state[0][i] is the number of mutations which individual i possesses
state[1:][i] are the identities of the mutations which individual i possesses
*/

int* mutant_identities;
int* mutant_counts;

int mut_counter;
int num_mutations, index_born, index_die; 
double t, t_store, r1;
double totalrate = replication_rate*N_POP; // time to next event is a simple constant

double mean_mutations_indiv, std_mutations_indiv;

//int nTimeSteps = (int)(TFinal/DELTAT_STORE)+1;
int nStore, nPrint, printCounter;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Simulate ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

for (int i=0; i< NREPEATS; i++){
	if(i%100 == 0) printf("i = %d\n", i);
	
	// reset state
	state = malloc((MAX_MUTNS+1) * sizeof(int*));
	for (int i = 0; i < (MAX_MUTNS+1); i++)	state[i] = calloc(N_POP, sizeof(int));

	mut_counter = 0; // reset identity counter of mutations
	t = 0.0; t_store = 0.0; nStore = 1; printCounter = 0;

	while (t < TFinal){		
		r1 = RND; 
		t += -log(r1)/totalrate; // update time with exponential increment

		//store state on time lattice
		if (t>(t_store - 1E-12)){
			nPrint = floor((t - t_store)/(double)DELTAT_STORE); 
			for (int k=0;k<(nPrint+1);k++){	
				if (printCounter >= (int)round(((double)(TFinal+DELTAT_STORE)/DELTAT_STORE))){
				 printf("error1, %d,%.4f, %d, %d\n", printCounter,t_store, k, nPrint); //to debug.
				 return 0;
				}

				mutant_identities = calloc(MAX_MUTNS, sizeof(int*));
				mutant_counts = calloc(MAX_MUTNS, sizeof(int*));

				Get_idtys_counts_mutants(state, mutant_identities, mutant_counts);

				fpStat = fopen(str, "a");
				fprintf(fpStat, "%f,%d", printCounter*DELTAT_STORE, i);

				//Print_summary_string(state, fpStat, mutant_identities, mutant_counts);
				Print_n_mut_n_het(state, fpStat, mutant_counts);				
				
				mean_std_mutants_per_indiv(state, &mean_mutations_indiv, &std_mutations_indiv);
				fprintf(fpStat, ",%f,%f", mean_mutations_indiv, std_mutations_indiv);				

				fprintf(fpStat, "\n");
				fclose(fpStat);

				free(mutant_identities);
				free(mutant_counts);

				printCounter ++; 
				t_store = (nStore*DELTAT_STORE);
				nStore += 1;

				if (t_store < TFinal + DELTAT_STORE + (double)DELTAT_STORE/10.0 && t_store > TFinal+DELTAT_STORE - (double)DELTAT_STORE/10.0) break;				

			}
		}

		// update state
		num_mutations = gsl_ran_binomial(r, mutation_rate, LEN_GENOME); // compute number of mutations in this replication event		
		index_born = (int)round(RND*(N_POP-1));
		index_die = (int)round(RND*(N_POP-1));

		if(num_mutations == 0) Moran_update_state(state,index_born,index_die);
		else{Moran_update_state_mutate(state, index_born, index_die, &mut_counter, num_mutations);}

	}

	// Clean up
	for (int i = 0; i < MAX_MUTNS; i++){
		free(state[i]);
	}
	free(state);
	
}

clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Total time taken = %.3f\n", time_spent);

return 0;
}