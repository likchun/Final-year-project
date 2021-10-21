/* main.c */

// Compile and run this file to start simulation for a network of N effective neuron nodes.
// Variables and settings are defined in HEADER FILE 'param.h'.


/* includes & defines */

// Specific variables & settings are included in the HEADER FILE (param.h)
#define HEADER_FILE "param.h"
#include HEADER_FILE

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
// #include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_DEPRECATE
#define PI 3.14159265358979323846


/* gloabl variables declarations */

float potential[N] = {0}, recovery[N] = {0};
float coupling[N][N] = {0}, spike_timestamp[N][SpikeBuffer] = {0}, **temporary;
int node_type[N] = {0}, spike_count[N] = {0}, exc_node_num[N] = {0}, inh_node_num[N] = {0};
int exc_node_idx[N][N] = {0}, inh_node_idx[N][N] = {0};


/* function prototypes */

void SuppressInhibitory();
void EnhanceInhibitory();
void SuppressExcitatory();
void EnhanceExcitatory();

int random_integer(int upper_limit, int lower_limit);
float gaussian_random();
float update_potential(float recovery_variable, float membrane_potential, float current);
float update_recovery(int node_idx, float recovery_variable, float membrane_potential);
float conductance_exc(int node_idx, float time);
float conductance_inh(int node_idx, float time);
void ErrorOccur(int error_code, char const * error_message, int error_number);
bool IsExactlyOneBooleanTrue(bool, bool, bool, bool);
void SuppressOrEnhanceNetwork();
void SetInitialValues();
void PopulateCouplingMatrix(FILE* fp_read);
void FastReferenceToNodeType();
void ClassifyNodeType();
void AllocateSpaceForTemporary();
void Export_INFO(FILE* fp_info);
void Export_CNFG(FILE* fp_cnfg);


int main(int argc, char *argv[])
{

    /* local variables */

    int idx = 0, jdx = 0, timestep = 0;
    float t = 0, sqrt_dt = sqrt(dt), current = 0, noise = 0;
    float potential_temp = 0;

    if (history_truncation == false) {
        trunc_time_exc = tn;
        trunc_time_inh = tn;
    }


    /* input & outputs */

    FILE* check_fout_exist_data = NULL;
    FILE* check_fout_exist_spik = NULL;
    char user_input[20];
    // Check if previous data file exists, prevent unintentional overwriting
    if ((check_fout_exist_data = fopen(output_data, "r")) || (check_fout_exist_spik = fopen(output_spik, "r"))) {
        printf("\nPrevious data file exists. If you choose to OVERWRITE it, input \"Overwrite\" (case sensitive), otherwise quit: ");
        scanf("%s", user_input);
        if (strcmp(user_input, "Overwrite") != 0) {
            fclose(check_fout_exist_data);
            ErrorOccur(6, NULL, 0);
        } else {
            fclose(check_fout_exist_data);
        }
    }

    FILE* fp_read = NULL;
	if (!(fp_read = fopen(input_data, "r"))) {
		ErrorOccur(0, input_data, 0);
	}


    /* function calls */

    // Seed for random function
    srand(seed_for_random);

    SetInitialValues();

    PopulateCouplingMatrix(fp_read);

    FastReferenceToNodeType();

    ClassifyNodeType();
	
	SuppressOrEnhanceNetwork();

    fclose(fp_read);

    // Initialization ends here
    printf("\n[notice] Initialization completed.\n\n");

    time_t start_date_time = time(NULL);
    struct tm start_dt = *localtime(&start_date_time);

    printf("[notice] Starting time: %02d:%02d:%02d.\n\n", start_dt.tm_hour, start_dt.tm_min, start_dt.tm_sec);
    printf("[info] Network: %s\n[info] T: %0.0f, dt: %f\n[info] Random number seed: %0.0f\n", input_data, tn, dt, seed_for_random);
	if (suppress_inhibitory) { printf("\n[info] Suppress Inhibitory: YES\n[info] Suppression level: %f\n[info] Inhibitory S.D.: %f", alter_inh_k, alter_inh_sd); }
    if (enhance_inhibitory)  { printf("\n[info] Enhance Inhibitory: YES\n[info] Enhancement level: %f\n[info] Inhibitory S.D.: %f", alter_inh_k, alter_inh_sd); }
    if (suppress_excitatory) { printf("\n[info] Suppress Excitatory: YES\n[info] Suppression level: %f\n[info] Excitatory S.D.: %f", alter_inh_k, alter_inh_sd); }
    if (enhance_excitatory)  { printf("\n[info] Enhance Excitatory: YES\n[info] Enhancement level: %f\n[info] Excitatory S.D.: %f", alter_inh_k, alter_inh_sd); }


    // Start timer
    clock_t tic = clock();

    // Export the initial data points for every nodes
    if (output_potential_enabled == true) {
        AllocateSpaceForTemporary();
        for (idx = 0; idx < N; idx++) {
            temporary[idx][0] = potential[idx];
        }
    }
    
    /* main loop */

    // Check 'output_potential_enabled' outside of the loop so as to avoid unnecessary branching within the loop
    if (output_potential_enabled == true) {
        while (t < tn) {
            t += dt;
            timestep = (int) (t / dt);
            for (idx = 0; idx < N; idx++) {
                
                // Generating noise for calculation in v(t), white gaussian noise is used
                noise = sigma * gaussian_random();

                // Main calculations
                potential_temp = potential[idx];

                current = beta*(conductance_exc(idx, t)*(thresh_v_exc-potential_temp) - conductance_inh(idx, t)*(potential_temp-thresh_v_inh));

                potential[idx] += (update_potential(potential_temp, recovery[idx], current)) * dt + noise * sqrt_dt;
                recovery[idx] += update_recovery(idx, potential_temp, recovery[idx]) * dt;
                
                // Temporarily store potential into temporary[][]
                temporary[idx][timestep] = potential[idx];
                
                // Reset the membrane potential if it is greater than 30 mV
                if (potential[idx] >= 30) {
                    if (node_type[idx] == -1) {
                        potential[idx] = c_inh;
                        recovery[idx] += d_inh;
                    } else {
                        potential[idx] = c_exc;
                        recovery[idx] += d_exc;
                    }
                    spike_timestamp[idx][spike_count[idx]] = t;
                    spike_count[idx]++;
                }
            }
        }

        // Export data for OUT_POTV.txt
        FILE* fp_out = fopen(output_data, "w");
        for (idx = 0; idx < N; idx++) {
            fprintf(fp_out, "%0.4f", temporary[idx][0]);
            for (jdx = 1; jdx < (int)(tn/dt); jdx++) {
                fprintf(fp_out, "\t%0.4f", temporary[idx][jdx]);
            }
            fprintf(fp_out, "\n");
        }
        fclose(fp_out);

    } else {
        while (t < tn)
        {
            t += dt;
            for (idx = 0; idx < N; idx++) {

                noise = sigma * gaussian_random();

                potential_temp = potential[idx];

                current = beta*(conductance_exc(idx, t)*(thresh_v_exc-potential_temp) - conductance_inh(idx, t)*(potential_temp-thresh_v_inh));

                potential[idx] += (update_potential(potential_temp, recovery[idx], current)) * dt + noise * sqrt_dt;
                recovery[idx] += update_recovery(idx, potential_temp, recovery[idx]) * dt;
                
                // No time series data output here
                
                if (potential[idx] >= 30) {
                    if (node_type[idx] == -1) {
                        potential[idx] = c_inh;
                        recovery[idx] += d_inh;
                    } else {
                        potential[idx] = c_exc;
                        recovery[idx] += d_exc;
                    }
                    spike_timestamp[idx][spike_count[idx]] = t;
                    spike_count[idx]++;
                }
            }
        }
    }

    // Stop timer
    clock_t toc = clock();
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;

    time_t end_date_time = time(NULL);
    struct tm end_dt = *localtime(&end_date_time);


    /* outputs */
    
    // Export data for OUT_SPIK.txt
    // Format: node number, number of spikes, timestamps of spikes
    if (output_spikeinfo_enabled == true) {
        FILE* fp_spik = fopen(output_spik, "w");
        for (idx = 0; idx < N; idx++) {
            fprintf(fp_spik, "%d\t%d", idx+1, spike_count[idx]);
            for (jdx = 0; jdx < spike_count[idx]; jdx++) {
                fprintf(fp_spik, "\t%0.2f", spike_timestamp[idx][jdx]);
            }
            fprintf(fp_spik, "\n");
        }
        fclose(fp_spik);
    }
    
    // Export data for OUT_INFO.txt
    FILE* fp_info = fopen(output_info, "w");
    fprintf(fp_info, "Completed on %d-%02d-%02d %02d:%02d:%02d\n", end_dt.tm_year + 1900, end_dt.tm_mon + 1, end_dt.tm_mday, end_dt.tm_hour, end_dt.tm_min, end_dt.tm_sec);
    fprintf(fp_info, "Execution time: %f seconds\n\n", time_spent);
    Export_INFO(fp_info);
    fclose(fp_info);

    // Export data for INI_CNFG.txt
    FILE* fp_cnfg = fopen(output_ini_cnfg, "w");
    Export_CNFG(fp_cnfg);
    fclose(fp_cnfg);

    printf("\n[notice] Completed. Time taken: %f seconds.\n\n", time_spent);

    if (output_potential_enabled == true) {
        free(temporary);
    }

    return 0;
}


/* functions (suppression / enhancement) */

void SuppressInhibitory() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (coupling[i][j] < 0) {
                coupling[i][j] += alter_inh_k * alter_inh_sd;
                if (coupling[i][j] > 0) {
                    coupling[i][j] = 0;
                }
            }
        }
    }
}

void EnhanceInhibitory() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (coupling[i][j] < 0) {
                coupling[i][j] -= alter_inh_k * alter_inh_sd;
            }
        }
    }
}

void SuppressExcitatory() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
			if (coupling[i][j] > 0) {
				coupling[i][j] -=  alter_exc_k * alter_exc_sd;
				if (coupling[i][j] < 0) {
					coupling[i][j] = 0;
				}
			}
        }
    }
}

void EnhanceExcitatory() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
			if (coupling[i][j] > 0) {
				coupling[i][j] +=  alter_exc_k * alter_exc_sd;
			}
        }
    }
}


/* functions (computational) */

// Random number generator, draw integers from an interval uniformly
// includes lower bound, but exclude upper bound
int random_integer(int upper, int lower) {
    return (rand() % (upper - lower + 1)) + lower;
}

// Random Gaussian number generator, using Box-Muller method
float gaussian_random() {
    // '(float)rand() / (float)RAND_MAX' generate a random number within [0,1), i.e. includes 0,
    // which will produce erroneous result in Box-Muller method; x, y should be random number within (0,1)
    // Redraw if '(float)rand() / (float)RAND_MAX' gives 0
    float x = (float)rand() / (float)RAND_MAX;
    while (x == 0) { x = (float)rand() / (float)RAND_MAX; }
    float y = (float)rand() / (float)RAND_MAX;
    while (y == 0) { y = (float)rand() / (float)RAND_MAX; }
    return sqrt(-2 * log(x)) * cos(2 * PI * y);
}

float update_potential(float potential, float recovery, float current) {
	return ((0.04 * potential * potential) + (5 * potential) + 140 - recovery + current);
}

float update_recovery(int idx, float potential, float recovery) {
    if (node_type[idx] == -1) {
        return (a_inh * (b_inh * potential - recovery));
    } else {
        return (a_exc * (b_exc * potential - recovery));
    }
}

float conductance_exc(int idx, float time) {
    int j = 0, k = 0;
    float spike_summ = 0, node_summ = 0, time_diff = 0;
    for (j = 0; j < exc_node_num[idx]; j++) {
        for (k = 0; k < spike_count[exc_node_idx[idx][j]]; k++) {
            time_diff = (time - spike_timestamp[exc_node_idx[idx][j]][k]);
            if (time_diff < trunc_time_exc) {
                spike_summ += exp(-1 * time_diff / tau_exc);
            }
        }
        node_summ += (coupling[idx][exc_node_idx[idx][j]] * spike_summ);
        spike_summ = 0;
    }
    return node_summ;
}

float conductance_inh(int idx, float time) {
    int j = 0, k = 0;
    float spike_summ = 0, node_summ = 0, time_diff = 0;
    for (j = 0; j < inh_node_num[idx]; j++) {
        for (k = 0; k < spike_count[inh_node_idx[idx][j]]; k++) {
            time_diff = (time - spike_timestamp[inh_node_idx[idx][j]][k]);
            if (time_diff < trunc_time_inh) {
                spike_summ += exp(-1 * time_diff / tau_inh);
            }
        }
        node_summ += (-1 * coupling[idx][inh_node_idx[idx][j]] * spike_summ);
        spike_summ = 0;
    }
    return node_summ;
}


/* functions (utility) */

// Handle runtime errors; exit and terminate all process
void ErrorOccur(int error_code, char const * error_message, int error_number) {
    printf("\n<warning> Fatal runtime error encountered.\n");
    switch (error_code) {
        case 0:
            printf("<warning> File %s cannot be opened.\n", error_message);
            break;
        case 1:
            printf("<warning> Memory allocation error in %s.\n", error_message);
            break;
        case 2:
            printf("<warning> Funtion read error: 'fgets'; 'token=NULL' at row %d.\n", error_number);
            break;
        case 3:
            printf("<warning> Function conversion error: 'strtod(%s)' at row %d.\n          Possible solution: increase value of InputBuffer (in HEADER FILE).\n", error_message, error_number);
            break;
        case 4:
            printf("<warning> Classification error: the outgoing links of each node, when exist, are assumed be all excitatory or all inhibitory, inconsistency is detected.\n");
            break;
        case 5:
            printf("<warning> Setting error in: %s.\n Current setting is %d. Make changes in the HEADER FILE 'param.h'.\n", error_message, error_number);
            break;
        case 6:
            printf("[notice] Terminated by user.\n");
            break;
		case 7:
			printf("<warning> You can only suppress or enhance one type of node.\n          In 'param.h', '/* synaptic weights suppression, enhancement */' section, you can only enable one of them.\n");
			break;
		default:
            printf("<warning> Unknown error.\n");
            break;
    }
    printf("\n");
    exit(1);
}

bool IsExactlyOneBooleanTrue(bool b1, bool b2, bool b3, bool b4) {
	bool areAnyTrue = false;
	bool areTwoTrue = false;
	bool bool_arr[4] = {b1, b2, b3, b4};
	for (int i = 0; (!areTwoTrue) && (i < 4); i++) {
		areTwoTrue = (areAnyTrue && bool_arr[i]);
		areAnyTrue |= bool_arr[i];
	}
	return ((areAnyTrue) && (!areTwoTrue));
}

void SuppressOrEnhanceNetwork() {
	if (suppress_inhibitory || enhance_inhibitory || suppress_excitatory || enhance_excitatory) {
		if (!IsExactlyOneBooleanTrue(suppress_inhibitory,enhance_inhibitory,suppress_excitatory,enhance_excitatory)) {
			ErrorOccur(7, NULL, 0);
		}
	}

    if (suppress_inhibitory) {
        SuppressInhibitory();
    } else if (enhance_inhibitory) {
        EnhanceInhibitory();
	} else if (suppress_excitatory) {
        SuppressExcitatory();
    } else if (enhance_excitatory) {
        EnhanceExcitatory();
    }
}

// Setting initial values v0 & u0 for each node, using 3 methos, refer to HEADER FILE
void SetInitialValues() {
    int i = 0;
    for (i = 0; i < N; i++) {
        if (init_val_rnd_method == 0) {
            potential[i] = v0;
            recovery[i] = u0;
        } else if (init_val_rnd_method == 1) {
            potential[i] = v0 + gaussian_random() * v0_rnd_sd;
            recovery[i] = u0 + gaussian_random() * u0_rnd_sd;
        } else if (init_val_rnd_method == 2) {
            potential[i] = random_integer(v0_rnd_max, v0_rnd_min);
            recovery[i] = random_integer(u0_rnd_max, u0_rnd_min);
        } else {
            ErrorOccur(5, "randomize initial values generating method; 0: gaussian, 1: draw from an interval evenly", init_val_rnd_method);
        }   
    }
}

void PopulateCouplingMatrix(FILE* fp_read) {
    int i = 0, j = 0;
    int current_row = 0;
    double delimit_value;
    char* token, line[InputBuffer], * endptr;

	for (i = 0; i < N; i++) {
		++current_row;
		token = fgets(line, InputBuffer, fp_read);
		if (token == NULL) {
			ErrorOccur(2, NULL, current_row);
		}
		token = strtok(line, Delimiter);
		for (j = 0; j < N; j++) {
			delimit_value = strtod(token, &endptr);
			if (token == endptr) {
				ErrorOccur(3, token, current_row);
			}
			coupling[i][j] = delimit_value;
			token = strtok(NULL, Delimiter);
		}
	}
}

// In-memory reference to index of exc and inh nodes; prepare for optimized loops for faster experience
void FastReferenceToNodeType() {
    int i = 0, j = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (coupling[i][j] > 0) {
                exc_node_idx[i][exc_node_num[i]] = j;
                exc_node_num[i]++;
            } else if (coupling[i][j] < 0) {
                inh_node_idx[i][inh_node_num[i]] = j;
                inh_node_num[i]++;
            }
        }
    }
}

// Classify each node into appropriate type: 1: EXCI, -1: INHI, 0: UNCL
void ClassifyNodeType() {
    int i = 0, j = 0;;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			// Classified as E-nodes
			if (coupling[i][j] > 0) {
				if (node_type[j] == -1) {
                    printf("\nAt coupling matrix [%d, %d]", i+1, j+1);
					ErrorOccur(4, NULL, 0);
                }
				node_type[j] = 1;
			}
			// Classified as I-nodes
			else if (coupling[i][j] < 0) {
				if (node_type[j] == 1) {
					printf("\nAt coupling matrix [%d, %d]", i+1, j+1);
					ErrorOccur(4, NULL, 0);
                }
                node_type[j] = -1;
			}
            // Otherwise, keep it zero
		}
	}
}

void AllocateSpaceForTemporary() {
    int i = 0;
    temporary = (float **) calloc(N, sizeof(float*));
    if (!temporary) {
        ErrorOccur(1, "temporary storage matrix temporary[]", 0);
    }
    for (i = 0; i < N; i++) {
        temporary[i] = (float*) calloc((int)(tn/dt), sizeof(float));
        if (!temporary[i]) {
            ErrorOccur(1, "temporary storage matrix temporary[][]", 0);
        }
    }
}

void Export_INFO(FILE* fp_info) {
    fprintf(fp_info, ">> Simulation Variables <<\n\n");
	fprintf(fp_info, "Input file (coupling matrix): %s\n", input_data);
    fprintf(fp_info, "Number of nodes used, N: %d\n", N);
    fprintf(fp_info, "Time step, dt: %f\n", dt);
    fprintf(fp_info, "Total simulation time, tn: %f\n", tn);
    fprintf(fp_info, "Strength of noise, sigma: %f\n", sigma);
    fprintf(fp_info, "Type of noise: Gaussian white\n");
	fprintf(fp_info, "\n>> Suppress INH, Enhance EXC <<\n\n");
	fprintf(fp_info, "Suppress inhibitory: %s\n", suppress_inhibitory ? "true" : "false");
	fprintf(fp_info, "Enhance inhibitory: %s\n", enhance_inhibitory ? "true" : "false");
	fprintf(fp_info, "  suppression / enhancement level, k: %f\n", alter_inh_k);
	fprintf(fp_info, "  suppression / enhancement sd: %f\n", alter_inh_sd);
	fprintf(fp_info, "Suppress excitatory: %s\n", suppress_excitatory ? "true" : "false");
	fprintf(fp_info, "Enhance excitatory: %s\n", enhance_excitatory ? "true" : "false");
	fprintf(fp_info, "  suppression / enhancement level, k: %f\n", alter_exc_k);
	fprintf(fp_info, "  suppression / enhancement sd: %f\n", alter_exc_sd);
    fprintf(fp_info, "\n>> Parameters for the Izhikevich's Neuron Spiking Model <<\n\n");
    fprintf(fp_info, "Method to generate randomize initial values: ");
    if (init_val_rnd_method == 0) {
        fprintf(fp_info, "none, fixed values\n");
        fprintf(fp_info, "  initial values, v0: %f, u0: %f\n", v0, u0);
    } else if (init_val_rnd_method == 1) {
        fprintf(fp_info, "Gaussian, random values follow normal distribution\n");
        fprintf(fp_info, "  mean of v0: %f, u0: %f\n", v0, u0);
        fprintf(fp_info, "  s.d. of v0: %f, u0: %f\n", v0_rnd_sd, u0_rnd_sd);
    } else if (init_val_rnd_method == 2) {
        fprintf(fp_info, "draw from intervals uniformly\n");
        fprintf(fp_info, "  v0 in range (%d, %d), u0 in range (%d, %d)\n", v0_rnd_min, v0_rnd_max, u0_rnd_min, u0_rnd_max);\
    }
    fprintf(fp_info, "\nParameters for excitatory node:\n  a: %f, b: %f, c: %f, d: %f\n", a_exc, b_exc, c_exc, d_exc);
    fprintf(fp_info, "Parameters for inhibitory node:\n  a: %f, b: %f, c: %f, d: %f\n", a_inh, b_inh, c_inh, d_inh);
    fprintf(fp_info, "\nOther parameters for model:\n  threshold potential (excitatory): %f\n  threshold potential (inhibitory): %f\n", thresh_v_exc, thresh_v_inh);
    fprintf(fp_info, "  tau_exc: %f, tau_inh: %f\n", tau_exc, tau_inh);
    fprintf(fp_info, "  conductance G amplifying constant, beta: %f\n", beta);
	fprintf(fp_info, "\n>> Simulation Settings <<\n\n");
    fprintf(fp_info, "Seed for generating random numbers: %f\n", seed_for_random);
    fprintf(fp_info, "Spiking history truncation: %s\n", history_truncation ? "true" : "false");
    if (history_truncation == true) {
    fprintf(fp_info, "  truncation time (excitatory): %f\n", trunc_time_exc);
    fprintf(fp_info, "  truncation time (inhibitory): %f\n", trunc_time_inh);
    }
    fprintf(fp_info, "\n>> Other Settings <<\n\n");
    fprintf(fp_info, "Spike buffer: %d\n", SpikeBuffer);
    fprintf(fp_info, "Input buffer: %d\n", InputBuffer);
    if ("\t" == Delimiter) {
        fprintf(fp_info, "Delimiter: \"\\t\" (tab)\n");
    } else {
        fprintf(fp_info, "Delimiter: \"%s\"\n", Delimiter);
    }
    fprintf(fp_info, "\n\nNote: all units are in ms or mV");
}

void Export_CNFG(FILE *fp_cnfg) {
    // line 1: simualtion parameters
    fprintf(fp_cnfg, "%d\t%f\t%f\t%f\n", N, dt, tn, sigma);
	// line 2: synaptic weights suppression, enhancement
	fprintf(fp_cnfg, "%d\t%d\t%f\t%f\t%d\t%d\t%f\t%f\n", suppress_inhibitory, enhance_inhibitory, alter_inh_k, alter_inh_sd, suppress_excitatory, enhance_excitatory, alter_exc_k, alter_exc_sd);
    // line 3: model parameters
    fprintf(fp_cnfg, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", a_exc, b_exc, c_exc, d_exc, a_inh, b_inh, c_inh, d_inh, thresh_v_exc, thresh_v_inh, tau_exc, tau_inh, beta);
    fprintf(fp_cnfg, "%f\t%f\t%d\t%f\t%f\t%d\t%d\t%d\t%d\n", v0, u0, init_val_rnd_method, v0_rnd_sd, u0_rnd_sd, v0_rnd_max, v0_rnd_min, u0_rnd_max, u0_rnd_min);
    // line 4: simulation settings
    fprintf(fp_cnfg, "%f\t%d\t%f\t%f\n", seed_for_random, history_truncation ? 1 : 0, trunc_time_exc, trunc_time_inh);
	// line 5: output settings
    fprintf(fp_cnfg, "%d\t%d\n", output_spikeinfo_enabled ? 1 : 0, output_potential_enabled ? 1 : 0);
    // line 6: output paths
    fprintf(fp_cnfg, "%s\t%s\t%s\t%s\t%s\n", input_data, output_data, output_info, output_spik, output_ini_cnfg);
    // line 7: buffers and delimiter
    fprintf(fp_cnfg, "%d\t%d\t%s", SpikeBuffer, InputBuffer, Delimiter);
}
