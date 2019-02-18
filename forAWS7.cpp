/*
*  Find how much one sweep shifts the ancestry at another locus in space
* This code is to find the survival probability in the branching regime
*
* Modification of Mike McLaren's meta_hap_stepping_stone.c
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include<time.h>
#include<vector>
#include<map>

// Turn off/on logging of genotype counts in alldata.dat file with 0/1
#define LOGGING 0
// Haplotype indices
#define ab 0
#define Ab 1
#define aB 2
#define AB 3
#define Lsweep 200
// Number of haplotypes (4) and number of genotypes (10)
#define HTYPES 4

// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(unsigned int **n, double**x, double**xmig, double w[4], double r,
	double mig, unsigned int N, unsigned int L, int wrap, int *midnose1) {
	double D; // coefficient of linkage disequilibrium
	int i, j, k, flag = 1,flag2=1;
	// Recombination
	for (i = 0; i < L; i++) {
		// Calculate haplotype frequencies
		for (j = 0; j < 4; j++) {
			x[i][j] = ((float)n[i][j]) / N;
		}

		// Recombination
		D = x[i][0] * x[i][3] - x[i][1] * x[i][2];
		x[i][0] -= r * D;
		x[i][1] += r * D;
		x[i][2] += r * D;
		x[i][3] -= r * D;
	}

	// Migration
	if (wrap) {
		// Wraparound (no boundaries), so that first and last demes are neighbors
		for (j = 0; j < 4; j++) {
			// leftmost deme:
			xmig[0][j] = (1 - mig) * x[0][j] + 0.5 * mig * (x[L - 1][j] + x[1][j]);
			//rightmost deme:
			xmig[L - 1][j] = (1 - mig) * x[L - 1][j] + 0.5 * mig * (x[L - 2][j] + x[0][j]);
		}
	}
	else {
		// 0 and L are really boundaries
		for (j = 0; j < 4; j++) {
			// leftmost deme:
			xmig[0][j] = (1 - 0.5 * mig) * x[0][j] + 0.5 * mig * x[1][j];
			//rightmost deme:
			xmig[L - 1][j] = (1 - 0.5*mig) * x[L - 1][j] + 0.5 * mig * x[L - 2][j];
		}
	}
	// All interior demes
	for (i = 1; i < L - 1; i++) {
		for (j = 0; j < 4; j++) {
			xmig[i][j] = (1 - mig) * x[i][j] + 0.5 * mig * (x[i - 1][j] + x[i + 1][j]);
		}
	}

	// Selection and sampling within demes

	for (i = 0; i < L; i++) {
		// Selection (within deme)
		for (j = 0; j < 4; j++) {
			if (xmig[i][j] == 1.0) { // a genotype is fixed in the deme, so no need to do sampling
				for (k = 0; k < 4; k++) n[i][k] = 0;
				n[i][j] = N;
				break;
			}
			x[i][j] = xmig[i][j] * w[j];
		}
		// Normalization of x[i] done in multinomial function, so dividing by
		// average fitness is not necessary
		gsl_ran_multinomial(R, 4, N, x[i], n[i]);
  	
		if (n[i][1]>=N/2.)
                {
			*midnose1 = i;
		}
	}
}

/*function roundint to round a double to an integer*/
int roundint(double x) {
	int result;

	result = floor(x);
	if (x - result >= 0.5) result++;

	return result;
} /*roundint*/


void mutate(unsigned int **n, int *focal, int midnose1, int chi,int N)
{
   *focal = midnose1+chi;
    if(gsl_ran_flat(R,0,1)<double(n[*focal][0])/N)
    {
       n[*focal][0]-=1;
       n[*focal][2]+=1;
    }
    else
    {
       n[*focal][1]-=1;
       n[*focal][3]+=1;
    }
}

int main(int argc, char *argv[]) {
	int trails;
	long SEED;
	int wrap; //whether the range wraps around
	double r, mig, s1, s2; // recomb rate, migration rate, selection
	double w[4]; // fitnesses of the four genotypes
	unsigned int N, t, tfinal; // deme size, number of demes, deme to paint with neutral allele, time, max number of generations
	unsigned int ntot[4]; // numbers of the four genotypes in the total population
	FILE *datafile;
	FILE *paramfile;
	char *outfile;
	char *outfile1 = (char*)malloc(1000);
	int i, j, k;
	int midnose1;
        int tappear,focal;
        unsigned int L;
        int chi;
	j = 1;
	trails = atof(argv[j++]);
	chi = atof(argv[j++]);
	wrap = roundint(atof(argv[j++]));
	mig = atof(argv[j++]);
        tappear = atof(argv[j++]);
	N = (unsigned int)roundint(atof(argv[j++]));
        L = (unsigned int)roundint(atof(argv[j++]));
	s1 = 0.05;
	s2 = atof(argv[j++]);
	r = atof(argv[j++]);
	outfile = argv[j++];
	w[0] = 1;
	w[1] = 1 + s1;
	w[2] = 1 + s2;
	w[3] = 1 + s1 + s2;


	tfinal = 500000;// (unsigned int)roundint(atof(argv[j++]));
	SEED = (-1)*time(NULL);// atof(argv[j++]);
	strcpy(outfile1, outfile);

	datafile = fopen(strcat(outfile, ".txt"), "w");
	strcpy(outfile, outfile1);

	paramfile = fopen(strcat(outfile, "_params.txt"), "w");
	strcpy(outfile, outfile1);
	fprintf(paramfile, " trails = %d\n chi = %d\n bc = %d\n mig = %f\n N = %d\n L=%d\n s2 = %f\n r = %lf\n", trails,chi,wrap, mig, N,L, s2, r);
	fclose(paramfile);


	gsl_rng_env_setup();
	T = gsl_rng_default;
	R = gsl_rng_alloc(T);
	gsl_rng_set(R, SEED);
 
	// Initialize population:
	unsigned int **n = (unsigned int**)malloc(L * sizeof(unsigned int*)); // numbers of the four genotypes in each deme
	double **x = (double**)malloc(L * sizeof(double*));
	double **xmig = (double**)malloc(L * sizeof(double*));
	for (i = 0; i < L; i++)
	{
		n[i] = (unsigned int*)malloc(4 * sizeof(unsigned int));
		x[i] = (double*)malloc(4 * sizeof(double));
		xmig[i] = (double*)malloc(4 * sizeof(double));
	}
        
	for (k = 0; k < trails; k++)
	{
                focal=0;
		// Left-most deme fixed for sweeping allele:
		for (i = 0; i < 10; i++)
		{
			n[i][0] = 0;
			n[i][1] = N;
			n[i][2] = 0;
			n[i][3] = 0;
		}
		for (i = 10;i < L; i++)
		{
			n[i][0] = N;
			n[i][1] = 0;
			n[i][2] = 0;
			n[i][3] = 0;
		}
                midnose1 = 9;
		// Run until one of the alleles fixes (or goes extinct):
		for (t = 0; t < tfinal; t++) {
	
   			// Evolve for one generation
			next_gen(n, x, xmig, w, r, mig, N,L, wrap, &midnose1);
                        //mutate a random member in deme chi from the midpoint
                        if (t==tappear)
                            mutate(n,&focal, midnose1,chi,N);	
			// Add up genotype counts within each deme
			for (j = 0; j < 4; j++) {
				ntot[j] = 0;
			}
			for (i = 0; i < L; i++) {
				for (j = 0; j < 4; j++) {
						ntot[j] += n[i][j];
				}
			}
                        if (t>tappear)
                        {
                           if (!((t-tappear)%20))
                           {
                               fprintf(datafile,"%d %d %d\n",ntot[0],ntot[1],ntot[2]);
                           }
                        }
			// Stop if one of the alleles is fixed or extinct:
			if (ntot[0] == N*L || ntot[1] == N*L || (ntot[2] +ntot[3])== N*L )
                        {
                                fprintf(datafile,"%d\n",(ntot[2] +ntot[3])== N*L);
				break;
                        }
		}

	}
	fclose(datafile);
        for (i = 0; i < L; i++)
	{
            free(n[i]);
            free(x[i]);
            free(xmig[i]);
	}
        free(n);
        free(x);
        free(xmig);
        free(outfile1);
        free(datafile);
        free(paramfile);
        printf("done");
	return 0;
}
