# two-alleles-branching-regime
Interference between a sweeping allele that propagates as a wave and an allele who is in the branching regime 
# for the c++ file:
## input:
trials-number of trials to run the simulation, chi - the distance of the focal allele from the sweeping allele,
wrap - boolean, if wrap=1 then the system is on a ring else it is has 'walls' on both sides, tappear - the time the focal allele appears at; 
prior to that time the sweeping allele is the only mutation, N - number of agents in a deme, L - size of the system (measured in demes),
s2 - the fitness of the focal allele, r - recombination rate, outfile - name of the output file

## output:
a file where each row is total number of wildtype, sweeping allele and focal allele.  Recorded every 20 generations. 
When simulation is done writes if the focal allele fixed or went extinct. 

## comments:
Needs gsl random libary to run. The initial distance between the alleles is defined as the distnace between the focal allele and the 
point on the wave where the occupancy is half the population in the deme.   

# Jupyter notebook:
Read files generated from forAWS.cpp. It assumes that the file name has a structure: string$i$x.txt. Where $i is an index corresponding to the recombination rate, $x is the initial seperation between the alleles measured in lattice sites.
