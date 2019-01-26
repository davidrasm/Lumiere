package beast.evolution.fitnessmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.OptionalDouble;
import java.util.stream.DoubleStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.SmallNumber;

/**
 * @author David Rasmussen
 */
@Description("Multi-site fitness model" +
             "used for flu H3N2 analysis tracking genotypes.")
public class FluGenotypeFitnessModel extends MultiSiteFitnessModel {
	
	// Input for site-wise fitness effects
	public Input<RealParameter> fitEffectsInput =
			new Input<>("mutFitEffects", "Site specific fit effects");
	
	// Mutation rates
	public Input<RealParameter> mutationMatrix =
			new Input<>("mutationMatrix", "Flattened mutation matrix, can be asymmetric, diagnonal entries omitted");
	
	protected Double[] fitEffects;
	
	@Override
	public void update() {

		fitEffects = fitEffectsInput.get().getValues();
		double mutRate = mutationMatrix.get().getValues()[0];
		//double backMutRate = mutationMatrix.get().getValues()[1];
		
        // Get genotype fitness values
        fitSpaceValues = new Double[genotypes.length];
        for (int i=0; i<genotypes.length; i++) {
        	fitSpaceValues[i] = fitEffects[i];
        }
        
		// Set up fitness space
		fitSpaceStates = genotypes.length;
		
		// Set up fitness space gammas
		fitSpaceGammas = new Double[fitSpaceStates*(fitSpaceStates-1)];
        
        // Compute transition rates between genotypes
        for (int u=0; u<genotypes.length; u++) {
        	for (int v=0; v<genotypes.length; v++) {
        		if (u != v) {
        			int index = u * (fitSpaceStates - 1) + ((u > v) ? v : v-1);
        			int[] genotype_i = genotypes[u];
                    int[] genotype_j = genotypes[v];
                    double rate = 0;
                    if (hammingDistance(genotype_i, genotype_j) < 2) { // if genotype j can be reached by a single mutation
                    	int mutDirection = mutationDirection(genotype_i, genotype_j);
                        if (mutDirection < 0) {
                        	rate = mutRate;// forward mutation
                        } else {
                        	rate = mutRate; // backward mutation
                        }
                    }
                    fitSpaceGammas[index] = rate;
        		}	
        	}
		}  
		
		
	}
	
	@Override
	public void setUp(Integer[] states, int sites) {
		
		int maxState = states[0]; // assumes all sites have same # of states
		seqStates = states;
		seqSites = sites;
		siteFitnessEffects = new double[seqSites][maxState];
		
		fitEffects = fitEffectsInput.get().getValues();
		double mutRate = mutationMatrix.get().getValues()[0];
		//double backMutRate = mutationMatrix.get().getValues()[1];
		
		// Define which genotypes to track - with 60 genotypes
		genotypes = new int[53][seqSites];
		genotypes[0] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[1] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[2] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[3] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[4] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0};
		genotypes[5] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
		genotypes[6] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0};
		genotypes[7] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0};
		genotypes[8] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0};
		genotypes[9] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0};
		genotypes[10] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
		genotypes[11] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,1,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
		genotypes[12] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
		genotypes[13] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0};
		genotypes[14] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0};
		genotypes[15] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0};
		genotypes[16] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0};
		genotypes[17] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0};
		genotypes[18] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0};
		genotypes[19] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0};
		genotypes[20] = new int[] {0,0,0,0,0,0,0,0,0,0,1,0,1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0};
		genotypes[21] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
		genotypes[22] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
		genotypes[23] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[24] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[25] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[26] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[27] = new int[] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[28] = new int[] {0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[29] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[30] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[31] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[32] = new int[] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[33] = new int[] {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0};
		genotypes[34] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0};
		genotypes[35] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[36] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
		genotypes[37] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
		genotypes[38] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[39] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[40] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0};
		genotypes[41] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
		genotypes[42] = new int[] {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1};
		genotypes[43] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[44] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0};
		genotypes[45] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[46] = new int[] {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[47] = new int[] {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
		genotypes[48] = new int[] {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[49] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[50] = new int[] {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[51] = new int[] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		genotypes[52] = new int[] {0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
        
        // Get genotype fitness values
        fitSpaceValues = new Double[genotypes.length];
        for (int i=0; i<genotypes.length; i++) {
        	fitSpaceValues[i] = fitEffects[i];
        }
        
		// Set up fitness space 
		fitSpaceStates = genotypes.length;
		
		// Set up fitness space gammas
		fitSpaceGammas = new Double[fitSpaceStates*(fitSpaceStates-1)];
        
        // Compute transition rates between genotypes
        for (int u=0; u<genotypes.length; u++) {
        	for (int v=0; v<genotypes.length; v++) {
        		if (u != v) {
        			int index = u * (fitSpaceStates - 1) + ((u > v) ? v : v-1);
        			int[] genotype_i = genotypes[u];
                    int[] genotype_j = genotypes[v];
                    double rate = 0;
                    if (hammingDistance(genotype_i, genotype_j) < 2) { // if genotype j can be reached by a single mutation
                    	int mutDirection = mutationDirection(genotype_i, genotype_j);
                        if (mutDirection < 0) {
                        	rate = mutRate;// forward mutation
                        } else {
                        	rate = mutRate; // backward mutation
                        }
                    }
                    fitSpaceGammas[index] = rate;
        		}	
        	}
		}  
		
	}
	
	static ArrayList<Integer> indexOfAll(Integer element, ArrayList<Integer> list){
	    ArrayList<Integer> indexList = new ArrayList<Integer>();
	    for (int i = 0; i < list.size(); i++)
	        if(element.equals(list.get(i)))
	            indexList.add(i);
	    return indexList;
	}
	
	static int hammingDistance(int[] genotype1, int[] genotype2) {
		int dist = 0;
		for(int i = 0; i < genotype1.length; i++) {
			if (genotype1[i] != genotype2[i]) dist++;
		}
		return dist;
	}
	
	static int mutationDirection(int[] genotype1, int[] genotype2) {
		int dist = 0;
		for(int i = 0; i < genotype1.length; i++) {
			dist += genotype1[i] - genotype2[i];
		}
		return dist;
	}
	
	@Override
	public int getClosestFitSpaceState(double fit) {
		
        // Returns discrete point in fitness space closest to input fitValue
		double minAbsDistance = Double.POSITIVE_INFINITY;
		double absDiff = 0;
		int fitState = -1;
		for (int k=0; k<fitSpaceStates; k++) {
			absDiff = Math.abs(fitSpaceValues[k] - fit);
			if (absDiff < minAbsDistance) {
				fitState = k;
				minAbsDistance = absDiff;
			}
		}
		
		return fitState;
	}
	
	public double getMinDistance(ArrayList<Double> fitValues, double val) {
		
        // Returns minimum distance in fitValues from val
		double minAbsDistance = Double.POSITIVE_INFINITY;
		double absDiff = 0;
		for (int k=0; k<fitValues.size(); k++) {
			absDiff = Math.abs(fitValues.get(k) - val);
			if (absDiff < minAbsDistance) {
				minAbsDistance = absDiff;
			}
		}
		return minAbsDistance;
		
	}
	
	@Override
	public double[][] getMarginalFitnessEffects(double[][] g) {
		
		// Assumes each genotype has its own fitness value

		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] = g[site][state] / siteSum;
			}
		}
		
		// Compute marginalFitEffects conditional on state of each site
		double mf = 0; double[] gpr = new double[genotypes.length]; int type = 0; double gprSum = 0;
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		marginalFitEffects[seqSites] = new double[1];
		for (int geno=0; geno<genotypes.length; geno++) { // marginalize over all possible genotypes
			gpr[geno] = 1.0; // prob of being in this genotype
			for (int j=0; j<seqSites; j++) {
				type = genotypes[geno][j];
				gpr[geno] *= norm_g[j][type];
			}
			gprSum += gpr[geno];
		}
		
		// Normalize genotype probs and compute marginal fitness
		mf = 0; // marginal fitness
		if (gprSum <= 0) {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] / genotypes.length;
			}
		} else {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] * gpr[geno] / gprSum;
			}
		}
		
		marginalFitEffects[seqSites][0] = mf;
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
				marginalFitEffects[site][state] = mf;
			}
		}
		
		return marginalFitEffects;
	}
	
	// @Override
	public double[][] getMarginalFitnessEffectsDefunct(double[][] g) {
		
		// Assumes each genotype has its own fitness value

		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] = g[site][state] / siteSum;
			}
		}
		
		// Compute marginalFitEffects conditional on state of each site
		double mf = 0; double[] gpr = new double[genotypes.length]; int type = 0; double gprSum;
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
	
				gpr = new double[genotypes.length]; 
				gprSum = 0;
				for (int geno=0; geno<genotypes.length; geno++) { // marginalize over all possible genotypes
					if (genotypes[geno][site] == state) { // otherwise no prob of being in this genotype
						gpr[geno] = 1.0; // prob of being in this genotype
						for (int j=0; j<seqSites; j++) {
							if (site != j) { // already accounted for this site above
								type = genotypes[geno][j];
								gpr[geno] *= norm_g[j][type];
							}
						}
					gprSum += gpr[geno];
					}
				}
				
				// Normalize genotype probs and compute marginal fitness
				mf = 0; // marginal fitness
				if (gprSum <= 0) {
					for (int geno=0; geno<genotypes.length; geno++) {
						mf += fitEffects[geno] / genotypes.length;
					}
				} else {
					for (int geno=0; geno<genotypes.length; geno++) {
						mf += fitEffects[geno] * gpr[geno] / gprSum;
					}
				}
				marginalFitEffects[site][state] = mf;
				
			}
		}
		
		// Compute marginalFitness for full lineage
		marginalFitEffects[seqSites] = new double[1];
		mf = 0; // marginal fitness
		gpr = new double[genotypes.length];
		gprSum = 0;
		for (int geno=0; geno<genotypes.length; geno++) { // marginalize over all possible genotypes
			gpr[geno] = 1.0; // prob of being in this genotype
			for (int j=0; j<seqSites; j++) {
				type = genotypes[geno][j];
				gpr[geno] *= norm_g[j][type];
			}
			gprSum += gpr[geno];
		}
		
		// Normalize genotype probs and compute marginal fitness
		mf = 0; // marginal fitness
		if (gprSum <= 0) {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] / genotypes.length;
			}
		} else {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] * gpr[geno] / gprSum;
			}
		}
		marginalFitEffects[seqSites][0] = mf;
		
		return marginalFitEffects;
	}
	
	@Override
	public double[][] getMarginalFitnessEffectsSN(SmallNumber[][] g) {
		
		// Same as getMarginalFitnessEffects but works on SmallNumber input
		// SmallNumbers are only reverted back to doubles after normalizing marginal site densities
		// For this model can we compute marginalFitness once for all sites/states?
		
		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
			}
		}
		
		// Compute marginalFitness for full lineage
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites		
		marginalFitEffects[seqSites] = new double[1];
		double mf = 0; // marginal fitness
		double[] gpr = new double[genotypes.length];
		int type = 0;
		double gprSum = 0;
		for (int geno=0; geno<genotypes.length; geno++) { // marginalize over all possible genotypes
			gpr[geno] = 1.0; // prob of being in this genotype
			for (int j=0; j<seqSites; j++) {
				type = genotypes[geno][j];
				gpr[geno] *= norm_g[j][type];
			}
			gprSum += gpr[geno];
		}
		
		// Normalize genotype probs and compute marginal fitness
		if (gprSum <= 0) {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] / genotypes.length;
			}
		} else {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] * gpr[geno] / gprSum;
			}
		}
		
		marginalFitEffects[seqSites][0] = mf;
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
				marginalFitEffects[site][state] = mf;	
			}
		}
		
		return marginalFitEffects;
	}
	
	// @Override
	public double[][] getMarginalFitnessEffectsSNDefunct(SmallNumber[][] g) {
		
		// Same as getMarginalFitnessEffects but works on SmallNumber input
		// SmallNumbers are only reverted back to doubles after normalizing marginal site densities
		// For this model can we compute marginalFitness once for all sites/states?
		
		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
			}
		}
		
		// Compute marginalFitEffects conditional on state of each site
		double mf = 0; double[] gpr = new double[genotypes.length]; int type = 0; double gprSum;
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
	
				gpr = new double[genotypes.length]; 
				gprSum = 0;
				for (int geno=0; geno<genotypes.length; geno++) { // marginalize over all possible genotypes
					if (genotypes[geno][site] == state) { // otherwise no prob of being in this genotype
						gpr[geno] = 1.0; // prob of being in this genotype
						for (int j=0; j<seqSites; j++) {
							if (site != j) { // already accounted for this site above
								type = genotypes[geno][j];
								gpr[geno] *= norm_g[j][type];
							}
						}
					gprSum += gpr[geno];
					}
				}
				
				// Normalize genotype probs and compute marginal fitness
				mf = 0; // marginal fitness
				if (gprSum <= 0) {
					for (int geno=0; geno<genotypes.length; geno++) {
						mf += fitEffects[geno] / genotypes.length;
					}
				} else {
					for (int geno=0; geno<genotypes.length; geno++) {
						mf += fitEffects[geno] * gpr[geno] / gprSum;
					}
				}
				marginalFitEffects[site][state] = mf;
				
			}
		}
		
		// Compute marginalFitness for full lineage
		marginalFitEffects[seqSites] = new double[1];
		mf = 0; // marginal fitness
		gpr = new double[genotypes.length];
		gprSum = 0;
		for (int geno=0; geno<genotypes.length; geno++) { // marginalize over all possible genotypes
			gpr[geno] = 1.0; // prob of being in this genotype
			for (int j=0; j<seqSites; j++) {
				type = genotypes[geno][j];
				gpr[geno] *= norm_g[j][type];
			}
			gprSum += gpr[geno];
		}
		
		// Normalize genotype probs and compute marginal fitness
		mf = 0; // marginal fitness
		if (gprSum <= 0) {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] / genotypes.length;
			}
		} else {
			for (int geno=0; geno<genotypes.length; geno++) {
				mf += fitEffects[geno] * gpr[geno] / gprSum;
			}
		}
		marginalFitEffects[seqSites][0] = mf;
		
		return marginalFitEffects;
	}
	
	/**
	 * Methods for logging
	 */
	@Override
	public void init(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void log(int sample, PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
