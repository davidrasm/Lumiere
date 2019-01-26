package beast.evolution.fitnessmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * @author David Rasmussen
 */
@Description("Multi-site fitness model" +
        		"allowing for site-specific mutation effects. P0 values are computed exactly." +
					"Assumes seqs are binary and fitness effects are multiplicative across sites")
public class SiteSpecificFitnessExactP0Model extends MultiSiteFitnessModel {

	// Input for site-wise fitness effects
	public Input<RealParameter> fitEffectsInput =
			new Input<>("mutFitEffects", "Site specific fit effects");
	
	// Input for fitness space values
	//public Input<RealParameter> fitSpaceValuesInput =
			//new Input<>("fitSpaceValues", "Points in fitness space to tracks");
	
	// Mutation rates
	public Input<RealParameter> mutationMatrix =
			new Input<>("mutationMatrix", "Flattened mutation matrix, can be asymmetric, diagnonal entries omitted");
	
	protected Double[] fitEffects;
	
	@Override
	public void update() {

		fitEffects = fitEffectsInput.get().getValues();
		double forwardMutRate = mutationMatrix.get().getValues()[0];
		double backMutRate = mutationMatrix.get().getValues()[1];
		
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k][1] = fitEffects[k]; // No longer 1-fitEffects
		}
		
        // Get genotype fitness values - needs to be updated
        Double[] genotypeFitValues = new Double[genotypes.length];
        for (int i=0; i<genotypes.length; i++) {
        	double fit = 1.0;
        	for (int j=0; j<seqSites; j++) {
        		fit *= siteFitnessEffects[j][genotypes[i][j]];
        	}
        	genotypeFitValues[i] = fit;
        }
        
        // Assign genotypes to discrete points in fitness space - needs to be updated
        ArrayList<Integer> genotypeFitTypes = new ArrayList<Integer>();
        for (int i=0; i<genotypes.length; i++) {
        	genotypeFitTypes.add(getClosestFitSpaceState(genotypeFitValues[i]));
        }
        
		// Set up fitness space - values should be Double[]?
		fitSpaceValues = genotypeFitValues;
		fitSpaceStates = genotypes.length;
		
        // Compute transition rates between genotypes - rates will not change unless mutation rates change
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
                        	rate = forwardMutRate;// forward mutation
                        } else {
                        	rate = backMutRate; // backward mutation
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
		double forwardMutRate = mutationMatrix.get().getValues()[0];
		double backMutRate = mutationMatrix.get().getValues()[1];
		
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k][0] = 1.0;
			siteFitnessEffects[k][1] = fitEffects[k]; // No longer 1-fitEffects
		}
		
		// Enumerate genotypes - only need to do once
        genotypes = new int[(int) Math.pow(maxState,seqSites)][seqSites];
        for (int i=0; i<seqSites; i++) {     
            int tiles = (int) Math.pow(maxState,i+1);
            int size = genotypes.length / tiles;
            for (int j=0; j<(tiles/2); j++) {
                int startLoc = (2 * j * size);
                int endLoc = startLoc + size;
                for (int k=startLoc; k<endLoc; k++) {
                	genotypes[k][i] = 1;
                }
            }  
        }
        
        // Get genotype fitness values - needs to be updated
        Double[] genotypeFitValues = new Double[genotypes.length];
        for (int i=0; i<genotypes.length; i++) {
        	double fit = 1.0;
        	for (int j=0; j<seqSites; j++) {
        		fit *= siteFitnessEffects[j][genotypes[i][j]];
        	}
        	genotypeFitValues[i] = fit;
        }
        
        // Assign genotypes to discrete points in fitness space - needs to be updated
        //ArrayList<Integer> genotypeFitTypes = new ArrayList<Integer>();
        //for (int i=0; i<genotypes.length; i++) {
        	//genotypeFitTypes.add(getClosestFitSpaceState(genotypeFitValues[i]));
        //}
		
		// Set up fitness space - values should be Double[]?
		fitSpaceValues = genotypeFitValues;
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
                        	rate = forwardMutRate;// forward mutation
                        } else {
                        	rate = backMutRate; // backward mutation
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
