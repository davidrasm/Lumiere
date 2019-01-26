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
@Description("4 Genotype fitness model" +
             "where fitness depends on genotype of lineage -- for testing approximations of MFBD model.")
public class FourGenotypeFitnessModel extends MultiSiteFitnessModel {

	// Input for site-wise fitness effects
	public Input<RealParameter> fitEffectsInput =
			new Input<>("mutFitEffects", "Site specific fit effects");
	
	// Input for fitness space values
	public Input<RealParameter> fitSpaceValuesInput =
			new Input<>("fitSpaceValues", "Points in fitness space to tracks");
	
	// Mutation rates
	public Input<RealParameter> mutationMatrix =
			new Input<>("mutationMatrix", "Flattened mutation matrix, can be asymmetric, diagnonal entries omitted");
	
	protected Double[] fitEffects;
	protected int[][] genotypes;
	
	@Override
	public void update() {

		fitEffects = fitEffectsInput.get().getValues();
		
		// Only one 'site' under this model since we're tracking all genotypes
		siteFitnessEffects[0][0] = 1.0; // genotype 00
		siteFitnessEffects[0][1] = fitEffects[0]; // genotype 10
		siteFitnessEffects[0][2] = fitEffects[1]; // genotype 01
		siteFitnessEffects[0][3] = fitEffects[0] * fitEffects[1]; // genotype 11 -- assuming multiplicative fitness
		
		// Set up fitness space
		fitSpaceValues = new Double[seqStates[0]];
		fitSpaceValues[0] = 1.0;
		fitSpaceValues[1] = fitEffects[0];
		fitSpaceValues[2] = fitEffects[1];
		fitSpaceValues[3] = fitEffects[0] * fitEffects[1]; 
		fitSpaceStates = fitSpaceValues.length;
		
		// Set up fitness space gammas
		fitSpaceGammas = mutationMatrix.get().getValues();
		
	}
	
	@Override
	public void setUp(Integer[] states, int sites) {
		
		int maxState = states[0]; // assumes all sites have same # of states
		seqStates = states;
		seqSites = sites;
		siteFitnessEffects = new double[seqSites][maxState];
		
		fitEffects = fitEffectsInput.get().getValues();
		
		// Only one 'site' under this model since we're tracking all genotypes
		siteFitnessEffects[0][0] = 1.0; // genotype 00
		siteFitnessEffects[0][1] = fitEffects[0]; // genotype 10
		siteFitnessEffects[0][2] = fitEffects[1]; // genotype 01
		siteFitnessEffects[0][3] = fitEffects[0] * fitEffects[1]; // genotype 11 -- assuming multiplicative fitness
		
		// Set up fitness space - values should be Double[]?
		fitSpaceValues = new Double[maxState];
		fitSpaceValues[0] = 1.0;
		fitSpaceValues[1] = fitEffects[0];
		fitSpaceValues[2] = fitEffects[1];
		fitSpaceValues[3] = fitEffects[0] * fitEffects[1]; 
		fitSpaceStates = fitSpaceValues.length;
		
		// Set up fitness space gammas
		fitSpaceGammas = mutationMatrix.get().getValues();        
		
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
