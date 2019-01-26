package beast.evolution.fitnessmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.OptionalDouble;
import java.util.stream.DoubleStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * @author David Rasmussen
 */
@Description("Multi-site fitness model" +
             "where all sites are neutral and have no effect on fitness.")
public class FluSiteSpecificFitnessModel extends MultiSiteFitnessModel {
	
	

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
		
		//for (int k=0; k<seqSites; k++) {
			//siteFitnessEffects[k][1] = fitEffects[k]; // No longer 1-fitEffects
		//}
		
		siteFitnessEffects = new double[seqSites][];
		int fitIndex = 0;
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k] = new double[seqStates[k]];
			siteFitnessEffects[k][0] = 1.0;
			for (int l=1; l<seqStates[k]; l++) {
				siteFitnessEffects[k][l] = fitEffects[fitIndex];
				fitIndex++;
			}
		}
		
        // Get genotype fitness values - needs to be updated
        double[] genotypeFitValues = new double[genotypes.length];
        for (int i=0; i<genotypes.length; i++) {
        	double fit = 1.0;
        	for (int j=0; j<seqSites; j++) {
        		fit *= siteFitnessEffects[j][genotypes[i][j]];
        	}
        	genotypeFitValues[i] = fit;
        }
        
        // Find min and max fitness values of all genotypes
        Double minFitVal = Arrays.stream(genotypeFitValues).min().getAsDouble();
        Double maxFitVal = Arrays.stream(genotypeFitValues).max().getAsDouble();
        
        // Dynamic grid between min and max genotype fitness values
        double dFitStep = 0.05;
        if (maxFitVal - minFitVal < 0.0000001) {
        	fitSpaceStates = 1;
        	fitSpaceValues = new Double[] {minFitVal};
        } else {
        	fitSpaceStates = (int) ((maxFitVal - minFitVal) / dFitStep) + 2;
        	fitSpaceValues = new Double[fitSpaceStates];
        	double fitVal = maxFitVal;
        	int cntr = 0;
        	while (fitVal > minFitVal){
        		fitSpaceValues[cntr] = fitVal;
        		fitVal -= dFitStep;
        		cntr++;
        	}
        	fitSpaceValues[fitSpaceStates-1] = minFitVal;
        }
		fitSpaceGammas = new Double[fitSpaceStates*(fitSpaceStates-1)];
        
        // Assign genotypes to discrete points in fitness space - needs to be updated
        ArrayList<Integer> genotypeFitTypes = new ArrayList<Integer>();
        for (int i=0; i<genotypes.length; i++) {
        	genotypeFitTypes.add(getClosestFitSpaceState(genotypeFitValues[i]));
        }
        
        // Compute transition rates in fitness space - needs to be updated
        for (int u=0; u<fitSpaceStates; u++) {
        	for (int v=0; v<fitSpaceStates; v++) {
        		if (u != v) {
        			int index = u * (fitSpaceStates - 1) + ((u > v) ? v : v-1);
        			ArrayList<Integer> types_in_u = indexOfAll(u, genotypeFitTypes);
        			ArrayList<Integer> types_in_v = indexOfAll(v, genotypeFitTypes);
        			if (types_in_u.isEmpty() || types_in_v.isEmpty()) {
        				fitSpaceGammas[index] = 0.0;
        				continue;
        			}
                    double sumOfRates = 0;
        			for (int i=0; i<types_in_u.size(); i++) {
        				for (int j=0; j<types_in_v.size(); j++) {
                            int[] genotype_i = genotypes[types_in_u.get(i)];
                            int[] genotype_j = genotypes[types_in_v.get(j)];
                            if (hammingDistance(genotype_i, genotype_j) < 2) { // if genotype j can be reached by a single mutation
                            	int mutDirection = mutationDirection(genotype_i, genotype_j);
                                if (mutDirection < 0) {
                                    sumOfRates += forwardMutRate;// forward mutation
                                } else {
                                    sumOfRates += backMutRate; // backward mutation
                                }
                            }
        				}
        			}
        			fitSpaceGammas[index] = sumOfRates / types_in_u.size();
        		}
        	}
        }
		
		
	}
	
	@Override
	public void setUp(Integer[] states, int sites) {
		
		seqStates = states;
		seqSites = sites;
		
		//siteFitnessEffects = new double[seqSites][seqStates];
		fitEffects = fitEffectsInput.get().getValues();
		siteFitnessEffects = new double[seqSites][];
		int fitIndex = 0;
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k] = new double[seqStates[k]];
			siteFitnessEffects[k][0] = 1.0;
			for (int l=1; l<seqStates[k]; l++) {
				siteFitnessEffects[k][l] = fitEffects[fitIndex];
				fitIndex++;
			}
		}
		
		double forwardMutRate = mutationMatrix.get().getValues()[0];
		double backMutRate = mutationMatrix.get().getValues()[1];
		
		// Enumerate all single and double mutant genotypes
		ArrayList<int[]> genotypeList = new ArrayList<>();
		for (int i=0; i<seqSites; i++) { 
			for (int j=1; j<seqStates[i]; j++) {
				int[] geno = new int[seqSites];
				geno[i] = j; // single mutant at this site
				genotypeList.add(geno);
				for (int k=0; k<seqSites; k++) {
					if (i!=k) {
						for (int l=1; l<seqStates[k]; l++) {
							geno = new int[seqSites];
							geno[i] = j;
							geno[k] = l;
							genotypeList.add(geno);
						}
					}
				}
			}
		}
		
		// Convert back to int[][]
		genotypes = new int[genotypeList.size()][seqSites];
		for (int i=0; i<genotypeList.size(); i++) {
			genotypes[i] = genotypeList.get(i);
		}
		
        // Get genotype fitness values - needs to be updated
        double[] genotypeFitValues = new double[genotypes.length];
        for (int i=0; i<genotypes.length; i++) {
        	double fit = 1.0;
        	for (int j=0; j<seqSites; j++) {
        		fit *= siteFitnessEffects[j][genotypes[i][j]];
        	}
        	genotypeFitValues[i] = fit;
        }
      
        // Find min and max fitness values of all genotypes
        Double minFitVal = Arrays.stream(genotypeFitValues).min().getAsDouble();
        Double maxFitVal = Arrays.stream(genotypeFitValues).max().getAsDouble();
        
        // Dynamic grid between min and max genotype fitness values
        double dFitStep = 0.05;
        if (maxFitVal - minFitVal < 0.0000001) {
        	fitSpaceStates = 1;
        	fitSpaceValues = new Double[] {minFitVal};
        } else {
        	fitSpaceStates = (int) ((maxFitVal - minFitVal) / dFitStep) + 2;
        	fitSpaceValues = new Double[fitSpaceStates];
        	double fitVal = maxFitVal;
        	int cntr = 0;
        	while (fitVal > minFitVal){
        		fitSpaceValues[cntr] = fitVal;
        		fitVal -= dFitStep;
        		cntr++;
        	}
        	fitSpaceValues[fitSpaceStates-1] = minFitVal;
        }
		fitSpaceGammas = new Double[fitSpaceStates*(fitSpaceStates-1)];
		//System.out.println(fitSpaceValues);
        
        // Assign genotypes to discrete points in fitness space - needs to be updated
        ArrayList<Integer> genotypeFitTypes = new ArrayList<Integer>();
        for (int i=0; i<genotypes.length; i++) {
        	genotypeFitTypes.add(getClosestFitSpaceState(genotypeFitValues[i]));
        }
        
        // Compute transition rates in fitness space - needs to be updated
        for (int u=0; u<fitSpaceStates; u++) {
        	for (int v=0; v<fitSpaceStates; v++) {
        		if (u != v) {
        			int index = u * (fitSpaceStates - 1) + ((u > v) ? v : v-1);
        			ArrayList<Integer> types_in_u = indexOfAll(u, genotypeFitTypes);
        			ArrayList<Integer> types_in_v = indexOfAll(v, genotypeFitTypes);
        			if (types_in_u.isEmpty() || types_in_v.isEmpty()) {
        				fitSpaceGammas[index] = 0.0;
        				continue;
        			}
                    double sumOfRates = 0;
        			for (int i=0; i<types_in_u.size(); i++) {
        				for (int j=0; j<types_in_v.size(); j++) {
                            int[] genotype_i = genotypes[types_in_u.get(i)];
                            int[] genotype_j = genotypes[types_in_v.get(j)];
                            if (hammingDistance(genotype_i, genotype_j) < 2) { // if genotype j can be reached by a single mutation
                            	int mutDirection = mutationDirection(genotype_i, genotype_j);
                                if (mutDirection < 0) {
                                    sumOfRates += forwardMutRate;// forward mutation
                                } else {
                                    sumOfRates += backMutRate; // backward mutation
                                }
                            }
        				}
        			}
        			fitSpaceGammas[index] = sumOfRates / types_in_u.size();
        		}
        	}
        }
        
        //System.out.println(fitSpaceGammas);
		
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
