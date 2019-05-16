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
             "used for Ebola analysis tracking 9 genotypes + geographic effects. Includes one additional site representing geo location. Uses exact P0s for each genotype")
public class EbolaGeoFitnessModel extends MultiSiteFitnessModel {
	
	// Input for site-wise fitness effects
	public Input<RealParameter> fitEffectsInput =
			new Input<>("mutFitEffects", "Site specific fit effects");
	
	public Input<RealParameter> geoEffectsInput =
			new Input<>("geoFitEffects", "Geo location fit effects");
	
	// Mutation rates
	public Input<RealParameter> mutationMatrix =
			new Input<>("mutationMatrix", "Flattened mutation matrix, can be asymmetric, diagnonal entries omitted");
	
	protected Double[] fitEffects;
	protected Double[] geoEffects;
	
	@Override
	public void update() {
		
		geoEffects = geoEffectsInput.get().getValues();
		geoEffects[0] = 1.0; // hard coded this so effect of geo locs is always relative to Guinea [0]
		//fitEffects = fitEffectsInput.get().getValues();
		
		int totalGenos = fitEffectsInput.get().getValues().length * geoEffects.length;
		fitEffects = new Double[totalGenos];
		for (int fe=0; fe<fitEffectsInput.get().getValues().length; fe++) {
			for (int gl=0; gl<geoEffects.length; gl++) {
				fitEffects[(fe*geoEffects.length) + gl] = fitEffectsInput.get().getValues()[fe] * geoEffects[gl];
			}
		}
		
		double forwardMutRate = mutationMatrix.get().getValues()[0]; // 0 -> 1 mutations
		double backMutRate = mutationMatrix.get().getValues()[2]; // 1 -> 0 mutations
		
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
		
		seqStates = states;
		seqSites = sites;
		
		geoEffects = geoEffectsInput.get().getValues();
		geoEffects[0] = 1.0; // hard coded this so effect of geo locs is always relative to Guinea [0]
		int totalGenos = fitEffectsInput.get().getValues().length * geoEffects.length;
		fitEffects = new Double[totalGenos];
		for (int fe=0; fe<fitEffectsInput.get().getValues().length; fe++) {
			for (int gl=0; gl<geoEffects.length; gl++) {
				fitEffects[(fe*geoEffects.length) + gl] = fitEffectsInput.get().getValues()[fe] * geoEffects[gl];
			}
		}
		
		double forwardMutRate = mutationMatrix.get().getValues()[0]; // 0 -> 1 mutations
		double backMutRate = mutationMatrix.get().getValues()[2]; // 1 -> 0 mutations
		
		// Define which genotypes to track - with 9 genotypes x 3 geo locations
		genotypes = new int[totalGenos][seqSites];
		genotypes[0] = new int[] {0,0,0,0,0,0,0,0,0,0}; // wt
		genotypes[1] = new int[] {0,0,0,0,0,0,0,0,0,1}; // wt
		genotypes[2] = new int[] {0,0,0,0,0,0,0,0,0,2}; // wt
		genotypes[3] = new int[] {0,1,0,0,0,0,0,0,0,0}; // A82V
		genotypes[4] = new int[] {0,1,0,0,0,0,0,0,0,1}; // A82V
		genotypes[5] = new int[] {0,1,0,0,0,0,0,0,0,2}; // A82V
		genotypes[6] = new int[] {0,0,0,0,1,0,0,0,0,0}; // P330S
		genotypes[7] = new int[] {0,0,0,0,1,0,0,0,0,1}; // P330S
		genotypes[8] = new int[] {0,0,0,0,1,0,0,0,0,2}; // P330S
		genotypes[9] = new int[] {0,0,1,0,1,0,0,0,1,0}; // P330S+N107D+G480D
		genotypes[10] = new int[] {0,0,1,0,1,0,0,0,1,1}; // P330S+N107D+G480D
		genotypes[11] = new int[] {0,0,1,0,1,0,0,0,1,2}; // P330S+N107D+G480D
		genotypes[12] = new int[] {0,1,0,0,0,0,1,0,0,0}; // A82V+R410S
		genotypes[13] = new int[] {0,1,0,0,0,0,1,0,0,1}; // A82V+R410S
		genotypes[14] = new int[] {0,1,0,0,0,0,1,0,0,2}; // A82V+R410S
		genotypes[15] = new int[] {0,1,0,0,0,0,1,1,0,0}; // A82V+R410S+K439E
		genotypes[16] = new int[] {0,1,0,0,0,0,1,1,0,1}; // A82V+R410S+K439E
		genotypes[17] = new int[] {0,1,0,0,0,0,1,1,0,2}; // A82V+R410S+K439E
		genotypes[18] = new int[] {1,1,0,0,0,0,0,0,0,0}; // A82V+R29K
		genotypes[19] = new int[] {1,1,0,0,0,0,0,0,0,1}; // A82V+R29K
		genotypes[20] = new int[] {1,1,0,0,0,0,0,0,0,2}; // A82V+R29K
		genotypes[21] = new int[] {0,1,0,1,0,0,0,0,0,0}; // A82V+T230A
		genotypes[22] = new int[] {0,1,0,1,0,0,0,0,0,1}; // A82V+T230A
		genotypes[23] = new int[] {0,1,0,1,0,0,0,0,0,2}; // A82V+T230A
		genotypes[24] = new int[] {0,1,0,0,0,1,0,0,0,0}; // A82V+I371V
		genotypes[25] = new int[] {0,1,0,0,0,1,0,0,0,1}; // A82V+I371V
		genotypes[26] = new int[] {0,1,0,0,0,1,0,0,0,2}; // A82V+I371V
        
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
						for (int j=0; j<seqSites; j++) { // was seqSites - 1 because last site represents geo loc
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
						for (int j=0; j<seqSites; j++) { // was seqSites - 1 because last site represents geo loc
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
