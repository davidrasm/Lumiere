package beast.evolution.fitnessmodel;

import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.SmallNumber;

/**
 * @author David Rasmussen
 */
@Description("Multi-site fitness model" +
				"This is the Abstract class for all multi-site models")
public abstract class MultiSiteFitnessModel extends CalculationNode implements Loggable {
    
    //public boolean dirty;
    //public boolean reject = false; 
    
	// Sequence space variables
	protected Integer[] seqStates;
    protected int seqSites; // number of sites in seq align
    
    // Fitness space variables
    protected int fitSpaceStates; // nr of discrete states/classes in fitness space
    protected Double[] fitSpaceValues; // fitness values corresponding to the states in fitness space
    protected Double[] fitSpaceGammas; // transition rates between fitness classes (linearly indexed)
    
    protected double[][] siteFitnessEffects;
    
    protected int[][] genotypes;
    
	@Override
	public void initAndValidate() {
	}
	
	public void setUp(Integer[] states, int sites) {
		seqStates = states;
		seqSites = sites;
		siteFitnessEffects = new double[seqSites][];
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k] = new double[states[k]];
		}
	}
	
	public abstract void update();
	
	public abstract int getClosestFitSpaceState(double fit);
	
	public double[][] getSiteFitnessEffects() {
		return siteFitnessEffects;
	}
	
	public int getFitSpaceStates() {
		return fitSpaceStates;
	}
	
	public Double[] getFitnessSpaceValues() {
		return fitSpaceValues;
	}
	
	public Double[] getFitnessSpaceGammas() {
		return fitSpaceGammas;
	}
	
	public int[][] getGenotypes() {
		return genotypes;
	}
	
	public double[][] getMarginalFitnessEffects(double[][] g) {
		
		// Base implementation: Assumes multiplicative fitness effects across sites
		
		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
		double prodOfWeightedFitEffects = 1.0;
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
			double w = 0; // weighted fitness effect of site
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] = g[site][state] / siteSum;
				w += norm_g[site][state] * siteFitnessEffects[site][state];
			}
			weightedFitEffects[site] = w;
			prodOfWeightedFitEffects *= w;
		}
		
		// Marginal fitness at each site
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
				marginalFitEffects[site][state] = siteFitnessEffects[site][state] * prodOfWeightedFitEffects / weightedFitEffects[site]; 
			}
		}
		
		// Marginal fitness for entire lineage
		marginalFitEffects[seqSites] = new double[1];
		marginalFitEffects[seqSites][0] = prodOfWeightedFitEffects;
		
		return marginalFitEffects;
	}
	
	public double[][] getMarginalFitnessEffectsSN(SmallNumber[][] g) {
		
		// Same as getMarginalFitnessEffects but works on SmallNumber input
		// SmallNumbers are only reverted back to doubles after normalizing marginal site densities
		
		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
		double prodOfWeightedFitEffects = 1.0;
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
			double w = 0; // weighted fitness effect of site
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
				w += norm_g[site][state] * siteFitnessEffects[site][state];
			}
			weightedFitEffects[site] = w;
			prodOfWeightedFitEffects *= w;
		}
		
		// Marginal fitness at each site
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
				marginalFitEffects[site][state] = siteFitnessEffects[site][state] * prodOfWeightedFitEffects / weightedFitEffects[site]; 
			}
		}
		marginalFitEffects[seqSites] = new double[1];
		marginalFitEffects[seqSites][0] = prodOfWeightedFitEffects;
		
		return marginalFitEffects;
	}
	

}
