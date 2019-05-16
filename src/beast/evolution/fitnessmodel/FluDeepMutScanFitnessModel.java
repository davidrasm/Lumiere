package beast.evolution.fitnessmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections; 
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
             "used for flu H3N2 analysis with mutational fitness effects predicted by deep mutational scanning.")
public class FluDeepMutScanFitnessModel extends MultiSiteFitnessModel {

	// Input for site-wise fitness effects
	public Input<RealParameter> fitEffectsInput =
			new Input<>("mutFitEffects", "Site specific fit effects");
	
	// Input for fitness space values
	//public Input<RealParameter> fitSpaceValuesInput =
			//new Input<>("fitSpaceValues", "Points in fitness space to tracks");
	
	// Mutation rates
	public Input<RealParameter> mutationMatrix =
			new Input<>("mutationMatrix", "Flattened mutation matrix, can be asymmetric, diagnonal entries omitted");
	
	// Scaling coefficient for mutFitEffects
	public Input<RealParameter> fitEffectsScalarInput =
			new Input<>("mutFitEffectsScalar", "Scaling coefficient for mutFitEffects");
	
	// Scaling exponent for mutFitEffects
	public Input<RealParameter> fitEffectsExponentInput =
			new Input<>("mutFitEffectsExponent", "Scaling exponent for mutFitEffects");
	
	protected Double[] fitEffects;
	protected double fitEffectsScalar;
	protected double fitEffectsExponent;
	
	@Override
	public void update() {

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
		
		fitEffectsScalar = fitEffectsScalarInput.get().getValue();
		fitEffectsExponent = fitEffectsExponentInput.get().getValue();
		
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
		
		fitEffectsScalar = fitEffectsScalarInput.get().getValue();
		fitEffectsExponent = fitEffectsExponentInput.get().getValue();
		
		/**
		 * Find min and max fitness values of all genotypes
		 * Could just pick reasonable min and maxFitVals like 0 and 3.0
		 */
        Double minFitVal = 0.1;
        Double maxFitVal = 2.0;
        //for (int k=0; k<seqSites; k++) {
        	//double min = Arrays.stream(siteFitnessEffects[k]).min().getAsDouble();
        	//double max = Arrays.stream(siteFitnessEffects[k]).max().getAsDouble();
        	//minFitVal += min - 1.0;
        	//maxFitVal += max - 1.0;
        //}
        //minFitVal =  1.0 + Math.pow(fitEffectsScalar * minFitVal, fitEffectsExponent);
        //minFitVal = Math.max(0.0, minFitVal);
        //maxFitVal =  1.0 + Math.pow(fitEffectsScalar * maxFitVal, fitEffectsExponent);
        
        // Dynamic grid between min and max genotype fitness values
        double dFitStep = 0.1;
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
		Arrays.fill(fitSpaceGammas, 0.0);
		
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
	
	public double log2(double x) {
		return Math.log(x) / Math.log(2.0);
	}
	
	@Override
	public double[][] getMarginalFitnessEffects(double[][] g) {
		
		// Additive fitness effects across sites with fitness remapping 
		
		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
		double sumOfWeightedFitEffects = 0.0;
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
			double w = 0; // weighted fitness effect of site
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] = g[site][state] / siteSum;
				w += norm_g[site][state] * siteFitnessEffects[site][state];
			}
			weightedFitEffects[site] = log2(w);
			sumOfWeightedFitEffects += log2(w);
		}
		
		// Marginal fitness at each site
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		double fv;
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
		        fv = log2(siteFitnessEffects[site][state]) + (sumOfWeightedFitEffects - weightedFitEffects[site]);
				fv = 1.0 + (fitEffectsScalar * fv);
		        fv = Math.max(0.01, fv);
		        fv = Math.pow(fv, fitEffectsExponent);
				marginalFitEffects[site][state] = fv; 
			}
		}
		
		// Marginal fitness for entire lineage
		marginalFitEffects[seqSites] = new double[1];
		fv = sumOfWeightedFitEffects;
		fv = 1.0 + (fitEffectsScalar * fv);
        fv = Math.max(0.01, fv);
        fv = Math.pow(fv, fitEffectsExponent);
		marginalFitEffects[seqSites][0] = fv; 

		return marginalFitEffects;
	}
	
	@Override
	public double[][] getMarginalFitnessEffectsSN(SmallNumber[][] g) {
		
		// Same as getMarginalFitnessEffects but works on SmallNumber input
		// SmallNumbers are only reverted back to doubles after normalizing marginal site densities
		// Additive fitness effects across sites with fitness remapping 
		
		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
		double sumOfWeightedFitEffects = 0.0;
		for (int site=0; site<seqSites; site++) {
			norm_g[site] = new double[seqStates[site]];
			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
			double w = 0; // weighted fitness effect of site
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
				w += norm_g[site][state] * siteFitnessEffects[site][state];
			}
			weightedFitEffects[site] = log2(w);
			sumOfWeightedFitEffects += log2(w);
		}
		
		// Marginal fitness at each site
		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
		double fv;
		for (int site=0; site<seqSites; site++) {
			marginalFitEffects[site] = new double[seqStates[site]];
			for (int state=0; state<seqStates[site]; state++) {
		        fv = log2(siteFitnessEffects[site][state]) + (sumOfWeightedFitEffects - weightedFitEffects[site]);
				fv = 1.0 + (fitEffectsScalar * fv);
		        fv = Math.max(0.01, fv);
		        fv = Math.pow(fv, fitEffectsExponent);
				marginalFitEffects[site][state] = fv; 
			}
		}
		marginalFitEffects[seqSites] = new double[1];
		fv = sumOfWeightedFitEffects;
		fv = 1.0 + (fitEffectsScalar * fv);
        fv = Math.max(0.01, fv);
        fv = Math.pow(fv, fitEffectsExponent);
		marginalFitEffects[seqSites][0] = fv; 
		
		return marginalFitEffects;
	}
	
//	@Override
//	public double[][] getMarginalFitnessEffects(double[][] g) {
//		
//		// Multiplicative fitness effects across sites with scaling exponents
//		
//		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
//		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
//		double prodOfWeightedFitEffects = 1.0;
//		for (int site=0; site<seqSites; site++) {
//			norm_g[site] = new double[seqStates[site]];
//			double siteSum = 0;
//			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
//			double w = 0; // weighted fitness effect of site
//			for (int state=0; state<seqStates[site]; state++) {
//				norm_g[site][state] = g[site][state] / siteSum;
//				w += norm_g[site][state] * siteFitnessEffects[site][state];
//			}
//			weightedFitEffects[site] = w;
//			prodOfWeightedFitEffects *= w;
//		}
//		
//		// Marginal fitness at each site
//		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
//		double fv;
//		for (int site=0; site<seqSites; site++) {
//			marginalFitEffects[site] = new double[seqStates[site]];
//			for (int state=0; state<seqStates[site]; state++) {
//				fv = siteFitnessEffects[site][state] * prodOfWeightedFitEffects / weightedFitEffects[site];
//				marginalFitEffects[site][state] = Math.pow(fv, fitEffectsExponent);
//			}
//		}
//		
//		// Marginal fitness for entire lineage
//		marginalFitEffects[seqSites] = new double[1];
//		marginalFitEffects[seqSites][0] = Math.pow(prodOfWeightedFitEffects, fitEffectsExponent);
//		
//		return marginalFitEffects;
//	}
	
//	@Override
//	public double[][] getMarginalFitnessEffectsSN(SmallNumber[][] g) {
//		
//		// Multiplicative fitness effects across sites with scaling exponents
//		
//		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
//		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
//		double prodOfWeightedFitEffects = 0.0;
//		for (int site=0; site<seqSites; site++) {
//			norm_g[site] = new double[seqStates[site]];
//			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
//			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
//			double w = 0; // weighted fitness effect of site
//			for (int state=0; state<seqStates[site]; state++) {
//				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
//				w += norm_g[site][state] * siteFitnessEffects[site][state];
//			}
//			weightedFitEffects[site] = w;
//			prodOfWeightedFitEffects *= w; // -1.0 so neutral mutations have no effect on sum over sites
//		}
//		
//		// Marginal fitness at each site
//		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
//		double fv;
//		for (int site=0; site<seqSites; site++) {
//			marginalFitEffects[site] = new double[seqStates[site]];
//			for (int state=0; state<seqStates[site]; state++) {
//				fv = siteFitnessEffects[site][state] * prodOfWeightedFitEffects / weightedFitEffects[site];
//				marginalFitEffects[site][state] = Math.pow(fv, fitEffectsExponent);
//			}
//		}
//		marginalFitEffects[seqSites] = new double[1];
//		marginalFitEffects[seqSites][0] = Math.pow(prodOfWeightedFitEffects, fitEffectsExponent); 
//		
//		return marginalFitEffects;
//	}
	
//	@Override
//	public double[][] getMarginalFitnessEffects(double[][] g) {
//		
//		// Additive fitness effects across sites with fitness remapping 
//		
//		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
//		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
//		double sumOfWeightedFitEffects = 0.0;
//		for (int site=0; site<seqSites; site++) {
//			norm_g[site] = new double[seqStates[site]];
//			double siteSum = 0;
//			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
//			double w = 0; // weighted fitness effect of site
//			for (int state=0; state<seqStates[site]; state++) {
//				norm_g[site][state] = g[site][state] / siteSum;
//				w += norm_g[site][state] * siteFitnessEffects[site][state];
//			}
//			weightedFitEffects[site] = w - 1.0;
//			sumOfWeightedFitEffects += w - 1.0; // -1.0 so neutral mutations have no effect on sum over sites
//		}
//		
//		// Marginal fitness at each site
//		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
//		double fv;
//		for (int site=0; site<seqSites; site++) {
//			marginalFitEffects[site] = new double[seqStates[site]];
//			for (int state=0; state<seqStates[site]; state++) {
//		        fv = (siteFitnessEffects[site][state] - 1.0) + (sumOfWeightedFitEffects - weightedFitEffects[site]);
//				fv = 1.0 + Math.pow(fitEffectsScalar * fv, fitEffectsExponent);
//		        fv = Math.max(0.01, fv);
//				marginalFitEffects[site][state] = fv; 
//			}
//		}
//		
//		// Marginal fitness for entire lineage
//		marginalFitEffects[seqSites] = new double[1];
//		fv = sumOfWeightedFitEffects;
//		fv = 1.0 + Math.pow(fitEffectsScalar * fv, fitEffectsExponent);
//        fv = Math.max(0.01, fv);
//		marginalFitEffects[seqSites][0] = fv; 
//
//		return marginalFitEffects;
//	}
	
//	@Override
//	public double[][] getMarginalFitnessEffectsSN(SmallNumber[][] g) {
//		
//		// Same as getMarginalFitnessEffects but works on SmallNumber input
//		// SmallNumbers are only reverted back to doubles after normalizing marginal site densities
//		// Additive fitness effects across sites with fitness remapping 
//		
//		double[][] norm_g = new double[seqSites][]; // g normalized by sum at each site
//		double [] weightedFitEffects = new double[seqSites]; // fitness effect of each site weighted by norm_g
//		double sumOfWeightedFitEffects = 0.0;
//		for (int site=0; site<seqSites; site++) {
//			norm_g[site] = new double[seqStates[site]];
//			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
//			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
//			double w = 0; // weighted fitness effect of site
//			for (int state=0; state<seqStates[site]; state++) {
//				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
//				w += norm_g[site][state] * siteFitnessEffects[site][state];
//			}
//			weightedFitEffects[site] = w - 1.0;
//			sumOfWeightedFitEffects += w - 1.0; // -1.0 so neutral mutations have no effect on sum over sites
//		}
//		
//		// Marginal fitness at each site
//		double[][] marginalFitEffects = new double[seqSites+1][]; // site fitness effects marginalized over all other sites
//		double fv;
//		for (int site=0; site<seqSites; site++) {
//			marginalFitEffects[site] = new double[seqStates[site]];
//			for (int state=0; state<seqStates[site]; state++) {
//		        fv = (siteFitnessEffects[site][state] - 1.0) + (sumOfWeightedFitEffects - weightedFitEffects[site]);
//				fv =  1.0 + Math.pow(fitEffectsScalar * fv, fitEffectsExponent);
//		        fv = Math.max(0.01, fv);
//				marginalFitEffects[site][state] = fv;
//			}
//		}
//		marginalFitEffects[seqSites] = new double[1];
//		fv = sumOfWeightedFitEffects;
//		fv = 1.0 + Math.pow(fitEffectsScalar * fv, fitEffectsExponent);
//        fv = Math.max(0.01, fv);
//		marginalFitEffects[seqSites][0] = fv; 
//		
//		return marginalFitEffects;
//	}
	
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
