package beast.evolution.fitnessmodel;

import java.io.PrintStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * @author David Rasmussen
 */
@Description("Multi-site fitness model" +
             "where mutations at all sites have same deleterious effect 'mutFitCost' on fitness.")
public class MutLoadFitnessModel extends MultiSiteFitnessModel {

	// Input for mutant fit cost
	public Input<RealParameter> fitCostInput =
			new Input<>("mutFitCost", "Selective cost of mutations, assumed invariant across sites");
	
	// Mutation rates
	public Input<RealParameter> mutationMatrix =
			new Input<>("mutationMatrix", "Flattened mutation matrix, can be asymmetric, diagnonal entries omitted");
	
	protected double fitCost;
	
	@Override
	public void update() {

		fitCost = fitCostInput.get().getValue();
		double forwardMutRate = mutationMatrix.get().getValues()[0];
		double backMutRate = mutationMatrix.get().getValues()[1];
		
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k][1] = 1.0 - fitCost;
		}
		
		// Update fitSpaceValues
		for (int k=0; k<fitSpaceStates; k++) {
			fitSpaceValues[k] = Math.pow(1-fitCost,k);
		}
		
		// Update fitSpaceGammas (not necessary here)
		for (int i = 0; i < fitSpaceStates; i++) {
			for (int j = 0; j < fitSpaceStates; j++) {
				if (j == (i-1)) {
					int index = i * (fitSpaceStates - 1) + j; // Indexing is broken here
					fitSpaceGammas[index] = backMutRate * i; // back mutation from i --> j
				} else if (j == (i+1)) {
					int index = i * (fitSpaceStates - 1) + (j-1); // Indexing is broken here
					fitSpaceGammas[index] = forwardMutRate * (seqSites - i); // forward mutation from i -- j
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
		
		fitCost = fitCostInput.get().getValue();
		double forwardMutRate = mutationMatrix.get().getValues()[0];
		double backMutRate = mutationMatrix.get().getValues()[1];
		
		for (int k=0; k<seqSites; k++) {
			siteFitnessEffects[k][0] = 1.0;
			siteFitnessEffects[k][1] = 1.0 - fitCost;
		}
		
		// Set up fitness space
		fitSpaceStates = seqSites + 1; // plus one for wildtype w/ no mutations
		fitSpaceValues = new Double[fitSpaceStates];
		for (int k=0; k<fitSpaceStates; k++) {
			fitSpaceValues[k] = Math.pow(1-fitCost,k);
		}
		
		// Set up fitness space gammas
		fitSpaceGammas = new Double[fitSpaceStates*(fitSpaceStates-1)];
		for (int i = 0; i < fitSpaceStates; i++) {
			for (int j = 0; j < fitSpaceStates; j++) {
				if (j == (i-1)) {
					int index = i * (fitSpaceStates - 1) + j; // Indexing is broken here
					fitSpaceGammas[index] = backMutRate * i; // back mutation from i --> j
				} else if (j == (i+1)) {
					int index = i * (fitSpaceStates - 1) + (j-1); // Indexing is broken here
					fitSpaceGammas[index] = forwardMutRate * (seqSites - i); // forward mutation from i -- j
				}
			}
		}
		
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
