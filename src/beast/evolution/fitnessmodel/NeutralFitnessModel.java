package beast.evolution.fitnessmodel;

import java.io.PrintStream;

import beast.core.Description;

/**
 * @author David Rasmussen
 */
@Description("Multi-site fitness model" +
             "where all sites are neutral and have no effects on fitness.")

public class NeutralFitnessModel extends MultiSiteFitnessModel {
	
	@Override
	public void update() {
		
		// Nothing to update
		
	}
	
	@Override
	public void setUp(Integer[] states, int sites) {
		
		int maxState = states[0]; // assumes all sites have same # of states
		seqStates = states;
		seqSites = sites;
		siteFitnessEffects = new double[seqSites][maxState];
		
		for (int k=0; k<seqSites; k++) {
			for (int i=0; i<maxState; i++) {
				siteFitnessEffects[k][i] = 1.0; // all sites are neutral and have fitness = 1.0 under this model
			}
		}
		
		// Set up fitness space -- only one state with fitness 1.0
		fitSpaceStates = 1;
		fitSpaceValues = new Double[1];
		fitSpaceValues[0] = 1.0;
		fitSpaceGammas = new Double[1];
		fitSpaceGammas[0] = 0.0; // no transitions
		
	}
	
	@Override
	public int getClosestFitSpaceState(double fit) {
		return 0; // always zero for this model
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
