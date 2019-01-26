package beast.evolution.speciation;

import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.*;
import beast.core.Input;
import beast.core.Description;
import beast.core.util.Utils;

import beast.math.SmallNumber;
import beast.math.p0ge_InitialConditions;
import beast.math.p0ge_SiteConditions;
import beast.util.HeapSort;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * @author Denise Kuehnert
 * Date: Jul 2, 2013
 * Time: 10:28:16 AM
 * 
 * @author David Rasmussen
 * modified for the marginal fitness birth-death model (Dec., 2017)
 *
 */

@Description("This class implements the marginal fitness birth-death model with muliple evolving sites. " +
		"It is based on the original multi-deme birth-death model implemented by Denise Kuehnert." +
		"This version may not be compataible with some of the original features of the MTBD model like sampled ancestor trees." +
		"This version is implemented to prevent numerical underflowing, using so-called 'SmallNumbers', with the cost of additional computational complexity")

public class FitnessBirthDeathModel extends PiecewiseFitnessBirthDeathDistribution implements Loggable {


	public Input<TraitSet> tiptypes = new Input<>("tiptypes", "trait information for initializing traits (like node types/locations) in the tree");
	public Input<String> typeLabel = new Input<>("typeLabel", "type label in tree for initializing traits (like node types/locations) in the tree");
	public Input<IntegerParameter> tipTypeArray = new Input<IntegerParameter>("tipTypeArray", "integer array of traits (like node types/locations) in the tree, index corresponds to node number in tree");
	
	public Input<Boolean> storeNodeTypes = new Input<>("storeNodeTypes", "store tip node types? this assumes that tip types cannot change (default false)", false);

	private int[] nodeStates; // MFBD model currently does not use nodeStates - but these could hold an additional state not given by the sequence data
	protected int[][] tipSeqStates; // Added to hold tip sequence data

	Boolean print = false;

	@Override
	public void initAndValidate() {

		if ((tiptypes.get()==null?0:1) + (typeLabel.get()==null?0:1) + (tipTypeArray.get()==null?0:1) != 1 )
			throw new RuntimeException("Tip types need to be specified exactly once using either tiptypes OR typeLabel OR tipTypeArray.");

		TreeInterface tree = treeInput.get();

		checkOrigin(tree);

		ntaxa = tree.getLeafNodeCount();
		
		super.initAndValidate();

		if (storeNodeTypes.get()) {

			nodeStates = new int[ntaxa];

			for (Node node : tree.getExternalNodes()){
				nodeStates[node.getNr()] = getNodeState(node, true);
			}
		}

		int contempCount = 0;
		for (Node node : tree.getExternalNodes())
			if (node.getHeight()==0.)
				contempCount++;
		
		// Store tip sequence data from alignment
		tipSeqStates = new int[ntaxa][siteCount];
		for (Node node : tree.getExternalNodes()) {
	        int taxonIndex = getTaxonIndex(node.getID(), seqData); // taxonIndex should be same as nodeNr, but just to be safe
	        for (int i = 0; i < siteCount; i++) {
	            tipSeqStates[node.getNr()][i] = seqData.getPattern(taxonIndex, i); // Ambiguous states are ignored.
	        }
	        //System.out.println(Arrays.toString(tipSeqStates[node.getNr()]));
		}

		if (checkRho.get() && contempCount>1 && rho==null)
			throw new RuntimeException("Error: multiple tips given at present, but sampling probability \'rho\' is not specified.");

		collectTimes(T);
		setRho();
	}
	
    /**
    *
    * Copied method from TreeLikelihood
    *
    * @param taxon the taxon name as a string
    * @param data the alignment
    * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
    *         or -1 if the taxon is not in the alignment.
    */
   private int getTaxonIndex(String taxon, Alignment data) {
       int taxonIndex = data.getTaxonIndex(taxon);
       if (taxonIndex == -1) {
       	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
               taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
           }
           if (taxonIndex == -1) {
           	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
           }
       }
       return taxonIndex;
	}

	void computeRhoTips(){

		double tipTime;

		for (Node tip : treeInput.get().getExternalNodes()) {

			tipTime = T-tip.getHeight();
			isRhoTip[tip.getNr()] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TO DO: make a warning that rho sampling precision is with 1e-10. Maybe do a threshold to the type of dating associated with the data?
				if (Math.abs(time-tipTime) <  globalPrecisionThreshold && rho[getNodeState(tip,false)*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

			}
		}
	}

	/**
	 * Implementation of getG with Small Number structure for ge equations. Avoids underflowing of integration results.
	 * WARNING: getG and getGSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param t
	 * @param PG0
	 * @param t0
	 * @param node
	 * @return
	 */
	public p0ge_SiteConditions getG(double t, p0ge_SiteConditions PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)


		if (node.isLeaf()) {
			
			// pInitialConditions are no longer indexed by node numbers, but by times. But these initial conditions are never used anyways, so set to time 0
			System.arraycopy(pInitialConditions[0], 0, PG0.conditionsOnP, 0, fitnessModel.getFitSpaceStates()); // copy pInitialConditions to PG0
			
		}

		//return getG(t,  PG0,  t0, pg_integrator, PG, T, maxEvalsUsed);
		return getGFastEuler(t,  PG0,  t0, pg_integrator, PG, T, maxEvalsUsed); // with fast Euler integration

	}

	/**
	 * WARNING: calculateTreeLogLikelihood allows use of both classic and non-underflowing methods. Some chunks of code are therefore present in two similar versions in this method.
	 * When modifying one of the versions, one should check if the other version also needs the corresponding changes.
	 */
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {

		Node root = tree.getRoot();

		if (origin.get()==null)
			T = root.getHeight();
		else
			updateOrigin(root);


		collectTimes(T);
		setRho();

		if ((orig < 0) || updateRates() < 0 ||  (times[totalIntervals-1] > T)) {
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		double[] noSampleExistsProp = new double[]{0.0};
		SmallNumber[] PrSN = new SmallNumber[siteCount]; // PrSN for each site
		SmallNumber full_tree_pD = new SmallNumber(0); // pD for entire tree
		for (int k=0; k<siteCount; k++) PrSN[k] = new SmallNumber();

		try{  // start calculation

			// Compute initial conditions at all tips in tree
			//pInitialConditions = getAllInitialConditionsForP(tree); // P0 conditions are now computed on a grid of times for all lineages 
			
			// Compute P0 conditions at all time points
			pInitialConditions = getConditionsForPAllTimes(tree);

			if (conditionOnSurvival.get()) {

				noSampleExistsProp = pInitialConditions[pInitialConditions.length-1];

				if (print) System.out.println("\nnoSampleExistsProp = " + noSampleExistsProp[0] + ", " + noSampleExistsProp[1]);

			}

			p0ge_SiteConditions pSN = new p0ge_SiteConditions();

			if ( orig > 0 ) { // If origin is before root(?)
					pSN = calculateSubtreeLikelihood(root,0,orig);
			} else {
				int childIndex = 0;
				if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

					pSN = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					childIndex = Math.abs(childIndex-1);

					p0ge_SiteConditions p1SN = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());
					
					for (int k=0; k<siteCount; k++) { // Now over all sites
						for (int i=0; i<seqStates[k]; i++) {
							pSN.conditionsOnG[k][i] = SmallNumber.multiply(pSN.conditionsOnG[k][i], p1SN.conditionsOnG[k][i]);
						}
					}

			}

			if (print) System.out.print("final p per state = ");
				
				// Compute final prob density at root for each site
				double[][] marginalFitEffects = fitnessModel.getMarginalFitnessEffectsSN(pSN.conditionsOnG);
				
				// Only need these if approxEProbs == false
				double[][] marginalSiteProbs = getMarginalSiteProbsSN(pSN.conditionsOnG);
				int[][] genotypes = fitnessModel.getGenotypes();
				
				for (int k=0; k<siteCount; k++) {
					for (int root_state=0; root_state<seqStates[k]; root_state++){
	
						if (pSN.conditionsOnG[k][root_state].getMantissa()>0) {
							if (!conditionOnSurvival.get()){
								PrSN[k] = SmallNumber.add(PrSN[k], pSN.conditionsOnG[k][root_state].scalarMultiply(freq[root_state]));
							} else {
								double expected_p0;
								if (approxEProbs) {
									
									int fitClass = fitnessModel.getClosestFitSpaceState(marginalFitEffects[k][root_state]);
									expected_p0 = noSampleExistsProp[fitClass];
									
								} else {
									
									expected_p0 = 0;
									double[] gpr = new double[genotypes.length];
									double gprSum = 0;
									for (int geno=0; geno<genotypes.length; geno++) {
										if (genotypes[geno][k] == root_state) { // otherwise no prob of being in this genotype
											gpr[geno] = 1.0;
											for (int j=0; j<siteCount; j++) {
												if (k != j) { // already accounted for this site above
													int type = genotypes[geno][j];
													gpr[geno] *= marginalSiteProbs[j][type];
												}
											}
											gprSum += gpr[geno];
										}
									}
									
									if (gprSum <= 0) {
										for (int geno=0; geno<genotypes.length; geno++) {
											expected_p0 += noSampleExistsProp[geno] / genotypes.length;
										}
									} else {
										for (int geno=0; geno<genotypes.length; geno++) {
											expected_p0 += noSampleExistsProp[geno] * gpr[geno] / gprSum;
										}
									}
									
								}
								PrSN[k] = SmallNumber.add(PrSN[k], pSN.conditionsOnG[k][root_state].scalarMultiply(freq[root_state]).scalarMultiply(1/(1-expected_p0)));
							}
						}
						//if (print) System.out.print(pSN.conditionsOnP[root_state] + "\t" + pSN.conditionsOnG[root_state] + "\t");
					}
				}
				full_tree_pD = pSN.conditionsOnG[siteCount][0]; // should full_tree_pD be conditioned on survival as well?


		}catch(Exception e){

			if (e instanceof ConstraintViolatedException){throw e;}

			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		maxEvalsUsed = Math.max(maxEvalsUsed, PG.maxEvalsUsed);

		// Need take product over all sites and divide by total tree density
		logP = 0;
		for (int k=0; k<siteCount; k++) {
			logP += PrSN[k].log();
		}
		logP -= (siteCount-1) * full_tree_pD.log(); // 'divide' by full prob of tree

		if (print) System.out.println("\nlogP = " + logP);

		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		if (SAModel && !(removalProbability.get().getDimension()==n && removalProbability.get().getValue()==1.)) {
			int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
			logP +=  Math.log(2)*internalNodeCount;
		}
		return logP;
	}

	private int getNodeState(Node node, Boolean init){

		try {

			if (!storeNodeTypes.get() || init){

				int nodestate = tiptypes.get() != null ?
						(int) tiptypes.get().getValue((node.getID())) :
						typeLabel.get()!=null ?
								((node instanceof MultiTypeNode) ? ((MultiTypeNode) node).getNodeType() : -2) :
								(int) tipTypeArray.get().getArrayValue((node.getNr()));

				if (nodestate == -2) {
					Object d = node.getMetaData(typeLabel.get());

					if (d instanceof Integer) nodestate = (Integer) node.getMetaData(typeLabel.get());
					else if
							(d instanceof Double) nodestate = (((Double) node.getMetaData(typeLabel.get())).intValue());
					else if
							(d instanceof int[]) nodestate = (((int[]) node.getMetaData(typeLabel.get()))[0]);
				}

				return nodestate;

			}
			else return nodeStates[node.getNr()];

		}catch(Exception e){
			throw new ConstraintViolatedException("Something went wrong with the assignment of types to the nodes (node ID="+node.getID()+"). Please check your XML file!");
		}
	}

	/**
	 * Implementation of calculateSubtreeLikelihood with Small Number structure. Avoids underflowing of integration results.
	 * @param node
	 * @param from
	 * @param to
	 * @return
	 */
	p0ge_SiteConditions calculateSubtreeLikelihood(Node node, double from, double to) {

		double[] pconditions = new double[fitnessModel.getFitSpaceStates()]; // E prob of no sampled descendants
		SmallNumber[][] gconditions = new SmallNumber[siteCount+1][]; // now holds gconditions for all sites
		for (int k=0; k<siteCount; k++) {
			gconditions[k] = new SmallNumber[seqStates[k]];
			for (int i=0; i<seqStates[k]; i++) {
				gconditions[k][i] = new SmallNumber();
			}
		}
		gconditions[siteCount] = new SmallNumber[1]; // final 'site' holds gconditions for full tree
		gconditions[siteCount][0] = new SmallNumber();
		p0ge_SiteConditions init = new p0ge_SiteConditions(pconditions, gconditions); // with gconditions for all sites

		int index = Utils.index(to,times, totalIntervals); // time interval index for params with change points

		if (node.isLeaf()){ // sampling event
			
			int nodestate;

			if (!isRhoTip[node.getNr()]) {
	
				for (int k=0; k<siteCount; k++) {
					nodestate = tipSeqStates[node.getNr()][k];
					if (samplingCoupledToRemoval) {
						init.conditionsOnG[k][nodestate] = new SmallNumber(death[index] * psi[index]); // removal with sampling
					} else {
						init.conditionsOnG[k][nodestate] = SAModel?
							new SmallNumber((r[index] + pInitialConditions[node.getNr()][nodestate]*(1-r[index])) // Handling of sampling probs under SAModel unchecked here and probably not correct
									*psi[index]) // with SA: ψ_i(r + (1 − r)p_i(τ))
							: new SmallNumber(psi[index]); // psi is sampling rate
					}
				}
				
				// Need to include sampling probs for full_tree_pD
				if (samplingCoupledToRemoval) {
					init.conditionsOnG[siteCount][0] = new SmallNumber(death[index] * psi[index]); // removal with sampling
				} else {
					init.conditionsOnG[siteCount][0] = new SmallNumber(psi[index]); // psi is sampling rate
				}

			}	else {
				
				for (int k=0; k<siteCount; k++) {
					nodestate = tipSeqStates[node.getNr()][k];
					init.conditionsOnG[k][nodestate] = SAModel? 
							new SmallNumber((r[index] + pInitialConditions[node.getNr()][nodestate]/(1-rho[index])*(1-r[index]))
									*rho[index])  :
										new SmallNumber(rho[index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ+δ)) //the +δ is translated by dividing p_i with 1-ρ_i (otherwise there's one too many "*ρ_i" )
					
				}
				
				// Need to include sampling probs for full_tree_pD
				init.conditionsOnG[siteCount][0] = new SmallNumber(rho[index]); // rho is the sampling fraction
				
			}

			if (print) System.out.println("Sampling at time " + (T-to));

			return getG(from, init, to, node);
		}


		else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

			if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor
				
				throw new ConstraintViolatedException("Error: Sampled ancestors not yet implemented in multi-site fitness models!");

			} else {   // birth / infection event

				int childIndex = 0;
				if (node.getChild(1).getNr() > node.getChild(0).getNr())
					childIndex = 1; // always start with the same child to avoid numerical differences

				p0ge_SiteConditions g0 = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				childIndex = Math.abs(childIndex - 1);

				p0ge_SiteConditions g1 = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				if (print)
					System.out.println("Infection at time " + (T - to));//+ " with p = " + p + "\tg0 = " + g0 + "\tg1 = " + g1);

				// Need to compute marginal fitness for each child lineage from SmallNumbers
				double[][] marginalFitEffectsChild0 = fitnessModel.getMarginalFitnessEffectsSN(g0.conditionsOnG);
				double[][] marginalFitEffectsChild1 = fitnessModel.getMarginalFitnessEffectsSN(g1.conditionsOnG);
				
				// Don't forget to copy over init.conditionsOnP;
				System.arraycopy(g0.conditionsOnP, 0, init.conditionsOnP, 0, fitnessModel.getFitSpaceStates()); // no longer necessary
				
				for (int k=0; k<siteCount; k++) {	
					for (int i=0; i<seqStates[k]; i++) {
						final double child0Lambda = birth[index] * marginalFitEffectsChild0[k][i];
						final double child1Lambda = birth[index] * marginalFitEffectsChild1[k][i];
						init.conditionsOnG[k][i] = SmallNumber.multiply(g0.conditionsOnG[k][i], g1.conditionsOnG[k][i]).scalarMultiply(child0Lambda + child1Lambda); // Mutation can't occur at birth event so both child lineages need to be in state i BUT either child could have been parent
					}
				}
				
		        // Update G for full_tree_pD
		        final double child0Lambda = birth[index] * marginalFitEffectsChild0[siteCount][0];
		        final double child1Lambda = birth[index] * marginalFitEffectsChild1[siteCount][0];
		        init.conditionsOnG[siteCount][0] = SmallNumber.multiply(g0.conditionsOnG[siteCount][0], g1.conditionsOnG[siteCount][0]).scalarMultiply(child0Lambda + child1Lambda);
				
				//if (init.isSumGZero(siteCount)) {
					//System.out.println("Prob density G evaluated as zero at a branching event at time : " + to + " for node: " + node.getNr());
				//}
				//System.out.println(node.getNr());
				
			}
		}


		else {// found a single child node

			throw new RuntimeException("Error: Single child-nodes found (although not using sampled ancestors)");
		}

		if (print){
			System.out.print("p after subtree merge = ");
			for (int i=0;i<n;i++) System.out.print(init.conditionsOnP[i] + "\t");
			for (int i=0;i<n;i++) System.out.print(init.conditionsOnG[i] + "\t");
			System.out.println();
		}

		return getG(from, init, to, node);
	}

	// used to indicate that the state assignment went wrong
	protected class ConstraintViolatedException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public ConstraintViolatedException(String s) {
			super(s);
		}

	}

	@Override
	public void init(PrintStream out){

		super.init(out);

		if (tipTypeArray.get()!=null) {
			IntegerParameter types = tipTypeArray.get();
			for (int i = 0; i < types.getDimension(); i++) {
				out.print(tipTypeArray.get().getID() + ("_node") + (i+1) + "\t");
			}

		}

	}

	@Override
	public void log(int sampleNr, PrintStream out) {

		super.log(sampleNr, out);

		if (tipTypeArray.get()!=null) {
			for (int i = 0; i < tipTypeArray.get().getDimension(); i++) {
				out.print(treeInput.get().getNode(i).getID() + "\t");

			}
		}
	}


}

