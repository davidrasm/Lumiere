package beast.evolution.speciation;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Utils;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.UnsortedAlignment;
import beast.evolution.alignment.UnsortedFilteredAlignment;
import beast.evolution.fitnessmodel.MultiSiteFitnessModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.math.MultiSiteScaledNumbers;
import beast.math.ScaledNumbers;
import beast.math.SmallNumber;
import beast.math.SmallNumberScaler;
import beast.math.p0_ODE;
import beast.math.p0ge_InitialConditions;
import beast.math.p0ge_ODE;
import beast.math.p0ge_SiteConditions;
import beast.util.HeapSort;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created with IntelliJ IDEA.
 * User: Denise
 * Date: 22.08.14
 * Time: 14:05
 * 
 * @author David Rasmussen
 * modified for the marginal fitness birth-death model (Dec., 2017)
 * 
 */
@Citation("Rasmussen D.A. & Stadler T. 2019. " +
		"Coupling adaptive molecular evolution to phylodynamics using fitness-dependent birth-death models" +
		"Add ref.")

@Description("Piece-wise constant rates are assumed to be ordered by state and time. First k entries of an array give " +
		"values belonging to type 1, for intervals 1 to k, second k intervals for type 2 etc.")
public abstract class PiecewiseFitnessBirthDeathDistribution extends SpeciesTreeDistribution {


	public Input<RealParameter> frequencies =
			new Input<>("frequencies", "The frequencies for each type",  Input.Validate.REQUIRED);

	public Input<RealParameter> origin =
			new Input<>("origin", "The origin of infection x1");

	public Input<Boolean> originIsRootEdge =
			new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);

	public Input<Integer> maxEvaluations =
			new Input<>("maxEvaluations", "The maximum number of evaluations for ODE solver", 1000000);

	public Input<Boolean> conditionOnSurvival =
			new Input<>("conditionOnSurvival", "condition on at least one survival? Default true.", true);

	public Input<Double> relativeTolerance =
			new Input<>("relTolerance", "relative tolerance for numerical integration", 1e-7);

	public Input<Double> absoluteTolerance =
			new Input<>("absTolerance", "absolute tolerance for numerical integration", 1e-100 /*Double.MIN_VALUE*/);

	// the interval times for the migration rates
	public Input<RealParameter> migChangeTimesInput =
			new Input<>("migChangeTimes", "The times t_i specifying when migration rate changes occur", (RealParameter) null);

	// the interval times for the birth rate
	public Input<RealParameter> birthRateChangeTimesInput =
			new Input<>("birthRateChangeTimes", "The times t_i specifying when birth/R rate changes occur", (RealParameter) null);

	// the interval times for the birth rate among demes
	public Input<RealParameter> b_ijChangeTimesInput =
			new Input<>("birthRateAmongDemesChangeTimes", "The times t_i specifying when birth/R among demes changes occur", (RealParameter) null);

	// the interval times for the death rate
	public Input<RealParameter> deathRateChangeTimesInput =
			new Input<>("deathRateChangeTimes", "The times t_i specifying when death/becomeUninfectious rate changes occur", (RealParameter) null);

	// the interval times for sampling rate
	public Input<RealParameter> samplingRateChangeTimesInput =
			new Input<>("samplingRateChangeTimes", "The times t_i specifying when sampling rate or sampling proportion changes occur", (RealParameter) null);

	// the interval times for removal probability
	public Input<RealParameter> removalProbabilityChangeTimesInput =
			new Input<RealParameter>("removalProbabilityChangeTimes", "The times t_i specifying when removal probability changes occur", (RealParameter) null);

	public Input<RealParameter> intervalTimes =
			new Input<>("intervalTimes", "The time t_i for all parameters if they are the same", (RealParameter) null);

	public Input<Boolean> migTimesRelativeInput =
			new Input<>("migTimesRelative", "True if migration rate change times specified relative to tree height? Default false", false);

	public Input<Boolean> b_ijChangeTimesRelativeInput =
			new Input<>("birthRateAmongDemesTimesRelative", "True if birth rate change times specified relative to tree height? Default false", false);

	public Input<Boolean> birthRateChangeTimesRelativeInput =
			new Input<>("birthRateTimesRelative", "True if birth rate change times specified relative to tree height? Default false", false);

	public Input<Boolean> deathRateChangeTimesRelativeInput =
			new Input<>("deathRateTimesRelative", "True if death rate change times specified relative to tree height? Default false", false);

	public Input<Boolean> samplingRateChangeTimesRelativeInput =
			new Input<>("samplingRateTimesRelative", "True if sampling rate times specified relative to tree height? Default false", false);

	Input<Boolean> removalProbabilityChangeTimesRelativeInput =
			new Input<Boolean>("removalProbabilityTimesRelative", "True if removal probability change times specified relative to tree height? Default false", false);

	public Input<BooleanParameter> reverseTimeArraysInput =
			new Input<>("reverseTimeArrays", "True if the time arrays are given in backwards time (from the present back to root). Order: 1) birth 2) death 3) sampling 4) rho 5) r 6) migration. Default false." +
					"Careful, rate array must still be given in FORWARD time (root to tips).");

	// the times for rho sampling
	public Input<RealParameter> rhoSamplingTimes =
			new Input<>("rhoSamplingTimes", "The times t_i specifying when rho-sampling occurs", (RealParameter) null);
	public Input<Boolean> contemp =
			new Input<>("contemp", "Only contemporaneous sampling (i.e. all tips are from same sampling time, default false)", false);

	public Input<RealParameter> birthRate =
			new Input<>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time");
	public Input<RealParameter> deathRate =
			new Input<>("deathRate", "The deathRate vector with birthRates between times");
	public Input<RealParameter> samplingRate =
			new Input<>("samplingRate", "The sampling rate per individual");      // psi

	public Input<RealParameter> m_rho =
			new Input<>("rho", "The proportion of lineages sampled at rho-sampling times (default 0.)");


	public Input<RealParameter> R0 =
			new Input<>("R0", "The basic reproduction number");
	public Input<RealParameter> becomeUninfectiousRate =
			new Input<>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (through recovery or sampling)", Input.Validate.XOR, deathRate);
	public Input<RealParameter> samplingProportion =
			new Input<>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRate);

	public Input<BooleanParameter> identicalRatesForAllTypesInput =
			new Input<>("identicalRatesForAllTypes", "True if all types should have the same 1) birth 2) death 3) sampling 4) rho 5) r 6) migration rate. Default false.");

	public Input<RealParameter> R0_base =
			new Input<>("R0_base",
					"The basic reproduction number for the base pathogen class, should have the same dimension as " +
					"the number of time intervals.");
	public Input<RealParameter> lambda_ratio =
			new Input<>("lambda_ratio",
					"The ratio of basic infection rates of all other classes when compared to the base lambda, " +
					"should have the dimension of the number of pathogens - 1, as it is kept constant over intervals.");

	public Input<RealParameter> migrationMatrix =
			new Input<>("migrationMatrix", "Flattened migration matrix, can be asymmetric, diagnonal entries omitted");
	public Input<RealParameter> migrationMatrixScaleFactor =
			new Input<>("migrationMatrixScaleFactor", "A real number with which each migration rate entry is scaled.");


	public Input<RealParameter> birthRateAmongDemes =
			new Input<>("birthRateAmongDemes", "birth rate vector with rate at which transmissions occur among locations");

	public Input<RealParameter> R0AmongDemes =
			new Input<>("R0AmongDemes", "The basic reproduction number determining transmissions occur among locations");


	public Input<RealParameter> removalProbability =
			new Input<RealParameter>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");


	public Input<Integer> stateNumber =
			new Input<>("stateNumber", "The number of states or locations", Input.Validate.REQUIRED);

	public Input<RealParameter> adjustTimesInput =
			new Input<>("adjustTimes", "Origin of MASTER sims which has to be deducted from the change time arrays");
	// <!-- HACK ALERT for reestimation from MASTER sims: adjustTimes is used to correct the forward changetimes such that they don't include orig-root (when we're not estimating the origin) -->

	public Input<Boolean> useRKInput =
			new Input<>("useRK", "Use fixed step size Runge-Kutta integrator with 1000 steps. Default false", false);

	public Input<Boolean> checkRho = new Input<>("checkRho", "check if rho is set if multiple tips are given at present (default true)", true);
	
	public Input<Double> dtTimeStepInput =
			new Input<>("dtTimeStep", "dt time step to use for numerical integration of likelihood", 1e-6);
	
	public Input<Boolean> approxEProbsInput = new Input<>("approxEProbs", "approximate E (p0) probs in fitness space (default true)", true);

	//  TO DO CHECKER QUE C'EST PAS POSSIBLE DE REMETTRE 1e-20
	public final static double globalPrecisionThreshold = 1e-10;

	double T;
	double orig;
	int ntaxa;

	p0_ODE P;
	p0ge_ODE PG;

	FirstOrderIntegrator pg_integrator;
	public int maxEvalsUsed;
	public Double minstep;
	public Double maxstep;

	// these four arrays are totalIntervals in length
	protected Double[] birth;
	Double[] death;
	Double[] psi;
	Double[] rho;
	Double[] r;

	/**
	 * The number of change points in the birth rate, b_ij, death rate, sampling rate, rho, r
	 */
	int migChanges;
	int birthChanges;
	int b_ij_Changes;
	int deathChanges;
	int samplingChanges;
	int rhoChanges;
	int rChanges;


	public boolean SAModel;

	/**
	 * The number of times rho-sampling occurs
	 */
	int rhoSamplingCount;
	Boolean constantRho;
	Boolean[] isRhoTip;
	Boolean[] isRhoInternalNode;

	/**
	 * Total interval count
	 */
	int totalIntervals;
	int n;  // number of states / locations

	protected List<Double> migChangeTimes = new ArrayList<>();
	protected List<Double> birthRateChangeTimes = new ArrayList<>();
	protected List<Double> b_ijChangeTimes = new ArrayList<>();
	protected List<Double> deathRateChangeTimes = new ArrayList<>();
	protected List<Double> samplingRateChangeTimes = new ArrayList<>();
	protected List<Double> rhoSamplingChangeTimes = new ArrayList<>();
	protected List<Double> rChangeTimes = new ArrayList<Double>();

	Boolean contempData;
	SortedSet<Double> timesSet = new TreeSet<>();

	protected Double[] times = new Double[]{0.};

	protected Boolean transform;

	Boolean migTimesRelative = false;
	Boolean birthRateTimesRelative = false;
	Boolean b_ijTimesRelative = false;
	Boolean deathRateTimesRelative = false;
	Boolean samplingRateTimesRelative = false;
	Boolean rTimesRelative = false;
	Boolean[] reverseTimeArrays;

	Double[] M;
	Double[] b_ij;
	Boolean birthAmongDemes = false;

	Double[] freq;

	double[][] pInitialConditions;
	double[] sortedNodes;
	
	// For tracking p0 as a time series
	Double[] integrationTimes;
	double[][] pConditions;
	
	// Fitness model specific params
	public Input<MultiSiteFitnessModel> fitnessModelInput =
			new Input<>("fitnessModel", "The multi-site fitness model",  Input.Validate.REQUIRED);
	public MultiSiteFitnessModel fitnessModel;
	public Double dtTimeStep; // added for fast Euler integration
	public Boolean samplingCoupledToRemoval = false; // is sampling coupled to removal?
	public Boolean approxEProbs = true;
	
	// Alignment input and params
	public Input<UnsortedAlignment> alignInput = new Input<>("alignment", "sequence alignment");
	public UnsortedAlignment seqData;
	public int siteCount;
	public Integer[] seqStates; // can now be variable across sites
	public int maxState; 

	@Override
	public void initAndValidate() {
		
		seqData = alignInput.get(); // Get seq data from align for each tip
		
		siteCount = seqData.getPatternCount(); // Number of evolving sites
		maxState = seqData.getMaxDataState(); // Maximum # of states at each site
		seqStates = seqData.getMaxDataStates(); // # of states at each site
		
		// Set up fitness model
		fitnessModel = fitnessModelInput.get();
		fitnessModel.setUp(seqStates, siteCount);
		
		// Set dt time step for numerical integration
		dtTimeStep = dtTimeStepInput.get();

		identicalRatesForAllTypes = new Boolean[]{false, false, false, false, false, false};
		if (identicalRatesForAllTypesInput.get()!=null)
			identicalRatesForAllTypes = identicalRatesForAllTypesInput.get().getValues();

		if (removalProbability.get() != null) SAModel = true;
		
		if (approxEProbsInput.get()!=null) {
			approxEProbs = approxEProbsInput.get();
		}

		birth = null;
		b_ij = null;
		death = null;
		psi = null;
		rho = null;
		r = null;
		birthRateChangeTimes.clear();
		deathRateChangeTimes.clear();
		samplingRateChangeTimes.clear();
		if (SAModel) rChangeTimes.clear();
		totalIntervals = 0;
		n = stateNumber.get();

		birthAmongDemes = (birthRateAmongDemes.get() !=null || R0AmongDemes.get()!=null);

		Double factor;
		if (migrationMatrix.get()!=null) {
			M = migrationMatrix.get().getValues();

			if (migrationMatrixScaleFactor.get()!=null) {
				factor = migrationMatrixScaleFactor.get().getValue();
				for (int i = 0; i < M.length; i++) M[i] *= factor;
			}

			if (n>1 && M.length != n*(n-1)) {
				double timeChanges = 0;
				if (migChangeTimesInput.get()!=null) {
					timeChanges = migChangeTimesInput.get().getDimension();
				} else if(intervalTimes.get() != null){
					timeChanges = intervalTimes.get().getDimension();
				}
				if (timeChanges == 0 || M.length != n*(n-1)*timeChanges ) 
					throw new RuntimeException("Migration matrix dimension is incorrect!");
			}
			//migChanges = migrationMatrix.get().getDimension()/Math.max(1,(n*(n-1))) - 1;  // <-- migrationMatrix has different dimension than n = 1
			migChanges = migrationMatrix.get().getDimension()/Math.max(1,(maxState*(maxState-1))) - 1;  // <-- migrationMatrix has dimension n = seqStates

		}


		else if (!birthAmongDemes) throw new RuntimeException("Error in BDMM setup: need to specify at least one of the following: migrationMatrix, R0AmongDemes, birthRateAmongDemes");

		birthRateTimesRelative = birthRateChangeTimesRelativeInput.get();
		b_ijTimesRelative = b_ijChangeTimesRelativeInput.get();
		migTimesRelative = migTimesRelativeInput.get();
		deathRateTimesRelative = deathRateChangeTimesRelativeInput.get();
		samplingRateTimesRelative = samplingRateChangeTimesRelativeInput.get();
		if (SAModel) rTimesRelative = removalProbabilityChangeTimesRelativeInput.get();

		reverseTimeArrays = new Boolean[]{false, false, false, false, false, false};
		if (reverseTimeArraysInput.get()!= null )  {
			Boolean[] r = reverseTimeArraysInput.get().getValues();
			for (int i=0; i<r.length; i++)
				reverseTimeArrays[i] = r[i];
		}

		rhoSamplingCount = 0;
		contempData = contemp.get();

		if (birthRate.get() == null && R0.get() == null && R0_base.get() == null && lambda_ratio.get() == null) {
			throw new RuntimeException("Either birthRate, R0, or R0_base and R0_ratio need to be specified!");
		} else if ((birthRate.get() != null && R0.get() != null)
				|| (R0.get() != null && (R0_base.get() != null || lambda_ratio.get() != null))
				|| (birthRate.get() != null && (R0_base.get() != null || lambda_ratio.get() != null))) {
			throw new RuntimeException("Only one of birthRate, or R0, or R0_base and lambda_ratio need to be specified!");
		} else if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null) {
			
			// For condition where birthRate, deathRate and samplingProportion is defined
			// In this case psi = samplingRate
			
			// Using untransformed parameterization
			transform = false;
			death = deathRate.get().getValues();
			psi = samplingRate.get().getValues();
			birth = birthRate.get().getValues();
			if (SAModel) r = removalProbability.get().getValues();

			if (birthRateAmongDemes.get()!=null ){

				birthAmongDemes = true;
				b_ij=birthRateAmongDemes.get().getValues();
			}
		} else if (birthRate.get() != null && deathRate.get() != null && samplingProportion.get() != null) {
			
			// For condition where birthRate, deathRate and samplingProportion is defined
			// In this case psi = samplingProportion
			
			// Using untransformed parameterization
			transform = false;
			samplingCoupledToRemoval = true; // sampling must be coupled to removal under this parameterization
			death = deathRate.get().getValues();
			psi = samplingProportion.get().getValues();
			birth = birthRate.get().getValues();
			if (SAModel) r = removalProbability.get().getValues();

			if (birthRateAmongDemes.get()!=null ){

				birthAmongDemes = true;
				b_ij=birthRateAmongDemes.get().getValues();
			}
		} else if ((R0.get() != null || (R0_base.get() != null && lambda_ratio.get() != null)) && becomeUninfectiousRate.get() != null && samplingProportion.get() != null) {
			transform = true;
		} else {
			throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0 (or R0_base AND R0_ratio), becomeUninfectiousRate and samplingProportion!");
		}

		if (transform) {

			if (R0AmongDemes.get()!=null) {
				birthAmongDemes = true;
				b_ij_Changes = R0AmongDemes.get().getDimension()/Math.max(1,(n*(n-1))) - 1;
			}

			if (birthChanges < 1) {
				if (R0.get()!=null) {
					birthChanges = R0.get().getDimension() / n - 1;
				} else {
					birthChanges = R0_base.get().getDimension() - 1;
				}
			}
			samplingChanges = samplingProportion.get().getDimension()/n - 1;
			deathChanges = becomeUninfectiousRate.get().getDimension()/n - 1;

		} else {    //todo: b d s param doesn't work yet with rate changes (unless all parameters have equally many)

			if (birthChanges < 1) birthChanges = birthRate.get().getDimension()/n - 1;
			if (birthAmongDemes) b_ij_Changes = birthRateAmongDemes.get().getDimension()/(n*(n-1)) - 1;
			deathChanges = deathRate.get().getDimension()/n - 1;
			if (samplingCoupledToRemoval) {
				samplingChanges = samplingProportion.get().getDimension()/n - 1;
			} else {
				samplingChanges = samplingRate.get().getDimension()/n - 1;
			}
		}


		if (SAModel) rChanges = removalProbability.get().getDimension()/n -1;

		if (m_rho.get()!=null) {
			rho = m_rho.get().getValues();
			rhoChanges = m_rho.get().getDimension()/n - 1;
		}

		freq = frequencies.get().getValues();

		double freqSum = 0;
		for (double f : freq) freqSum+= f;
		if (freqSum!=1.)
			throw new RuntimeException("Error: frequencies must add up to 1 but currently add to " + freqSum + ".");

	}
	
		
	public double[][] fastEulerIntegrateGxP0(double[][] g, double to, double from, int index) {

		// Solve ge equations - these are the D_{N} probs in the notation of Stadler and Bonhoeffer, 2013
		// Assumes entire time series of p0 is known - these are the E probs in the notation of Stadler and Bonhoeffer, 2013
		// Signs (+/-) on derivatives are reversed
		
		int k, l;
		
		//final double dt = from - to; // should be negative to reverse flows
		//if (dt > 0) {
			//System.out.println("WARNING: Found non-positive integration time step!!!");
		//}
		
		double[][] gDot = new double[siteCount+1][]; // holds derivatives
		
		// Only need these if approxEProbs = false
		//if (approxEProbs == false) {
			int[][] genotypes = fitnessModel.getGenotypes(); 
			double[][] marginalSiteProbs = getMarginalSiteProbs(g);	
		//}
		
		// Notes on integration times
		// We are always integrating starting at 'to'/'nextTo' back to 'from'
		// dt will always be negative	
			
		double nextTo = to - dtTimeStep;
		double dt = nextTo - to; // should be negative
		to = nextTo; // do we need this other 'to'?
		boolean stepBack = true;
		while (stepBack) {
		
			if (nextTo < from) {
				dt = nextTo - from;
				stepBack = false;
			}
			
			// Get p0 conditions for this time point
			int p0Index = Utils.index(to,integrationTimes, integrationTimes.length);
			// fixed integration p0 prob, previously was: double[] p0 = pConditions[p0Index]
			double[] p0 = pConditions[integrationTimes.length - p0Index - 1]; 
		
			// Check dt updates
			//System.out.println("dt = " + dt);
			//System.out.println(nextTo);
			
			// Update marginal fit effects:
			double[][] marginalFitEffects = fitnessModel.getMarginalFitnessEffects(g);
			if (approxEProbs == false) {
				marginalSiteProbs = getMarginalSiteProbs(g); // Only need this if approxEProbs == false
			}
			
			for (int site=0; site<=siteCount; site++) {
				
				if (site == siteCount) { // final site represents entire lineage for computing p(T|S,\theta)
					
					gDot[site] = new double[1]; // only one "state"
				
					k = index;
					double lambda = birth[k] * marginalFitEffects[siteCount][0]; // lambda = expected birth rate based on marginal fitness
					
					// Compute expected p0 for entire lineage
					double expected_p0;
					if (approxEProbs) { // default is true
					
						int fitClass = fitnessModel.getClosestFitSpaceState(marginalFitEffects[site][0]);
						expected_p0 = p0[fitClass];
						
					} else {
						
						expected_p0 = 0;
						double[] gpr = new double[genotypes.length];
						double gprSum = 0;
						for (int geno=0; geno<genotypes.length; geno++) {
							gpr[geno] = 1.0;
							for (int j=0; j<siteCount; j++) {
								int type = genotypes[geno][j];
								gpr[geno] *= marginalSiteProbs[j][type];
							}
							gprSum += gpr[geno];
						}
						
						if (gprSum <= 0) {
							for (int geno=0; geno<genotypes.length; geno++) {
								expected_p0 += p0[geno] / genotypes.length;
							}
						} else {
							for (int geno=0; geno<genotypes.length; geno++) {
								expected_p0 += p0[geno] * gpr[geno] / gprSum;
							}
						}
						
					}
					
					if (samplingCoupledToRemoval) {
						gDot[site][0] = + (lambda+death[k]
								- 2*lambda*expected_p0)*g[site][0];
					} else {
						gDot[site][0] = + (lambda+death[k]+psi[k]
								- 2*lambda*expected_p0)*g[site][0];
					}
					
					g[site][0] += gDot[site][0] * dt; // Update g according to derivatives in gDot
					
				} else { // for a given site
					
					gDot[site] = new double[seqStates[site]];
				
					// Loop over all seqState
					for (int i=0; i<seqStates[site]; i++){
						
						//k = i*totalIntervals + index; <-- should just be k = index since i iterates over seq states
						k = index;
						
						double lambda = birth[k] * marginalFitEffects[site][i]; // lambda = expected birth rate based on marginal fitness
						
						double expected_p0;
						if (approxEProbs) { // default is true so old code will run as before
						
							int fitClass = fitnessModel.getClosestFitSpaceState(marginalFitEffects[site][i]);
							expected_p0 = p0[fitClass];
							
						} else {
							
							expected_p0 = 0;
							double[] gpr = new double[genotypes.length];
							double gprSum = 0;
							for (int geno=0; geno<genotypes.length; geno++) {
								if (genotypes[geno][site] == i) { // otherwise no prob of being in this genotype
									gpr[geno] = 1.0;
									for (int j=0; j<siteCount; j++) {
										if (site != j) { // already accounted for this site above
											int type = genotypes[geno][j];
											gpr[geno] *= marginalSiteProbs[j][type];
										}
									}
									gprSum += gpr[geno];
								}
							}
							
							if (gprSum <= 0) {
								for (int geno=0; geno<genotypes.length; geno++) {
									expected_p0 += p0[geno] / genotypes.length;
								}
							} else {
								for (int geno=0; geno<genotypes.length; geno++) {
									expected_p0 += p0[geno] * gpr[geno] / gprSum;
								}
							}
							
						}
						
						if (samplingCoupledToRemoval) {
							gDot[site][i] = + (lambda+death[k]
									- 2*lambda*expected_p0)*g[site][i];
						} else {
							gDot[site][i] = + (lambda+death[k]+psi[k]
									- 2*lambda*expected_p0)*g[site][i];
						}
			
						// Loop over all other states j
						for (int j=0; j<seqStates[site]; j++){
			
							l = (i*(maxState-1)+(j<i?j:j-1))*totalIntervals + index;
							//l = (i*(n-1)+(j<i?j:j-1))*totalIntervals + index;
							//l = i*(seqStates[site]-1)+(j<i?j:j-1); // linear indexing 
							//l = i*(maxState-1)+(j<i?j:j-1); // made this max state so it works regardless of seqStates[site]
							
							if (i!=j){
							
								//No birthAmongDemes in multi-site fitness model
								//if (b_ij!=null){     // infection among demes
									//gDot[n+i] += b_ij[l]*g[n+i]; // prob of no transmission from i to j
									//if (!augmented) {
										//gDot[n+i] -= b_ij[l]* ( g[i]*g[n+j] + g[j]*g[n+i]);
									//}
								//}
								
								if (M[0]!=null) { // mutation
									gDot[site][i] += M[l] * g[site][i]; // prob of no mutation from i to j
									gDot[site][i] -= M[l] * g[site][j]; // mutation from i to j
								}
								
							}
						}
					}
				
					// Update g according to derivatives in gDot
					for (int i=0; i<seqStates[site]; i++) {
						g[site][i] += gDot[site][i] * dt;
					}
				
				}
				
			}
			
			// Update trackers
			nextTo = to - dtTimeStep;
			dt = nextTo - to;
			to = nextTo;
	
		}
		
		return g;
	}
	
	public double[][] getMarginalSiteProbs(double[][] g) {
		
		// Get (normalized) marginal site probabilities for all sites
		double[][] norm_g = new double[siteCount][]; // g normalized by sum at each site
		for (int site=0; site<siteCount; site++) {
			norm_g[site] = new double[seqStates[site]];			
			double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum += g[site][state];
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] = g[site][state] / siteSum;
			}
		}
		
		return norm_g;
	}
	
	public double[][] getMarginalSiteProbsSN(SmallNumber[][] g) {
		
		// Same as getMarginalSiteProbs but works on SmallNumber input
		// SmallNumbers are only reverted back to doubles after normalizing marginal site densities
		
		double[][] norm_g = new double[siteCount][]; // g normalized by sum at each site
		for (int site=0; site<siteCount; site++) {
			norm_g[site] = new double[seqStates[site]];	
			SmallNumber siteSum = new SmallNumber(); //double siteSum = 0;
			for (int state=0; state<seqStates[site]; state++) siteSum = SmallNumber.add(siteSum, g[site][state]); //siteSum += g[site][state];
			for (int state=0; state<seqStates[site]; state++) {
				norm_g[site][state] =  SmallNumber.divide(g[site][state], siteSum).revert(); // only revert once normalized
			}
		}

		return norm_g;
	}
	
	/**
	 * Fast implementation of getG using Euler integration
	 * Implementation of getG with Small Number structure for the ge equations. Avoids underflowing of integration results.
	 * WARNING: getG and getGSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param t
	 * @param PG0
	 * @param t0
	 * @return
	 */
	public p0ge_SiteConditions getGFastEuler(double t, p0ge_SiteConditions PG0, double t0,
			FirstOrderIntegrator pg_integrator, p0ge_ODE PG, Double T, int maxEvalsUsed){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)


		try {

			if (Math.abs(T-t) < globalPrecisionThreshold|| Math.abs(t0-t) < globalPrecisionThreshold ||  T < t) {
				return PG0;
			}

			double from = t;
			double to = t0;
			double oneMinusRho;

			double threshold  = T/10;

			int indexFrom = Utils.index(from, times, times.length);
			int index = Utils.index(to, times, times.length);

			int steps = index - indexFrom;
			if (Math.abs(from-times[indexFrom]) < globalPrecisionThreshold ) steps--;
			if (index>0 && Math.abs(to-times[index-1]) < globalPrecisionThreshold ) {
				steps--;
				index--;
			}
			index--;

			// pgScaled contains the set of initial conditions scaled made to fit the requirements on the values 'double' can represent. It also contains the factor by which the numbers were multiplied
			MultiSiteScaledNumbers pgScaled = SmallNumberScaler.scaleAllSites(PG0);
			
			// integrationResults will temporarily store the results of	 each integration step as 'doubles', before converting them back to 'SmallNumbers'
			double[] p0IntegrationResults = new double[fitnessModel.getFitSpaceStates()];
			double[][] gIntegrationResults = new double[siteCount+1][];
			for (int site=0; site<siteCount; site++) {
				gIntegrationResults[site] = new double[seqStates[site]];
			}
			gIntegrationResults[siteCount] = new double[1];

			while (steps > 0){

				from = times[index];
								
				
				gIntegrationResults = fastEulerIntegrateGxP0(pgScaled.getEquation(),to,from,index);
				PG0 = SmallNumberScaler.unscaleAllSites(p0IntegrationResults, gIntegrationResults, pgScaled.getScalingFactor());

				// Needs to be done for each site as well
				if (rhoChanges>0){
					oneMinusRho = 1-rho[index]; //rho[i*totalIntervals + index];
					//for (int i=0; i<fitnessModel.getFitSpaceStates(); i++){
						//PG0.conditionsOnP[i] *= oneMinusRho; // no longer need to do this here since pre-solving for P0 at all times
					//}
					for (int site=0; site<=siteCount; site++) {
						if (site == siteCount) {
							PG0.conditionsOnG[site][0] = PG0.conditionsOnG[site][0].scalarMultiply(oneMinusRho);
						} else {
							for (int i=0; i<seqStates[site]; i++) {	
								PG0.conditionsOnG[site][i] = PG0.conditionsOnG[site][i].scalarMultiply(oneMinusRho);
							}
						}
					}
				}

				to = times[index];

				steps--;
				index--;

				// 'rescale' the results of the last integration to prepare for the next integration step
				pgScaled = SmallNumberScaler.scaleAllSites(PG0);
			}
			
			// Integrate backwards starting w/ 'to' back to final 'from'
			gIntegrationResults = fastEulerIntegrateGxP0(pgScaled.getEquation(),to,from,index+1);
			PG0 = SmallNumberScaler.unscaleAllSites(p0IntegrationResults, gIntegrationResults, pgScaled.getScalingFactor());
			
			//if (PG0.isSumGZero(siteCount)) {
				//System.out.println("Prob density G evaluated as zero along a branch at time : " + to);
			//}
			//System.out.println();
			
		}catch(Exception e){

			throw new RuntimeException("couldn't calculate g");
		}

		if (pg_integrator.getEvaluations() > maxEvalsUsed) maxEvalsUsed = pg_integrator.getEvaluations();

		return PG0;
	}

	void setRho(){

		isRhoTip = new Boolean[ treeInput.get().getLeafNodeCount()];
		Arrays.fill(isRhoTip,false);

		isRhoInternalNode = new Boolean[ treeInput.get().getInternalNodeCount()];
		Arrays.fill(isRhoInternalNode,false);

		if (m_rho.get() != null) {

			constantRho = !(m_rho.get().getDimension() > n);

			if (m_rho.get().getDimension() <= n && (rhoSamplingTimes.get()==null || rhoSamplingTimes.get().getDimension() < 2)) {
				if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() <= n && samplingProportion.get().getValue() == 0.) || // todo:  instead of samplingProportion.get().getValue() == 0. need checked that samplingProportion[i]==0 for all i=0..n-1
						(samplingRate.get() != null && samplingRate.get().getDimension() <= 2 && samplingRate.get().getValue() == 0.))) {                              // todo:  instead of samplingRate.get().getValue() == 0. need checked that samplingRate[i]==0 for all i=0..n-1

					// check if data set is contemp!
					for (Node node : treeInput.get().getExternalNodes()){
						if (node.getHeight()>0.) throw new RuntimeException("Error in analysis setup: Parameters set for entirely contemporaneously sampled data, but some nodeheights are > 0!");
					}

					contempData = true;
					System.out.println("BDMM: setting contemp=true.");
				}
			}

			if (contempData) {
				if (m_rho.get().getDimension() != 1 && m_rho.get().getDimension() != n)
					throw new RuntimeException("when contemp=true, rho must have dimension 1 (or equal to the stateNumber)");

				else {
					rho = new Double[n*totalIntervals];
					Arrays.fill(rho, 0.);
					Arrays.fill(isRhoTip, true);
					for (int i=1; i<=n; i++)  rho[i*totalIntervals - 1] = m_rho.get().getValue(i-1);

					rhoSamplingCount = 1;
				}
			}
			else {
				Double[] rhos = m_rho.get().getValues();
				rho = new Double[n*totalIntervals];
				Arrays.fill(rho, 0.);
				for (int i = 0; i < totalIntervals; i++) {
					for (int j=0;j<n;j++){
						rho[j*totalIntervals+i]= rhoSamplingChangeTimes.contains(times[i]) ? (rhos[constantRho? j : j*(1+rhoChanges)+rhoSamplingChangeTimes.indexOf(times[i])]) : 0.;
					}
				}
				computeRhoTips();
			}


		} else {
			rho = new Double[n*totalIntervals];
			Arrays.fill(rho, 0.);
		}

	}

	abstract void computeRhoTips();


	/**
	 * Collect all the times of parameter value changes and rho-sampling events
	 */
	void collectTimes(double maxTime) {

		timesSet.clear();

		getChangeTimes(maxTime, migChangeTimes,
				migChangeTimesInput.get() != null ? migChangeTimesInput.get() : intervalTimes.get(),
						migChanges, migTimesRelative, reverseTimeArrays[5]);

		getChangeTimes(maxTime, birthRateChangeTimes,
				birthRateChangeTimesInput.get() != null ? birthRateChangeTimesInput.get() : intervalTimes.get(),
						birthChanges, birthRateTimesRelative, reverseTimeArrays[0]);

		getChangeTimes(maxTime, b_ijChangeTimes,
				b_ijChangeTimesInput.get() != null ? b_ijChangeTimesInput.get() : intervalTimes.get(),
						b_ij_Changes, b_ijTimesRelative, reverseTimeArrays[0]);

		getChangeTimes(maxTime, deathRateChangeTimes,
				deathRateChangeTimesInput.get() != null ? deathRateChangeTimesInput.get() : intervalTimes.get(),
						deathChanges, deathRateTimesRelative, reverseTimeArrays[1]);

		getChangeTimes(maxTime, samplingRateChangeTimes,
				samplingRateChangeTimesInput.get() != null ? samplingRateChangeTimesInput.get() : intervalTimes.get(),
						samplingChanges, samplingRateTimesRelative, reverseTimeArrays[2]);

		getChangeTimes(maxTime, rhoSamplingChangeTimes,
				rhoSamplingTimes.get()!=null ? rhoSamplingTimes.get() : intervalTimes.get(),
						rhoChanges, false, reverseTimeArrays[3]);

		if (SAModel) getChangeTimes(maxTime, rChangeTimes,
				removalProbabilityChangeTimesInput.get() != null ? removalProbabilityChangeTimesInput.get() : intervalTimes.get(),
						rChanges, rTimesRelative, reverseTimeArrays[4]);

		for (Double time : migChangeTimes) {
			timesSet.add(time);
		}

		for (Double time : birthRateChangeTimes) {
			timesSet.add(time);
		}

		for (Double time : b_ijChangeTimes) {
			timesSet.add(time);
		}

		for (Double time : deathRateChangeTimes) {
			timesSet.add(time);
		}

		for (Double time : samplingRateChangeTimes) {
			timesSet.add(time);
		}

		for (Double time : rhoSamplingChangeTimes) {
			timesSet.add(time);
		}

		if (SAModel) {
			for (Double time : rChangeTimes) {
				timesSet.add(time);
			}
		}


		times = timesSet.toArray(new Double[timesSet.size()]);
		totalIntervals = times.length;

	}

	/**
	 * set change times
	 */
	public void getChangeTimes(double maxTime, List<Double> changeTimes, RealParameter intervalTimes, int numChanges, boolean relative, boolean reverse) {
		changeTimes.clear();

		if (intervalTimes == null) { //equidistant

			double intervalWidth = maxTime / (numChanges + 1);

			double end;
			for (int i = 1; i <= numChanges; i++) {
				end = (intervalWidth) * i;
				changeTimes.add(end);
			}
			end = maxTime;
			changeTimes.add(end);

		} else {

			if (!reverse && intervalTimes.getValue(0) != 0.0) {
				throw new RuntimeException("First time in interval times parameter should always be zero.");
			}

			if (numChanges > 0 && intervalTimes.getDimension() != numChanges + 1) {
				throw new RuntimeException("The time interval parameter should be numChanges + 1 long (" + (numChanges + 1) + ").");
			}

			int dim = intervalTimes.getDimension();

			double end;
			for (int i = (reverse?0:1); i < dim; i++) {
				end = reverse ? (maxTime - intervalTimes.getValue(dim - i - 1)) : intervalTimes.getValue(i);
				if (relative) end *= maxTime;
				if (end < maxTime) changeTimes.add(end);
			}

			if (adjustTimesInput.get()!=null){

				double iTime;
				double aTime = adjustTimesInput.get().getValue();

				for (int i = 0 ; i < numChanges; i++){

					iTime = intervalTimes.getArrayValue(i+1);

					if (aTime<iTime) {
						end = iTime - aTime;
						if
						(changeTimes.size() > i) changeTimes.set(i, end);
						else
							if (end < maxTime)
								changeTimes.add(end);
					}
				}
			}
			end = maxTime;

			changeTimes.add(end);
		}
	}


	void updateBirthDeathPsiParams(){

		Double[] birthRates = birthRate.get().getValues();
		Double[] deathRates = deathRate.get().getValues();
		Double[] samplingRates = new Double[0];
		if (samplingCoupledToRemoval) {
			samplingRates = samplingProportion.get().getValues(); // psi will be samplingProportions
		}
		else {
			samplingRates = samplingRate.get().getValues(); // psi will be samplingRates
		}
		Double[] removalProbabilities = new Double[1];

		if (SAModel) {
			removalProbabilities = removalProbability.get().getValues();
			r =  new Double[n*totalIntervals];
		}

		int state;

		for (int i = 0; i < n*totalIntervals; i++) {

			state =  i/totalIntervals;

			birth[i] = (identicalRatesForAllTypes[0]) ? birthRates[index(times[i%totalIntervals], birthRateChangeTimes)] :
				birthRates[birthRates.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state];
			death[i] = (identicalRatesForAllTypes[1]) ? deathRates[index(times[i%totalIntervals], deathRateChangeTimes)] :
				deathRates[deathRates.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state];
			psi[i] = (identicalRatesForAllTypes[2]) ? samplingRates[index(times[i%totalIntervals], samplingRateChangeTimes)] :
				samplingRates[samplingRates.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state];
			if (SAModel) r[i] = (identicalRatesForAllTypes[4]) ? removalProbabilities[index(times[i%totalIntervals], rChangeTimes)] :
				removalProbabilities[removalProbabilities.length > n ? (rChanges+1)*state+index(times[i%totalIntervals], rChangeTimes) : state];

		}

	}


	void updateAmongParameter(Double[] param, Double[] paramFrom, int nrChanges, List<Double> changeTimes){

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int dt = 0; dt < totalIntervals; dt++) {
					if (i != j) {
						param[(i * (n - 1) + (j < i ? j : j - 1)) * totalIntervals + dt]
								= paramFrom[(paramFrom.length > (n * (n - 1)))
								            ? (nrChanges + 1) * (n - 1) * i + index(times[dt], changeTimes)
								            : (i * (n - 1) + (j < i ? j : j - 1))];
					}
				}
			}
		}

	}
	
	// Same as updateAmongParameter, but uses maxState instead of n since # of states might be variable across sites
	void updateAmongMigrationMatrix(Double[] param, Double[] paramFrom, int nrChanges, List<Double> changeTimes){

		for (int i = 0; i < maxState; i++) {
			for (int j = 0; j < maxState; j++) {
				for (int dt = 0; dt < totalIntervals; dt++) {
					if (i != j) {
						param[(i * (maxState - 1) + (j < i ? j : j - 1)) * totalIntervals + dt]
								= paramFrom[(paramFrom.length > (maxState * (maxState - 1)))
								            ? (nrChanges + 1) * (maxState - 1) * i + index(times[dt], changeTimes)
								            : (i * (maxState - 1) + (j < i ? j : j - 1))];
					}
				}
			}
		}

	}
	
//	void updateAmongMigrationMatrixCorrected(Double[] param, Double[] paramFrom, int nrChanges, List<Double> changeTimes){
//
//		for (int i = 0; i < maxState; i++) {
//			for (int j = 0; j < maxState; j++) {
//				for (int dt = 0; dt < totalIntervals; dt++) {
//					if (i != j) {
//						param[(i * (maxState - 1) + (j < i ? j : j - 1)) * totalIntervals + dt]
//								= paramFrom[(paramFrom.length > (maxState * (maxState - 1)))
//								            ? (nrChanges + 1) * ((maxState - 1) * i + (j < i ? j : j - 1)) + index(times[dt], changeTimes) // <-- edited to skip diagonal
//								            : (i * (maxState - 1) + (j < i ? j : j - 1))];
//					}
//				}
//			}
//		}
//
//	}



	void updateRho(){
		if (m_rho.get() != null && (m_rho.get().getDimension()==1 ||  rhoSamplingTimes.get() != null)) {

			Double[] rhos = m_rho.get().getValues();
			rho = new Double[n*totalIntervals];
			int state;

			for (int i = 0; i < totalIntervals*n; i++) {

				state =  i/totalIntervals;

				rho[i]= rhoChanges>0?
						rhoSamplingChangeTimes.contains(times[i]) ? rhos[rhos.length > n ? (rhoChanges+1)*state+index(times[i%totalIntervals], rhoSamplingChangeTimes) : state] : 0.
								: rhos[0];
			}
		}
	}


	/**
	 * @param t the time in question
	 * @return the index of the given time in the list of times, or if the time is not in the list, the index of the
	 *         next smallest time
	 *         This index function should only be used in transformParameters(), for likelihood calculations the times List needs to be used (e.g. with Untils.index(...))
	 */
	public int index(double t, List<Double> times) {

		int epoch = Collections.binarySearch(times, t);

		if (epoch < 0) {
			epoch = -epoch - 1;
		}

		return epoch;
	}


	public void transformWithinParameters(){

		Double[] p = samplingProportion.get().getValues();
		Double[] ds = becomeUninfectiousRate.get().getValues();
		Double[] R;
		if (R0.get() != null) {
			R = R0.get().getValues();
		} else {
			Double[] l_ratio = lambda_ratio.get().getValues();
			Double[] R_sens = R0_base.get().getValues();

			int totalIntervals = R_sens.length;
			int totalTypes = l_ratio.length + 1;
			R = new Double[totalIntervals * totalTypes];
			for (int i=0; i < totalIntervals; i++) {
				R[i] = R_sens[i];
				for (int j=1; j < totalTypes; j++) {
					double lambda = R_sens[i] * ds[ds.length > totalTypes ? index(times[i%totalIntervals], deathRateChangeTimes) : 0];
					R[i + totalIntervals * j] = (lambda * l_ratio[j - 1]) / ds[ds.length > totalTypes ? (deathChanges+1)*j+index(times[i%totalIntervals], deathRateChangeTimes) : j];
				}
			}
		}

		Double[] removalProbabilities = new Double[1];
		if (SAModel) removalProbabilities = removalProbability.get().getValues();

		int state;

		for (int i = 0; i < totalIntervals*n; i++){

			state =  i/totalIntervals;

			birth[i] = ((identicalRatesForAllTypes[0]) ? R[index(times[i%totalIntervals], birthRateChangeTimes)] :
				R[R.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state])
					* ((identicalRatesForAllTypes[1]) ? ds[index(times[i%totalIntervals], deathRateChangeTimes)] :
						ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state]);

			if (!SAModel) {
				psi[i] = ((identicalRatesForAllTypes[2]) ? p[index(times[i%totalIntervals], samplingRateChangeTimes)] :
					p[p.length > n ? (samplingChanges + 1) * state + index(times[i % totalIntervals], samplingRateChangeTimes) : state])
						* ((identicalRatesForAllTypes[1]) ? ds[index(times[i%totalIntervals], deathRateChangeTimes)] :
							ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state]);

				death[i] = ((identicalRatesForAllTypes[1]) ? ds[index(times[i%totalIntervals], deathRateChangeTimes)] :
					ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state])
						- psi[i];
			}

			else {
				r[i] = (identicalRatesForAllTypes[4]) ? removalProbabilities[index(times[i%totalIntervals], rChangeTimes)] :
					removalProbabilities[removalProbabilities.length > n ? (rChanges+1)*state+index(times[i%totalIntervals], rChangeTimes) : state];

				psi[i] = ((identicalRatesForAllTypes[2]) ? p[index(times[i%totalIntervals], samplingRateChangeTimes)] :
					p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state])
						* ((identicalRatesForAllTypes[1]) ? ds[index(times[i%totalIntervals], deathRateChangeTimes)] :
							ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state])
						/ (1+(r[i]-1)*
								((identicalRatesForAllTypes[2]) ? p[index(times[i%totalIntervals], samplingRateChangeTimes)] :
									p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]));


				death[i] = ((identicalRatesForAllTypes[1]) ? ds[index(times[i%totalIntervals], deathRateChangeTimes)] :
					ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state])
						- psi[i]*r[i];
			}
		}

	}


	public void transformAmongParameters(){

		Double[] RaD = (birthAmongDemes) ? R0AmongDemes.get().getValues() : new Double[1];
		Double[] ds = becomeUninfectiousRate.get().getValues();

		if (birthAmongDemes)    {

			for (int i = 0; i < n; i++){

				for (int j=0; j<n ; j++){

					for (int dt=0; dt<totalIntervals; dt++){

						if (i!=j){
							b_ij[(i*(n-1)+(j<i?j:j-1))*totalIntervals+dt]
									= RaD[(RaD.length>(n*(n-1)))
									      ?  (b_ij_Changes+1)*(n-1)*i + index(times[dt], b_ijChangeTimes)
									      : (i*(n-1)+(j<i?j:j-1))]
									    		  * ds[ds.length > n ? (deathChanges+1)*i+index(times[dt], deathRateChangeTimes) : i];
						}
					}
				}

			}
		}
	}


	void checkOrigin(TreeInterface tree){

		if (origin.get()==null){
			T = tree.getRoot().getHeight();
		}
		else {

			updateOrigin(tree.getRoot());

			if (!Boolean.valueOf(System.getProperty("beast.resume")) && orig < 0)
				throw new RuntimeException("Error: origin("+T+") must be larger than tree height("+tree.getRoot().getHeight()+")!");
		}

	}


	void updateOrigin(Node root){

		T = origin.get().getValue();
		orig = T - root.getHeight();

		if (originIsRootEdge.get()) {

			orig = origin.get().getValue();
			T = orig + root.getHeight();
		}

	}


	void setupIntegrators(){   // set up ODE's and integrators

		if (minstep == null) minstep = T*1e-100;
		if (maxstep == null) maxstep = T/10;

		Boolean augmented = this instanceof FitnessBirthDeathModel; //BirthDeathMigrationModel;
		
		// Set up integrators to solve in dimension of p0, not g equations
		int fitStates = fitnessModel.getFitSpaceStates();
		Double[] fitValues = fitnessModel.getFitnessSpaceValues();
		Double[] gammas = fitnessModel.getFitnessSpaceGammas();		
		
		// Need to pass fitValues into p0_ODE here!!
		P = new p0_ODE(birth, ((birthAmongDemes) ? b_ij : null), death,psi, gammas, fitStates, fitValues, totalIntervals, times, samplingCoupledToRemoval); // P integrator operates only in fitness space
		
		PG = new p0ge_ODE(birth, ((birthAmongDemes) ? b_ij : null), death,psi,M, fitStates, totalIntervals, T, times, P, maxEvaluations.get(), augmented); // set up for fitness space, not seq space

		p0ge_ODE.globalPrecisionThreshold = globalPrecisionThreshold;

		//		// TO DO CLEAN UP IF DOES NOT WORK
		//		double[] absoluteToleranceVector = new double [2*n];
		//		double[] relativeToleranceVector = new double [2*n];
		//		for(int i = 0; i<n; i++) {
		//			absoluteToleranceVector[i] = 1e-50;
		//			absoluteToleranceVector[i+n] = 1e-180;
		//			relativeToleranceVector[i] = 1e-7;
		//			relativeToleranceVector[i+n] = 1e-7;
		//		}


		if (!useRKInput.get()) {
			pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get()); //new HighamHall54Integrator(minstep, maxstep, absolutePrecision.get(), tolerance.get()); //new DormandPrince853Integrator(minstep, maxstep, absolutePrecision.get(), tolerance.get()); //new DormandPrince54Integrator(minstep, maxstep, absolutePrecision.get(), tolerance.get()); // 

			//			pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteToleranceVector, relativeToleranceVector);

			pg_integrator.setMaxEvaluations(maxEvaluations.get());

			PG.p_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
			PG.p_integrator.setMaxEvaluations(maxEvaluations.get());
		} else {
			pg_integrator = new ClassicalRungeKuttaIntegrator(T / 1000);
			PG.p_integrator = new ClassicalRungeKuttaIntegrator(T / 1000);

		}
	}



	/**
	 * Obtain element of rate matrix for migration model for use in likelihood
	 * calculation.
	 *
	 * @param i
	 * @param j
	 * @return Rate matrix element.
	 */
	public double getNbyNRate(int i, int j) {
		if (i==j)
			return 0;

		int offset = getArrayOffset(i, j);

		if (migrationMatrixScaleFactor.get()==null)
			return migrationMatrix.get().getValue(offset);
		else
			return migrationMatrixScaleFactor.get().getValue()*migrationMatrix.get().getValue(offset);

	}


	/**
	 * Obtain offset into "rate matrix" and associated flag arrays.
	 *
	 * @param i
	 * @param j
	 * @return Offset (or -1 if i==j)
	 */
	protected int getArrayOffset(int i, int j) {

		if (i==j)
			throw new RuntimeException("Programmer error: requested migration "
					+ "rate array offset for diagonal element of "
					+ "migration rate matrix.");


		if (j>i)
			j -= 1;
		return i*(n-1)+j;   // todo: check if this is correct!!!
	}



	// Interface requirements:

	@Override
	public List<String> getArguments() {
		return null;
	}


	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

	@Override
	public boolean requiresRecalculation(){
		return true;
	}

	/**
	 *
	 * @param integrator
	 * @param PG
	 * @param to
	 * @param pgScaled
	 * @param from
	 * @return
	 */
	public ScaledNumbers safeIntegrate(FirstOrderIntegrator integrator, p0ge_ODE PG, double to, ScaledNumbers pgScaled, double from){

		// if the integration interval is too small, nothing is done (to prevent infinite looping)
		if(Math.abs(from-to) < globalPrecisionThreshold /*(T * 1e-20)*/) return pgScaled;

		if(T>0 && Math.abs(from-to)>T/2 ) {
			pgScaled = safeIntegrate(integrator, PG, to, pgScaled, from + (to-from)/2);
			pgScaled = safeIntegrate(integrator, PG, from + (to-from)/2, pgScaled, from);
		} else {

			double relativeToleranceConstant = 1e-7;
			double absoluteToleranceConstant = 1e-100;
			double[] absoluteToleranceVector = new double [2*n];
			double[] relativeToleranceVector = new double [2*n];
			for(int i = 0; i<n; i++) {
				absoluteToleranceVector[i] = absoluteToleranceConstant;
				if(pgScaled.getEquation()[i+n] > 0) {
					absoluteToleranceVector[i+n] = Math.max(1e-310, pgScaled.getEquation()[i+n]*absoluteToleranceConstant);
				} else {
					absoluteToleranceVector[i+n] = absoluteToleranceConstant;
				}
				relativeToleranceVector[i] = relativeToleranceConstant;
				relativeToleranceVector[i+n] = relativeToleranceConstant;
			}

			double[] integrationResults = new double[pgScaled.getEquation().length];
			int a = pgScaled.getScalingFactor();
			int n = pgScaled.getEquation().length/2;


			integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteToleranceVector, relativeToleranceVector);
			integrator.integrate(PG, to, pgScaled.getEquation(), from, integrationResults);


			double[] pConditions = new double[n];
			SmallNumber[] geConditions = new SmallNumber[n];
			for (int i = 0; i < n; i++) {
				pConditions[i] = integrationResults[i];
				geConditions[i] = new SmallNumber(integrationResults[i+n]);
			}
			pgScaled = SmallNumberScaler.scale(new p0ge_InitialConditions(pConditions, geConditions));
			pgScaled.augmentFactor(a);

		}

		return pgScaled;
	}

	/**
	 * Find the lowest non-zero value among the integration results, only on the ge part.
	 * @param values
	 * @return
	 */
	public double getMinIntegrationOnGe(double[] values) {
		double min=Double.MAX_VALUE;
		if (values.length < 2)
			throw new RuntimeException("Invalid inital-conditions array size");
		for (int i = (values.length -1); i> (values.length/2 - 1) ; i-- ) {
			if (values[i] < min && values[i]!=0 ) min=values[i];
		}
		return min;
	}

	public p0ge_ODE getPG() {
		return PG;
	}

	/*
	 * Find all initial conditions for all future integrations on p0 equations
	 * Not currently used with multi-site fitness models b/c we find p0 conditions once at all times for all lineages
	 * @param tree
	 * @return an array of arrays storing the initial conditions values
	 */
	public double[][] getAllInitialConditionsForP(TreeInterface tree){
		
		int leafCount = tree.getLeafNodeCount();
		double[] leafHeights = new double[leafCount];
		int[] indicesSortedByLeafHeight  =new int[leafCount];
		for (int i=0; i<leafCount; i++){
			leafHeights[i] = T - tree.getNode(i).getHeight();
			// System.out.println(nodeHeight[i]);
			indicesSortedByLeafHeight[i] = i;
		}

		HeapSort.sort(leafHeights, indicesSortedByLeafHeight);
		//"sort" sorts in ascending order, so we have to be careful since the integration starts from the leaves at height T and goes up to the root at height 0 (or >0)

		double[][] pInitialCondsAtLeaves = new double[leafCount + 1][fitnessModel.getFitSpaceStates()];

		double t = leafHeights[indicesSortedByLeafHeight[leafCount-1]];

		boolean rhoSampling =  (m_rho.get()!=null);

		pInitialCondsAtLeaves[indicesSortedByLeafHeight[leafCount-1]] = PG.getP(t, rhoSampling, rho);
		double t0 = t;

		if (leafCount >1 ){
			for (int i = leafCount-2; i>-1; i--){
				t = leafHeights[indicesSortedByLeafHeight[i]];

				//If the next higher leaf is actually at the same height, store previous results and skip iteration
				if (Math.abs(t-t0) < globalPrecisionThreshold) {
					t0 = t;
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] = pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]];
					continue;
				} else {
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] = PG.getP(t, pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]], t0, rhoSampling, rho);
					t0 = t;
				}

			}
		}


		pInitialCondsAtLeaves[leafCount] = PG.getP(0, pInitialCondsAtLeaves[indicesSortedByLeafHeight[0]], t0, rhoSampling, rho);

		return pInitialCondsAtLeaves;
	}
	
	/*
	 * Compute p0 conditions for a grid of integration times
	 * @param tree
	 * @return an array of arrays storing the initial conditions values
	 */
	public double[][] getConditionsForPAllTimes(TreeInterface tree){
		
		// Get all times				
		double endTime = T;
		int dtSteps = (int) (T / dtTimeStep) + 1;
		integrationTimes = new Double[dtSteps];
		
		double time = 0.0;
		int cntr = 0;
		while (time < endTime){
			integrationTimes[cntr] = time;
			time += dtTimeStep;
			cntr++;
		}
		integrationTimes[dtSteps-1] = T; // ensure final time is included

		pConditions = new double[dtSteps][fitnessModel.getFitSpaceStates()];

		double t = T - integrationTimes[0];

		boolean rhoSampling =  (m_rho.get()!=null);

		pConditions[0] = PG.getP(t, rhoSampling, rho);
		double t0 = t;

		for (int i = 1; i<dtSteps; i++){
			t = T - integrationTimes[i];
			if (Math.abs(t-t0) < globalPrecisionThreshold) {
				t0 = t;
				pConditions[i] = pConditions[i-1];
				continue;
			} else {
				pConditions[i] = PG.getP(t, pConditions[i-1], t0, rhoSampling, rho);
				t0 = t;
			}
		}

		return pConditions;
	}

	protected Double updateRates() {

		birth = new Double[n*totalIntervals];
		death = new Double[n*totalIntervals];
		psi = new Double[n*totalIntervals];
		b_ij = new Double[totalIntervals*(n*(n-1))];
		M = new Double[totalIntervals*(maxState*(maxState-1))]; // <-- M has dimension n = seqStates;
		if (SAModel) r =  new Double[n * totalIntervals];

		if (transform) {
			transformParameters();
		}
		else {
			
			// birthAmongDemes should not be used in multi-site fitness models
			
			Double[] birthAmongDemesRates = new Double[1];

			if (birthAmongDemes) birthAmongDemesRates = birthRateAmongDemes.get().getValues();

			updateBirthDeathPsiParams();

			if (birthAmongDemes) {

				updateAmongParameter(b_ij, birthAmongDemesRates, b_ij_Changes, b_ijChangeTimes);
			}
		}

		if (migrationMatrix.get()!=null) {
			
			Double[] migRates = migrationMatrix.get().getValues();
			
			// Added this just for FluGenotypeFitnessModel
			// Can remove if no longer needed
			int migRatesSize = maxState * (maxState-1); 
			if (migRates.length != migRatesSize)  {
				double mRate = migRates[0];
				migRates = new Double[migRatesSize];
				for (int state=0; state<migRatesSize; state++) migRates[state] = mRate;
			}

			Double factor;
			if (migrationMatrixScaleFactor.get()!=null) {
				factor = migrationMatrixScaleFactor.get().getValue();
				for (int i = 0; i < migRates.length; i++) migRates[i] *= factor;
			}
			
			updateAmongMigrationMatrix(M, migRates, migChanges, migChangeTimes); // seems to work even without corrected method

		}

		updateRho();

		freq = frequencies.get().getValues();
		
		// Update fitness model
		fitnessModel.update();

		setupIntegrators();

		return 0.;
	}

	public void transformParameters(){

		transformWithinParameters();
		transformAmongParameters();
	}

	Boolean[] identicalRatesForAllTypes;

}
