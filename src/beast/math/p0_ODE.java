package beast.math;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.core.util.Utils;


/**
 * @author dkuh004
 *         Date: May 24, 2012
 *         Time: 6:42:00 PM
 */

public class p0_ODE implements FirstOrderDifferentialEquations {

	Double[] b;
	Double[] b_ij;
	Double[] d;
	Double[] s;

	Double[] M;
	
	Double[] fitValues;

	int dimension;
	int intervals;
	Double[] times;
	int index;
	
	Boolean samplingCoupledToRemoval;

	public p0_ODE(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M, int dimension, Double[] fitValues, int intervals, Double[] times, Boolean samplingCoupledToRemoval) {

		this.b = b;
		this.b_ij = b_ij;
		this.d = d;
		this.s = s;
		this.M = M;
		this.dimension = dimension;
		this.fitValues = fitValues;
		this.intervals = intervals;

		this.times = times;
		
		this.samplingCoupledToRemoval = samplingCoupledToRemoval;

	}

	public void updateRates(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M,  Double[] fitValues, Double[] times){

		this.b = b;
		this.b_ij = b_ij;
		this.d = d;
		this.s = s;
		this.M = M;
		this.fitValues = fitValues;
		this.times = times;

	}

	public int getDimension() {
		return this.dimension;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {

		index = Utils.index(t, times, intervals); //finds the index of the time interval t lies in 
		int k, l;

		for (int i = 0; i<dimension; i++){

			k = index;
		
			double lambda = b[k] * fitValues[i]; // need to pass in fitValues
			if (samplingCoupledToRemoval) {
				yDot[i] = + (lambda+d[k])*y[i] - (d[k] * (1-s[k])) - lambda*y[i]*y[i] ;
			} else {
				yDot[i] = + (lambda+d[k]+s[k])*y[i] - d[k] - lambda*y[i]*y[i] ;
			}

			for (int j=0; j<dimension; j++){

				//l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;
				l = i*(dimension-1)+(j<i?j:j-1); // linear indexing

				if (i!=j){

					// No birthAmongDemes in multi-site fitness model
					//if (b_ij!=null){     // infection among demes
						//yDot[i] += b_ij[l]*y[i]; 
						//yDot[i] -= b_ij[l]*y[i]*y[j];
					//}

					if (M[0]!=null) {// migration:
						yDot[i] += M[l] * y[i];
						yDot[i] -= M[l] * y[j];
					}
				}
			}
		}

	}

	/**
	 * Some basic tests for numerical integration with this class
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception{
		
		// 2d test
		Double[] b = {1.03,1.06};
		Double[] d = {1.,1.};
		Double[] s = {0.02,0.04};
		Double[] M = new Double[]{3.,4.};
		Double[] fitValues = new Double[]{1.,1.};

		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-20, 1.0e-9);//new ClassicalRungeKuttaIntegrator(.01); //
		FirstOrderDifferentialEquations ode = new p0_ODE(b,null,d,s,M, 2, fitValues, 1, new Double[]{0.}, false);
		double[] y0 = new double[]{1.,1.};
		double[] y = new double[2];

		//        StepHandler stepHandler = new StepHandler() {
		//            public void init(double t0, double[] y0, double t) {
		//            }
		//
		//            public void handleStep(StepInterpolator interpolator, boolean isLast) {
		//                double   t = interpolator.getCurrentTime();
		//                double[] y = interpolator.getInterpolatedState();
		//                System.out.println(t + " " + y[0] + " " + y[1]);
		//            }
		//        };
		//        integrator.addStepHandler(stepHandler);
		//
		//
		integrator.integrate(ode, 10, y0, 1, y);

		System.out.println("Solution: " +y[0]+" "+y[1]);
		//
		//
		//        // 3d test
		//        Double[] b = {1.03,1.06, 1.5};
		//        Double[] d = {1.,1., 1.2};
		//        Double[] s = {0.02,0.04, 0.1};
		//        Double[] M = {3., 1.,4.,1.,2., 2.};
		//
		//        FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);//new ClassicalRungeKuttaIntegrator(.01); //
		//        FirstOrderDifferentialEquations ode = new p0_ODE(b,d,s,M, 3);
		//        double[] y0 = new double[]{1.,1.,1.};
		//        double[] y = new double[3];
		//
		//        integrator.integrate(ode, 10, y0, 1, y);
		//
		//        System.out.println("Solution: " +y[0]+" "+y[1]+" "+y[2]);
	}

}


