package beast.math;

import beast.core.Description;

/**
 * Created by Jeremie Scire (jscire) on 24.06.16.
 */

@Description("Set of initial conditions for a system of ODEs, all increased by the same factor to prevent underflowing (which appears with numbers too small to be distinguished from zero, when using the type 'double')")
public class MultiSiteScaledNumbers {

	// scale factor(s)
	private int factor;

	// set(s) of initial conditions on which the scale factor(s) was/were applied
	private double[][] equation;
	
	private double[] p0equation;

	public MultiSiteScaledNumbers(int f, double[][] e, double[] p){
		this.factor = f;
		this.equation = e;
		this.p0equation = p;
	}

	public MultiSiteScaledNumbers(){
		this.factor = 0;
		this.equation = new double[0][0];
		this.p0equation = new double[0];
	}

	public double[][] getEquation(){
		return this.equation;
	}
	
	public double[] getP0Equation(){
		return this.p0equation;
	}
	
	public void setEquation(double[][] newEquation){
		this.equation = newEquation;
	}

	public int getScalingFactor(){
		return this.factor;
	}
	
	public void setScalingFactor(int newFactor){
		this.factor = newFactor; 
	}
	
	public void augmentFactor(int increase) {
		factor += increase;
	}
}
