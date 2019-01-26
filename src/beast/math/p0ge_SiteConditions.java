package beast.math;

import beast.core.Description;

/**
 * Created by Jeremie Scire (jscire)
 * Modified by David Rasmussen to hold g conditions for every site
 */

@Description("This type contains both sets of initial conditions for both equation types: p and ge")
public class p0ge_SiteConditions {
	
	int dimension;
	public SmallNumber[][] conditionsOnG;
	public double[] conditionsOnP;
	
	public p0ge_SiteConditions(double[] pcond, SmallNumber[][] gcond) {
		
		// This will no longer be true since pcond will be in fitness space and gcond will be in seq space
		//if(pcond.length != gcond.length) {
			//throw new RuntimeException("Incorrect initialization: difference of size between conditionsOnG and conditionsOnP");
		//}
		
		dimension = gcond.length; // or should this be pcond.length?
		conditionsOnP = pcond;
		conditionsOnG = gcond;
	}
	
	public p0ge_SiteConditions() {		
		dimension = 1;
		conditionsOnP = new double[] {0};
		conditionsOnG = new SmallNumber[][] {{new SmallNumber()}}; // proper init for a 2D array?
	}
	
	public double[] getConditionsOnP(){
		return this.conditionsOnP;
	}
	
	public SmallNumber[][] getConditionsOnG(){
		return this.conditionsOnG;
	}
	
	public boolean isSumGZero(int site) {
		SmallNumber sum = new SmallNumber(0.0);
		for (int i = 0; i < dimension; i++) {
			sum = SmallNumber.add(sum, conditionsOnG[site][i]);
		}
		double convertedSum = sum.revert();
		return (convertedSum <= 0);
	}

}
