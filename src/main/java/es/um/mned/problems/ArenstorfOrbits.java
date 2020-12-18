/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.problems;

import java.util.Arrays;
import es.um.mned.interpolation.StateFunction;
import es.um.mned.interpolation.*;
import es.um.mned.methods.*;
import es.um.mned.ode.*;
import es.um.mned.utils.*;

/**
 * https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/
 * @author paco
 */
public class ArenstorfOrbits extends InitialValueProblem {
	
    static public double PERIOD = 17.0652165601579625588917206249;

	static private double sMu = 0.012277471;
    static private double sMuPrime = 1-sMu;
    
    public ArenstorfOrbits(double t0, double[] x0) {
		super(t0, x0);
	}
    
    public double[] getDerivative(double t, double[] x) {
    	super.addToEvaluationCounter();
        double D1 = Math.pow((x[0]+sMu)*(x[0]+sMu) + x[2]*x[2],1.5);
        double D2 = Math.pow((x[0]-sMuPrime)*(x[0]-sMuPrime) + x[2]*x[2],1.5);
        return new double[] { 
            x[1], 
            x[0] + 2*x[3] - sMuPrime*(x[0]+sMu)/D1 - sMu*(x[0]-sMuPrime)/D2,
            x[3], 
            x[2] - 2*x[1] - sMuPrime*x[2]/D1 - sMu*x[2]/D2,
        };
    }
    
    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------

}
