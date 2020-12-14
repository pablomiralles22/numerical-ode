/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.problems;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.ode.*;

/**
 *
 * @author paco
 */
public class Rigid1D extends ExtendedInitialValueProblem {

    // ------------------
    // Implementation of ExtendedInitialValueProblem
    // ------------------

    public Rigid1D(double t0, double[] x0) {
		super(t0, x0);
	}
    
    public double[] getDerivative(double t, double[] x) {
    	super.addToEvaluationCounter();
        return new double[] { 
        	5 * Math.exp(5*t) * (x[0]-t) * (x[0]-t) + 1
        };
    }
    
    public double[] getDerivativeDY(double t, double[] x) {
    	return new double[] {
    		10 * Math.exp(5*t) * (x[0] - t)
    	};
    }

    // ------------------
    // End of implementation of ExtendedInitialValueProblem
    // ------------------

    static public class TrueSol implements StateFunction {
            
        public double[] getState(double time) {
            return new double[] { time - Math.exp(-5*time) };
        }
            
        public double getState(double time, int index) {
            return time - Math.exp(-5*time);
        }
                
    }
    
}
