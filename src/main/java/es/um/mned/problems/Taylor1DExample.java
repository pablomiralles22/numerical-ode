/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.problems;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.*;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.DisplaySolution;

/**
 *
 * @author paco
 */


public class Taylor1DExample extends InitialValueProblem {
    
    public Taylor1DExample(double t0, double[] x0) {
		super(t0, x0);
	}
        
    // ------------------
    // Implementation of InitialValueProblem
    // ------------------
    
    public double[] getDerivative(double t, double[] y) {
        return new double[] { y[0]-t*t+1 };
    }
    
    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------
    
    
}
