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
import es.um.mned.tools.*;

/**
 *
 * @author paco
 */
public class TwoBodyProblem extends InitialValueProblem {

	static private double sG = 8.6498928e-4;

    private double mM1 = 1988.5, mM2 = 0.0059724;
    private double mConstant = sG * (mM1 + mM2);
    
    // ------------------
    // Implementation of InitialValueProblem
    // ------------------

    public TwoBodyProblem(double t0, double[] x0) {
		super(t0, x0);
	}
    
    public double[] getDerivative(double t, double[] x) {
    	super.addToEvaluationCounter();
        double div  = Math.pow(x[0]*x[0]+x[2]*x[2],1.5);
        return new double[] { x[1], -mConstant * x[0] / div, 
                              x[3], -mConstant * x[2] / div };
    }

    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------    

}
