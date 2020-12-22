/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.problems;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.ode.*;

public class ParabolicThrowWithFriction extends InitialValueProblem {
    static private double mGravity = 9.8;

    private final double mFrictionCoefficient;
    private final double mMass = 1;
    
    private final double constant;
    
    public ParabolicThrowWithFriction(double t0, double[] x0, double coefficient) {
    	super(t0, x0);
        mFrictionCoefficient = coefficient;
        constant = mFrictionCoefficient/mMass;
    }
    // ------------------
    // Implementation of InitialValueProblem
    // ------------------

    public double[] getDerivative(double t, double[] x) {
    	super.addToEvaluationCounter();
        double speed = Math.sqrt(x[1]*x[1]+x[3]*x[3]);
        return new double[] { 
        		x[1],
        		-constant * x[1] * speed, 
                x[3],
                -constant * x[3] * speed - mGravity
                };
    }

    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------

}
