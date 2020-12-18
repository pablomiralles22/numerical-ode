/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.problems;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.interpolation.*;
import es.um.mned.methods.*;
import es.um.mned.ode.*;
import es.um.mned.utils.*;

/**
 *
 * @author paco
 */
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

    public double getInitialTime() { 
        return 0; 
    }
    
    public double[] getInitialState() { // x,vx, y,vy (m/s)
        return new double[] { 0, 100, 300, 0 };
    } 
    
    public double[] getDerivative(double t, double[] x) {
    	super.addToEvaluationCounter();
        double speed = Math.sqrt(x[1]*x[1]+x[3]*x[3]);
        return new double[] { x[1], -constant * x[1] * speed, 
                              x[3], -constant * x[3] * speed - mGravity };
    }

    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------

    
    static public class TrueSol implements StateFunction {
            public double[] getState(double time) {
                return new double[] { 100*time, 100, 300 - 0.5*mGravity*time*time, -mGravity};
            }
            public double getState(double time, int index) {
                switch (index) {
                    case 0 : return 100*time;
                    case 1 : return 100;
                    case 2 : return 300 - 0.5*mGravity*time*time;
                    case 3 : return -mGravity;
                    default : return Double.NaN;
                }
            }
                
            public double yZeroAt() {
                return Math.sqrt(2*300/mGravity);
            }

    }
}
