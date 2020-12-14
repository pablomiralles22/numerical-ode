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
import es.um.mned.tools.*;

/**
 *
 * @author paco
 */
public class SimpleHarmonicOscillator extends InitialValueProblem {
    static private double l = 0.7;
    static private double m = 1.0;
    static private double k = 1.5;

    static private double b = 0.; // 0.3
    static private double amp = 0.; // 0.4
    static private double freq = 1.3; // 2.4
    
    static private double Xo = 1.5;
    static private double Vo = 0;
    
    public SimpleHarmonicOscillator(double t0, double[] x0) {
    	super(t0, x0);
    }
    
    private double force(double time) {
        return amp*Math.sin(freq*time);
    }
    // ------------------
    // Implementation of InitialValueProblem
    // ------------------

    
    public double[] getDerivative(double t, double[] x) {
        super.addToEvaluationCounter();
        return new double[] { x[1], 
                              -k/m * (x[0]-l) - b/m*x[1] + force(t)/m, 
                              };
    }

    
    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------

    
    static public class TrueSol implements StateFunction {
            static double Fo = Math.sqrt(k/m);
            
            public double[] getState(double time) {
                return new double[] { (Xo-l)*Math.cos(Fo*time) + l, 
                                      -Fo*(Xo-l)*Math.sin(Fo*time) };
            }
            public double getState(double time, int index) {
                switch (index) {
                    case 0 : return (Xo-l)*Math.cos(Fo*time)+l;
                    case 1 : return -Fo*(Xo-l)*Math.sin(Fo*time);
                    default : return Double.NaN;
                }
            }
                

    }
    

}
