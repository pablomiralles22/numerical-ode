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
public class SimpleHarmonicOscillator extends InitialValueProblem {
    private double l = 0.7;
    private double m = 1.0;
    private double k = 1.5;

    private double b = 0.; // 0.3
    private double amp = 0.; // 0.4
    private double freq = 1.3; // 2.4
    
    private double Xo = 1.5;
    
    public SimpleHarmonicOscillator(
            double t0, double[] x0,
            double l, double m, double k,
            double b, double amp, double freq
            ) {
    	super(t0, x0);
        Xo = x0[0];
        this.l = l; this.m = m; this.k = k;
        this.b = b; this.amp = amp; this. freq = freq;
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

    
    private class TrueSol implements StateFunction {
            double Fo = Math.sqrt(k/m);
            
            public double[] getState(double time) {
                return new double[] { 
                    (Xo-l)*Math.cos(Fo*time) + l, 
                    -Fo*(Xo-l)*Math.sin(Fo*time)
                };
            }
            public double getState(double time, int index) {
                switch (index) {
                    case 0 : return (Xo-l)*Math.cos(Fo*time)+l;
                    case 1 : return -Fo*(Xo-l)*Math.sin(Fo*time);
                    default : return Double.NaN;
                }
            }

    }

    public StateFunction getTrueSol() {
        return new TrueSol();
    }

}
