/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.ode;

import java.util.Arrays;

/**
 * Interface for an InitialValueProblem of Ordinary Differential Equations
 * 
 * @author F. Esquembre
 * @version September 2020
 */
public abstract class InitialValueProblem {
	double t0;
	double[] x0;
	int evaluationCounter;
    
    public InitialValueProblem(double t0, double[] x0) {
		this.t0 = t0;
		this.x0 = x0;
		evaluationCounter = 0;
	}

	/**
     * The initial value of time (independent variable)
     * @return 
     */
    public double getInitialTime() {
    	return t0;
    }
    
    /**
     * The initial value of the state
     * @return a newly created array with a copy of the initial state
     */
    public double[] getInitialState() {
    	return Arrays.copyOf(x0, x0.length);
    }
    
    /**
     * Computes the derivative f(t,Y(t)) that defines the ODE
     * @param time the given time
     * @param state the given state
     * @return a newly created array with the value of the derivative
     */
    public abstract double[] getDerivative(double time, double[] state);
    
    public void addToEvaluationCounter() {
    	evaluationCounter++;
    }
    
    public int getEvaluationCounter() {
    	return evaluationCounter;
    }
    
    public void resetEvaluationCounter() {
    	evaluationCounter = 0;
    }
    
}
