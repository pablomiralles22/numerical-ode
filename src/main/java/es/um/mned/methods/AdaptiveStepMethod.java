/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.methods;

import java.util.ArrayList;
import es.um.mned.ode.InitialValueProblem;

/**
 * Abstract class for a Fixed Step Method to solve an InitialValueProblem
 * 
 * @author F. Esquembre
 * @version September 2020
 */
abstract public class AdaptiveStepMethod extends FixedStepMethod {
 
    protected double mTolerance = 1.0e-4;
    protected ArrayList<Double> mStepList = new ArrayList<Double>();

    public AdaptiveStepMethod(InitialValueProblem problem, double step) {
        super(problem,step);
    }
    
    public void setTolerance(double tolerance) {
        mTolerance = tolerance;
    }

    public double getTolerance() {
        return mTolerance;
    }
    
    public ArrayList<Double> getStepList() {
        return mStepList;
    }
    
}
