package es.um.mned.methods;

import es.um.mned.ode.ConvergenceException;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolutionPoint;

abstract public class AdaptiveStepMethod extends FixedStepMethod {
	
 
    protected double mTolerance = 1.0e-4;
    // protected ArrayList<Double> mStepList = new ArrayList<Double>();

    public AdaptiveStepMethod(InitialValueProblem problem, double step) {
        super(problem,step);
    }
    
    public void setTolerance(double tolerance) {
        mTolerance = tolerance;
    }

    public double getTolerance() {
        return mTolerance;
    }
    
    @Override
    public NumericalSolutionPoint step() throws ConvergenceException {
        currentUserTime += getStep();
        double time = solveUpTo(currentUserTime);
        if(time == Double.NaN) return null;
        return new NumericalSolutionPoint(
                currentUserTime,
                getSolution().getState(currentUserTime)
                );
    }
    
}
