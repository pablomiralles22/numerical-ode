package es.um.mned.methods;

import java.util.Optional;

import es.um.mned.ode.ConvergenceException;
import es.um.mned.ode.Event;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolutionPoint;

abstract public class AdaptiveStepMethod extends FixedStepMethod {
	
	private final double DEFAULT_TOL = 1e-4;
 
    protected double mTolerance;
    protected double mCurrentStep;
    protected double mMinimumStepAllowed; // Non-convergence minimum

    public AdaptiveStepMethod(
    		InitialValueProblem problem,
    		double step,
    		Optional<Double> tol,
    		Optional<Double> minStep,
    		Optional<Event> event
    		) {
        super(problem,step, event);
        mCurrentStep = step;
        mMinimumStepAllowed = Math.abs(minStep.orElse(step / 1.0e6));
        mTolerance = tol.orElse(DEFAULT_TOL);
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
