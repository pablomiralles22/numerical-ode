package es.um.mned.methods;

import java.util.Optional;

import es.um.mned.ode.ConvergenceException;
import es.um.mned.ode.Event;
import es.um.mned.ode.InitialValueProblem;

public class AdaptiveStepEulerMethod extends AdaptiveStepMethod {
	
    private double[] mHalfStepState;
    private double[] mFullStepState; 
    
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public AdaptiveStepEulerMethod(
    		InitialValueProblem problem,
    		double step,
    		Optional<Double> tolerance,
    		Optional<Double> minStep,
    		Optional<Event> event
    		) {
        super(problem, step, tolerance, minStep, event);
        
        mHalfStepState = problem.getInitialState();
        mFullStepState = problem.getInitialState();
    }
    
    @Override
    public int getOrder() {
    	return 1;
    }

    /**
     * Extrapolated Euler method implementation
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     * @throws ConvergenceException 
     */
    public double doStep(double deltaTime, double time, double[] state) throws ConvergenceException {
        while (Math.abs(mCurrentStep)>=mMinimumStepAllowed) {
            double[] derivative = mProblem.getDerivative(time, state);
            double halfStep = mCurrentStep/2;
            for (int i=0; i<state.length; i++) {
                mHalfStepState[i] = state[i] + halfStep     * derivative[i];
                mFullStepState[i] = state[i] + mCurrentStep * derivative[i];
            }
            derivative = mProblem.getDerivative(time+halfStep, mHalfStepState);
            double error = 0;
            for (int i=0; i<state.length; i++) {
                mHalfStepState[i] += halfStep * derivative[i];
                double errorInIndex = mHalfStepState[i]-mFullStepState[i];
                error += errorInIndex*errorInIndex;
            }
            error = Math.sqrt(error);
            if (error<mTolerance*Math.abs(mCurrentStep)) {
                for (int i=0; i<state.length; i++) {
                    //state[i] = mHalfStepState[i]; 
                    state[i] = 2*mHalfStepState[i] - mFullStepState[i];
                }
                time += mCurrentStep;
                // Adapt step
                if (error<1.0e-10) mCurrentStep = 2*mCurrentStep;
                else {
                    double q = 0.84*(mTolerance*Math.abs(mCurrentStep))/error;
                    mCurrentStep *= q;
                }
                //System.out.println ("ACCEPTED: t = "+time+ " New step is "+mCurrentStep+ " error = "+error);
                return time;
            }
            // Try a new smaller step
            double q = 0.84*(mTolerance*Math.abs(mCurrentStep))/error;
            mCurrentStep *= q;
//            System.out.println ("REJECTED: t = "+time+ " New step is "+mCurrentStep+ " error = "+error);
        }
        throw new ConvergenceException("Adaptative Euler Method did not converge.");
//        // Was not able to reach tolerance before going below mMinimumStepAllowed
//        return Double.NaN; 
    }
    
}
