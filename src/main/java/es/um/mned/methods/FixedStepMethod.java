package es.um.mned.methods;

import java.util.Iterator;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.ode.NumericalSolutionPoint;
import es.um.mned.ode.Event;
import es.um.mned.interpolation.StateFunction;
import es.um.mned.utils.BisectionMethod;
import es.um.mned.utils.ConvergenceException;

abstract public class FixedStepMethod {
 
    /**
     * Computes the max error of two solutions taken one at half the step of the other
     * @param fullStep
     * @param halfStep
     * @return 
     */
    static public double maxHalfStepError (NumericalSolution fullStep, NumericalSolution halfStep) {
        Iterator<NumericalSolutionPoint> iteratorFull = fullStep.iterator();
        Iterator<NumericalSolutionPoint> iteratorHalf = halfStep.iterator();
        double maxError = 0;
        while (iteratorFull.hasNext() && iteratorHalf.hasNext()) {
            double[] stateFull = iteratorFull.next().getState();
            double[] stateHalf = iteratorHalf.next().getState();
            double estimatedError = 0;
            for (int i=0; i<stateFull.length; i++) {
                double estimatedErrorInI = Math.abs(stateHalf[i]-stateFull[i]); 
                estimatedError = Math.max(estimatedError,estimatedErrorInI);
            }
            maxError = Math.max(maxError, estimatedError);
            if (!iteratorHalf.hasNext()) return maxError;
            iteratorHalf.next();
        }
        return maxError;
    }
    
    private double mStep;
    private NumericalSolution mSolution;
    private Event event = null;
	protected double currentUserTime;
    protected InitialValueProblem mProblem;
    
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    protected FixedStepMethod(InitialValueProblem problem, double step) {
        mProblem = problem;
        mStep = step;
        currentUserTime = problem.getInitialTime();
        mSolution = new NumericalSolution(problem, getOrder());
    }

    protected FixedStepMethod(InitialValueProblem problem, double step, Event event) {
        this(problem, step);
        this.event = event;
    }
    
    abstract public int getOrder();
    
    /**
     * Particular method implementation
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     * @throws ConvergenceException 
     */
    abstract public double doStep(double deltaTime, double time, double[] state) throws ConvergenceException;
    
    /**
     * Get the step
     * @return the initial step given
     */
    public double getStep() {
        return mStep;
    }
    
    public Event getEvent() {
		return event;
	}

	public void setEvent(Event event) {
		this.event = event;
	}


	// Converts event to state function
    private static class EventStateFunction implements StateFunction {
        private StateFunction interpolator;
        private Event event;

        public EventStateFunction(StateFunction interpolator, Event event) {
            this.interpolator = interpolator;
            this.event = event;
        }
        
        public double[] getState(double t) {
            return new double[] {
                event.crossFunction(t, interpolator.getState(t))
            };
        }

        public double getState(double t, int index) {
            return event.crossFunction(t, interpolator.getState(t));
        }
    }

    // Checks if the event occurred between the last 2 points
    private boolean checkEvent() throws ConvergenceException {
        if(event == null) return false;
        int size = mSolution.getSize();
        if(size <= 1) return false;

        NumericalSolutionPoint p1 = mSolution.get(size-2);
        NumericalSolutionPoint p2 = mSolution.get(size-1);
        
        

        if(event.crossFunction(p1.getTime(), p1.getState())
                * event.crossFunction(p2.getTime(), p2.getState()) <= 0) {
        	
        	StateFunction interpolator = mSolution.getInterpolator(size-1);

            double zero = BisectionMethod.findZero(
                    new EventStateFunction(interpolator, event),
                    p1.getTime(),
                    p2.getTime(),
                    event.getTolerance(),
                    0
                    );
            event.crossAction(zero, interpolator.getState(zero));
            return true;
        }
        return false;
    }

    protected double solveUpTo(double maxTime) throws ConvergenceException {
        NumericalSolutionPoint lastPoint = mSolution.getLastPoint();
        double time = lastPoint.getTime();
        double[] state = lastPoint.getState();
        if (mStep > 0) {
            while (time < maxTime) {
                time = doStep(mStep,time,state);
                if (Double.isNaN(time)) return Double.NaN;
                mSolution.add(time, state);
                if(checkEvent() && event.stopCondition())
                	break;
            }
        } 
        else if (mStep < 0) {
            while (time > maxTime) {
                time = doStep(mStep,time,state);
                if (Double.isNaN(time)) return Double.NaN;
                mSolution.add(time, state);
                if(checkEvent() && event.stopCondition())
                	break;
            }
        }
        return time;
    }
    
    /**
     * Steps the problem once
     * @return the newly computed solution point, null if there was any error
     * @throws ConvergenceException 
     */
    public NumericalSolutionPoint step() throws ConvergenceException {
    	currentUserTime += mStep;
    	
        NumericalSolutionPoint lastPoint = mSolution.getLastPoint();
        double time = lastPoint.getTime();
        double[] state = lastPoint.getState();
        time = doStep(mStep,time,state);
        
        if (Double.isNaN(time)) return null;
        
        
        NumericalSolutionPoint point = mSolution.add(time, state);
        checkEvent(); // No need to check if blocking
        
        return point;
    }
    
    /**
     * Iteratively steps the problem until time equals or exceeds finalTime
     * @param finalTime the time which we want to reach or exceed
     * @return the actual time of the last computed solution point (may differ -exceed- the requested finalTime).
     * returns NaN if there was any error in the solving
     * @throws ConvergenceException 
     */
    public NumericalSolution solve(double finalTime) throws ConvergenceException {
        currentUserTime = finalTime;
        solveUpTo(finalTime);
        return mSolution;
    }
    
    /**
     * Gets the solution computed so far
     * @return an instance of NumericalSolution
     */
    public NumericalSolution getSolution() { return mSolution; }

}
