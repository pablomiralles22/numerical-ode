package es.um.mned.methods;

import es.um.mned.ode.Event;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.utils.ConvergenceException;

public class AdaptiveStepPredictorCorrector4Method extends AdaptiveStepMethod {
    static public final int sSTEPS = 4;

    private double mCurrentStep;
    private double mMinimumStepAllowed; // Non-convergence minimum

    protected boolean mMustRestart=true;
    protected double[] mPredictorState, mCorrectorState;
    protected double[] mAuxState; // Required by the RK starter
    protected double[]   mTimes       = new double[sSTEPS-1];   // Times taken at restart
    protected double[][] mStates      = new double[sSTEPS-1][]; // ordered 2 = (i-2), 1 = (i-1) , 0 = i
    protected double[][] mDerivatives = new double[sSTEPS][];   // ordered 3 = (i-3), 2 = (i-2), 1 = (i-1) , 0 = i

    private int dim;
    private int queueTop = -1;
    private double[][] queue;

    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public AdaptiveStepPredictorCorrector4Method(InitialValueProblem problem, double step, double tolerance) {
        super(problem,step);
        super.setTolerance(tolerance);
        mCurrentStep = step;
        mMinimumStepAllowed = Math.abs(step)/1.0e6;
        mPredictorState = problem.getInitialState();
        mCorrectorState = problem.getInitialState();
        for (int i=0; i<mStates.length; i++) mStates[i] = problem.getInitialState();
        mAuxState = problem.getInitialState();
        // ------
        dim = problem.getInitialState().length;
        queue = new double[sSTEPS-1][dim];
    }
    
    public AdaptiveStepPredictorCorrector4Method(InitialValueProblem problem, double step, double tolerance, Event event) {
        this(problem, step, tolerance);
        super.setEvent(event);
    }
    
    
    @Override
    public int getOrder() {
    	return 4;
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
        if(queueTop >= 0) { // QUEUED STEPS
            System.arraycopy(queue[queueTop--], 0, state, 0, state.length);
            return time + deltaTime;
        }
        // GENERATE NEW POINTS
        while (Math.abs(mCurrentStep)>=mMinimumStepAllowed) {
            double h24 = mCurrentStep/24.0;
            double currentTime=time;
            double[] currentState=state;
            if (mMustRestart) {
                restartMethod(time, state);
                currentTime  = mTimes[0];
                currentState = mStates[0];
            }
            // Predictor: 4-steps Adams-Bashford
            mDerivatives[0] = mProblem.getDerivative(currentTime, currentState);
            for  (int i=0; i<state.length; i++) {
                mPredictorState[i] = currentState[i] + h24 * ( 55*mDerivatives[0][i] - 59*mDerivatives[1][i] + 37*mDerivatives[2][i] -9*mDerivatives[3][i]);
            }
            // Corrector: 3-steps Adams-Moulton 
            double[] derivativeIp1 = mProblem.getDerivative(currentTime+mCurrentStep, mPredictorState);
            for (int i=0; i<state.length; i++) {
                mCorrectorState[i] = currentState[i] + h24 * ( 9*derivativeIp1[i] + 19*mDerivatives[0][i] -5*mDerivatives[1][i] + mDerivatives[2][i]);
            }
            // CALCULATE ERROR
            double norm = 0;
            for (int i=0; i<state.length; i++) {
                double diffInIndex = mCorrectorState[i]-mPredictorState[i];
                norm += diffInIndex*diffInIndex;
            }
            norm = Math.sqrt(norm);
            double error = 19.0 * norm / 270.0;
            double maxErrorAllowed = mTolerance * Math.abs(mCurrentStep);
            if (error < maxErrorAllowed) { // ACCEPT
                time = currentTime + mCurrentStep;
                
                if (mMustRestart) { // return the first one of the 4, queue the rest.
                    System.arraycopy(mStates[sSTEPS - 2], 0, state, 0, dim);
                    for (int i = sSTEPS - 2; i >= 1; i--)
                        System.arraycopy(mStates[i-1], 0, queue[i], 0, dim);
                    System.arraycopy(mCorrectorState, 0, queue[0], 0, dim);
                    queueTop = sSTEPS-2;
                } else { // just return current point
                    System.arraycopy(mCorrectorState, 0, state, 0, state.length);
                }
                
                if (error < maxErrorAllowed*0.1) { // error is really small --> adapt step
                    if (norm < 1.0e-16) { // Prevent division by zero
                        mCurrentStep = 2 * mCurrentStep;
                    } 
                    else {
                        double q = 1.5*Math.pow(maxErrorAllowed/norm, 0.25);
                        q = Math.min(4, q); // Do not grow too much
                        //System.out.print ("ACCEPTED: t = "+time+ " Old step is "+mCurrentStep+ " error = "+error);
                        mCurrentStep *= q;
                        //System.out.println ("  New step is "+mCurrentStep+" state= "+state[0]);
                    }   
                    mMustRestart = true;
                } else {
                    //System.out.println ("ACCEPTED: t = "+time+ " with step "+mCurrentStep+ " error = "+error);
                    for (int i=mDerivatives.length-1; i>0; i--) // Prepare next step
                        System.arraycopy(mDerivatives[i-1],0,mDerivatives[i],0,state.length);
                    mMustRestart = false;
                }
                return time;
            }
            // Try a new smaller step
            double q = 1.5*Math.pow(maxErrorAllowed/norm, 0.25);
            q = Math.max(q, 0.1); // Do not shrink too much
            mCurrentStep *= q;
            //System.out.println ("REJECTED: t = "+time+ " New step is "+mCurrentStep+ " error = "+error);
            mMustRestart = true;
        }
        throw new ConvergenceException("Adaptative Predictor Corrector Method did not converge.");
//        // Was not able to reach tolerance before going below mMinimumStepAllowed
//        return Double.NaN; 
    }

    protected void restartMethod(double time, double[] state) {
        //System.out.println ("Restarting RK: t = "+time+ " with step "+mCurrentStep+" state= "+state[0]);
        mDerivatives[3] = mProblem.getDerivative(time, state);
        mTimes[2] = rungeKuttaStep(mCurrentStep, time, state, mStates[2], mDerivatives[3]);
        
        mDerivatives[2] = mProblem.getDerivative(mTimes[2], mStates[2]);
        mTimes[1] = rungeKuttaStep(mCurrentStep, mTimes[2], mStates[2], mStates[1], mDerivatives[2]);
        
        mDerivatives[1] = mProblem.getDerivative(mTimes[1], mStates[1]);
        mTimes[0] = rungeKuttaStep(mCurrentStep, mTimes[1], mStates[1], mStates[0], mDerivatives[1]);
        // Yes, one could write this using a for loop...
    }

    
    protected double rungeKuttaStep(double deltaTime, double time, double[] state, double[] newState, double[] k1) {
        double h2 = deltaTime/2.0;
        for (int i=0; i<state.length; i++) {
            mAuxState[i] = state[i] + h2 * k1[i];
        }
        double[] k2 = mProblem.getDerivative(time+h2, mAuxState);
        for (int i=0; i<state.length; i++) {
            mAuxState[i] = state[i] + h2 * k2[i];
        }
        double[] k3 = mProblem.getDerivative(time+h2, mAuxState);
        for (int i=0; i<state.length; i++) {
            mAuxState[i] = state[i] + deltaTime * k3[i];
        }
        double[] k4 = mProblem.getDerivative(time+deltaTime, mAuxState);
        double h6 = deltaTime/6;
        for (int i=0; i<state.length; i++) {
            newState[i] = state[i] + h6 * (k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }
        return time+deltaTime;
    }
    
}
