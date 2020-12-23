package es.um.mned.scripts;

import es.um.mned.methods.*;
import es.um.mned.ode.*;
import es.um.mned.problems.Rigid1D;

public class ImplicitMethodsRigid1D {

	public static void main(String[] args) throws Exception {
		// Params
		double maxTime = 8;
		double tolerance = 1e-10;
		double hStep = 1e-3;
		
		double t0 = 0.;
		double[] x0 = new double[] {-1};
		
		// Problem
        ExtendedInitialValueProblem problem = new Rigid1D(t0, x0);
        
        // Methods
//        FixedStepBackwardsEulerNewton1DMethod method = new FixedStepBackwardsEulerNewton1DMethod(problem,hStep);
//        FixedStepTrapezoidalNewton1DMethod method = new FixedStepTrapezoidalNewton1DMethod(problem, hStep);
        FixedStepBDFNewton1DMethod method = new FixedStepBDFNewton1DMethod(problem, 3, hStep);
//        FixedStepMethod method = new FixedStepPredictorCorrector4Method(problem,10);
//        FixedStepMethod method = new AdaptiveStepRKFehlbergMethod(problem,hStep, tolerance);
        
        method.setNewtonTolerance(tolerance);
        
        try {
			method.solve(maxTime);
		} catch (ConvergenceException e) {
			e.printStackTrace();
			System.exit(0);
		}

        // Get max error using analytical solution
        double err = method.getSolution().getMaxError(new Rigid1D.TrueSol());
        System.out.println("Error = " + err);
        // Plot steps if adaptative
//        if (method instanceof AdaptiveStepMethod) DisplaySequence.plot(((AdaptiveStepMethod) method).getSolution().getStepList());
        // Print number of evaluations
        System.out.println ("Evaluations = "+problem.getEvaluationCounter());
	}
}
