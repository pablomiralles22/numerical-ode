package es.um.mned.scripts;

import es.um.mned.methods.*;
import es.um.mned.ode.*;
import es.um.mned.tools.*;
import es.um.mned.problems.Rigid1D;

public class ImplicitMethodsRigid1D {

	public static void main(String[] args) {
		double maxTime = 8;
		double tolerance = 1e-10;
        ExtendedInitialValueProblem problem = new Rigid1D();

        double hStep = 1e-5;
        FixedStepMethod method = new BackwardsEulerNewton1DMethod(problem,hStep, tolerance);
        // FixedStepMethod method = new TrapezoidalNewton1DMethod(problem, hStep, tolerance);
        //method = new FixedStepPredictorCorrector4Method(problem,10);
        //method = new AdaptiveStepRKFehlbergMethod(problem,hStep, tolerance);
        
        double lastTime = method.solve(maxTime);
        if (Double.isNaN(lastTime)) {
            System.out.println ("Method broke!");
        }
        // DisplaySolution.listError(method.getSolution(), new Rigid1D.TrueSol(), new int[]{0});
        double err = method.getSolution().getMaxError(new Rigid1D.TrueSol());
        System.out.println("Error = " + err);
        if (method instanceof AdaptiveStepMethod) DisplaySequence.plot(((AdaptiveStepMethod) method).getStepList());
        System.out.println ("Evaluations = "+method.getEvaluationCounter());

	}
}
