package es.um.mned.scripts;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.problems.Taylor1DExample;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.DisplaySolution;

public class Taylor1DExampleScript {
	
	static private double sTo = 0;
    static private double sYo = 0.5;
	
	static public class Taylor1DExample2ndOrderSolver extends FixedStepMethod {

        public Taylor1DExample2ndOrderSolver(InitialValueProblem problem, double step) {
            super(problem, step);
        }
        
        
        @Override
        public int getOrder() {
        	return 2;
        }
    
        @Override
        public double doStep(double h, double t, double[] state) {
            double wi = state[0];
            state[0] = wi +  h * ( (1+h/2.0) * (wi-t*t+1) - h*t );
            return t + h;
        }
    }

    static public class Taylor1DExample4thdOrderSolver extends FixedStepMethod {

        public Taylor1DExample4thdOrderSolver(InitialValueProblem problem, double step) {
            super(problem, step);
        }
        
        
        @Override
        public int getOrder() {
        	return 4;
        }
    
        @Override
        public double doStep(double h, double t, double[] state) {
            double wi = state[0];
            state[0] = wi +  h * ( (1 + h/2. + h*h/6. + h*h*h/24.) * (wi-t*t) - h*t * (1 + h/3. + h*h/24.) 
                                 + (1 + h/2. - h*h/6. - h*h*h/24.)) ;
            return t + h;
        }
    }

    
    static public class TrueSol implements StateFunction {    
        public double[] getState(double t) {
            return new double[] { (t+1)*(t+1) - 0.5*Math.exp(t) };
        }
        public double getState(double t, int index) {
            switch (index) {
                case 0 : return (t+1)*(t+1) - 0.5*Math.exp(t);
                default : return Double.NaN;
            }
        }         
    }
	
	public static void main(String[] args) throws ConvergenceException {
        double maxTime = 2;
        double h = 0.1;
        
        Taylor1DExample problem = new Taylor1DExample(sTo, new double[] {sYo});
        TrueSol sol = new TrueSol();

        FixedStepMethod taylor2Method = new Taylor1DExample2ndOrderSolver(problem,h);
        FixedStepMethod taylor4Method = new Taylor1DExample4thdOrderSolver(problem,h);
        FixedStepMethod eulerMethod     = new FixedStepEulerMethod (problem,h);
        FixedStepMethod eulerHalfMethod = new FixedStepEulerMethod (problem,h/2);

        eulerMethod.solve(maxTime);
        eulerHalfMethod.solve(maxTime);
        taylor2Method.solve(maxTime);
        taylor4Method.solve(maxTime);
        NumericalSolution extrapolated = FixedStepEulerMethod.extrapolate(eulerMethod.getSolution(),eulerHalfMethod.getSolution());

        int[] indexes = new int[]{0};
        System.out.println ("Euler h        : Max Error = " + eulerMethod.getSolution().getMaxError(sol, indexes));
        System.out.println ("Euler h/2      : Max Error = " + eulerHalfMethod.getSolution().getMaxError(sol, indexes));
        System.out.println ("Extrapolation  : Max Error = " + extrapolated.getMaxError(sol, indexes));
        System.out.println ("Taylor order 2 : Max Error = " + taylor2Method.getSolution().getMaxError(sol, indexes));
        System.out.println ("Taylor order 4 : Max Error = " + taylor4Method.getSolution().getMaxError(sol, indexes));
        //DisplaySolution.listError(taylorMethod4.getSolution(), new TrueSol(), indexes,10);

        
        if (true) {
            DisplaySolution.timePlot(taylor4Method.getSolution());
        }
    }

}
