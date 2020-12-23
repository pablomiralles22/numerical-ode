package es.um.mned.ode;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.function.Consumer;
import java.util.stream.IntStream;
import java.util.HashMap;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.interpolation.Interpolator;

public class NumericalSolution implements StateFunction{
	private InitialValueProblem ivp;
    private ArrayList<NumericalSolutionPoint> pointList;
    private HashMap<Integer, StateFunction> interpolators;
    private int order;
    
    /*
     * ==================================================
     * Constructors
     * ==================================================
     */
    
    /**
     * Creates another solution for extrapolation purposes.
     * The new solution uses the same problem, mind it when
     * checking the number of evaluations.
     * @param other another numerical solution
     * @return new empty solution ready to fill
     */
    public static NumericalSolution createExtrapolationSol(NumericalSolution other) {
    	NumericalSolution extrapolatedSol = new NumericalSolution();
    	extrapolatedSol.ivp = other.ivp;
    	extrapolatedSol.order = other.order + 1;
    	return extrapolatedSol;
    }
   
	/**
	 * Private constructor for the previous one
	 */
	private NumericalSolution() {
		pointList = new ArrayList<>();
		interpolators = new HashMap<>();
	}
    
    /**
     * Creates a NumericalSolution with the initial condition as first point
     * @param problem the InitialValueProblem being solved
     */
    public NumericalSolution(InitialValueProblem problem, int order) {
    	this();
    	ivp = problem;
    	interpolators = new HashMap<>();
        pointList.add(new NumericalSolutionPoint(problem.getInitialTime(),problem.getInitialState()));
        this.order = order;
    }

    /*
     * ==================================================
     * Methods to access and modify
     * ==================================================
     */
    
	/**
     * Adds a solution point
     * 
     * @param time the time of the point to add: t
     * @param state the state: Y(t)
     * @return  the newly added NumericalSolutionPoint if successful, null if failed
     */
    public NumericalSolutionPoint add(double time, double[] state) {
        NumericalSolutionPoint point = new NumericalSolutionPoint(time,state);
        if (pointList.add(point)) return point;
        return null;
    }

    /**
     * @param index
     * @return point at position index
     */
    public NumericalSolutionPoint get(int index) {
    	return pointList.get(index);
    }
    
    /**
     * Gets the last point of the solution
     * @return 
     */
    public NumericalSolutionPoint getLastPoint() { 
        return pointList.get(pointList.size()-1);
    }
    
    /**
     * Removes last howMany elements from the solution.
     * I've commented it because I don't use it and it does not
     * play well with how I handle interpolation
     * @param howMany number of elements to remove
     */
//    public void removeLast(int howMany) {
//        int size = pointList.size();
//        howMany = Math.min(howMany, size);
//        for (int i=0,last=size-1; i<howMany; i++,last--) {
//        	interpolators.remove(last);
//            pointList.remove(last);
//        }
//    }
    
    /**
     * Returns an iterator to iterate through the whole list of points in the solution
     * @return 
     */
    public Iterator<NumericalSolutionPoint> iterator() {
        return pointList.iterator();
    }
    
    /**
     * Exports the ability to do a forEach
     * @param c
     */
    public void forEach(Consumer<NumericalSolutionPoint> c) {
    	pointList.forEach(c);
    }
    
    /**
     * Returns an iterator to iterate through the last few points in the solution
     * @param numberOfPoints int, number of points at the end
     * @return 
     */
    public Iterator<NumericalSolutionPoint> iterator(int numberOfPoints) {
        int size = pointList.size();
        return pointList.subList(size-numberOfPoints, size).iterator();
    }
    
    /**
     * @return size of the list of points calculated
     */
    public int getSize() {
    	return pointList.size();
    }
    
    /**
     * Creates an interpolator for the range [t_{index-1},t_index]
     * @param index
     */
    private void putInterpolator(int index) {
        // in case we don't have enough points. The worst case scenerio
        // in this library is having 1 interpolator with worse precision
        // in the very first interval, and only in methods with order>=4
    	int nPoints = Math.min((order+2)/2, pointList.size());
        int start;
        if(index - (nPoints+1)/2 < 0) 
            start = 0;
        else if(index + nPoints/2 >= pointList.size())
            start = pointList.size() - nPoints;
        else
            start = index - (nPoints+1)/2;
        
        HashMap<Double, double[][]> m = new HashMap<>();
        for(int i=start, currentOrder=-1; i<start+nPoints; ++i) {
            NumericalSolutionPoint p = pointList.get(i);
            double time = p.getTime();
            double[] state = p.getState();
            if(currentOrder + 1 < order) {
            	m.put(time, new double[][] {state, ivp.getDerivative(time,state)});
            	currentOrder += 2;
            } else {
            	m.put(time, new double[][] {state});
            	currentOrder += 1;
            }
        }
        interpolators.put(index, new Interpolator(m));
    }
    
    /**
     * Gets the interpolator for range [t_{index-1},t_index], and creates
     * it if it does not exist.
     * @param index
     * @return
     */
    public StateFunction getInterpolator(int index) {
    	if(index < 0 || index > pointList.size())
    		return null;
    	
    	if(!interpolators.containsKey(index))
    		putInterpolator(index);
    	
    	return interpolators.get(index);
    }
    
    /*
     * ============================================================
     * State Function methods
     * ============================================================
     */
    
	/**
	 * Interpolates to get the value for t \in [a,b]
	 */
	public double[] getState(double t) {
		if (pointList.isEmpty() || t < pointList.get(0).getTime() || t > pointList.get(pointList.size() - 1).getTime())
			System.err.println("Evaluation out of bounds, precision not guaranteed.");

		int l = 1, r = pointList.size() - 1;
		while (l <= r) {
			int m = (l + r) / 2;
			if (pointList.get(m).getTime() <= t)
				l = m + 1;
			else
				r = m - 1;
		}

		l = Math.min(l, pointList.size() - 1);

		return getInterpolator(l).getState(t);
	}

	/**
	 * Interpolates to get the value of position index for t \in [a,b]
	 */
	public double getState(double t, int index) {
		if (pointList.isEmpty() || t < pointList.get(0).getTime() || t > pointList.get(pointList.size() - 1).getTime())
			throw new IllegalArgumentException("t is not the defined domain for this solution");

		int l = 1, r = pointList.size() - 1;
		while (l <= r) {
			int m = (l + r) / 2;
			if (pointList.get(m).getTime() <= t)
				l = m + 1;
			else
				r = m - 1;
		}

		l = Math.min(l, pointList.size() - 1);

		return getInterpolator(l).getState(t, index);
	}
   
	/*
	* ============================================================
	* Different utils for a numerical solution
	* ============================================================
	*/
   
	/**
	 * Returns the maximum error in the list of points comparing only
	 * certain indexes
	 * @param analyticalSolution real solution of the problem
	 * @param indexList lists of indexes to compare
	 * @return
	 */
	public double getMaxError(StateFunction analyticalSolution, int[] indexList) {
		double err = pointList.stream()
				.map(p -> {
					double aux = 0.0;
					double[] state = p.getState();
					double t = p.getTime();
					for (int i = 0; i < indexList.length; ++i) {
						double diff = (state[indexList[i]] - analyticalSolution.getState(t, indexList[i]));
						aux += diff * diff;
					}
					return aux;
				})
				.max(Double::compare)
				.orElse(0.0);
		return Math.sqrt(err);
	}
   
    /**
     * Returns the maximum error in the list of points
     * @param analyticalSolution
     * @return
     */
    public double getMaxError(StateFunction analyticalSolution) {
    	if(pointList.isEmpty()) return 0.;
    	
    	int dim = pointList.get(0).getState().length;
    	int[] indexes = IntStream.range(0, dim).toArray();
    	
    	return getMaxError(analyticalSolution, indexes);
    }
    
    /**
     * @return list of steps used
     */
    public ArrayList<Double> getStepList() {
    	ArrayList<Double> stepList = new ArrayList<>(pointList.size()-1);
    	for(int i=0; i+1 < pointList.size(); ++i)
    		stepList.add(pointList.get(i+1).getTime() - pointList.get(i).getTime());
    	return stepList;
    }
}
