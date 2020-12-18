/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.ode;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.interpolation.Interpolator;

/**
 * Numerical Solution is a list of NumericalSolutionPoints (t,Y(t)), 
 * t is time (the independent variable), Y(t) is state at t.
 * 
 * @author F. Esquembre
 * @version September 2020
 */
public class NumericalSolution implements StateFunction{
	private InitialValueProblem ivp;
    private ArrayList<NumericalSolutionPoint> pointList;
    private HashMap<Integer, StateFunction> interpolators;
    private int order;
    
    /*
     * Creates another solution for extrapolation purposes.
     * The new solution uses the same problem, mind it when
     * checking the number of evaluations.
     */
    public static NumericalSolution createExtrapolationSol(NumericalSolution other) {
    	return new NumericalSolution(other.ivp, other.order + 1);
    }
   

    /**
     * Creates a NumericalSolution with the initial condition as first point
     * @param problem the InitialValueProblem being solved
     */
    public NumericalSolution(InitialValueProblem problem, int order) {
    	ivp = problem;
    	interpolators = new HashMap<>();
        pointList = new ArrayList<>();
        pointList.add(new NumericalSolutionPoint(problem.getInitialTime(),problem.getInitialState()));
        this.order = order;
    }

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
    
    public void removeLast(int howMany) {
        int size = pointList.size();
        howMany = Math.min(howMany, size);
        for (int i=0,last=size-1; i<howMany; i++,last--) {
            pointList.remove(last);
        }
    }
    
    /**
     * Returns an iterator to iterate through the whole list of points in the solution
     * @return 
     */
    public Iterator<NumericalSolutionPoint> iterator() {
        return pointList.iterator();
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
    
    public int getSize() {
    	return pointList.size();
    }
    
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
    
    public StateFunction getInterpolator(int index) {
    	if(index < 0 || index > pointList.size())
    		return null;
    	
    	if(!interpolators.containsKey(index))
    		putInterpolator(index);
    	
    	return interpolators.get(index);
    }
    
    public double[] getState(double t) {
    	if(pointList.isEmpty() || t<pointList.get(0).getTime()
    			|| t>pointList.get(pointList.size()-1).getTime())
    		System.err.println("Evaluation out of bounds, precision not guaranteed.");
    	
    	int l = 1, r = pointList.size()-1;
    	while(l <= r) {
    		int m = (l+r) / 2;
    		if(pointList.get(m).getTime() <= t) l = m+1;
    		else r = m-1;
    	}

        l = Math.min(l, pointList.size()-1);

		return getInterpolator(l).getState(t);
    }

   public double getState(double t, int index) {
	   if(pointList.isEmpty() || t<pointList.get(0).getTime()
   			|| t>pointList.get(pointList.size()-1).getTime())
    		System.err.println("Evaluation out of bounds, precision not guaranteed.");
   	
	   	int l = 1, r = pointList.size()-1;
	   	while(l <= r) {
	   		int m = (l+r) / 2;
	   		if(pointList.get(m).getTime() <= t) l = m+1;
	   		else r = m-1;
	   	}

        l = Math.min(l, pointList.size()-1);
	   	
	   	return getInterpolator(l).getState(t, index);
   }
   
    /*
     * Returns max error given analytical solution
     */
    public double getMaxError(StateFunction analyticalSolution) {
    	double err = 0.0;
    	for(NumericalSolutionPoint p : pointList) {
    		double currentErr = 0.0;
    		for(int i=0; i < p.getState().length; ++i) {
    			double diff = (p.getState(i) - analyticalSolution.getState(p.getTime(), i));
    			currentErr += diff*diff;
    		}
    		err = Math.max(err, currentErr);
    	}
    	return Math.sqrt(err);
    }
    
    public ArrayList<Double> getStepList() {
    	ArrayList<Double> stepList = new ArrayList<>(pointList.size()-1);
    	for(int i=0; i+1 < pointList.size(); ++i)
    		stepList.add(pointList.get(i+1).getTime() - pointList.get(i).getTime());
    	return stepList;
    }
}
