/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.ode;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.HashMap;

import es.um.mned.interpolation.HermiteInterpolator;
import es.um.mned.interpolation.StateFunction;
import es.um.mned.interpolation.Interpolator;

/**
 * Numerical Solution is a list of NumericalSolutionPoints (t,Y(t)), 
 * t is time (the independent variable), Y(t) is state at t.
 * 
 * @author F. Esquembre
 * @version September 2020
 */
public class NumericalSolution {
	private InitialValueProblem ivp;
    private ArrayList<NumericalSolutionPoint> pointList;
    private HashMap<Integer, StateFunction> interpolators;
    private int order;
   
//    /**
//     * Creates an empty NumericalSolution
//     * @param problem the InitialValueProblem being solved
//     */
//    public NumericalSolution() {
//        pointList = new ArrayList<>();
//    }
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
        for (int i=0,last=size-1; i<howMany; i++,last--) pointList.remove(last);
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
    
    public double[] getState(double t) throws Exception {
    	if(pointList.isEmpty() || t<pointList.get(0).getTime()
    			|| t>pointList.get(pointList.size()-1).getTime())
    		throw new Exception("Evaluation out of bounds.");
    	
    	int l = 0, r = pointList.size()-1;
    	while(l <= r) {
    		int m = (l+r) / 2;
    		if(pointList.get(m).getTime() <= t) l = m+1;
    		else r = m-1;
    	} // doesn't matter if t == r due to the way I select points
    	
    	if(!interpolators.containsKey(l))
            putInterpolator(l);
		return interpolators.get(l).getState(t);
    }

   public double getState(double t, int index) throws Exception {
	   if(pointList.isEmpty() || t<pointList.get(0).getTime()
   			|| t>pointList.get(pointList.size()-1).getTime())
   		throw new Exception("Evaluation out of bounds.");
   	
	   	int l = 0, r = pointList.size()-1;
	   	while(l <= r) {
	   		int m = (l+r) / 2;
	   		if(pointList.get(m).getTime() <= t) l = m+1;
	   		else r = m-1;
	   	} // doesn't matter if t == r due to the way I select points
	   	
	   	if(!interpolators.containsKey(l))
	           putInterpolator(l);
			return interpolators.get(l).getState(t)[index];
   }
    
}
