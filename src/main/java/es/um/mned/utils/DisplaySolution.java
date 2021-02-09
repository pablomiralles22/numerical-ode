/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.utils;

import java.util.Iterator;
import java.util.stream.IntStream;

import org.opensourcephysics.frames.PlotFrame; 
import org.opensourcephysics.controls.AbstractAnimation;
import org.opensourcephysics.controls.AnimationControl;
import javax.swing.JFrame;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.ode.NumericalSolutionPoint;
import org.opensourcephysics.display.Dataset;


public class DisplaySolution {
	private static final int MAX_POINTS = 10000; // Estimated empirically, might change

	
	/*
	 * ===================================================
	 * Lists
	 * ===================================================
	 */

    static public void list(NumericalSolution solution, int[] indexes) {
        Iterator<NumericalSolutionPoint> iterator = solution.iterator();
        while (iterator.hasNext()) {
            NumericalSolutionPoint point = iterator.next();
            System.out.print("time="+point.getTime());
            for (int i=0; i<indexes.length;i++) {
                System.out.println ("x["+i+"] = "+point.getState(indexes[i]));
            }
            System.out.println();
        }
    }

    /*
	 * ===================================================
	 * Time plotting
	 * ===================================================
	 */
    
    static public void timePlot(NumericalSolution solution) {
    	int skip = Math.max(solution.getSize() / MAX_POINTS - 1, 0);
        timePlot(solution, skip);
    }
    
    static public void timePlot(NumericalSolution solution, int skip) {
        int dimension = solution.get(0).getState().length;
        int[] indexes = IntStream.range(0, dimension).toArray();
        timePlot(solution, indexes,skip);
    }
    
    static public void timePlot(NumericalSolution solution, int[] indexes) {
    	int skip = Math.max(solution.getSize() / MAX_POINTS - 1, 0);
        timePlot(solution, indexes, skip);
    }
    
    static public void timePlot(NumericalSolution solution, int[] indexes, int skip) {
        PlotFrame frame = new PlotFrame ("time" , "x[*]" , "Time plot frame") ;
        frame.setConnected(true); // sets default to connect dataset points
        frame.setSize(800,600);
        for (int i=0; i<indexes.length;i++) {
            frame.setMarkerShape(i, Dataset.NO_MARKER);
            //frame.setMarkerColor(0,java.awt.Color.RED);
            frame.setXYColumnNames (i , "time" ,"x["+i+"]") ; // sets names for each dataset
            
            for(int j=0; j < solution.getSize(); j += skip+1) {
            	NumericalSolutionPoint point = solution.get(j);
            	frame.append(i, point.getTime(), point.getState(indexes[i]));
            }
        }
        frame.setVisible ( true ) ;
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;
    }
    
    /*
	 * ===================================================
	 * State plotting
	 * ===================================================
	 */
    
    static public void statePlot(NumericalSolution solution, int index1, int index2) {
    	int skip = Math.max(solution.getSize() / MAX_POINTS - 1, 0);
        statePlot(solution, index1, index2, skip);
    }
    
    static public void statePlot(NumericalSolution solution, int index1, int index2, int skip) {
        PlotFrame frame = new PlotFrame ("x["+index1+"]" , "x["+index2+"]"  , "State plot frame") ;
        frame.setConnected(true); // sets default to connect dataset points
        frame.setSize(800,600);
        frame.setMarkerShape(0, Dataset.NO_MARKER);
        frame.setMarkerColor(0,java.awt.Color.BLUE);
        frame.setXYColumnNames (0 , "x["+index1+"]" , "x["+index2+"]") ; // sets names for each dataset
        
        for(int j=0; j < solution.getSize(); j += skip+1) {
        	NumericalSolutionPoint point = solution.get(j);
        	frame.append(0, point.getState(index1), point.getState(index2));
        }
        
        frame.setVisible ( true ) ;
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;
    }
    
    /*
	 * ===================================================
	 * Animations
	 * ===================================================
	 */
    
    static public void animateTimePlot(NumericalSolution solution, int[] indexes) {
        AnimationControl.createApp(new AnimateSequence(solution,indexes));
    }
    
    static private class AnimateSequence extends AbstractAnimation {
        NumericalSolution mSolution;
        int[] mIndexes;
        PlotFrame mFrame;
        Iterator<NumericalSolutionPoint> mIterator;

        
        public AnimateSequence(NumericalSolution solution, int[] indexes) {
            mSolution = solution;
            mIndexes = indexes;
            mFrame = new PlotFrame ("time" , "x[*]" , "Time plot animation") ;
            mFrame.setConnected(true); // sets default to connect dataset points
            mFrame.setSize(800,600);
            for (int i=0; i<indexes.length;i++) {
                mFrame.setXYColumnNames (i , "time" ,"x["+i+"]") ; // sets names for each dataset
            }
            mFrame.setVisible ( true ) ;
            mFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;            
        }
        public void initializeAnimation () {
//            mIndex = 0;
//            mFrame.append(0, mIndex, mSequence[mIndex]);
            mIterator = mSolution.iterator();
            doStep();

        }
        
        public void resetAnimation () {
            mIterator = mSolution.iterator();
            mFrame.clearData();
            super.resetAnimation();
        }

        public void doStep () { 
            if (!mIterator.hasNext()) super.stopAnimation();
            else {
                NumericalSolutionPoint point = mIterator.next();
                for (int i=0; i<mIndexes.length;i++) {
                    mFrame.append(i, point.getTime(), point.getState(mIndexes[i]));
                }
                mFrame.repaint();
                try { Thread.sleep (500) ; } catch(InterruptedException ie) {}
            }
        }

    }
    
    
//  static public void listError(NumericalSolution solution, 
//                               StateFunction trueSolution,
//                               int[] indexes) {
//      listError(solution.iterator(),trueSolution, indexes);
//  }
//  
//  static public void listError(NumericalSolution solution, 
//                               StateFunction trueSolution,
//                               int[] indexes, int numberOfPoints) {
//      listError(solution.iterator(numberOfPoints),trueSolution, indexes);
//  }
//
//  static public double maximumError(NumericalSolution solution, 
//                                  StateFunction trueSolution,
//                                  int[] indexes) {
//      return maximumError(solution.iterator(),trueSolution, indexes);
//  }
//  
//  static private double maximumError(Iterator<NumericalSolutionPoint> iterator, StateFunction trueSolution, int[] indexes) {
//      double maxError = 0.0;
//      while (iterator.hasNext()) {
//          NumericalSolutionPoint point = iterator.next();
//          double[] trueState = trueSolution.getState(point.getTime());
//          double error = 0;
//          for (int i=0; i<indexes.length;i++) {
//              error += Math.abs(point.getState(indexes[i])-trueState[indexes[i]]);
//          }
//          maxError = Math.max(maxError, error);
//      }
//      return maxError;
//  }
//
//  static private double listError(Iterator<NumericalSolutionPoint> iterator, StateFunction trueSolution, int[] indexes) {
//      double maxError = 0.0;
//      while (iterator.hasNext()) {
//          NumericalSolutionPoint point = iterator.next();
//          double[] trueState = trueSolution.getState(point.getTime());
//          System.out.print("time="+point.getTime());
//          double error = 0;
//          for (int i=0; i<indexes.length;i++) {
//              System.out.print(", x["+indexes[i]+"]="+point.getState(indexes[i])+", vs "+trueState[indexes[i]]);
//              error += Math.abs(point.getState(indexes[i])-trueState[indexes[i]]);
//          }
//          maxError = Math.max(maxError, error);
//          System.out.println(", error = "+error);
//      }
//      return maxError;
//  }

}
