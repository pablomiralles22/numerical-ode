/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.utils;

import java.util.ArrayList;
import java.util.Iterator;
import org.opensourcephysics.frames.PlotFrame; 
import org.opensourcephysics.controls.AbstractAnimation;
import javax.swing.JFrame;
import org.opensourcephysics.controls.AnimationControl;
import org.opensourcephysics.display.Dataset;

/**
 *
 * @author F. Esquembre
 * @version September 2020
 */
public class DisplaySequence {
    
    static public void list(double[] x) {
        for (int  i=0; i<x.length; i++) System.out.println ("x["+i+"] = "+x[i]);
    }

    static public void plot(ArrayList<Double> list) {
        PlotFrame frame = new PlotFrame ("n" , "x[n]" , "Sequence plot frame") ;
        frame.setSize(800,400);
        frame.setConnected(true); // sets default to connect dataset points
        frame.setMarkerShape(0, Dataset.NO_MARKER);        frame.setMarkerSize(0, 5);
        frame.setMarkerColor(0,java.awt.Color.RED);
        frame.setXYColumnNames (0 , "n" ,"x[n]") ; // sets nammes for first dataset
        Iterator<Double> iterator = list.iterator();
        int i=0;
        while (iterator.hasNext()) {
            frame.append(0, i, iterator.next().doubleValue());
            i++;
        }
        frame.setVisible ( true ) ;
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;
    }
    
    static public void plot(double[] x) {
        PlotFrame frame = new PlotFrame ("n" , "x[n]" , "Sequence plot frame") ;
        frame.setSize(800,400);
        frame.setConnected(false); // sets default to connect dataset points
        frame.setMarkerSize(0, 5);
        frame.setMarkerColor(0,java.awt.Color.RED);
        frame.setXYColumnNames (0 , "n" ,"x[n]") ; // sets nammes for first dataset
        for (int i=0; i<x.length; i++) {
            frame.append(0, i, x[i]);
        }
        frame.setVisible ( true ) ;
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;
    }
    
    static public void animate(double[] seq) {
        AnimationControl.createApp(new AnimateSequence(seq));
    }
    
    static private class AnimateSequence extends AbstractAnimation {
        double[] mSequence;
        int mIndex;
        PlotFrame mFrame;
        
        public AnimateSequence(double[] x) {
            mSequence = x;
            mIndex = 0;
            mFrame = new PlotFrame ("n" , "x[n]" , "Sequence animation") ;
            mFrame.setConnected(false); // sets default to connect dataset points
            mFrame.setSize(800,400);
            mFrame.setMarkerSize(0, 5);
            mFrame.setMarkerColor(0,java.awt.Color.RED);
            mFrame.setXYColumnNames (0 , "n" ,"x[n]") ; // sets nammes for first dataset
            mFrame.setVisible ( true ) ;
            mFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;            
        }
        public void initializeAnimation () {
            mIndex = 0;
            mFrame.append(0, mIndex, mSequence[mIndex]);

        }
        
        public void resetAnimation () {
            mIndex = 0;
            mFrame.clearData();
            super.resetAnimation();
        }

        public void doStep () { 
            mIndex++;
            control.println("Index = "+mIndex);
            if (mIndex>=mSequence.length) super.stopAnimation();
            else {
                mFrame.append(0, mIndex, mSequence[mIndex]);
                mFrame.repaint();
                try { Thread.sleep (500) ; } catch(InterruptedException ie) {}

            }
        }
    }

}
