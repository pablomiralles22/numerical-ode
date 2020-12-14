/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.interpolation;

/**
 *
 * @author paco
 */
public interface ExtendedStateFunction extends StateFunction {
    
    public double[] getDerivative(double time);

    public double getDerivative(double time, int index);
    
}
