package es.um.mned.interpolation;

/**
 *
 * @author paco
 */
public interface StateFunction {
    
    public double[] getState(double time);

    public double getState(double time, int index);
    
}
