package es.um.mned.interpolation;

public interface ExtendedStateFunction extends StateFunction {
    
    public double[] getDerivative(double time);

    public double getDerivative(double time, int index);
    
}
