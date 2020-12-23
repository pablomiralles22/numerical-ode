package es.um.mned.interpolation;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

/**
 *
 * @author paco
 */
public class Interpolator implements StateFunction {

    private double[] r;
    private double[][] coeffs;
    private int dim;
    private int n;

    public Interpolator(Map<Double, double[][]> m) {

        dim = m.values().iterator().next()[0].length;
        n = 0;
        for(double[][] al : m.values()) {
            n += al.length;
        }

        r = new double[n];
        {
            int i=0;
            for(Double point : m.keySet()) {
                for(int k=0; k < m.get(point).length; ++k)
                    r[i++] = point.doubleValue();
            }
        }

        double[][][] dd = new double[n][n][dim];

        for(int i=0; i<n; ++i) {
            dd[0][i] = m.get(r[i])[0];
        }

        coeffs = new double[n][dim];
        double factorial = 1;
        for(int i=1; i<n; ++i) {
            factorial *= i;
            for(int j=i; j<n; ++j) {
                if(r[j] == r[j-i])
                    for(int h=0; h<dim; ++h)
                        dd[i][j][h] = m.get(r[j])[i][h] / factorial;
                else
                    for(int h=0; h<dim; ++h)
                        dd[i][j][h] = (dd[i-1][j][h]-dd[i-1][j-1][h]) / (r[j]-r[j-i]);
            }
        }

        for(int i=0; i<n; ++i)
        	System.arraycopy(dd[i][i], 0, coeffs[i], 0, dim);
//            coeffs[i] = Arrays.copyOf(dd[i][i], dim);
    }

    public double getState(double t, int index) {
        double ans = coeffs[0][index];
        double aux = 1.;
        for(int i=1; i<n; ++i) {
            aux *= (t - r[i-1]);
            ans += aux * coeffs[i][index];
        }
        return ans;
    }

    public double[] getState(double t) {
        double[] ans = Arrays.copyOf(coeffs[0], dim);
        double aux = 1.;
        
        for(int i=1; i<n; ++i) {
        	aux *= (t - r[i-1]);
        	for(int h=0; h<dim; ++h)
                ans[h] += aux * coeffs[i][h];
        }
        return ans;
    }
    
    public static void main(String[] args) {
    	double[][] l1, l2, l3;
    	l1 = new double[][] {new double[] {2}, new double[] {-8}, new double[] {56}};
    	l2 = new double[][] {new double[] {1}};
    	l3 = new double[][] {new double[] {2}, new double[] {8}};
    	HashMap<Double, double[][]> hm = new HashMap<>();
    	hm.put(-1., l1);
    	hm.put(0., l2);
    	hm.put(1., l3);
    	Interpolator inter = new Interpolator(hm);
    	System.out.println(inter.getState(3, 0));
    	System.out.println(inter.n);
    }

}
