package es.um.mned.ode;


public abstract class Event {
	
	private double tolerance;
	
	public Event(double tolerance) {
		this.tolerance = tolerance;
	}
	
	public Event() {
		this(0.);
	}

	public double getTolerance() {
		return tolerance;
	}

	public abstract double crossFunction(double time, double[] state);
	
	public abstract void crossAction(double time, double[] state);
	
	public abstract boolean stopCondition();

}
