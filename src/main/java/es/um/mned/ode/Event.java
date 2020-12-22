package es.um.mned.ode;


public abstract class Event {
	
	private boolean blocking;
	private double tolerance;
	
	public Event(boolean blocking, double tolerance) {
		this.blocking = blocking;
		this.tolerance = tolerance;
	}

	public boolean isBlocking() {
		return blocking;
	}

	public double getTolerance() {
		return tolerance;
	}

	public abstract double crossFunction(double time, double[] state);
	
	public abstract void crossAction(double time, double[] state);

}
