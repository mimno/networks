package edu.cornell.infosci.mimno;

public class UnicodeBarplot {

	public static final String[] BARS = { " ", "\u2581", "\u2582", "\u2583", "\u2584", "\u2585", "\u2586", "\u2587", "\u2588" };

	public static String getBar(double x, double min, double max) {
		if (x > max) { x = max; }
		if (x < min) { x = min; }
		return BARS[ (int) Math.round(8.0 * (x - min) / (max - min)) ];
	}

	public static String getBar(boolean b) {
		if (b) { return "\u2580"; }
		else { return "\u2584"; }
	}
	
	public static String getBars(double[] sequence) {
		double max = Double.NEGATIVE_INFINITY;
		double min = Double.POSITIVE_INFINITY;
		
		for (double x : sequence) {
			if (x > max) { max = x; }
			if (x < min) { min = x; }
		}

		return getBars(sequence, min, max);
	}

	public static String getBars(double[] sequence, double min, double max) {

		StringBuilder out = new StringBuilder();
		for (double x :sequence) {
			out.append(getBar(x, min, max));
		}
		return out.toString();
	}

	public static String getBars(int[] sequence) {
		double max = Double.NEGATIVE_INFINITY;
		double min = 0; //double min = Double.POSITIVE_INFINITY;
		
		for (int x : sequence) {
			if (x > max) { max = x; }
			//if (x < min) { min = x; }
		}

		StringBuilder out = new StringBuilder();
		for (int x :sequence) {
			out.append(getBar((double) x, min, max));
		}
		return out.toString();
	}

	public static void main (String[] args) {
		StringBuilder builder = new StringBuilder();
		for (int i = 0; i < 100; i++) {
			builder.append(UnicodeBarplot.getBar(i % 4 == 0));
		}
		System.out.println(builder);

		builder = new StringBuilder();
		for (int i = 0; i < 100; i++) {
			builder.append(UnicodeBarplot.getBar(i, 0, 100));
		}
		System.out.println(builder);
	}
}