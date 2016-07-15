package util;

import java.util.ArrayList;

public class MatchedFilter {

	public static double[][] getAvg(ArrayList<double[]> observations){
		if(observations.isEmpty()) return null;
		double[] avg = new double[observations.get(0).length];
		for(double[] o : observations){
			for(int i=0;i<avg.length;i++){
				avg[i] += o[i];
			}
		}
		float sum = 0;
		for(double v : avg){
			sum += v;
		}
		for(int i=0;i<avg.length;i++){
			avg[i] /= Math.abs(sum);
		}
		
		double[][] ret = new double[avg.length][1];
		for(int i=0; i<avg.length;i++){
			ret[i][0] = avg[i];
		}
		
		return ret;
	}
	
	
	public static double[] getFilter(double[][] avgedSignal, ArrayList<double[]> noise){	
		double[][] R = getCovMatrix(noise);
		double[][] Rinv = MC.invert(R);
		double[][] Rinvs = MC.multiply(Rinv, avgedSignal);
		double div = Math.sqrt(MC.multiply(MC.transpose(avgedSignal),Rinvs)[0][0]);
		//System.out.println(div);
		double[][] h = MC.multiply(Rinvs, 1/div);
		double[] filter = new double[h.length];
		for(int i=0;i<filter.length;i++)
			filter[i] = h[i][0];
		return filter;
	}
	
	
	private static double[][] getCovMatrix(ArrayList<double[]> noise){
		double num = 0;
		double[][] cov = new double[noise.get(0).length][noise.get(0).length];
		//double[][] avg = getAvg(noise);
		for(double[] v : noise){
			double[][] arg = new double[v.length][1];
			for(int i=0;i<v.length;i++){
				arg[i][0] = v[i];
			}
			num ++;
			cov = MC.sum(cov, MC.multiply(arg, MC.transpose(arg)),1);
		}
		cov = MC.multiply(cov, 1/(num-1));
		return cov;
	}
	
}
