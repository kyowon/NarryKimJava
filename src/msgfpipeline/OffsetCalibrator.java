package msgfpipeline;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction.Parametric;
import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;

import msgf.Tolerance;
import msgfpipeline.parser.MSGFPlusParser.MSGFPlusPSM;

public class OffsetCalibrator { // just recalibrate precursor mz
	private ArrayList<MSGFPlusPSM> allPSMs = null;//new ArrayList<MSGFPlusPSM>();
	private ArrayList<MSGFPlusPSM> filteredPSMs = null;
	private ArrayList<Integer> hqPSMIndices = null;
	private HashMap<Integer, Double> offsetMap = null;
	private Tolerance tolerance = null;
	
	public OffsetCalibrator(ArrayList<MSGFPlusPSM> allPSMs, Tolerance tolerance, double qvalueThreshold){
		this.allPSMs = allPSMs;
		Collections.sort(this.allPSMs);
		this.tolerance = tolerance;
		this.hqPSMIndices = new ArrayList<Integer>();
		int n = 0;
		for(int i=0; i<allPSMs.size();i++){
			MSGFPlusPSM psm = allPSMs.get(i);
			if(psm.getQvalue() > qvalueThreshold) continue;
			hqPSMIndices.add(i);
			//if(psm.getPrecursorErr() > tolerance.getToleranceAsDa(psm.getCharge() * psm.getPrecursor())) continue;
			n++;
		}
		System.out.println("Qualified psms before calibration : " + n);
		if(n < 5) return;
		calculateOffsets();
		System.out.println("Calculating offsets done");
		updatePrecursorMzs();
		System.out.println("Precursor recalibration done");
	}
	
	private void calculateOffsets(){
		offsetMap = new HashMap<Integer, Double>();
		offsetMap.put(allPSMs.get(0).getScanNumber(), .0);
		offsetMap.put(allPSMs.get(allPSMs.size()-1).getScanNumber(), .0);
		double[] x = new double[hqPSMIndices.size()]; 
		double[] y = new double[x.length];
		double preMax = Double.NEGATIVE_INFINITY;
		for(int i=0;i<hqPSMIndices.size();i++){
			MSGFPlusPSM c = allPSMs.get(hqPSMIndices.get(i));
			x[i] = Math.max(preMax + .0001, c.getScanNumber());
			preMax = preMax < x[i]? x[i] : preMax;
			
			//if(i>0 && x[i] == preMax) x[i] = x[i-1] + .00001;
			y[i] = c.getPrecursorErr();
		}
		LoessInterpolator interpolator = new LoessInterpolator(.1, 10);
		double[] z = interpolator.smooth(x, y);
		for(int i=0;i<x.length;i++){
			//System.out.println(x[i] +  " " + y[i] + " " + z[i] );
			offsetMap.put((int)x[i], Double.isNaN(z[i]) ? y[i] : z[i]);
		}
		
		
	}
	
	private void updatePrecursorMzs(){
		double[] coeff = new double[2];
		filteredPSMs = new ArrayList<MSGFPlusPSM>();
		
		int currentIndex = 0;
		
		for(int i=0;i<allPSMs.size();i++){
			//System.out.println(i + " " + allPSMs.size());
			MSGFPlusPSM psm = allPSMs.get(i);
			CurveFitter<Parametric> fitter = new CurveFitter<Parametric>(new LevenbergMarquardtOptimizer());
			if(currentIndex < hqPSMIndices.size() && i > hqPSMIndices.get(currentIndex)) currentIndex ++;
			for(int j=-1;j<=0;j++){
				int k = currentIndex + j;
				int c;
				if(k<0){
					c = 0;
					j = 100;
				}
				else if(k>=hqPSMIndices.size()){
					c = allPSMs.size()-1;
					j = 100;
				}
				else{
					c = hqPSMIndices.get(k);
				}
				//System.out.println(c + " " + offsetMap.get(allPSMs.get(c).getScanNumber()));
				//System.out.println(n++);
				fitter.addObservedPoint(allPSMs.get(c).getScanNumber(), offsetMap.get(allPSMs.get(c).getScanNumber()));
			}
			//System.out.println("hh");
			coeff = fitter.fit(new PolynomialFunction.Parametric(), new double[coeff.length]);
			//System.out.println(coeff[0] + " " + coeff[1]);
			PolynomialFunction fitted = new PolynomialFunction(coeff);
			
			double offset = fitted.value(psm.getScanNumber());
			psm.setRecalibratedPrecursor((float) (psm.getPrecursor() - offset));
			if(Math.abs(psm.getRecalibratedPrecursorErr()) > tolerance.getToleranceAsDa(psm.getPrecursor() * psm.getCharge())) continue;
			filteredPSMs.add(psm);
		}
	}
	
	public ArrayList<MSGFPlusPSM> getFilteredPSMs(){
		return filteredPSMs;
	}
	
	
}
