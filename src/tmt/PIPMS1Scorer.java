package tmt;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import tmt.FilterWIthPIP.TmtPSM;
import msgf.Tolerance;
import msutil.Composition;
import msutil.Peak;
import msutil.Spectrum;

public class PIPMS1Scorer {
	private Spectrum spectrum; // ms 1 spectrum
	//private Composition composition;  // composition should include protons due to charge state..
	private Tolerance tolerance;
	private IsotopeDistributionGenerator isotope;
	private int charge;
	private TmtPSM psm;
	
	public PIPMS1Scorer(Spectrum spectrum, TmtPSM psm, Tolerance tolerance, int maxBin){
		this.spectrum = spectrum;
		this.psm = psm;
		this.tolerance = tolerance;
		this.charge = psm.getCharge();
		this.isotope = new IsotopeDistributionGenerator(psm.getComposition(), psm.getMz() * psm.getCharge(), maxBin, psm.getTmtTagCount());
	}
	
	
	public double getScore(){
		float minMass = (float) ((isotope.getMonoIsotopeMass()-Composition.ISOTOPE)/charge);
		float maxMass = isotope.getMaxMass()/charge;
		ArrayList<Peak> allPeaks = spectrum.getPeakListByMassRange(minMass, maxMass);
		double norm = 0;
		//System.out.println("Charge : " + charge);
		for(Peak peak : allPeaks){
			norm += peak.getIntensity() * peak.getIntensity();
		//	System.out.println("all " + peak);
		}
		
		norm = Math.sqrt(norm);
		
		//System.out.println(" im " + isotope.getMonoIsotopeMass()/ charge);
	//	System.out.println(minMass + " " + maxMass);
		
		ArrayList<Peak> selectedPeaks = new ArrayList<Peak>();
		double sum = 0;		
		HashMap<Float, Double> isotopeMap = isotope.getDistMap(charge);
	//	System.out.println(isotopeMap.keySet());
		ArrayList<Float> mzs = new ArrayList<Float>(isotopeMap.keySet()); 
		Collections.sort(mzs);
		for(float mz : mzs){
			Peak peak = spectrum.getPeakByMass(mz, tolerance);
			//System.out.print("\tselected " + mz);
			if(peak != null){
				selectedPeaks.add(peak);
				sum += peak.getIntensity() * isotopeMap.get(mz);
				//System.out.println("\t" + peak);
			}//else System.out.println("\tMissing");
		}
		
		//System.out.println("_____________________________________________");
		//
		//if((sum) / (norm+.1) == 0 ){
		//	System.out.println(allPeaks + "\n\t" + mzs + "\n\t" + selectedPeaks);
		//	System.out.println(psm.getComposition() + " " + psm.getMz() * psm.getCharge()+ " " +  psm.getTmtTagCount());
		//	for(float mz : mzs){ 
		//		System.out.println(isotopeMap.get(mz));
		//	}
		//}
		return (sum) / (norm+.1);		
	}
}
