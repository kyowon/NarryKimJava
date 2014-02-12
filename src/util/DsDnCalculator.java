package util;

import java.util.ArrayList;
import java.util.HashMap;

import net.sf.samtools.util.Tuple;

public class DsDnCalculator {
	
	static private HashMap<String, ArrayList<Double>> synSiteMap = null;
	static private HashMap<Tuple<String, String>, Double> synSubMap = null;
	
	static{
		synSiteMap = new HashMap<String, ArrayList<Double>>();
		
	}
	
	static private double getNonSyn(String codonSeq, int index){
		if(!codonSeq.contains("-")){
			Codon codon = Codon.getStandardCodon(codonSeq);
			String aa = codon.getCodingAA();
			double i = 0;
			
			for(Codon mcodon : codon.getMutatedStandardCodons(index)){
				if(mcodon.getCodingAA().equals(aa)) continue;
				i++;
			}
			return i/3;		
		}else{
			ArrayList<Codon> mcodons = Codon.getMatchingStandardCodons(codonSeq.replace('-', '.'));
			double i = 0;
			for(Codon mc : mcodons){
				System.out.println(mc.getNucleotideSeq() + " " + getNonSyn(mc.getNucleotideSeq(), index));
				i += getNonSyn(mc.getNucleotideSeq(), index);
			}
			return i/mcodons.size();
		}
	}
	
	
	public double getDsDnRatio(String[] sequences){
		return 0;
	}
	
	public static void main(String[] args){
		System.out.println(DsDnCalculator.getNonSyn("ACT", 2));
	}
	
}
