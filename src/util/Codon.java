package util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

public class Codon extends ArrayList<Nucleotide>{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	static private HashMap<String, Codon> standardCodons;
	static private HashMap<String, ArrayList<Codon>> aaCodonstable;
	static{		
		standardCodons = new HashMap<String, Codon>();
		aaCodonstable = new HashMap<String, ArrayList<Codon>>();
		
		HashMap<String, String> codonAAtable = new HashMap<String, String>();
		codonAAtable.put ("TTT", "F");
		codonAAtable.put ("TTC", "F");
		codonAAtable.put ("TTA", "L");
		codonAAtable.put ("TTG", "L");
		codonAAtable.put ("TCT", "S");
		codonAAtable.put ("TCC", "S");
		codonAAtable.put ("TCA", "S");
		codonAAtable.put ("TCG", "S");
		codonAAtable.put ("TAT", "Y");
		codonAAtable.put ("TAC", "Y");
		codonAAtable.put ("TAA", "X");
		codonAAtable.put ("TAG", "X");
		codonAAtable.put ("TGA", "X");
		codonAAtable.put ("TGT", "C");
		codonAAtable.put ("TGC", "C");
		codonAAtable.put ("TGG", "W");
		codonAAtable.put ("CTT", "L");
		codonAAtable.put ("CTC", "L");
		codonAAtable.put ("CTA", "L");
		codonAAtable.put ("CTG", "L");
		codonAAtable.put ("CCT", "P");
		codonAAtable.put ("CCC", "P");
		codonAAtable.put ("CCA", "P");
		codonAAtable.put ("CCG", "P");
		codonAAtable.put ("CAT", "H");
		codonAAtable.put ("CAC", "H");
		codonAAtable.put ("CAA", "Q");
		codonAAtable.put ("CAG", "Q");
		codonAAtable.put ("CGT", "R");
		codonAAtable.put ("CGC", "R");
		codonAAtable.put ("CGA", "R");
		codonAAtable.put ("CGG", "R");
		codonAAtable.put ("ATT", "I");
		codonAAtable.put ("ATC", "I");
		codonAAtable.put ("ATA", "I");
		codonAAtable.put ("ATG", "M");
		codonAAtable.put ("ACT", "T");
		codonAAtable.put ("ACC", "T");
		codonAAtable.put ("ACA", "T");
		codonAAtable.put ("ACG", "T");
		codonAAtable.put ("AAT", "N");
		codonAAtable.put ("AAC", "N");
		codonAAtable.put ("AAA", "K");
		codonAAtable.put ("AAG", "K");
		codonAAtable.put ("AGT", "S");
		codonAAtable.put ("AGC", "S");
		codonAAtable.put ("AGA", "R");
		codonAAtable.put ("AGG", "R");
		codonAAtable.put ("GTT", "V");
		codonAAtable.put ("GTC", "V");
		codonAAtable.put ("GTA", "V");
		codonAAtable.put ("GTG", "V");
		codonAAtable.put ("GCT", "A");
		codonAAtable.put ("GCC", "A");
		codonAAtable.put ("GCA", "A");
		codonAAtable.put ("GCG", "A");
		codonAAtable.put ("GAT", "D");
		codonAAtable.put ("GAC", "D");
		codonAAtable.put ("GAA", "E");
		codonAAtable.put ("GAG", "E");
		codonAAtable.put ("GGT", "G");
		codonAAtable.put ("GGC", "G");
		codonAAtable.put ("GGA", "G");
		codonAAtable.put ("GGG", "G");
		
		for(String seq : codonAAtable.keySet()){
			standardCodons.put(seq, new Codon(seq).setCodingAA(codonAAtable.get(seq)));
		}
		
		for(Codon codon : standardCodons.values()){
			String aa = codon.getCodingAA();
			if(!aaCodonstable.containsKey(aa)) aaCodonstable.put(aa, new ArrayList<Codon>());
			aaCodonstable.get(aa).add(codon);
		}
	}
	
	private String nucleotideSeq;
	private String codingAA;
	
	private Codon(String seq){
		nucleotideSeq = seq;
		for(int i=0;i<3;i++){
			this.add(Nucleotide.getStandardNucleotide((seq.substring(i, i+1))));
		}
	}
	private Codon setCodingAA(String aa){
		codingAA = aa;
		return this;
	}
	
	public String getNucleotideSeq() { return nucleotideSeq; }
	public String getCodingAA(){return codingAA;}
	
	static public Codon getStandardCodon(String seq){		
		return standardCodons.get(seq.toUpperCase());
	}
	
	static public ArrayList<Codon> getStandardCodons(){
		return new ArrayList<Codon>(standardCodons.values());
	}
	
	static public ArrayList<Codon> getCodingCodons (String aa){
		return aaCodonstable.get(aa.toUpperCase());
	}
	
	static public String getAminoAcids(String seq){
		StringBuffer aas = new StringBuffer();
		for(int i=0;i<seq.length()/3;i++){
			String codon = seq.substring(i*3, i*3+3);
			Codon scodon = getStandardCodon(codon);
			if(scodon == null) break;
			String aa = scodon.getCodingAA();
			if(aa.equals("X")) break;
			aas.append(aa);
		}
		
		return aas.toString();
	}
	
	static public ArrayList<Codon> getMatchingStandardCodons(String regex){
		ArrayList<Codon> matchingCodons = new ArrayList<Codon>();
		Pattern pattern = Pattern.compile(regex);
		for(Codon codon : getStandardCodons()){
			if(!pattern.matcher(codon.getNucleotideSeq()).find()) continue;
			matchingCodons.add(codon);
		}
		return matchingCodons;
	}
	
	static public HashSet<String> getStopCodonStringSet(){
		HashSet<String> stopCodons = new HashSet<String>();
		stopCodons.add ("TAA");
		stopCodons.add ("TAG");
		stopCodons.add ("TGA");
		return stopCodons;
	}
	
	
	public ArrayList<Codon> getMutatedStandardCodons(int index){
		String seq = getNucleotideSeq();
		StringBuffer regex = new StringBuffer();
		for(int i=0;i<seq.length();i++){
			if(i == index){
				regex.append("[^");
				regex.append(seq.charAt(i));
				regex.append("]");
			}else regex.append(seq.charAt(i));
		}		
		return Codon.getMatchingStandardCodons(regex.toString());
	}
	
	public static void main(String[] args){
		System.out.println(Codon.getStandardCodon("GGG").getCodingAA());
		for(Codon codon : Codon.getStandardCodon("GGG").getMutatedStandardCodons(2)){
			System.out.println(codon.getNucleotideSeq() + " " + codon.getCodingAA());
		}
	}
	
}
