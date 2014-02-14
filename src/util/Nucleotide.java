package util;

import java.util.ArrayList;
import java.util.HashMap;

public class Nucleotide {
	private String name; // add composition and etc in the future...
//	private double mass;
	private static HashMap<String, Nucleotide> standardNucleotides;
	private Nucleotide(String name){
		this.name = name;
	//	this.mass = mass;
	}
	
	static{
		standardNucleotides = new HashMap<String, Nucleotide>();
		standardNucleotides.put("A", new Nucleotide("A"));
		standardNucleotides.put("C", new Nucleotide("C"));
		standardNucleotides.put("G", new Nucleotide("G"));
		standardNucleotides.put("T", new Nucleotide("T"));
	}

	static public ArrayList<Nucleotide> getStandardNucleotides(){
		return new ArrayList<Nucleotide>(standardNucleotides.values());
	}
	
	static public Nucleotide getStandardNucleotide(String name){
		return standardNucleotides.get(name);
	}
	
	
	public String getName() { return name; }

	@Override
	public int hashCode() {
		return this.name.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(!(obj instanceof Nucleotide)) return false;
		return ((Nucleotide)obj).getName().equals(this.getName());
	}
	
	public static String getComplementarySeq(String seq){
		StringBuffer cc = new StringBuffer();
		char[] nas = seq.toCharArray();
		for(int i = nas.length-1;i>=0;i--){
			char na = nas[i];
			cc.append(getComplementaryNA(na));
		}
		return cc.toString();
	}
	
	private static char getComplementaryNA(char na){
		if(na == '-') return '-';
		if(na == 'A') return 'T';
		if(na == 'T') return 'A';
		if(na == 'C') return 'G';
		return 'C';
	}
	
}
