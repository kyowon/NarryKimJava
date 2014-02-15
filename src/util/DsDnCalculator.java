package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.MafParser;


public class DsDnCalculator {
	
	static private HashMap<String, ArrayList<Double>> nonSynSiteMap = null;
	static private HashMap<ArrayList<String>, Double> nonSynSubMap = null;
	
	static{
		nonSynSiteMap = new HashMap<String, ArrayList<Double>>();
		nonSynSubMap = new HashMap<ArrayList<String>, Double>();
		
		for(String codonSeq : getCodonSeqsWithWildCards()){
			ArrayList<Double> nc = new ArrayList<Double>();
			for(int i=0;i<3;i++)
				nc.add(getNonSynSite(codonSeq, i));
			nonSynSiteMap.put(codonSeq,nc);
			for(String codonSeq2 : getCodonSeqsWithWildCards()){
				//if(codonSeq.equals("ACT")) System.out.println(codonSeq2);
				ArrayList<String> key = new ArrayList<String>();
				key.add(codonSeq); key.add(codonSeq2);
				nonSynSubMap.put(key, getNonSynSub(codonSeq, codonSeq2));
			}			
		}
	}
	
	static private ArrayList<String> getCodonSeqsWithWildCards(){
		ArrayList<String> codons = new ArrayList<String>();
		for(Codon codon : Codon.getStandardCodons())
			codons.add(codon.getNucleotideSeq());
		
		codons.add("---");
		codons.add("A--");
		codons.add("C--");
		codons.add("G--");
		codons.add("T--");
		
		codons.add("-A-");
		codons.add("-C-");
		codons.add("-G-");
		codons.add("-T-");
		
		codons.add("--A");
		codons.add("--C");
		codons.add("--G");
		codons.add("--T");
		
		codons.add("AA-"); codons.add("A-A"); codons.add("-AA");
		codons.add("AC-"); codons.add("A-C"); codons.add("-AC");
		codons.add("AG-"); codons.add("A-G"); codons.add("-AG");
		codons.add("AT-"); codons.add("A-T"); codons.add("-AT");
		
		codons.add("CA-"); codons.add("C-A"); codons.add("-CA");
		codons.add("CC-"); codons.add("C-C"); codons.add("-CC");
		codons.add("CG-"); codons.add("C-G"); codons.add("-CG");
		codons.add("CT-"); codons.add("C-T"); codons.add("-CT");
		
		codons.add("GA-"); codons.add("G-A"); codons.add("-GA");
		codons.add("GC-"); codons.add("G-C"); codons.add("-GC");
		codons.add("GG-"); codons.add("G-G"); codons.add("-GG");
		codons.add("GT-"); codons.add("G-T"); codons.add("-GT");
		
		codons.add("TA-"); codons.add("T-A"); codons.add("-TA");
		codons.add("TC-"); codons.add("T-C"); codons.add("-TC");
		codons.add("TG-"); codons.add("T-G"); codons.add("-TG");
		codons.add("TT-"); codons.add("T-T"); codons.add("-TT");
		
		return codons;
	}
	
	static private double getNonSynSite(String codonSeq, int index){
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
			//	System.out.println(mc.getNucleotideSeq() + " " + getNonSynSite(mc.getNucleotideSeq(), index));
				i += getNonSynSite(mc.getNucleotideSeq(), index);
			}
			return i/mcodons.size();
		}
	}
	
	static private double getNonSynSub(String codonSeq1, String codonSeq2){
		if(!codonSeq1.contains("-") && !codonSeq2.contains("-")){
			if(Codon.getStandardCodon(codonSeq1).getCodingAA().equals(Codon.getStandardCodon(codonSeq2).getCodingAA())) return 0;
			return 1;
		}else{
			int j=0;
			double i = 0;
			ArrayList<Codon> mcodons1 = Codon.getMatchingStandardCodons(codonSeq1.replace('-', '.'));
			ArrayList<Codon> mcodons2 = Codon.getMatchingStandardCodons(codonSeq2.replace('-', '.'));
			for(Codon mc1 : mcodons1){
				for(Codon mc2 : mcodons2){
					j++;
					i += getNonSynSub(mc1.getNucleotideSeq(), mc2.getNucleotideSeq());
				}
			}
			return i/j;
		}
	}
	
	static public double getNonSynSiteSum(String sequence){
		double ret = 0;
		for(int i=0;i<sequence.length()/3;i++){
			for(double v : nonSynSiteMap.get(sequence.substring(i*3, i*3+3))){
				ret += v;
			}
		}
		return ret;
	}
	
	static public double getNonSynSubSum(String sequence1, String sequence2){
		double ret = 0;
		for(int i=0;i<sequence1.length()/3;i++){
			ArrayList<String> key = new ArrayList<String>();
			key.add(sequence1.substring(i*3, i*3+3));
			key.add(sequence2.substring(i*3, i*3+3));
			//System.out.println(key);
			ret += nonSynSubMap.get(key);
		}
		return ret;
	}
	
	static public double getDsDnRatio(String[] sequences){
		double nNSSubs = 0;
		double nNSSites = 0;
		int n1=0, n2=0;
		for(int i=0; i<sequences.length; i++){
			sequences[i] = sequences[i].toUpperCase();
			//if(!isPlusStrand) sequences[i] = Nucleotide.getComplementarySeq(sequences[i]);
			//System.out.println(sequences[i]);
		}
		
		/*String[] seqs = new String[sequences.length];
		for(int i=0; i<seqs.length; i++){
			seqs[i] = "";
		}
		for(int i=0;i<sequences[0].length();i++){
			if(sequences[0].charAt(i) != '-'){
				for(int k=0; k<sequences.length; k++){
					seqs[k] += sequences[k].charAt(i);
				}
			}
		}
		sequences = seqs;
		*/
		for(int i=0; i<sequences.length; i++){
			//sequences[i] = sequences[i].toUpperCase();
			n1++;
			nNSSites += getNonSynSiteSum(sequences[i]);
			for(int j=0;j<sequences.length; j++){
				if(i==j) continue;
				n2++;
				nNSSubs += getNonSynSubSum(sequences[i], sequences[j]);
			}
		}
		if(n1>0) nNSSites /= n1;
		if(n2>0) nNSSubs /= n2;
		
		double nSSites = sequences[0].length() - nNSSites;
		double nSSubs = sequences[0].length()/3 - nNSSubs;
	//	System.out.println(nNSSites + " " + n1);
		return (nNSSubs * nSSites) / (nNSSites * nSSubs);
	}
	
	public static void out(String mafFileDir, String outFile){		
		try {
			PrintStream out = new PrintStream(outFile);
			
			for(File mafFile : new File(mafFileDir).listFiles()){
				if(!mafFile.getName().endsWith("maf.gz"))continue;
				MafParser mp = new MafParser(mafFile.getAbsolutePath());
				
				for(String contig : mp.getContigs()){
					out.println("CONTIG\t"+contig);
					ArrayList<Integer> sPositions = mp.getStartPositions(contig);
					ArrayList<Integer> ePositions = mp.getEndPositions(contig);
					for(int i=0; i<sPositions.size();i++){
						int sp = sPositions.get(i);
						int ep = ePositions.get(i);
						for(int cp = sp-2;cp<ep+3;cp++){
							double dsdnp = getDsDnRatio(mp.getSeqs(contig, cp, true, 3));
							double dsdnm = getDsDnRatio(mp.getSeqs(contig, cp, true, 3));
							out.println(cp + "\t" + dsdnp + "\t" + dsdnm);
						}
						
					}
					
				}
			}
			
			
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	public static void main(String[] args){		
		String[] s = {
			//	"ACTCCGAACGGGGCGTTAGAGTTGAAACCCGTTAGA", 
				"----CAGACCATAACTAAAA-----ACAATAA------ATAACT-----------------GATGCACAA-TCTGAAC---TGTGCAGTGCTTATACAAACTTA--ATA",
				"----CAGACGGTACCTAAAA-----ACAGGAA------ATGACT-----------------GACGCACAA-TCCGGAC---TGTGC-CTGCTTGTGTGAAGTTA--ATA",
		//		"----TCTACTGTAACTAAAA-----TAGATAACCCCTCATGACC-----------------TAGTCACTT-TCTGAA--------------------------------",
		//		"----CC-----TGACTAAA-----------------------------------------------ATAAGTCCAGAA---TGTATAATACTTGTGCAAATTTA--ATA",
		//		"----tctactgtaattaaaatagctatagtagactgtagacacc-----------------aatgcataaatccagaa---gatataatgcttgtgcgaattgg--ata",
		//		"AAAAGCCACTGAAAC--AA------AAAGTGT------AACGCTGCTGCAGAAAGTAGTGCAATGCTTCAATCCAGAA---GACACATTGGCTGTGGGGG-TCA--GTA",
			//	"------------------------------------",
				//"ACTCCGAACGGGGCGTTAGAGTTGAAACCCGTTAGA",
				//"------------------------------------",
			//	"------------------------------------",
				
				//"------------------------------------",
				//"------------------------------------",
				//"------------------------------------",
				//"------------------------------------",
				//"------------------------------------",
				//"------------------------------------",
				//"------------------------------------",
			};
		
		
		//System.out.println(nonSynSubMap.get(new Tuple<String, String>("ACT", "ACG")));
		
		//System.out.println(getDsDnRatio(s));
		
		out("/media/kyowon/Data1/RPF_Project/genomes/mm9/maf", "/media/kyowon/Data1/RPF_Project/genomes/mm9/maf/out.txt");
		//System.out.println(nonSynSiteMap.get("TTA"));
		//System.out.println(nonSynSiteMap.get("ATA"));
		//System.out.println(getNonSynSiteRatioSum("ACTCCGAACGGGGCGTTAGAGTTGAAACCCGTTAGA"));
		//System.out.println(getNonSynSiteRatioSum("ACGCCGATCGGCGCGATAGGGTTCAAGCTCGTACGA"));
		//System.out.println(DsDnCalculator.getNonSyn("ACT", 2));
		//System.out.println(nonSynSiteMap);
		//
	}
	
}
