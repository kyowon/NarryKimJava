package rpf;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;

import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;
import parser.ZeroBasedFastaParser;

public class MakeProteinFasta {
	private static HashMap<String, String> codonTable;
	static {
		codonTable = new HashMap<String, String>();
		codonTable.put("ATT", "I");
		codonTable.put("ATC", "I");
		codonTable.put("ATA", "I");
		
		codonTable.put("CTT", "L");
		codonTable.put("CTC", "L");
		codonTable.put("CTA", "L");
		codonTable.put("CTG", "L");
		codonTable.put("TTA", "L");
		codonTable.put("TTG", "L");
		
		codonTable.put("GTT", "V");
		codonTable.put("GTC", "V");
		codonTable.put("GTA", "V");
		codonTable.put("GTG", "V");
		
		codonTable.put("TTT", "F");
		codonTable.put("TTC", "F");		

		codonTable.put("ATG", "M");
		
		codonTable.put("TGT", "C");
		codonTable.put("TGC", "C");
		
		codonTable.put("GCT", "A");
		codonTable.put("GCC", "A");
		codonTable.put("GCA", "A");
		codonTable.put("GCG", "A");
		
		codonTable.put("GGT", "G");
		codonTable.put("GGC", "G");
		codonTable.put("GGA", "G");
		codonTable.put("GGG", "G");
		
		codonTable.put("CCT", "P");
		codonTable.put("CCC", "P");
		codonTable.put("CCA", "P");
		codonTable.put("CCG", "P");

		codonTable.put("ACT", "T");
		codonTable.put("ACC", "T");
		codonTable.put("ACA", "T");
		codonTable.put("ACG", "T");

		codonTable.put("TCT", "S");
		codonTable.put("TCC", "S");
		codonTable.put("TCA", "S");
		codonTable.put("TCG", "S");
		codonTable.put("AGT", "S");
		codonTable.put("AGC", "S");
		
		codonTable.put("TAT", "Y");
		codonTable.put("TAC", "Y");
		
		codonTable.put("TGG", "W");
		
		codonTable.put("CAA", "Q");
		codonTable.put("CAG", "Q");
		
		codonTable.put("AAT", "N");
		codonTable.put("AAC", "N");
		
		codonTable.put("CAT", "H");
		codonTable.put("CAC", "H");

		codonTable.put("GAA", "E");
		codonTable.put("GAG", "E");
		
		codonTable.put("GAT", "D");
		codonTable.put("GAC", "D");
		
		codonTable.put("AAA", "K");
		codonTable.put("AAG", "K");
		
		codonTable.put("CGT", "R");
		codonTable.put("CGC", "R");
		codonTable.put("CGA", "R");
		codonTable.put("CGG", "R");
		codonTable.put("AGA", "R");
		codonTable.put("AGG", "R");
		
		codonTable.put("TAA", "X");
		codonTable.put("TAG", "X");
		codonTable.put("TGA", "X");
	}
	
	private int peptideLength = 50;
	private ZeroBasedFastaParser fasta = null;
	private ScoringOutputParser scoringOutputParser = null;
	private String outFastaFile;
	private double scoreThreshold;
	
	public MakeProteinFasta(String scoringFile, String fastaFile, String outFastaFile, double scoreThreshold){
		scoringOutputParser = new ScoringOutputParser(scoringFile);
		fasta = new ZeroBasedFastaParser(fastaFile);
		this.scoreThreshold = scoreThreshold;
		this.outFastaFile = outFastaFile;
		generate();
	}
	
	private void generate(){
		PrintStream out;
		try {
			out = new PrintStream(outFastaFile);
			for(ScoredPosition position : scoringOutputParser.getPositions()){
				StringBuffer peptide = new StringBuffer();
				if(position.getScore() < scoreThreshold) continue;
				if(position.isAnnotated()) continue;
				String contig = position.getContig();
				boolean isPlusStrand = position.isPlusStrand();
				int pos = position.getPosition(); 
				int start = isPlusStrand? pos : pos - peptideLength * 3 + 1;
				int end = isPlusStrand? pos + peptideLength * 3 : pos + 1;
						
				String nas = fasta.getSequence(contig, start, end);
				if(nas == null) continue;
				if(!isPlusStrand){
					StringBuffer tnas = new StringBuffer();
					for(int i=nas.length()-1;i>=0;i--){
						char na = nas.charAt(i);
						tnas.append(getComplementaryNA(na));
					}
					nas = tnas.toString();
				}
				
				for(int i=0;i<nas.length();i+=3){
					char[] dst = new char[3];
					nas.getChars(i, i+3, dst, 0);
					StringBuffer codon = new StringBuffer();
					for(char na : dst){
						codon.append(na);
					}
					String aa = codonTable.get(codon.toString());
					if(aa.equals("X")) break;
					peptide.append(aa);
				}
				if(peptide.length() < 7) continue;
				//System.out.println(position);
				String pepName = ">" + contig+"_" + pos + "_" + (isPlusStrand? '+' : '-');
				out.println(pepName);
				out.println(peptide);
				//if(!isPlusStrand)break;
			}
			out.close();
		} catch (FileNotFoundException e) {			
			e.printStackTrace();
		}
	}
	
	
	
	public static String getComplementaryCodon(String codon){
		StringBuffer cc = new StringBuffer();
		char[] nas = codon.toCharArray();
		for(int i = nas.length-1;i>=0;i--){
			char na = nas[i];
			cc.append(getComplementaryNA(na));
		}
		return cc.toString();
	}
	
	private static char getComplementaryNA(char na){
		if(na == 'A') return 'T';
		if(na == 'T') return 'A';
		if(na == 'C') return 'G';
		return 'C';
	}
	
	public static void main(String[] args){
		System.out.println(codonTable.size());
		MakeProteinFasta test = new MakeProteinFasta("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10mNew.sorted.plus.cov.score.tsv.windowed.tsv",
				"/media/kyowon/Data1/RPF_Project/data/hg19.fa",
				"/media/kyowon/Data1/RPF_Project/data/hg19_protein_Thy_2.1.fasta",
				2.1);
	}
}
