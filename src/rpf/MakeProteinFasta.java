package rpf;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import net.sf.samtools.util.BufferedLineReader;
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
	
	private int peptideLength = 33;
	private ZeroBasedFastaParser fasta = null;
	private ScoringOutputParser scoringOutputParser = null;
	private String outFastaFile;
	private double scoreThreshold;
	private HashSet<String> excludingKeywords = null;
	
	public MakeProteinFasta(String scoringFile, String fastaFile, String outFastaFile, double scoreThreshold){
		scoringOutputParser = new ScoringOutputParser(scoringFile);
		fasta = new ZeroBasedFastaParser(fastaFile);
		this.scoreThreshold = scoreThreshold;
		this.outFastaFile = outFastaFile;
	
	}
	
	private void getExcludingKeywords(String blastFile){
		excludingKeywords = new HashSet<String>();
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(blastFile));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("Query= ")){
					excludingKeywords.add(s.split(" ")[1]);
				}
			}
			
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	
	//
	private void generate(){
		PrintStream out;
		HashMap<String, Integer> suffixMap = new HashMap<String, Integer>();
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
				
				String suffix = "_etc";
				if(position.getGeneName() == null) suffix = "_IG";
				else{
					if(position.getGBGeneName().startsWith("LINC")) suffix = "_LINC";
					else{
						if(pos >= position.getCdsStart() && pos < position.getCdsEnd()) suffix = "_dORF";
						else{
							if(isPlusStrand && pos < position.getCdsStart() - peptideLength * 3) suffix = "_uORF";
							else if(!isPlusStrand && pos > position.getCdsEnd() + peptideLength * 3) suffix = "_uORF";
						}
					}				
				}
				
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
					String codonString = codon.toString();
					String aa = codonTable.get(codonString);
					if(i==0 && (!codonString.equals("ATG") && !codonString.equals("CTG"))){
						suffix+= "_NonStartCodon"; break;
					}
					if(aa.equals("X")) break;
					peptide.append(aa);
				}
				//System.out.println(position);
				String pepName = contig+":" + pos + "_" + (isPlusStrand? '+' : '-');
				if(excludingKeywords != null && excludingKeywords.contains(pepName)) continue;
				
				
				
				if(!suffixMap.containsKey(suffix)) suffixMap.put(suffix, 0);
				suffixMap.put(suffix, suffixMap.get(suffix)+1);
				//if(!isPlusStrand)break;
				if(suffix.contains("_etc")) continue;	
				if(peptide.length() < 7) continue;
				
				out.println(">" + pepName + suffix);
				out.println(peptide);
			}
			out.close();
			System.out.println(suffixMap);
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
		String key = "Noco";
		//System.out.println(codonTable.size());
		MakeProteinFasta test = new MakeProteinFasta("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/"+ key + "_Harr10mNew.sorted.plus.cov.score.tsv.windowed.tsv",
				"/media/kyowon/Data1/RPF_Project/data/hg19.fa",
				"/media/kyowon/Data1/RPF_Project/data/hg19_protein_" + key + "_2.1.fasta",
				2.1);
		
		//test.getExcludingKeywords("/media/kyowon/Data1/RPF_Project/data/hg19_protein_" + key + "_1.8.blast");
		test.generate();
	}
}
