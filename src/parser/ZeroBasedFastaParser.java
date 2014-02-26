package parser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

public class ZeroBasedFastaParser {
	private HashMap<String, byte[]> refSeqs;
	
	
	public ZeroBasedFastaParser(String fastaFile){
		read(fastaFile);		
	}
	
	public Set<String> getContigs(){
		return refSeqs.keySet();
	}
	
	public boolean containsContig(String contig){
		return refSeqs.containsKey(contig);
	}
	
	public String getSequence(String contig){
		return getSequence(contig, 0, getLength(contig));
	}
	
	
	/**
    * Get the sequence
    * @param contig the contig
    * @param start the zero-based start position, inclusive. Should be smaller than end
    * @param end the zero-based end position, exclusive.
    * @return the sequence as it is read from the file.
    */
	public String getSequence(String contig, int start, int end){
		if(!refSeqs.containsKey(contig)) return new String();
		byte[] refSeq = refSeqs.get(contig);
		StringBuffer seq = new StringBuffer();
		end = getLength(contig) < end? getLength(contig) : end;
		
		for(int i=start;i<end;i++){
			if(i<0) continue;
			if(i>=refSeq.length) break;
			seq.append((char)refSeq[i]);
		}	
		////System.out.println(seq.toString());
		//String ret = start>end? getComplementaryCodon(seq.toString()): seq.toString();		
		return seq.toString();
	}
	

	/**
	    * Get the sequence
	    * @param contig the contig
	    * @param positions the zero-based positions. All inclusive, does not need to be sorted
	    * @return the sequence as it is read from the file.
	    */
	public String getSequence(String contig, List<Integer> positions){ 
		byte[] refSeq = refSeqs.get(contig);
		StringBuffer seq = new StringBuffer();
		for(int i : positions){
			if(i<0) continue;
			if(i>=refSeq.length) break;
			seq.append((char)refSeq[i]);
		}		
		
		return seq.toString();
	}
	
	/**
	    * Get the complementary sequence (order reversed)
	    * @param seq the sequence to be complemented
	    * @param reverseOrder reverse order
	    * @return the complementary sequence
	    */
	public static String getComplementarySequence(String seq, boolean reverseOrder){
		StringBuffer cc = new StringBuffer();
		char[] nas = seq.toCharArray();
		if(reverseOrder){
			for(int i = nas.length-1;i>=0;i--){
				char na = nas[i];
				cc.append(getComplementaryNA(na));
			}
		}else{
			for(int i = 0;i<nas.length;i++){
				char na = nas[i];
				cc.append(getComplementaryNA(na));
			}
		}
		return cc.toString();
	}
	
	private static char getComplementaryNA(char na){
		if(na == 'A') return 'T';
		if(na == 'T') return 'A';
		if(na == 'C') return 'G';
		if(na == 'G') return 'C';
		if(na == 'a') return 't';
		if(na == 't') return 'a';
		if(na == 'c') return 'g';
		return 'c';
	}
	
	public int getLength(String contig){
		if(!refSeqs.containsKey(contig)) return 0;
		return refSeqs.get(contig).length;
	}
	
	private void read(String fastaFile){
		refSeqs = new HashMap<String, byte[]>();
		try {
			IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(fastaFile));
			ReferenceSequence ref = null;
			while((ref = fasta.nextSequence())!=null){
				refSeqs.put(ref.getName(), ref.getBases());
			}
			fasta.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
		
	}

	//chr10 s mm9.chr10 53411805 5000 - 129993255 GACCTCAGGGGCCCTCTGGGTTATACTGGCTACTCAGGCAACTCATGGGTAAGCTGG----GCAGTGCACAGCTCTGTAGAGGTCAGTATAGGTGTGTTACCAAAGTCCCTGTGGAATGGCTTGCCTTCCAGGTTCTGCACTATGGGAACCCTGCCAACAGGAGCTTTGGCTTCCAGAGGGTTAGATCTGCCAGATTACAAGGGCCACTCTGAACAAGGCAGGCACTGGTGACTTGCTTCCTCACAAGGTCATAAGGACAACAGAGCAGGAGAGCCAGCCCAAGGTCAGACTCAGGTGCAGTTGAAGAAAGCCTCAGGCATGAATGTCCTCAAATCACAGGTTGCCAGTGGACACTGACTTCTAGGAAAGCTGTTCTCACAACTGCCAAAGAGAAGGGAAGCTGCAGG-TGGGGGCTGGGAGTCAGGATCATGGCGGATGCCTGCTTC---------------CTGGGGGATGCCCTG---TGCTGGTCCCTGTCTCTGACCTGTCTCTGAATCGCTCCTGTTCTCTAGGAGCTGGGGGGTGGGGGTGGTCTAGACCATGGGAGGTGTAGTCTTTCAGTAGCTGCCTGTGTGACAGTGTGA-CTGAGCCCAGCTTGTGATTATGAGTTCCTTGAAAACTTAGATGTTTGCTCCAGGGCCACAATCCTCCCACCCCAACCCTCCTGTTCAGCACTCAGCGGGCACCTACCAGGTACTGAACCTGAGTCAGATGGGTAGGAACTGTATCCTGGGACTCTCTTCCATCTCACCTAG-TCC--TTAGCTCTTAGAATCTAAGGTGTAGACAGTGTCCCCT---CTCCACCTTCCTGTGGGTGCAGGGATCTCTTTCTCTTTACCTTGGGTTCTCTCC--------CCAGCTTGG------CGGCTGACATCTTCTGGACCAAGCAATCTCCATGCTTCCTCT------ACAGATCCTGTGTGACAACA----TACTCTCTGTCTCCAAAGGGTCTGAGGTTATCCTTG-----------ATCAGGCAGTGGCAGGGAAATGGGTTGGTGCTGTTCCTGGCTCCCCGGCTTCTATTTTCCCTC-----------------------------------------------------------------------TTATCTGAAATCCTGACTGTTGCTG----------CTTGTGGGCACACAGAGTGATGCCAGACTGAAGAGCTGATGGGACCCCTTACACAGGGCTGGG-AGAGACAGAGATGGTGGCCT------------------TTAG------TAGACCAGGACCTTAGTGTTCCCTGATAGTGCAAGT-----------------------------GTTTCTGTGACCgcgtgtgtgcgtgcgtgcg---------------------tgcgtgcgtgcgtgcgtgcg--------------------------tgtgtgcgtgtgtgtg---------------------------------------------------------tgtgtgtgtgtgtgtgtgtgtgtgtTTGAGTATGGAGTAGaca----------------------------tgtgtgtgtccatgaatg------tgtatgcatgtgcatatgt-a-------tatgtaGT----------ACAAGAC---TCCAGATGGGGTCCAAAGCTTAGACTCT--GAGCTGTGAATCCAT----------------------CCACC-----CAGTGAA--------TACTGTGTA---------------------------------------------TTCAA--CACTCCAAAGTGAGCAGCA------------------------------------------------------------GG-----GTGTGGGGGGAACTAAGCA-AGAGTGAGTGCTGG-------TGTGGGTACTCAGTGTGGGAGGGACTGTGTGTTCCTAC-----ATAGATAGGGCCCT--TGAAATAGGAGCTATGGCCACTCCAGTGCC--------TCCAGCCCTGGCTTTTGTTTTTGCCCAGGAATTTCATG------------GGTATAGAGCCAGCCTGCCAGTCCAGGGACCCT----AGAAATAAGCATCTGGA-TGGGGC--AGGGACCAGGATATGACATCACTTTAGGGTCTTTCCTCGACCCCCCCTCCTCCCTGCCA---------------TGGAGCTGAGCCC----CCACTACTCAAACATGCAGGA------GTAGGGGG-CAGGGTTCTGGACTTCGGCCCTGGGA---TTTTCCAGAGATACTCTGAGGCATGAAATCACATTAGGA------AATCTTCAAGCA-------CACTGTTCACTGGGCATTTCTCAAGCA--GACGCAATCATTGTGGGGACTAACTTGGAGATGTAGAGCTCA--CAGGCACAGATGTGGGC---CTATGG--------CTGCCTCACTCTTAGAGAGGG-ACTTG--GTATCCAGCTG-TTTGCAGCT--GGACTCTGAGCCAAAC-GGCACAGCCCTGCCCTTTAGACAGGCATTGATTGATG--------------TTGATGTGGG-ACAGTTGG---------GGACATGCCTT-TGTGGCTGTTTCATCACA-----CAGAAGCATTGCC--ACTGCAGCCTTGGGGGAGGGTGGCAGTAACTCT-GT---CCACACCACAGCT-TGGAGAAAGTGGAGCCAGGAGGGAC---AGCTGGCTGAG-----GCTGTAGCCTGGTAC------------------------ATGCTACATCCATGTGTACAAGCTGGTTACATGAGTGTGTG----------TGTATG------------------------------------AGTGAGTATGGTTTTTGTGCATGCGCAGGCTGCTAGC-----------TTCTGTTTCTTC---------------------------------------------TGCCTCTTTGTCTGT---------------GGAATGGGC-------------GCTCACAGGGC----------------CTCACCTAGCA--------TGTGCCTGCCTGGCTGGTGTTGGGGActctgc-----------tgggcttcatggctttg-------------------ggctgtgtgttc---------------------ctgtctatcctggttcttttgtcactggggaaactgagggaggaa--gccac------------------agagggctgctggagtggagcccagtcctGCTCTGAGTTATGAATGGAACTGAGGGTAGCCAGCCTTGGCAGCAGCT---TCCCCATAGTTGATCTCAGTGTTTCCACCTTCCTGCCTGGTGGCCATTGGACACTTCTTTTATTGTCTTGTCTTTtccctccctctctccctctctccctccctcccctccctctctccctccctctctccctctctccctccctcccctccctctctccctccctcccttcctctcttcctATGTGCATTGATGTTTCGTATgagggcgttgggtcctggagttacagatggttgtgagccaccgtgtggttgctgggaattgaactcaggacctctggaagaccagtcggtaaccttaaccactgagccacc----------CTGGACACTTGTTTTCATTGAGATGTTTTGGCCTCTCATTCTCCGGTCTCTGCCTGAACTCAGTCCCTGTCCCTCCTCCTATGTAGCTAACCTGTTTCTGCTGGCCCAGATCTGCACTGCCTTAGACTTCATACAGAGCTGGGGAGGAGCAGAGTTTTTTTttttttttttttttttagatttatttatttattatatgtaagtacactgtagctatcttcagacactccagaagagggcgtcagtcagatctcattacggatggttgtgagccaccatgtggttgctgggatttgaactctggaccttcagaagagcagtcgggtgctcttaaccactgagccatctcaccagcccAGGAGCAGAGTTTTTGAATAAATGCCAACCTTTCCTGGGTCCGGGGGTGGTGTGGGGTGGGGTCCCTG-GTGCTGTTCTGGGTGGAAGTGTTCAGATCTATCTGGGGGCACCTGTAACAGCCTATAGAGGCTGAGCTCTTCCCACTAGACCTTGTAAAAAAGGGGAACCAGGGCAACTTACAATTCCAGTGCTCCACACCCCAAAGTATTTATACCCCATTGGGCACCCTTGGGGACCATCTTATGAATTTCAATAGTGATGG-----AGGTGAGGCAGAGAGCTATGCCTTGTAGGGTCTTAGATGTAGGACCACAGAACGTGGGGGG-TTACAGGGTGTTCCATCTCCTTGGCTTCATTTCCGTGGATGTTGACTTTGATTAGATACAGAACGTCTCTGAAATTCGAAAACCCAGCT-ACACCATGGACATCTCTGTCCAAAGCCCTAGGACATGTGTGTCTCCTAATCCCCAAGTCTAACTATCCAGAAGGTCTTGATGTTGGAGGGGCTGGGAAGGCATGTGGTCAGGTGGCCAGGCAAGTGTGTCCCTCCCTTTCCCATCTACATCTTGCTGTTGCCTAGGGCTCACAGTCAGACACTGAGGCAGGAGGAGActcccttcccttccccccccccccccccgcccccGTGGGAGTTCATACAGAGTCAGCACAGGAGGACATCCCTGGTGAGCATCCAGCTATATCTGTGAAGGGAATGTTGGAGAAGAACCTACCTGATCTAATGGGGTCAGGGCGGCTATGTGGAGGGAGGTTGGCAGGTACCCATGTGCTGTAACCCAGAGTAGGCTGCAGGATGGACAAGCCACTTATAGCCCCTGGTTTATTATAGATAGAGTGAGAT-TTTTACCTGAAGCTGTGTGAAGGAAAGAAACAGCACGGGGAACATCCAGCTGCAGACAGACCTAGCAAAGCAAATGCAAGAGTAACGGCAGTGCCCTGACTGTCTTTTCCTGGGACCCAGGCTGCTAAGCCATGCCGTGCGCTCTGCACGCAGCTGGAATCCATA-GACATGAGGCGGAAGCCTGCTGTTGATGGCATATCCCTTGGATACATCCCACAC---------------------------CTTTTTTCCTTTACCTTCCTCTTGAGGCCTGGCTA---ACTGCTACAGTCTGTGTTACTTTGAAGCAAC--CTCCAGGAATTCCTTTCCAAACACTTGCTGGAC-------------CTCAGAGCAGGTCACTGCAGAGAGAGCACCATAGGCAGATTATATGAGAGGATTTAAGAACAAACCAGAGAAAGAAGTAATGTGGCCACTGTGGCCCACAGAGAAAACTAGTGGACTCTTGGAGGGTCTAAGGGACAGATC-AAGTTGATGGCTTTTAGGCTTTGGGCAGGTAGGAGGGCTGT---------------------GGCAGGGGAAAGGAACGGACA-CCCAGCCTTTCTCTAAATCAACAACTGTAACACAAGAATTGGGACTGCAGCAAGCTGAGCGTTCTTCAAACTGGGTT---TCATCATAGAATGGTGGCCAGC--AGCCAGGA-CTTCCGGGT-GTTCT---CCAATGGAGGCCCTAGGATCTGGGGCGCGGTTTCCATGTCCCCG--------AGGACCCTGGCCAGGC-----------CCCCACGCCCACACCCACTCTCAGCA----GA---GACT-CAGCAGCTCTAAGTTACACAAGCTGGAGTTCCGC-TGGGTGTC----TGTTACCCCATCTTTCTGCTCTGGCTGCTCCTGGAGGTGCA-GCCAGG-------------AAGGAGTTAAA-------------TGGAATCTGAAGGTACCGTTTGGTGGATCACCTTGATTCTCAGGCTGGGTTC--ATTCCCTTCTGCCAGAGGGAAG-CGGTGGAGGGAAG---------GAGACTGG----GAACCAG---GA----------ATGAGGG------------AGG-----GT---------------------GGGGACC--TGGCAGGCC--CCACCTCCAGGA--------------CCACAGGAGCC
	
	public static void main(String[] args){
		ZeroBasedFastaParser fasta = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.fa");
		//chr5	4086564

		int off = 0;
		System.out.println(fasta.getSequence("chr5", 108863380-3, 108863380));
	//	System.out.println(fasta.getSequence("chr10", 53411805+5000+1, 53411805+1));
	///	System.out.println(fasta.getSequence("chr10", 53411805+5000-1, 53411805-1));
	}
}
