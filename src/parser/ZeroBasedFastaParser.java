package parser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
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
	
	public String getSequence(String contig, int start, int end){
		if(!refSeqs.containsKey(contig)) return null;
		byte[] refSeq = refSeqs.get(contig);
		StringBuffer seq = new StringBuffer();
		end = getLength(contig) < end? getLength(contig) : end;
		
		int tst=start;
		int tend=end;
		
		if(start>end){
			tst = end;
			tend = start;
		}
		
		for(int i=tst;i<tend;i++){
			if(i<0) continue;
			if(i>=refSeq.length) break;
			seq.append((char)refSeq[i]);
		}	
		String ret = start>end? getComplementaryCodon(seq.toString()): seq.toString();		
		return ret.toUpperCase();
	}
	
	public String getSequence(String contig, ArrayList<Integer> positions){ // positions should be sorted in ascending order
		if(!refSeqs.containsKey(contig)) return null;
		byte[] refSeq = refSeqs.get(contig);
		StringBuffer seq = new StringBuffer();
		for(int i : positions){
			if(i<0) continue;
			if(i>=refSeq.length) break;
			seq.append((char)refSeq[i]);
		}				
		return seq.toString().toUpperCase();
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

	
	public static void main(String[] args){
		ZeroBasedFastaParser test = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/data/hg19.fa");
		//System.out.println(test.getContigs());
//		System.out.println(test.getLength("1"));
		System.out.println(test.getSequence("chr20", 62151688, 62151688+150));
		//chr20_62151688_+
		//chr18_72057531
		
	}
}
