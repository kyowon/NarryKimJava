package parser;

import java.io.File;
import java.io.IOException;
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
		for(int i=start;i<end;i++){
			seq.append((char)refSeq[i]);
		}
				
		return seq.toString();
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
		ZeroBasedFastaParser test = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/test.fa");
		System.out.println(test.getContigs());
		System.out.println(test.getLength("1"));
		System.out.println(test.getSequence("1", 0, 1000));
		
	}
}
