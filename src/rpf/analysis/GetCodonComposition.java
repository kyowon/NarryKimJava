package rpf.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import parser.BufferedLineReader;
import parser.ZeroBasedFastaParser;

public class GetCodonComposition {
	public static void main(String[] args) throws IOException{
		String merged ="/media/kyowon/Data1/RPF_Project/samples/sample1/results/RPF1_2.1.csv";
		ZeroBasedFastaParser fasta = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		PrintStream outFasta = new PrintStream("/media/kyowon/Data1/RPF_Project/genomes/hg19_ATG_Annotated.fa");
		BufferedLineReader in = new BufferedLineReader((merged));
		String s;
		HashMap<String, Double> codonCounter = new HashMap<String, Double>();
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			if(!token[3].equals("ATG")) continue;
			String contig = token[0];
			boolean isPlusStrand = token[2].equals("+");
			boolean isAnnotated = token[12].equals("T");
			if(!token[11].equals("0")) continue;
			if(isAnnotated) continue;
			
			//System.out.println(token[2] + " " + token[12]);
			outFasta.println(">"+token[13]);
			String codon = null;
			if(isPlusStrand){
				int start = Integer.parseInt(token[1])+3;
				int end = start + 21;
				for(int i=start; i<end;i+=3){
					codon = fasta.getSequence(contig, i, i+2);
				}
				outFasta.println(fasta.getSequence(contig, start, end));
			}else{
				int start = Integer.parseInt(token[1])-2;
				int end = start - 21;
				for(int i=start; i>end;i-=3){
					codon = fasta.getSequence(contig, i, i-2);
				}
				outFasta.println(fasta.getSequence(contig, start, end));
			}
			
			if(codon !=null){
				if(!codonCounter.containsKey(codon)) codonCounter.put(codon, 0.0);
				codonCounter.put(codon, codonCounter.get(codon)+1);
			}
			
		}
		
		
		
		double sum = 0;
		for(String key : codonCounter.keySet()){
			sum += codonCounter.get(key);
		}
		for(String key : codonCounter.keySet()){
			codonCounter.put(key, codonCounter.get(key)/sum);
		}
		ArrayList<String> codons = new ArrayList<String>(codonCounter.keySet());
		Collections.sort(codons);
		
		for(String codon : codons){
			System.out.println("%" + codon + " \n " + codonCounter.get(codon) * 1000);
		}
		outFasta.close();
		in.close();
	}
}
