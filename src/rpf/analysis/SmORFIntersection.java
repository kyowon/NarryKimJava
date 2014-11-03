package rpf.analysis;

import java.io.IOException;
import java.util.HashSet;

import parser.BufferedLineReader;

public class SmORFIntersection {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		BufferedLineReader in1 = new BufferedLineReader("/media/kyowon/Data1/RPF_Project/samples/sample3/results/out_0.3_smORF.csv");
		BufferedLineReader in2 = new BufferedLineReader("/media/kyowon/Data1/RPF_Project/samples/sample5/results/out_0.3_smORF.csv");
		String s;
		HashSet<String> seqs1 = new HashSet<String>();
		HashSet<String> seqs2 = new HashSet<String>();
		HashSet<String> intersection = new HashSet<String>();
		
		while((s = in1.readLine())!=null){
			if(s.startsWith("Contig")) continue;
			String[] token = s.split("\t");
			seqs1.add(token[token.length-1]);
		}
		
		while((s = in2.readLine())!=null){
			if(s.startsWith("Contig")) continue;
			String[] token = s.split("\t");
			seqs2.add(token[token.length-1]);
			if(seqs1.contains(token[token.length-1])) intersection.add(token[token.length-1]);
		}
		
		in1.close();
		in2.close();
		
		System.out.println(seqs1.size());
	//	System.out.println(seqs1);
		System.out.println(seqs2.size());
	//	System.out.println(seqs2);
		System.out.println(intersection.size());
	//	System.out.println(intersection);
		
	}

}
