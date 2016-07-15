package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class Temp {

	public static void main(String[] args) throws IOException {
		HashSet<String> genes = new HashSet<String>();
		HashSet<String> sps = new HashSet<String>();

		BufferedLineReader in = new BufferedLineReader("/media/kyowon/Data1/fCLIP/genes.remove");
		String s;
		while((s=in.readLine())!=null){
			if(s.isEmpty()) continue;
			genes.add(s);
		}
		in.close();
		System.out.println(genes);
		in = new BufferedLineReader("/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.MT.csv");
		while((s=in.readLine())!=null){
			for(String gene : genes){
				if(s.contains(gene)){
					String[] token = s.split("\t");
					sps.add(token[0] + " " + token[1] + " " + token[2]+" ");
					break;
				}
			}
		}
		
		in.close();
		
		in = new BufferedLineReader("/media/kyowon/Data1/fCLIP/samples/sample8/results/misc/bamMerged_cs.m");
		PrintStream out = new PrintStream("/media/kyowon/Data1/fCLIP/samples/sample8/results/misc/bamMerged_cs_selected.m");
		while((s=in.readLine())!=null){
			if(!s.contains(",")){
				out.println(s);
				continue;
			}
			for(String sp : sps){
				if(s.endsWith(sp)){
					out.println(s);
					break;
				}
			}
			
		}
		out.close();
		in.close();
		
	}

}
