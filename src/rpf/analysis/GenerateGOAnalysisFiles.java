package rpf.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class GenerateGOAnalysisFiles {
	static public void main(String[] args) throws IOException{
		String file = "/media/kyowon/Data1/RPF_Project/samples/sample4/results/outC_0.3.csv";
		BufferedLineReader in = new BufferedLineReader(file);
		PrintStream falseOut = new PrintStream(file + ".false.csv");
		PrintStream trueOut = new PrintStream(file + ".true.csv");
		String s;
		
		HashSet<String> falseSet = new HashSet<String>();
		HashSet<String> trueSet = new HashSet<String>();
		while((s=in.readLine())!=null){
			if(s.startsWith("Contig")) continue;
			String[] token = s.split("\t");
			if(!token[5].startsWith("NM_")) continue;
		//	if(!token[7].equals("T")) continue;
			if(token[9].equals("false")){
				String[] st = token[8].split(";");
				
				if(!st[0].toUpperCase().equals(st[1].toUpperCase())){
					boolean isconsistent = true;
					for(int j=2;j<st.length;j++){
						if(st[j].toUpperCase().equals(st[1].toUpperCase())) continue;
						isconsistent = false;
						break;
					}
					if(isconsistent) falseSet.add(token[5]);
					continue;
				}				
			}
			trueSet.add(token[5]);
			
		}
		
		for(String f : falseSet){
			falseOut.println(f);
		}
		for(String t : trueSet){
			if(falseSet.contains(t)) continue;
			trueOut.println(t);
		}
		
		falseOut.close();
		trueOut.close();
		in.close();
	}
}
