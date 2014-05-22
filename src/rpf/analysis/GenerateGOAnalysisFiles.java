package rpf.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class GenerateGOAnalysisFiles {
	static public void main(String[] args) throws IOException{
		String file = "/media/kyowon/Data1/RPF_Project/samples/sample4/results/out_0.3.csv";
		BufferedLineReader in = new BufferedLineReader(file);
		PrintStream falseOut = new PrintStream(file + ".target.csv");
		PrintStream trueOut = new PrintStream(file + ".background.csv");
		String s;
		
		int[] targetGroup = {0,1}; // NOCO
		int[] controlGroup = {2,3}; // THY
		
		HashSet<String> targets = new HashSet<String>();
		HashSet<String> controls = new HashSet<String>();
		int caseNumber = 1;
		if(caseNumber == 1){
			targets.add("5U");
			targets.add("NC");
			controls.add("T");
			controls.add("TS");
			//controls.add("T");
			//controls.add("TS");
		}else{
			controls.add("NC");
			//targets.add("5U");
			//targets.add("T");
			//targets.add("TS");
		}
		
		HashSet<String> falseSet = new HashSet<String>();
		HashSet<String> trueSet = new HashSet<String>();
		while((s=in.readLine())!=null){
			if(s.startsWith("Contig")) continue;
			String[] token = s.split("\t");
			//if(!token[6].startsWith("NM_5")) continue;
		//	if(!token[7].equals("T")) continue;
			if(Double.parseDouble(token[9]) > 0){
				/*String[] st = token[8].split(";");
				boolean write = false;
				for(int t : targetGroup){
					if(targets.contains(st[t].toUpperCase())){
						write = true;
			//			break;
					}
					else{
						write = false;
						break;
					}
				}
				boolean write2 = false;
				for(int c : controlGroup){
					if(controls.contains(st[c].toUpperCase())){
						write2 = true;
						//break;
					}
					else{
						write2 = false;
						break;
					}
				}
				if(write && write2)*/
					
					
					falseSet.add(token[4] + "\t" + token[5]);			
			}else if(Double.parseDouble(token[9]) == 0.0)
			trueSet.add(token[4] + "\t" + token[5]);
			
		}
		
		for(String f : falseSet){
			falseOut.println(f);
		}
		System.out.println(falseSet.size());
		for(String t : trueSet){
			if(falseSet.contains(t)) continue;
			trueOut.println(t);
		}
		
		falseOut.close();
		trueOut.close();
		in.close();
	}
}
