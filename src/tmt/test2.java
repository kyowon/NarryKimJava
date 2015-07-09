package tmt;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import parser.BufferedLineReader;

public class test2 {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		BufferedLineReader outFiltered = new BufferedLineReader("/media/kyowon/Data1/MassSpec/tmt2/originalMerged.filtered.tsv");
		PrintStream out = new PrintStream("/media/kyowon/Data1/MassSpec/tmt2/merged.filtered.tsv");
		BufferedLineReader merged = new BufferedLineReader("/media/kyowon/Data1/MassSpec/tmt2/merged.tsv");
			
		HashSet<String> specs = new HashSet<String>();
		String s;
		
		while((s=outFiltered.readLine())!=null){
			if(s.startsWith("#"))continue;
			String[] token = s.split("\t");
			boolean toWrite = true;
			for(int i=37;i<42;i++){
				if(!token[i].equals("0")){
					toWrite = false;
					break;
				}
			}
			if(!toWrite) continue;
			specs.add(token[0] + token[2]);
		}
		String header =  "";
		ArrayList<String> lines = new ArrayList<String>();
		double[] intensities = new double[6];
		
		while((s=merged.readLine())!=null){
			if(s.startsWith("ID")){
				header = s;
				continue;
			}
			String[] token = s.split("\t");
			if(specs.contains(token[1] + token[2])){
				for(int i=21;i<27;i++){
					intensities[i-21] += Double.parseDouble(token[i]);
				}
				lines.add(s);
			}
		}
		out.println(header);
		
		for(String l : lines){
			String[] t = l.split("\t");
			if(t[11].contains("XXX")) continue;
			for(int i=0;i<21;i++){
				if(i != 11)
					out.print(t[i]+"\t");
				else{
					t[i] = t[i].replace("\"", "");
					if(t[i].contains("("))
						out.print(t[i].substring(0, t[i].indexOf('('))+"\t");
					else out.print(t[i]+"\t");
				}
			}
			for(int i=21;i<27;i++){
				out.print(Double.parseDouble(t[i]) / intensities[i-21] * intensities[0] );
				out.print("\t");
			}
			for(int i=27;i<t.length-1;i++){
				out.print(t[i]+"\t");
			}
			out.println(t[t.length-1]);
		}
		
		out.close();
		merged.close();
		outFiltered.close();
	}

}
