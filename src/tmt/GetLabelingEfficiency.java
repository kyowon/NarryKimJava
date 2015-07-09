package tmt;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;


public class GetLabelingEfficiency {
	public static void main(String[] args) {
		try {
			BufferedLineReader in = new BufferedLineReader((args[0]));
			PrintStream out = new PrintStream(args[1]);
			String massShift = args[2];
			String s;
			int nLabeled = 0;
			int lysLabeled = 0;
			int n = 0;
			int lys = 0;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				if(Float.parseFloat(token[15]) > 0.01) continue;
				String pep = token[8];
				n++;
				if(pep.startsWith(massShift)) nLabeled++;
				
				for(int i=0;i<pep.length();i++){
					if(pep.charAt(i) == 'K') lys ++;
					if(i+massShift.length()+1<=pep.length() && pep.substring(i+1, i+massShift.length()+1).equals(massShift)) lysLabeled ++;
				}
				
			}
						
			out.println("Number of Peptides : " + n);
			out.println("Number of Labeled N-terms : " + nLabeled);
			out.println("N-term efficiency : " + (float)nLabeled/n);
			out.println("Number of Lys : " + lys);
			out.println("Number of Labeled Lys : " + lysLabeled);
			out.println("Lys efficiency : " + (float)lysLabeled/lys);
			out.close();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

}
