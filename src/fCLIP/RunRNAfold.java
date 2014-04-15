package fCLIP;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class RunRNAfold {
	
	static ArrayList<Double> run(String seq){
		try {
			String[] cmd = {
					"/bin/sh",
					"-c",
					"echo "+ seq + " | RNAfold --noPS"
					};
			ProcessBuilder pr = new ProcessBuilder(cmd);
					 
			Process p = pr.start();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			StringBuilder builder = new StringBuilder();
			String line = null;
			while ( (line = br.readLine()) != null) {
			   builder.append(line);
			   builder.append(System.getProperty("line.separator"));
			}
			String result = builder.toString();
			
			String st = result.substring(result.indexOf('('), result.lastIndexOf('('));
			double max = 0;
			int d = 0;
			for(int i=0;i<st.length();i++){
				char c = st.charAt(i);
				if(c == '(') d ++;
				else if(c == ')') d --;
				max = max > d ? max : d;
			}
			
			ArrayList<Double> ret = new ArrayList<Double>();
			//System.out.println(result + " " + max);
			ret.add(Double.parseDouble(result.substring(result.lastIndexOf('(') + 1 , result.lastIndexOf(')'))));
			ret.add(max);
			return ret;
			//
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	static public void main(String[] args){
		System.out.println(RunRNAfold.run("TGCTGCTCTTATGCATGCAGTGATTTGAGGATTATTGCTCACGGTAAGAAAAACATGAGTAATATATTGTCTTCAAAGAACCTATGCTCTAACAGGAGAT"));
	}
	
}
