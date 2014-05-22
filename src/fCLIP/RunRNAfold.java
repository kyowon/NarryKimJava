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
			
			int si = result.indexOf('(');
			int ei = result.lastIndexOf('(');
			double max = 0;
			
			if(si <0 || ei < 0){
				
			}else{
				String st = result.substring(si, ei);
				int d = 0;
				for(int i=0;i<st.length();i++){
					char c = st.charAt(i);
					if(c == '(') d ++;
					else if(c == ')') d --;
					max = max > d ? max : d;
				}
			}
			
			ArrayList<Double> ret = new ArrayList<Double>();
			//System.out.println(result + " " + max);
			int i = result.lastIndexOf(')');
			if(ei >= 0 && i >= 0)
				ret.add(Double.parseDouble(result.substring(ei + 1 , i)));
			else ret.add(0.0);
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
