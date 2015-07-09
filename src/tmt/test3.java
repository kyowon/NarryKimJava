package tmt;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import msutil.Peptide;
import parser.BufferedLineReader;

public class test3 {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		PrintStream out = new PrintStream("/media/kyowon/Data1/MassSpec/tmt2/protein.tsv");
		BufferedLineReader merged = new BufferedLineReader("/media/kyowon/Data1/MassSpec/tmt2/merged.filtered.tsv");
		String s;
		HashMap<String, Integer> npsm = new HashMap<String, Integer>();
		HashMap<String, HashSet<String>> npep = new HashMap<String, HashSet<String>>();
		HashMap<String, double[]> intensity = new HashMap<String, double[]>();
		
		while((s=merged.readLine())!=null){
			if(s.startsWith("ID")) continue;
			String[] token = s.split("\t");
			
			String pro = token[11];
			String pep = new Peptide(token[9]).toString().toUpperCase();
			double[] ints = new double[6];
			for(int i=21;i<27;i++){
				ints[i-21] = Double.parseDouble(token[i]);
			}
			
			if(!npsm.containsKey(pro)){
				npsm.put(pro, 0);
				npep.put(pro, new HashSet<String>());
				intensity.put(pro, new double[6]);
			}
			
			npsm.put(pro, npsm.get(pro) + 1);
			npep.get(pro).add(pep);
			double[] intsb = intensity.get(pro);
			for(int i=0;i<6;i++){
				intsb[i] += ints[i];
			}
		}
		
		out.println("Protein\t#Spectra\t#Peptides\tIon1\tIon2\tIon3\tIon4\tIon5\tIon6");
		
		for(String pro : npsm.keySet()){
			out.print(pro + "\t" + npsm.get(pro) + "\t" + npep.get(pro).size());
			double[] ints = intensity.get(pro);
			for(int i=0;i<6;i++){
				out.print("\t" + ints[i]);
			}
			out.println();
			
		}
		
		
		out.close();
		merged.close();
	}

}
