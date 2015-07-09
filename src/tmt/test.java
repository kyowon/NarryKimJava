package tmt;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class test {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			BufferedLineReader merged = new BufferedLineReader("/media/kyowon/Data1/MassSpec/tmt2/merged.tsv");
			File originalFolder = new File("/media/kyowon/Data1/MassSpec/tmt2/original");
			PrintStream out = new PrintStream("/media/kyowon/Data1/MassSpec/tmt2/originalMerged.tsv");
			
			HashSet<String> specID = new HashSet<String>();
			String s;
			while((s=merged.readLine())!=null){
				if(s.startsWith("ID")) continue;
				String[] token = s.split("\t");
				specID.add(token[1] + token[2]);
			}
			
			
			merged.close();
			out.println("#SpecFile	SpecID	ScanNum	FragMethod	Precursor	IsotopeError	PrecursorError(ppm)	Charge	Peptide	Protein	DeNovoScore	MSGFScore	SpecEValue	EValue	QValue	PepQValue	BasePeakIntensity	BasePeakMZ	ReporterIonIntensityMax	Ion_126	Ion_127	Ion_128	Ion_129	Ion_130	Ion_131	Weighted Avg Pct Intensity Correction	Ion_126_ObsMZ	Ion_127_ObsMZ	Ion_128_ObsMZ	Ion_129_ObsMZ	Ion_130_ObsMZ	Ion_131_ObsMZ");
			for(File tsv : originalFolder.listFiles()){
				if(!tsv.getName().endsWith(".tsv")) continue;
				BufferedLineReader r = new BufferedLineReader(tsv.getAbsolutePath());
				String t;
				while((t=r.readLine())!=null){
					if(t.startsWith("#")) continue;
					String[] token = t.split("\t");
					if(token[9].contains(";")) continue; // remove shared peptides
					if(!token[8].contains("+229")) continue;
					if(Double.parseDouble(token[14]) > .01) continue;
					if(!specID.contains(token[0]+token[2])) continue;
					out.println(t);
				}
				
				r.close();
			}
			
			out.close();
			
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
