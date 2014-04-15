package parser;

import java.io.FileNotFoundException;
import java.io.PrintStream;

import msutil.Enzyme;

public class GenerateModAParameterFile {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String specFile = args[0];
		String fasta = args[1];
		String mode = args[2];
		String min = args[3];
		String max = args[4];
		
		Enzyme enzyme = null;
		int ei = Integer.parseInt(args[5]);
		if(ei == 0) enzyme = null;
		else if(ei == 1) enzyme = Enzyme.TRYPSIN;
		else if(ei == 2) enzyme = Enzyme.CHYMOTRYPSIN;
		else if(ei == 3) enzyme = Enzyme.LysC;
		else if(ei == 4) enzyme = Enzyme.LysN;
		else if(ei == 5) enzyme = Enzyme.ArgC;
		else if(ei == 6) enzyme = Enzyme.AspN;
		//Enzyme.TRYPSIN.getResidues().get(0).getResidue() + " " + Enzyme.TRYPSIN.isCTerm()
		boolean isCPlus57 = args[6].equals("true");
		
		try {
			PrintStream out = new PrintStream(args[7]);
			out.println("Spectra=" + specFile);
			out.println("Fasta=" + fasta);
			out.println("BlindMode=" + mode);
			out.println("MinModSize=" + min);
			out.println("MaxModSize=" + max);
			if(enzyme != null){
				out.print("Enzyme=etc, ");
				for(msutil.AminoAcid aa : enzyme.getResidues()){
					out.print(aa.getResidue());
				}
				out.println("/" + (enzyme.isNTerm()? "N" : "C"));
			}
			if(isCPlus57) out.println("ADD=C, 57.021464");
			out.println("PeptTolerance=0.5\nAutoPMCorrection=0\nFragTolerance=0.6\nMissedCleavage=0\n");
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

}
