package msgfpipeline;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import msgf.Tolerance;
import msgfpipeline.parser.MSGFPlusParser;
import msgfpipeline.parser.MSGFPlusParser.MSGFPlusPSM;
import msutil.Enzyme;

public class MSGFPipeline {
	
	public MSGFPipeline(String inputFilename, String outputFilename, String outputProteinFilename, String fasta, Enzyme enzyme, double qvalueThreshold, Tolerance tolerance){
		MSGFPlusParser parser = new MSGFPlusParser(inputFilename, fasta, enzyme);
		System.out.println("All psms : " + parser.getPsms().size());
		OffsetCalibrator calibrator = new OffsetCalibrator(parser.getPsms(), tolerance, qvalueThreshold);
		ArrayList<MSGFPlusPSM> filteredPSMs = calibrator.getFilteredPSMs();
		//filteredPSMs = parser.getPsms();
		System.out.println("Precursor filtered psms : " + filteredPSMs.size());
		
		FdrCalculator fdr = new FdrCalculator(filteredPSMs, outputFilename + ".tmp");
		fdr.calculate();
		PrintStream out;
		//int n = 0;
		try {
			out = new PrintStream(outputFilename);
			out.println(MSGFPlusParser.getHeader());
			ArrayList<MSGFPlusPSM>  qualifiedPSMs = new ArrayList<MSGFPlusPSM>();
			for(MSGFPlusPSM psm : filteredPSMs){
				if(psm.getQvalue() < qvalueThreshold){
					qualifiedPSMs.add(psm);
					out.println(psm);
				//	n ++;
				}
				
			}
			System.out.println("Qualified psms : " + qualifiedPSMs.size());
			out.close();	
			//new ProteinListFormatter(qualifiedPSMs, outputProteinFilename);
			
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	

	public static void main(String[] args){	
		//boolean test = true;
		int i = 0;
		String input = args[i++];//"/media/kyowon/Data1/MSGFPipeline/out.tsv";
		String output = args[i++]; //"/media/kyowon/Data1/MSGFPipeline/out.reformatted.tsv";
		String outputProtein = args[i++]; //"/media/kyowon/Data1/MSGFPipeline/out.reformatted.tsv";
		String fasta = args[i++];// "/media/kyowon/Data1/MSGFPipeline/uniprot_sprot.revCat.fasta";

		Enzyme enzyme = Enzyme.TRYPSIN;
		int enzymeIndex = Integer.parseInt(args[i++]);
		if(enzymeIndex == 1){
			enzyme = Enzyme.TRYPSIN;
		}else if(enzymeIndex == 2){
			enzyme = Enzyme.CHYMOTRYPSIN;
		}else if(enzymeIndex == 3){
			enzyme = Enzyme.LysC;
		}else if(enzymeIndex == 4){
			enzyme = Enzyme.LysN;
		}else if(enzymeIndex == 5){
			enzyme = Enzyme.GluC;
		}else if(enzymeIndex == 6){
			enzyme = Enzyme.ArgC;
		}else if(enzymeIndex == 7){
			enzyme = Enzyme.AspN;
		}else enzyme = null;
		
		double qval = Double.parseDouble(args[i++]);
		
		Tolerance tol = Tolerance.parseToleranceStr(args[i++]);//new Tolerance(10f, true);
		
		MSGFPlusPSM.setListAllProteins(false);
		
		for(File tsv : new File(input).listFiles()){
			if(tsv.getName().endsWith(".tsv")){
				new MSGFPipeline(tsv.getAbsolutePath(), output+ System.getProperty("file.separator")+tsv.getName(), outputProtein+ System.getProperty("file.separator")+tsv.getName(), fasta, enzyme, qval, tol);
			}
		}
		
		
		
	}
	
}
