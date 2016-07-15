package msgfpipeline;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import msgf.Tolerance;
import msgfpipeline.parser.MSGFPlusParser;
import msgfpipeline.parser.MSGFPlusParser.MSGFPlusPSM;
import msutil.Enzyme;

public class MSGFPipeline {
	
	public MSGFPipeline(String inputFilename, String outputFilename, String outputProteinFilename,
			String outputProteinFastaNameWithSharedPeptides, String outputProteinFastaNameWithoutSharedPeptides,
			String fasta, Enzyme enzyme, double qvalueThreshold, Tolerance tolerance, String[] inclusionProteinPrefixSet){
		MSGFPlusParser parser = new MSGFPlusParser(inputFilename, fasta, enzyme);
		System.out.println("All psms : " + parser.getPsms().size());
		int decoyHit = 0;
		for(MSGFPlusPSM psm : parser.getPsms()){
			if(psm.isDecoyHit()) decoyHit++;
		}
		System.out.println("Decoy psms : " + decoyHit);
		
		OffsetCalibrator calibrator = new OffsetCalibrator(parser.getPsms(), tolerance, qvalueThreshold);
		ArrayList<MSGFPlusPSM> filteredPSMs = calibrator.getFilteredPSMs();
		//qvalueThreshold = 1;
		
		//filteredPSMs = parser.getPsms();
		if(filteredPSMs == null || filteredPSMs.isEmpty()){
			System.out.println("Precursor filtered psms : 0 - output not generated");
			return;
		}
		System.out.println("Precursor filtered psms : " + filteredPSMs.size());
		
		if(inclusionProteinPrefixSet != null){
			FdrCalculator fdr = new FdrCalculator(filteredPSMs, outputFilename + ".exclude.tmp", inclusionProteinPrefixSet);
			fdr.calculate();
			//new File(outputFilename + ".exclude.tmp").deleteOnExit();
			ArrayList<MSGFPlusPSM> excludedPSMs = new ArrayList<MSGFPlusPSM>();
			for(MSGFPlusPSM psm : filteredPSMs){
				boolean toput = true;
				for(String pid : psm.getProteinIDs()){
					boolean toput2 = false;					
					for(String in : inclusionProteinPrefixSet){
						if(pid.startsWith(in) || pid.startsWith("XXX_"+in)){
							toput2 = true;	
							break;
						}
					}	
					toput = toput && toput2;
					if(!toput) break;
				}
				if(!toput && psm.getQvalue() <= qvalueThreshold) continue;
				excludedPSMs.add(psm);
			}
			System.out.println("Protein filtered psms : " + excludedPSMs.size());
			filteredPSMs = excludedPSMs;
		}
				
		HashSet<String> peptides = new HashSet<String>();
		FdrCalculator fdr = new FdrCalculator(filteredPSMs, outputFilename + ".tmp");
		fdr.calculate();
		//new File(outputFilename + ".tmp").deleteOnExit();
		PrintStream out;
		//int n = 0;
		try {
			out = new PrintStream(outputFilename);
			out.println(MSGFPlusParser.getHeader());
			ArrayList<MSGFPlusPSM>  qualifiedPSMs = new ArrayList<MSGFPlusPSM>();
			for(MSGFPlusPSM psm : filteredPSMs){
				if(psm.getQvalue() <= qvalueThreshold){
					boolean toput = true;
					if(inclusionProteinPrefixSet != null){
						if(psm.getProteinIDs() != null && psm.getProteinIDs().length > 0){
							for(String pn : psm.getProteinIDs()){
								boolean toput2 = false;
								for(String in : inclusionProteinPrefixSet){		
									if(pn.startsWith(in)){
										toput2 = true;
										break;
									}//else toput = false;
								}
								toput = toput && toput2;
								if(!toput) break;
							}
						}else toput = false;
					}
					
					if(toput){
						qualifiedPSMs.add(psm);
						peptides.add(psm.getPeptide());
						out.println(psm);
					}
				//	n ++;
				}
				
			}
			System.out.println("Qualified psms : " + qualifiedPSMs.size() + " Uniq Pep count : " + peptides.size());
			out.close();				
			
			new ProteinListFormatter(qualifiedPSMs, fasta, outputProteinFilename, outputProteinFastaNameWithSharedPeptides, outputProteinFastaNameWithoutSharedPeptides);
			
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	private static void printUsageAndExit(){
		System.out.println("java -jar -Xmx5g MSGFPipeline.jar [input dir containing csv] [output dir] [output protein dir] [fasta used for MSGF search] [enzyme ID] [FDR threshold] [ms1 tolerance]");
		System.out.println("enzyme ID : 1 - trypsin, 2 - chymotrypsin, 3 - LysC, 4 - LysN, 5 - GluC, 6 - ArgC, 7 - AspN");
		System.exit(0);
	}
	
	

	public static void main(String[] args){
		//boolean test = true;
		if(args.length < 7) printUsageAndExit();
		int i = 0;
		String input = args[i++];//"/media/kyowon/Data1/MSGFPipeline/out.tsv";
		String output = args[i++]; //"/media/kyowon/Data1/MSGFPipeline/out.reformatted.tsv";
		String outputProtein = args[i++]; //"/media/kyowon/Data1/MSGFPipeline/out.reformatted.tsv";
		String fasta = args[i++];// "/media/kyowon/Data1/MSGFPipeline/uniprot_sprot.revCat.fasta";
//h_;e_;1;2;3;4;5;6;7;8;9
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
		String[] inclusionProteinPrefixSet = null;
		if(args.length >i && args[i] != null){
			inclusionProteinPrefixSet = args[i].split(",");
			i++;
		}
		
		if(!fasta.contains("revCat")){
			fasta = fasta.substring(0, fasta.lastIndexOf('.')) + ".revCat" + fasta.substring(fasta.lastIndexOf('.'));
			System.out.println("fasta used : " + fasta);
		}
		
		//MSGFPlusPSM.proteinWeightThrehsold = 1.5e5;
		MSGFPlusPSM.setListAllProteins(true);
		
		if(!new File(output).exists()) new File(output).mkdir();
		if(!new File(outputProtein).exists()) new File(outputProtein).mkdir();
		
		for(File tsv : new File(input).listFiles()){
			if(tsv.getName().endsWith(".tsv")){
			//	if(new File(output+ System.getProperty("file.separator")+tsv.getName()).exists()) continue;
				System.out.println("Processing " + tsv.getName());
				new MSGFPipeline(tsv.getAbsolutePath(), output+ System.getProperty("file.separator")+tsv.getName(), 
						outputProtein+ System.getProperty("file.separator")+tsv.getName(), 
						outputProtein+ System.getProperty("file.separator")+tsv.getName() + ".shared.fasta",
						outputProtein+ System.getProperty("file.separator")+tsv.getName() + ".unshared.fasta",
						fasta, enzyme, qval, tol, inclusionProteinPrefixSet);
			}
		}	
	}	
}
