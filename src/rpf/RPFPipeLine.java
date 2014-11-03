package rpf;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import launcher.PhyloPLauncher;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.MafParser;
import parser.ZeroBasedFastaParser;
import rpf.parser.ScoringOutputParser;
import rpf.parser.ScoringOutputParser.ScoredPosition;
import util.DnDsCalculator;

public class RPFPipeLine {
	//private static String[] harrScoreOutFiles;
	private static String[] harrParamFiles;
	private static String[] harrBedFiles; 
	
	private static String[] rpfScoreOutFiles;
	private static String[] rpfParamFiles;
	private static String[] rpfBedFiles; 
	private static String[] rnaBedFiles; 
	private static HashSet<String> allowedCodons;
	private static String mafFolder;
	private static ZeroBasedFastaParser fastaParser = null;
	private static AnnotationFileParser annotationFileParser = null;
	private static double scoreThreshold = .1;
	//private static double rpfRPKMThreshold = 0;
	//private static double rnaRPKMThreshold = 0;
	private static double rnaRPKMThreshold = 1;
	private static double rpfRPKMThreshold = 5;
	private static double harrRPKMThreshold = 5;
	private static String outFile = null;
	private static int[][] groups = null;
	//private static String outControlFile = null;// RPKM 6.8 ...
	private static int leftWindowSize = 90;
	private static int rightWindowSize = 300;
	
	private static int positionQuantityChangeLength = 51; // should be 3 mult
	private static int positionQuantityOffset = 13; //TODO automation
	private static int maxLengthUntilStopcodon = 13000;
	//TODO : memeory efficient + RNA coordinate calculation + union scoredPosition considering long splices...
	
	private static boolean parseSegment(String s, int mode){
		boolean ret = true;
		String[] token = s.split(" ");
		if(token.length == 0 || token[0].isEmpty()) return ret;
		if(mode == 1){
			harrBedFiles = new String[token.length];
			System.out.print("Harr bed files : ");
			for(int i=0; i<token.length; i++){
				harrBedFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 3){
			rpfBedFiles = new String[token.length];
			System.out.print("RPF bed files : ");		
			for(int i=0; i<token.length; i++){
				rpfBedFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 5){
			rnaBedFiles = new String[token.length];
			System.out.print("RNA bed files : ");					
			for(int i=0; i<token.length; i++){
				rnaBedFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 7){
			fastaParser = new ZeroBasedFastaParser(s);
			System.out.println("Fasta file: " + s);
		}else if(mode == 8){
			annotationFileParser = new AnnotationFileParser(s);
			System.out.println("Annotation file: " + s);
		}else if(mode == 9){
			scoreThreshold = Double.parseDouble(s);
			System.out.println("Score threshold: " + s);
		}else if(mode == 10){
			outFile = s;
			System.out.println("Output file: " + s);
		}else if(mode == 11){
		//	outControlFile = s;
		}else if(mode == 12){
			mafFolder = s;
			System.out.println("Maf containing folder: " + s);
		}else if(mode == 13){
			System.out.print("Start codons : ");		
			allowedCodons = new HashSet<String>();
			for(int i=0; i<token.length; i++){
				allowedCodons.add(token[i].toUpperCase());
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 14){
			harrRPKMThreshold = Double.parseDouble(s);
			System.out.println("Harr RPKM threshold (within " + positionQuantityChangeLength + " nts): " + s);
		}else if(mode == 15){
			rpfRPKMThreshold = Double.parseDouble(s);
			System.out.println("RPF RPKM threshold (within " + positionQuantityChangeLength + " nts): " + s);
		}else if(mode == 16){
			rnaRPKMThreshold = Double.parseDouble(s);
			System.out.println("RNA RPKM threshold (within " + positionQuantityChangeLength + " nts): " + s);
		}else if(mode == 17){
			DnDsCalculator.numSpecies = Integer.parseInt(s);
			System.out.println("Number of species in maf file: " + s);
		}else if(mode == 18){
			System.out.print("Groups : ");		
			groups = new int[token.length][];
			for(int i=0;i<groups.length;i++){
				String[] st = token[i].split(",");
				groups[i] = new int[st.length];
				for(int j=0;j<st.length;j++){
					groups[i][j] = Integer.parseInt(st[j]);
					System.out.print(groups[i][j] + " ");
				}
				System.out.print("; ");
			}
			System.out.println();
		}else if(mode == 19) {
			PhyloPLauncher.setModFile(s);
		} else ret = false;
		
		
		
		//rpfRPKMThreshold
		return ret;
	}
	
	private static boolean parseArgs(String[] args){
		if(args == null) return false;
		int mode = -1;
		String s = "";
		for(int i=0;i<args.length;i++){	
			if(args[i].startsWith("-")) {
				if(mode > 0) if (!parseSegment(s.trim(), mode)) 
					return false;
				s = "";
				mode = -1;
			}
			if(args[i].equals("-harr")){
				mode = 1;
			}else if(args[i].equals("-rpf")){
				mode = 3;
			}else if(args[i].equals("-rna")){
				mode = 5;
			}else if(args[i].equals("-fasta")){
				mode = 7;
			}else if(args[i].equals("-ref")){
				mode = 8;
			}else if(args[i].equals("-scoreThreshold")){
				mode = 9;
			}else if(args[i].equals("-outputFile")){				
				mode = 10;
			}else if(args[i].equals("-outputControlFile")){				
				mode = 11;
			}else if(args[i].equals("-maf")){				
				mode = 12;
			}else if(args[i].equals("-startCodons")){
				mode = 13;
			}else if(args[i].equals("-harrRPKMThreshold")){				
				mode = 14;
			}
				/*else if(args[i].equals("-rpfRPKMThreshold")){				
			}
				mode = 13;
			}else if(args[i].equals("-rnaRPKMThreshold")){				
				mode = 14;
			}*/else if(args[i].equals("-rpfRPKMThreshold")){				
				mode = 15;
			}else if(args[i].equals("-rnaRPKMThreshold")){				
				mode = 16;
			}else if(args[i].equals("-numSpeciesInMafFile")){				
				mode = 17;
			}else if(args[i].equals("-groups")){
				mode = 18;
			}else if(args[i].equals("-phyloPModFile")){
				mode = 19;
			}
			
			
			//
			else{
				s += args[i] + " ";
			}
		}
		if(mode > 0) if (!parseSegment(s.trim(), mode)) return false;
		if(rpfBedFiles == null || rpfBedFiles == null || fastaParser == null || annotationFileParser == null || outFile == null) return false;
		//if(harrCovPlusFiles.length != harrCovMinusFiles.length) return false;
		
		if(harrBedFiles != null){
			harrParamFiles = new String[harrBedFiles.length];
			//harrScoreOutFiles = new String[harrCovPlusFiles.length];
			for(int i=0;i<harrParamFiles.length;i++){
				harrParamFiles[i] = harrBedFiles[i] +".param";
				//harrScoreOutFiles[i] = harrCovPlusFiles[i]+".score"+ scoreThreshold +".tsv";
			}
		}
		
		rpfParamFiles = new String[rpfBedFiles.length];
		rpfScoreOutFiles = new String[rpfBedFiles.length];
		for(int i=0;i<rpfParamFiles.length;i++){
			rpfParamFiles[i] = rpfBedFiles[i] +".param";
			rpfScoreOutFiles[i] = rpfBedFiles[i] +".score"+ scoreThreshold +".tsv";
		}
		
		if(groups == null)
			groups = new int[rpfScoreOutFiles.length][1];
		
		return true;
	}
	
	public static void main(String[] args){
		//ClassLoader classLoader = ExcelFile.class.getClassLoader();
	    
		if(!parseArgs(args)){
			printUsage();
			System.exit(0);
		}
				
		if(allowedCodons == null){
			allowedCodons = new HashSet<String>();
			allowedCodons.add("ATG"); allowedCodons.add("CTG");			
		}
		if(harrBedFiles!=null){
			for(int i=0;i<harrBedFiles.length;i++){
				if(!new File(harrParamFiles[i]).exists()){
					ScorerTrainer trainer = new ScorerTrainer(harrBedFiles[i], annotationFileParser, harrParamFiles[i]);
					System.out.println("Training for " + harrBedFiles[i]);
					trainer.train(leftWindowSize, rightWindowSize, 30); // harr
					System.out.println("Training done..");
				}				
			}
		}
		for(int i=0;i<rpfBedFiles.length;i++){
			if(!new File(rpfParamFiles[i]).exists()){
				ScorerTrainer trainer = new ScorerTrainer(rpfBedFiles[i], annotationFileParser, rpfParamFiles[i]);
				System.out.println("Training for " + rpfBedFiles[i]);
				trainer.train(leftWindowSize, rightWindowSize, 30); // rpf
				System.out.println("Training done..");
			}
		}
		for(int i=0;i<rpfBedFiles.length;i++){
			if(!new File(rpfScoreOutFiles[i]).exists()){
				System.out.println("Scoring for " + rpfBedFiles[i]);
				PrintStream out;
				try {
					out = new PrintStream(rpfScoreOutFiles[i]);				
					for(String contig: annotationFileParser.getContigs()){
						ArrayList<ScoredPosition> scoredPositions = new ArrayList<ScoredPosition>();
						Bed12Parser rpfBedParser = new Bed12Parser(rpfBedFiles[i], annotationFileParser, contig);
						Bed12Parser rnaBedParser = new Bed12Parser(rnaBedFiles[i], annotationFileParser, contig);
						
						Scorer scorer = new Scorer(rpfBedParser, rpfParamFiles[i], annotationFileParser, fastaParser);
						scoredPositions.addAll(scorer.getScoredPositions(scoreThreshold, allowedCodons));
						
						Quantifier rpfQuantifier = new Quantifier(rpfBedParser, annotationFileParser, fastaParser);
						ArrayList<ScoredPosition> rpfFilteredScoredPositions = ScoringOutputParser.getQuantityFilteredPositions(scoredPositions, rpfQuantifier, rpfRPKMThreshold);
						
						Quantifier rnaQuantifier = new Quantifier(rnaBedParser, annotationFileParser, fastaParser);
						ArrayList<ScoredPosition> rnaFilteredScoredPositions = ScoringOutputParser.getQuantityFilteredPositions(rpfFilteredScoredPositions, rnaQuantifier, rnaRPKMThreshold);
						
						ArrayList<ScoredPosition> quantityFilteredScoredPositions = rnaFilteredScoredPositions;
						
						if(harrBedFiles!=null){
							if(harrBedFiles.length == rpfBedFiles.length){ // synced
								Bed12Parser harrBedParser = new Bed12Parser(harrBedFiles[i], annotationFileParser, contig);
								
								Quantifier harrQuantifier = new Quantifier(harrBedParser, annotationFileParser, fastaParser);
								quantityFilteredScoredPositions = ScoringOutputParser.getQuantityFilteredPositions(rnaFilteredScoredPositions, harrQuantifier, harrRPKMThreshold);
							
								
							}else{
								Quantifier[] harrQuantifiers = new Quantifier[harrBedFiles.length];
								for(int j=0;j<harrQuantifiers.length;j++){
									Bed12Parser harrBedParser = new Bed12Parser(harrBedFiles[j], annotationFileParser, contig);
									harrQuantifiers[j] = new Quantifier(harrBedParser, annotationFileParser, fastaParser);
								}	
								quantityFilteredScoredPositions = ScoringOutputParser.getQuantityFilteredPositions(rnaFilteredScoredPositions, harrQuantifiers, harrRPKMThreshold);
							}
						}
						
						for(ScoredPosition position : quantityFilteredScoredPositions){
							out.println(position);
						}						
					}
					out.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			
				//scorer.scoreNWrite(-0.1, fastaParser, scoreOutFiles[i], true);
				System.out.println("Scoring done..");
			}
		}
		MafParser mafParser = new MafParser(mafFolder);
		mafParser.generateIndexFile();
		mafParser.readIndexFile();
		//DsDnOutParser dsdnParser = new DsDnOutParser(dsdnFile, annotationFileParser);
		
		MergeResults merge = new MergeResults(rpfScoreOutFiles, harrBedFiles, rpfBedFiles, rnaBedFiles, harrParamFiles, rpfParamFiles, groups, annotationFileParser, fastaParser, mafParser); 

		System.out.println("Merging results");
		merge.merge(outFile, scoreThreshold, positionQuantityChangeLength, positionQuantityOffset, maxLengthUntilStopcodon, allowedCodons);

		System.out.println("Merging done..");
	}
	
	private static void printUsage(){
		System.out.println("Usage:\n"
				+ "java -jar -Xmx20g RPFPipeLine.jar "
				+ "\n-harr [harr bed files, seperated by space] "
				+ "\n-rpf [rpf bed files, seperated by space]"
				+ "\n-rna [rna bed files, seperated by space]"
				+ "\n-maf [maf file folder]"
				+ "\n-fasta [fasta file]"
				+ "\n-ref [refFlat file]"
				+ "\n-phyloPModFile [phyloP mod file]"
				+ "\n-groups [group info (0-based, sample seperated by , and group seperated by space  For example, 0,1,2,3 4,5,6 means 0 1 2 3 are one group and 4 5 6 are the other)]"
				+ "\n-scoreThreshold [scoreThreshold - default 0.3]"
				//+ "\n-rpfRPKMThreshold [rpf RPKM threshold - default 0]"
				//+ "\n-rnaRPKMThreshold [rna RPKM threshold - default 0]"
				+ "\n-rnaRPKMthreshold [rna RPKM threshold (" + rightWindowSize  +" nt window) - default 5]"
				+ "\n-rpfRPKMthreshold [rpf RPKM threshold (" + rightWindowSize + " nt window) - default 5]"
				+ "\n-harrRPKMThreshold [harr RPKM threshold (" + positionQuantityChangeLength + " nt window) - default 10]"
				+ "\n-startCodons [start codons to consider - default ATG and CTG, seperated by space]"
				+ "\n-numSpeciesInMafFile [number of species in maf file - default 46 (for hg19)]"
				+ "\n-outputFile [output file]"
				//+ "\n-outputControlFile [outputControl file] (results with poor harr score)"
				);
	}
}
