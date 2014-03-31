package rpf;

import java.io.File;
import java.util.HashSet;

import parser.AnnotationFileParser;
import parser.MafParser;
import parser.ZeroBasedFastaParser;

public class RPFPipeLine {
	//private static String[] harrScoreOutFiles;
	private static String[] harrParamFiles;
	private static String[] harrCovPlusFiles; 
	private static String[] harrCovMinusFiles;
	
	private static String[] rpfScoreOutFiles;
	private static String[] rpfParamFiles;
	private static String[] rpfCovPlusFiles; 
	private static String[] rpfCovMinusFiles;
	private static String[] rnaCovPlusFiles; 
	private static String[] rnaCovMinusFiles;
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
	//private static String outControlFile = null;// RPKM 6.8 ...
	
	private static int positionQuantityChangeLength = 50;
	private static int positionQuantityOffset = 13; //TODO automation
	private static int maxLengthUntilStopcodon = 13000;
	
	
	private static boolean parseSegment(String s, int mode){
		boolean ret = true;
		String[] token = s.split(" ");
		if(mode == 1){
			harrCovPlusFiles = new String[token.length];
			System.out.print("Harr cov plus files : ");
			for(int i=0; i<token.length; i++){
				harrCovPlusFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 2){
			harrCovMinusFiles = new String[token.length];
			System.out.print("Harr cov minus files : ");			
			for(int i=0; i<token.length; i++){
				harrCovMinusFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 3){
			rpfCovPlusFiles = new String[token.length];
			System.out.print("RPF cov plus files : ");		
			for(int i=0; i<token.length; i++){
				rpfCovPlusFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 4){
			rpfCovMinusFiles = new String[token.length];
			System.out.print("RPF cov minus minus : ");		
			for(int i=0; i<token.length; i++){
				rpfCovMinusFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 5){
			rnaCovPlusFiles = new String[token.length];
			System.out.print("RNA cov plus files : ");					
			for(int i=0; i<token.length; i++){
				rnaCovPlusFiles[i] = token[i];
				System.out.print(token[i] + " ");
			}
			System.out.println();
		}else if(mode == 6){
			rnaCovMinusFiles = new String[token.length];
			System.out.print("RNA cov minus files : ");		
			for(int i=0; i<token.length; i++){
				rnaCovMinusFiles[i] = token[i];
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
		}else ret = false;
		
		
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
			if(args[i].equals("-harrP")){
				mode = 1;
			}else if(args[i].equals("-harrM")){				
				mode = 2;
			}else if(args[i].equals("-rpfP")){
				mode = 3;
			}else if(args[i].equals("-rpfM")){
				mode = 4;
			}else if(args[i].equals("-rnaP")){
				mode = 5;
			}else if(args[i].equals("-rnaM")){
				mode = 6;
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
			}else{
				s += args[i] + " ";
			}
		}
		if(mode > 0) if (!parseSegment(s.trim(), mode)) return false;
		if(harrCovPlusFiles == null || harrCovMinusFiles == null || fastaParser == null || annotationFileParser == null || outFile == null) return false;
		//if(harrCovPlusFiles.length != harrCovMinusFiles.length) return false;
		
		harrParamFiles = new String[harrCovPlusFiles.length];
		//harrScoreOutFiles = new String[harrCovPlusFiles.length];
		for(int i=0;i<harrParamFiles.length;i++){
			harrParamFiles[i] = harrCovPlusFiles[i]+".param";
			//harrScoreOutFiles[i] = harrCovPlusFiles[i]+".score"+ scoreThreshold +".tsv";
		}
		
		rpfParamFiles = new String[rpfCovPlusFiles.length];
		rpfScoreOutFiles = new String[rpfCovPlusFiles.length];
		for(int i=0;i<rpfParamFiles.length;i++){
			rpfParamFiles[i] = rpfCovPlusFiles[i]+".param";
			rpfScoreOutFiles[i] = rpfCovPlusFiles[i]+".score"+ scoreThreshold +".tsv";
		}
		
		
		return true;
	}
	
	public static void main(String[] args){
		
		if(!parseArgs(args)){
			printUsage();
			System.exit(0);
		}
				
		if(allowedCodons == null){
			allowedCodons = new HashSet<String>();
			allowedCodons.add("ATG"); allowedCodons.add("CTG");			
		}
		for(int i=0;i<harrCovPlusFiles.length;i++){
			if(!new File(harrParamFiles[i]).exists()){
				ScorerTrainer trainer = new ScorerTrainer(harrCovPlusFiles[i], harrCovMinusFiles[i], annotationFileParser, harrParamFiles[i]);
				System.out.println("Training for " + harrCovPlusFiles[i] + " and " + harrCovMinusFiles[i]);
				trainer.train(60, 40, 6); // harr
				System.out.println("Training done..");
			}
			/*if(!new File(harrScoreOutFiles[i]).exists()){
				System.out.println("Scoring for " + harrCovPlusFiles[i] + " and " + harrCovMinusFiles[i]);
				Scorer scorer = new Scorer(harrCovPlusFiles[i], harrCovMinusFiles[i], harrParamFiles[i], annotationFileParser, fastaParser);
				scorer.scoreNWrite(scoreThreshold, harrScoreOutFiles[i]+".beforeWindowFilter.tsv", false);
				scorer.writeWindowFilteredOutput(harrScoreOutFiles[i]+".beforeWindowFilter.tsv", harrScoreOutFiles[i], 30);
				//scorer.scoreNWrite(-0.1, fastaParser, scoreOutFiles[i], true);
				System.out.println("Scoring done..");
			}*/
		}
		
		for(int i=0;i<rpfCovPlusFiles.length;i++){
			if(!new File(rpfParamFiles[i]).exists()){
				ScorerTrainer trainer = new ScorerTrainer(rpfCovPlusFiles[i], rpfCovMinusFiles[i], annotationFileParser, rpfParamFiles[i]);
				System.out.println("Training for " + rpfCovPlusFiles[i] + " and " + rpfCovMinusFiles[i]);
				trainer.train(150, 150, 10); // rpf
				System.out.println("Training done..");
			}
			if(!new File(rpfScoreOutFiles[i]).exists()){
				System.out.println("Scoring for " + rpfCovPlusFiles[i] + " and " + rpfCovMinusFiles[i]);
				Scorer scorer = new Scorer(rpfCovPlusFiles[i], rpfCovMinusFiles[i], rpfParamFiles[i], annotationFileParser, fastaParser);
				if(!new File(rpfScoreOutFiles[i]+".beforeFilter.tsv").exists())
					scorer.scoreNWrite(scoreThreshold, allowedCodons, rpfScoreOutFiles[i]+".beforeFilter.tsv", false);
				
				Quantifier rpfQuantifier = new Quantifier(rpfCovPlusFiles[i], rpfCovMinusFiles[i], annotationFileParser, fastaParser);
				scorer.writeQuantityFilteredOutput(rpfScoreOutFiles[i]+".beforeFilter.tsv", rpfScoreOutFiles[i]+".rpfQuantityFiltered.tsv", rpfQuantifier, rpfRPKMThreshold, positionQuantityChangeLength, positionQuantityOffset);
				
				Quantifier rnaQuantifier = new Quantifier(rnaCovPlusFiles[i], rnaCovMinusFiles[i], annotationFileParser, fastaParser);
				scorer.writeQuantityFilteredOutput(rpfScoreOutFiles[i]+".rpfQuantityFiltered.tsv", rpfScoreOutFiles[i]+".rnaQuantityFiltered.tsv", rnaQuantifier, rnaRPKMThreshold, positionQuantityChangeLength, 0);	
				
				if(harrCovPlusFiles.length == rpfCovPlusFiles.length){ // synced
					Quantifier harrQuantifier = new Quantifier(harrCovPlusFiles[i], harrCovMinusFiles[i], annotationFileParser, fastaParser);
					scorer.writeQuantityFilteredOutput(rpfScoreOutFiles[i]+".rnaQuantityFiltered.tsv", rpfScoreOutFiles[i]+".harrQuantityFiltered.tsv", harrQuantifier, harrRPKMThreshold, positionQuantityChangeLength, positionQuantityOffset);	
				}else{
					Quantifier[] harrQuantifiers = new Quantifier[harrCovPlusFiles.length];
					for(int j=0;j<harrQuantifiers.length;j++){
						harrQuantifiers[j] = new Quantifier(harrCovPlusFiles[j], harrCovMinusFiles[j], annotationFileParser, fastaParser);
					}
					scorer.writeQuantityFilteredOutput(rpfScoreOutFiles[i]+".rnaQuantityFiltered.tsv", rpfScoreOutFiles[i]+".harrQuantityFiltered.tsv", harrQuantifiers, harrRPKMThreshold, positionQuantityChangeLength, positionQuantityOffset);
				}
				
				scorer.writeWindowFilteredOutput(rpfScoreOutFiles[i]+".harrQuantityFiltered.tsv", rpfScoreOutFiles[i], positionQuantityChangeLength);
				//scorer.scoreNWrite(-0.1, fastaParser, scoreOutFiles[i], true);
				System.out.println("Scoring done..");
			}
		}
		MafParser mafParser = new MafParser(mafFolder, annotationFileParser);
		mafParser.generateIndexFile();
		mafParser.readIndexFile();
		//DsDnOutParser dsdnParser = new DsDnOutParser(dsdnFile, annotationFileParser);
		
		MergeResults merge = new MergeResults(rpfScoreOutFiles, harrCovPlusFiles, harrCovMinusFiles, rpfCovPlusFiles, rpfCovMinusFiles, rnaCovPlusFiles, rnaCovMinusFiles, harrParamFiles, rpfParamFiles, annotationFileParser, fastaParser, mafParser);
		System.out.println("Merging results");
		merge.merge(outFile, scoreThreshold, positionQuantityChangeLength, positionQuantityOffset, maxLengthUntilStopcodon, allowedCodons);

		System.out.println("Merging done..");
	}
	
	private static void printUsage(){
		System.out.println("Usage:\n"
				+ "java -jar -Xmx5g RPFPipeLine.jar "
				+ "\n-harrP [harr coverage files (plus strand), seperated by space] "
				+ "\n-harrM [harr coverage files (minus strand), seperated by space]"
				+ "\n-rpfP [rpf coverage files (plus strand), seperated by space]"
				+ "\n-rpfM [rpf coverage files(minus strand), seperated by space]"
				+ "\n-rnaP [rna coverage files (plus strand), seperated by space]"
				+ "\n-rnaM [rna coverage files(minus strand), seperated by space]"
				+ "\n-maf [maf file folder]"
				+ "\n-fasta [fasta file]"
				+ "\n-ref [refFlat file]"
				+ "\n-scoreThreshold [scoreThreshold - default 0.3]"
				//+ "\n-rpfRPKMThreshold [rpf RPKM threshold - default 0]"
				//+ "\n-rnaRPKMThreshold [rna RPKM threshold - default 0]"
				+ "\n-rnaRPKMthreshold [rna RPKM threshold (" + positionQuantityChangeLength +" nt window) - default 5]"
				+ "\n-rpfRPKMthreshold [rpf RPKM threshold (" + positionQuantityChangeLength + " nt window) - default 5]"
				+ "\n-harrRPKMThreshold [harr RPKM threshold (" + positionQuantityChangeLength + " nt window) - default 10]"
				+ "\n-startCodons [start codons to consider - default ATG and CTG, seperated by space]"
				+ "\n-outputFile [output file]"
				//+ "\n-outputControlFile [outputControl file] (results with poor harr score)"
				);
	}
}
