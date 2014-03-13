package rpf;

import java.io.File;

import parser.AnnotationFileParser;
import parser.MafParser;
import parser.ZeroBasedFastaParser;

public class RPFPipeLine {
	private static String[] harrScoreOutFiles;
	private static String[] harrParamFiles;
	private static String[] harrCovPlusFiles; 
	private static String[] harrCovMinusFiles;
	
	//private static String[] rpfScoreOutFiles;
	private static String[] rpfParamFiles;
	private static String[] rpfCovPlusFiles; 
	private static String[] rpfCovMinusFiles;
	private static String[] rnaCovPlusFiles; 
	private static String[] rnaCovMinusFiles;

	private static String mafFolder;
	private static ZeroBasedFastaParser fastaParser = null;
	private static AnnotationFileParser annotationFileParser = null;
	private static double scoreThreshold = .1;
	private static double rpfRPKMThreshold = 0;
	private static double rnaRPKMThreshold = 5;
	private static String outFile = null;
	//private static String outControlFile = null;
	
	private static int positionQuantityChangeLength = 600;
	private static int positionQuantityOffset = 15;
	private static int maxLengthUntilStopcodon = 13000;
	
	
	private static boolean parseSegment(String s, int mode){
		boolean ret = true;
		String[] token = s.split(" ");
		if(mode == 1){
			harrCovPlusFiles = new String[token.length];
			for(int i=0; i<token.length; i++){
				harrCovPlusFiles[i] = token[i];
			}
		}else if(mode == 2){
			harrCovMinusFiles = new String[token.length];
			for(int i=0; i<token.length; i++){
				harrCovMinusFiles[i] = token[i];
			}
		}else if(mode == 3){
			rpfCovPlusFiles = new String[token.length];
			for(int i=0; i<token.length; i++){
				rpfCovPlusFiles[i] = token[i];
			}
		}else if(mode == 4){
			rpfCovMinusFiles = new String[token.length];
			for(int i=0; i<token.length; i++){
				rpfCovMinusFiles[i] = token[i];
			}
		}else if(mode == 5){
			rnaCovPlusFiles = new String[token.length];
			for(int i=0; i<token.length; i++){
				rnaCovPlusFiles[i] = token[i];
			}
		}else if(mode == 6){
			rnaCovMinusFiles = new String[token.length];
			for(int i=0; i<token.length; i++){
				rnaCovMinusFiles[i] = token[i];
			}
		}else if(mode == 7){
			fastaParser = new ZeroBasedFastaParser(s);
		}else if(mode == 8){
			annotationFileParser = new AnnotationFileParser(s);
		}else if(mode == 9){
			scoreThreshold = Double.parseDouble(s);
		}else if(mode == 10){
			outFile = s;
		}else if(mode == 11){
		//	outControlFile = s;
		}else if(mode == 12){
			mafFolder = s;
		}else if(mode == 13){
			rpfRPKMThreshold = Double.parseDouble(s);
		}else if(mode == 14){
			rnaRPKMThreshold = Double.parseDouble(s);
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
			}else if(args[i].equals("-rpfRPKMThreshold")){				
				mode = 13;
			}else if(args[i].equals("-rnaRPKMThreshold")){				
				mode = 13;
			}else{
				s += args[i] + " ";
			}
		}
		if(mode > 0) if (!parseSegment(s.trim(), mode)) return false;
		if(harrCovPlusFiles == null || harrCovMinusFiles == null || fastaParser == null || annotationFileParser == null || outFile == null) return false;
		//if(harrCovPlusFiles.length != harrCovMinusFiles.length) return false;
		
		harrParamFiles = new String[harrCovPlusFiles.length];
		harrScoreOutFiles = new String[harrCovPlusFiles.length];
		for(int i=0;i<harrParamFiles.length;i++){
			harrParamFiles[i] = harrCovPlusFiles[i]+".param";
			harrScoreOutFiles[i] = harrCovPlusFiles[i]+".score"+ scoreThreshold +".tsv";
		}
		
		rpfParamFiles = new String[rpfCovPlusFiles.length];
		//rpfScoreOutFiles = new String[rpfCovPlusFiles.length];
		for(int i=0;i<rpfParamFiles.length;i++){
			rpfParamFiles[i] = rpfCovPlusFiles[i]+".param";
	//		rpfScoreOutFiles[i] = rpfCovPlusFiles[i]+".score"+ scoreThreshold +".tsv";
		}
		
		
		return true;
	}
	
	public static void main(String[] args){
		
		if(!parseArgs(args)){
			printUsage();
			System.exit(0);
		}
				
		for(int i=0;i<harrCovPlusFiles.length;i++){
			if(!new File(harrParamFiles[i]).exists()){
				Trainer trainer = new Trainer(harrCovPlusFiles[i], harrCovMinusFiles[i], annotationFileParser, harrParamFiles[i]);
				System.out.println("Training for " + harrCovPlusFiles[i] + " and " + harrCovMinusFiles[i]);
				trainer.train(60, 90, 10); // harr
				System.out.println("Training done..");
			}
			if(!new File(harrScoreOutFiles[i]).exists()){
				System.out.println("Scoring for " + harrCovPlusFiles[i] + " and " + harrCovMinusFiles[i]);
				Scorer scorer = new Scorer(harrCovPlusFiles[i], harrCovMinusFiles[i], harrParamFiles[i], annotationFileParser, fastaParser);
				scorer.scoreNWrite(scoreThreshold, harrScoreOutFiles[i]+".beforeWindowFilter.tsv", false);
				scorer.writeWindowFilteredOutput(harrScoreOutFiles[i]+".beforeWindowFilter.tsv", harrScoreOutFiles[i], 30);
				//scorer.scoreNWrite(-0.1, fastaParser, scoreOutFiles[i], true);
				System.out.println("Scoring done..");
			}
		}
		
		for(int i=0;i<rpfCovPlusFiles.length;i++){
			if(!new File(rpfParamFiles[i]).exists()){
				Trainer trainer = new Trainer(rpfCovPlusFiles[i], rpfCovMinusFiles[i], annotationFileParser, rpfParamFiles[i]);
				System.out.println("Training for " + rpfCovPlusFiles[i] + " and " + rpfCovMinusFiles[i]);
				trainer.train(60, 150, 20); // rpf
				System.out.println("Training done..");
			}
			/*if(!new File(rpfScoreOutFiles[i]).exists()){
				System.out.println("Scoring for " + rpfCovPlusFiles[i] + " and " + rpfCovMinusFiles[i]);
				Scorer scorer = new Scorer(rpfCovPlusFiles[i], rpfCovMinusFiles[i], rpfParamFiles[i], annotationFileParser);
				scorer.scoreNWrite(scoreThreshold, fastaParser, rpfScoreOutFiles[i]+".beforeWindowFilter.tsv", false);
				scorer.writeWindowFilteredOutput(rpfScoreOutFiles[i]+".beforeWindowFilter.tsv", rpfScoreOutFiles[i], 25);
				//scorer.scoreNWrite(-0.1, fastaParser, scoreOutFiles[i], true);
				System.out.println("Scoring done..");
			}*/
		}
		MafParser mafParser = new MafParser(mafFolder, annotationFileParser);
		mafParser.generateIndexFile();
		mafParser.readIndexFile();
		//DsDnOutParser dsdnParser = new DsDnOutParser(dsdnFile, annotationFileParser);
		
		MergeResults merge = new MergeResults(harrScoreOutFiles, harrCovPlusFiles, harrCovMinusFiles, rpfCovPlusFiles, rpfCovMinusFiles, rnaCovPlusFiles, rnaCovMinusFiles, harrParamFiles, rpfParamFiles, annotationFileParser, fastaParser, mafParser);
		System.out.println("Merging results");
		merge.merge(outFile, scoreThreshold, rpfRPKMThreshold, rnaRPKMThreshold, positionQuantityChangeLength, positionQuantityOffset, maxLengthUntilStopcodon);
	//	merge = null;
	//	MergeResults mergeControl = new MergeResults(scoreOutFiles, harrCovPlusFiles, harrCovMinusFiles, rpfCovPlusFiles, rpfCovMinusFiles, rnaCovPlusFiles, rnaCovMinusFiles, paramFiles, annotationFileParser, fastaParser, mafParser);
	//	System.out.println("Merging results - control");
	//	mergeControl.merge(outControlFile, -0.05);
		
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
				+ "\n-rpfRPKMThreshold [rpf RPKM threshold - default 20]"
				+ "\n-rnaRPKMThreshold [rna RPKM threshold - default 20]"
				+ "\n-outputFile [output file]"
				//+ "\n-outputControlFile [outputControl file] (results with poor harr score)"
				);
	}
}
