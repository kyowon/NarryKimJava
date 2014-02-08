package rpf;

import java.io.File;

import parser.AnnotationFileParser;
import parser.ZeroBasedFastaParser;

public class RPFPipeLine {
	private static String[] scoreOutFiles;
	private static String[] harrCovPlusFiles; 
	private static String[] harrCovMinusFiles;
	private static String[] rpfCovPlusFiles; 
	private static String[] rpfCovMinusFiles;
	private static String[] rnaCovPlusFiles; 
	private static String[] rnaCovMinusFiles;
	private static String[] paramFiles;
	private static ZeroBasedFastaParser fastaParser = null;
	private static AnnotationFileParser annotationFileParser = null;
	private static double scoreThreshold = 2;
	private static String outFile = null;
	private static String outControlFile = null;
	
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
			outControlFile = s;
		}else ret = false;
		
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
			}else{
				s += args[i] + " ";
			}
		}
		if(mode > 0) if (!parseSegment(s.trim(), mode)) return false;
		if(harrCovPlusFiles == null || harrCovMinusFiles == null || fastaParser == null || annotationFileParser == null || outFile == null) return false;
		//if(harrCovPlusFiles.length != harrCovMinusFiles.length) return false;
		
		paramFiles = new String[harrCovPlusFiles.length];
		scoreOutFiles = new String[harrCovPlusFiles.length];
		for(int i=0;i<paramFiles.length;i++){
			paramFiles[i] = harrCovPlusFiles[i]+".param";
			scoreOutFiles[i] = harrCovPlusFiles[i]+".score.tsv";
		}
		
		return true;
	}
	
	public static void main(String[] args){
		
		if(!parseArgs(args)){
			printUsage();
			System.exit(0);
		}
				
		for(int i=0;i<harrCovPlusFiles.length;i++){
			if(!new File(paramFiles[i]).exists()){
				MatchedFilterTrainier trainer = new MatchedFilterTrainier(harrCovPlusFiles[i], harrCovMinusFiles[i], annotationFileParser, paramFiles[i]);
				System.out.println("Training for " + harrCovPlusFiles[i] + " and " + harrCovMinusFiles[i]);
				trainer.train(30, 50, 5); // harr
				System.out.println("Training done..");
			}
			if(!new File(scoreOutFiles[i]).exists()){
				System.out.println("Scoring for " + harrCovPlusFiles[i] + " and " + harrCovMinusFiles[i]);
				Scorer scorer = new Scorer(harrCovPlusFiles[i], harrCovMinusFiles[i], paramFiles[i], annotationFileParser);
				scorer.scoreNWrite(scoreThreshold, fastaParser, scoreOutFiles[i], false);
				scorer.scoreNWrite(-0.3, fastaParser, scoreOutFiles[i], true);
				System.out.println("Scoring done..");
			}
		}
		MergeResults merge = new MergeResults(scoreOutFiles, harrCovPlusFiles, harrCovMinusFiles, rpfCovPlusFiles, rpfCovMinusFiles, rnaCovPlusFiles, rnaCovMinusFiles, paramFiles, annotationFileParser, fastaParser);
		System.out.println("Merging results");
		merge.merge(outFile, scoreThreshold);
		merge = null;
		MergeResults mergeControl = new MergeResults(scoreOutFiles, harrCovPlusFiles, harrCovMinusFiles, rpfCovPlusFiles, rpfCovMinusFiles, rnaCovPlusFiles, rnaCovMinusFiles, paramFiles, annotationFileParser, fastaParser);
		System.out.println("Merging results - control");
		mergeControl.merge(outControlFile, -0.3);
		
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
				+ "\n-fasta [fasta file]"
				+ "\n-ref [refFlat file]"
				+ "\n-scoreThreshold [scoreThreshold - default 2]"
				+ "\n-outputFile [output file]"
				+ "\n-outputControlFile [outputControl file] (results with poor harr score)");
	}
}
