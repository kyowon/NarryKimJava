package rpf;

public class RPFPipeLine {
	private static String[][] plusCovFiles = null;
	private static String[][] minusCovFiles = null;
	private static String fasta = null;
	private static String refFlat = null;
	private static double scoreThreshold = 2;
	private static String outFile = null;
	private static String[][] paramFiles = null;
	private static String[] scoreOutFiles = null;
		
	
	private static boolean parseSegment(String s, int mode){
		String[] token = s.split(" ");
		if(mode == 1){
			if(plusCovFiles==null) plusCovFiles = new String[token.length][2];
			else if(plusCovFiles.length != token.length) return false;
			for(int i=0; i<token.length; i++){
				plusCovFiles[i][0] = token[i];
			}
		}else if(mode == 2){
			if(minusCovFiles==null) minusCovFiles = new String[token.length][2];
			else if(minusCovFiles.length != token.length) return false;
			for(int i=0; i<token.length; i++){
				minusCovFiles[i][0] = token[i];
			}
		}else if(mode == 3){
			if(plusCovFiles==null) plusCovFiles = new String[token.length][2];
			else if(plusCovFiles.length != token.length) return false;
			for(int i=0; i<token.length; i++){
				plusCovFiles[i][1] = token[i];
			}
		}else if(mode == 4){
			if(minusCovFiles==null) minusCovFiles = new String[token.length][2];
			else if(minusCovFiles.length != token.length) return false;
			for(int i=0; i<token.length; i++){
				minusCovFiles[i][1] = token[i];
			}
		}else if(mode == 5){
			fasta = s;
		}else if(mode == 6){
			refFlat = s;
		}else if(mode == 7){
			//annotatedSitesOnly = Boolean.parseBoolean(s);
		}else if(mode == 8){
			scoreThreshold = Double.parseDouble(s);
		}else if(mode == 9){
			outFile = s;
		}
		//System.out.println(mode + " " + s);
		return true;
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
			}else if(args[i].equals("-fasta")){
				mode = 5;
			}else if(args[i].equals("-ref")){
				mode = 6;
		//	}else if(args[i].equals("-annotatedSitesOnly")){
			//	mode = 7;
			}else if(args[i].equals("-scoreThreshold")){
				mode = 8;
			}else if(args[i].equals("-outputFile")){				
				mode = 9;
			}else{
				s += args[i] + " ";
			}
		}
		if(mode > 0) if (!parseSegment(s.trim(), mode)) return false;
		if(plusCovFiles == null || minusCovFiles == null || fasta == null || refFlat == null || outFile == null) return false;
		if(plusCovFiles.length != minusCovFiles.length) return false;
		
		paramFiles = new String[plusCovFiles.length][2];
		scoreOutFiles = new String[plusCovFiles.length];
		for(int i=0;i<paramFiles.length;i++){
			for(int j=0;j<paramFiles[i].length;j++){
				paramFiles[i][j] = plusCovFiles[i][j]+".param";
			}
			scoreOutFiles[i] = plusCovFiles[i][0]+".score.tsv";
		}
		
		return true;
	}
	
	public static void main(String[] args){
		
		if(!parseArgs(args)){
			printUsage();
			System.exit(0);
		}
				
		for(int i=0;i<plusCovFiles.length;i++){
			for(int j=0;j<plusCovFiles[i].length;j++){
				MatchedFilterTrainier trainer = new MatchedFilterTrainier(plusCovFiles[i][j], minusCovFiles[i][j], refFlat, paramFiles[i][j]);
				System.out.println("Training for " + plusCovFiles[i][j] + " and " + minusCovFiles[i][j]);
				if(j == 0) trainer.train(30, 50, 3); // harr
				else trainer.train(30, 200, 7); // rpf
				System.out.println("Training done..");
				if(j==0){
					System.out.println("Scoring for " + plusCovFiles[i][j] + " and " + minusCovFiles[i][j]);
					Scorer scorer = new Scorer(plusCovFiles[i][j], minusCovFiles[i][j], paramFiles[i][j]);
					scorer.setAnnotationFileFile(refFlat);
					scorer.scoreNWrite(0, fasta, scoreOutFiles[i]);
					scorer.writeWindowFilteredOutput(scoreOutFiles[i], scoreOutFiles[i] + ".windowed.tsv", 50);
					System.out.println("Scoring done..");
				}
			}
		}
		System.out.println("Merging results");
		MergeResults merge = new MergeResults(scoreOutFiles, plusCovFiles, minusCovFiles, paramFiles);
		merge.merge(outFile, scoreThreshold);
		System.out.println("Merging done..");
	}
	
	private static void printUsage(){
		System.out.println("Usage:\n"
				+ "java -jar -Xmx5g RPFPipeLine.jar "
				+ "\n-harrP [harr coverage files (plus strand), seperated by space] "
				+ "\n-harrM [harr coverage files (minus strand), seperated by space]"
				+ "\n-rpfP [rpf coverage files (plus strand), seperated by space]"
				+ "\n-rpfM [rpf coverage files(minus strand), seperated by space]"
				+ "\n-fasta [fasta file]"
				+ "\n-ref [refFlat file]"
				+ "\n-scoreThreshold [scoreThreshold - default 2]"
				+ "\n-outputFile [output file]");
	}
}
