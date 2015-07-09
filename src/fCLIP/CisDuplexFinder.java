package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import launcher.RNAcofoldLauncher;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.MirGff3FileParser;
import parser.ZeroBasedFastaParser;
import fCLIP.parser.BlatParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.BlatParser.BlatResult;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class CisDuplexFinder {
	
	public static class ScorePositionRunner implements Runnable { // TODO move this to FCLIP_Scorer.java
		private Bed12Parser bedParser;
		private FCLIP_Scorer scorer;
		private boolean isPlusStrand;
		//private ArrayList<ScoredPosition> positions;
		private Classifier classifier;
		double unpairedScoreThreshold, pairedScoreThreshold;
		private AnnotationFileParser annotationParser;
		private String outFileName;
		//= new FCLIP_Scorer(bedParser, fastaParser, mirParser, parameterFileName);
		
		ScorePositionRunner(String outFileName, String contig, String bedFileName, AnnotationFileParser annotationParser, ZeroBasedFastaParser fastaParser
				,MirGff3FileParser mirParser, String parameterFileName, Classifier classifier, double unpairedScoreThreshold, 
				double pairedScoreThreshold, boolean isPlusStrand){
			this.outFileName = outFileName;
			this.bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
			this.scorer = new FCLIP_Scorer(bedParser, fastaParser, mirParser, parameterFileName);
			this.isPlusStrand = isPlusStrand;
			this.classifier = classifier;
			this.unpairedScoreThreshold = unpairedScoreThreshold;
			this.pairedScoreThreshold = pairedScoreThreshold;
			this.annotationParser = annotationParser;
		}
		
		public void run() {
			PrintStream out;
			try {
				out = new PrintStream(outFileName);
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, classifier, unpairedScoreThreshold, pairedScoreThreshold, isPlusStrand)){					
					out.println(sp);					
				}
				out.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}			
		}
		
		//public ArrayList<ScoredPosition> getPositions() {return positions; }
		
	}
	
	public static void main(String[] args) {
		if(args.length == 13){
			String outFileName = args[0];
			String arffTrainOutFileName = args[1];
			String bedFileName = args[2];
			ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser(args[3]);
			AnnotationFileParser annotationParser = new AnnotationFileParser(args[4]);
			MirGff3FileParser mirParser = new MirGff3FileParser(args[5]);
			String parameterFileName = args[6];
			double unpairedScoreThreshold = Double.parseDouble(args[7]);
			double pairedScoreThreshold = Double.parseDouble(args[8]);
			int maxThreads = Integer.parseInt(args[9]);
			int flankingLength = Integer.parseInt(args[10]);
			int minpre = Integer.parseInt(args[11]);
			int maxpre = Integer.parseInt(args[12]);
			RNAcofoldLauncher.setSeqLength(minpre);
			
			FCLIP_Scorer.setFlankingNTNumber(flankingLength);
			FCLIP_Scorer.setMaxReadDiff(maxpre);
			FCLIP_Scorer.setMinReadDiff(minpre);
			//String pslOutFileName = outFileName + ".psl";
			//
			Classifier classifier = new Classifier(arffTrainOutFileName);
			try {
				ArrayList<String> tmps = new ArrayList<String>();
				ArrayList<Thread> threads = new ArrayList<Thread>();
				
				for(String contig : fastaParser.getContigs()){
					//if(!contig.equals("chr17")) continue; // TODO remove
					if(contig.length() > 5 || contig.equals("chrM")) continue;
					for(int i=0;i<2;i++){
						String t = outFileName + contig + i;
						ScorePositionRunner runner = new ScorePositionRunner(t, contig, bedFileName, 
								annotationParser, fastaParser, mirParser, parameterFileName, 
								classifier, unpairedScoreThreshold,
								pairedScoreThreshold, i == 0);
						
						Thread thread = new Thread(runner, contig + " " + (i == 0)); //Thread created       
						thread.start();
						tmps.add(t);
						threads.add(thread);
						if(threads.size() >= maxThreads){
							for(Thread rt : threads){
								try {
									rt.join();
								} catch (InterruptedException e) {						
									e.printStackTrace();
								}
							}
							threads.clear();
						}
						
					}				
				}
				
				for(Thread rt : threads){
					try {
						rt.join();
					} catch (InterruptedException e) {						
						e.printStackTrace();
					}
				}
				threads.clear();
				
				
				PrintStream out = new PrintStream(outFileName);
				out.println(ScoredPosition.getHeader());
				
				for(String t : tmps){
					BufferedLineReader in = new BufferedLineReader(t);
					String s;
					while((s=in.readLine())!=null){
						out.println(s);
					}
					in.close();
					new File(t).delete();
				}
				
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}			
			//new File(fastaOutFileName).delete();
		}
		if(args.length == 3){
			String outCsv = args[0];
			String csv = args[1];
			String psl = args[2];
			//BlatLauncher.run(fastaParser.getFastaFile(), fastaOutFileName, pslOutFileName, 100);
			ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
			BlatParser blatParser = new BlatParser(psl);
			try {
				PrintStream out = new PrintStream(outCsv);
				out.println(parser.getHeader());
				for(int i=0;i<parser.getPositions().size();i++){
					ScoredPosition position = parser.getPositions().get(i);
					ArrayList<BlatResult> r = blatParser.getResults(Integer.toString(i));
					if(r == null) position.setBlatHits(1);
					else position.setBlatHits(r.size());
					out.println(position);
				}
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		if(args.length == 2){
			String fastaOutFileName = args[0];
			String csv = args[1];
			int fastaIndex = 0;
			ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
			PrintStream outFasta;
			try {
				outFasta = new PrintStream(fastaOutFileName);
				for(ScoredPosition position : parser.getPositions()){
					outFasta.println(">"+ (fastaIndex++));
					outFasta.println(position.getSeq());
				}				
				outFasta.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			
		}
	}

}
