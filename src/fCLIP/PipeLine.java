package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import launcher.BlatLauncher;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.MirGff3FileParser;
import parser.ZeroBasedFastaParser;
import fCLIP.analysis.CheckRepeat;
import fCLIP.parser.BlatParser;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.BlatParser.BlatResult;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class PipeLine {

	
	static private void train(String outFileName, String arffTrainOutFileName, String bedFileName, ZeroBasedFastaParser fastaParser, AnnotationFileParser annotationParser, MirGff3FileParser mirParser, String parameterFileName
			, double unpairedScoreThreshold, double pairedScoreThreshold){
		if(new File(arffTrainOutFileName).exists())
			return;
		
		System.out.println("Training Classifier for " + bedFileName);
		HashSet<ScoredPosition> trainingPositions = new HashSet<ScoredPosition>();
		int annotated = 0;
		for(String contig : fastaParser.getContigs()){
			if(contig.length() > 5 || contig.equals("chrM")) continue;
			Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
			FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, fastaParser, mirParser, parameterFileName);
			for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, null, unpairedScoreThreshold, pairedScoreThreshold, true)){
				trainingPositions.add(sp);
			}
			for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, null, unpairedScoreThreshold, pairedScoreThreshold, false)){
				trainingPositions.add(sp);
			}		
		}		
		try {
			PrintStream outTrainArff = new PrintStream(arffTrainOutFileName);
			PrintStream outTrain = new PrintStream(outFileName);
			outTrainArff.println(Classifier.getArffHeader());
			outTrain.println(ScoredPositionOutputParser.ScoredPosition.getHeader());
			for(ScoredPosition sp : trainingPositions){
				if(sp.getMiRNAs() != null) 	annotated++;
			}
			
			System.out.println("Training size : " + trainingPositions.size() + " Annotated : " + annotated);
			int numUnannotated = 300;
			for(ScoredPosition sp : trainingPositions){
				if(sp.getMiRNAs() == null) numUnannotated--;
				if(numUnannotated < 0 && sp.getMiRNAs() == null) continue;
			
				outTrainArff.println(Classifier.toTrainingArffString(sp));
				outTrain.println(sp);
			}
			outTrainArff.close(); // training done		
			outTrain.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	static private void runCis(String outFileName, String arffTrainOutFileName, String bedFileName, ZeroBasedFastaParser fastaParser, AnnotationFileParser annotationParser, MirGff3FileParser mirParser, String parameterFileName
			, double unpairedScoreThreshold, double pairedScoreThreshold){
		if(new File(outFileName).exists()) return;
	
		String fastaOutFileName = outFileName + ".fa";
		String pslOutFileName = outFileName + ".psl";
		
		Classifier classifier = new Classifier(arffTrainOutFileName);
		try {
			PrintStream outFasta = new PrintStream(fastaOutFileName);
			
			int fastaIndex = 0;
			ArrayList<ScoredPosition> toWritePositions = new ArrayList<ScoredPositionOutputParser.ScoredPosition>();
			for(String contig : fastaParser.getContigs()){
				if(contig.length() > 5 || contig.equals("chrM")) continue;
				Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
				FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, fastaParser, mirParser, parameterFileName);
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, classifier, unpairedScoreThreshold, pairedScoreThreshold, true)){					
					toWritePositions.add(sp);
					outFasta.println(">"+ (fastaIndex++));
					outFasta.println(sp.getSeq());
				}
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, classifier, unpairedScoreThreshold, pairedScoreThreshold, false)){					
					toWritePositions.add(sp);
					outFasta.println(">"+ (fastaIndex++));
					outFasta.println(sp.getSeq());
				}
			}
			
			outFasta.close();
			BlatLauncher.run(fastaParser.getFastaFile(), fastaOutFileName, pslOutFileName, 100);
			BlatParser blatParser = new BlatParser(pslOutFileName);
			
			for(int i=0;i<toWritePositions.size();i++){
				ArrayList<BlatResult> r = blatParser.getResults(Integer.toString(i));
				if(r == null) toWritePositions.get(i).setBlatHits(1);
				else toWritePositions.get(i).setBlatHits(r.size());
			}
		
			PrintStream out = new PrintStream(outFileName);
			out.println(ScoredPosition.getHeader());

			for(ScoredPosition sp : toWritePositions){
				out.println(sp);
			}			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		new File(fastaOutFileName).delete();
		new File(pslOutFileName).delete();
	}
	
	static private void runTrans(String outFileNameM, String outFileNameU, String cisOutFileName, String arffTrainFileName, int blatHitThreshold, int seqLength, boolean sameDirection){
		if(new File(outFileNameM).exists()) return;
		ScoredPositionOutputParser sparser = new ScoredPositionOutputParser(cisOutFileName);
		System.out.println("Scoring " + outFileNameM + " pairs");	
		Classifier classifier = new Classifier(arffTrainFileName);
		String tmpM = outFileNameM + ".tmp.csv";
	//	String tmpU = outFileNameU + ".tmp.csv";
		if(!new File(tmpM).exists() )
			FCLIP_ScorePairs.getScoredPairs(sparser, classifier, blatHitThreshold, seqLength, sameDirection, tmpM, -1, 7);
		FCLIP_ScorePairs.setMatchedNumPositions(tmpM, outFileNameM);
	//	FCLIP_ScorePairs.setMatchedNumPositions(tmpU, outFileNameU);
		
		new File(tmpM).delete();
		//new File(tmpU).delete();
	}
	
	static public void filterTransPairs(String pairCsv, String outCsv, int num3pPaired, int num5pPaired){
		ScoredPairOutputParser parser = new ScoredPairOutputParser(pairCsv);
		PrintStream out;
		int sum = 0;
		int o2 = 0;
		try {
			out = new PrintStream(outCsv);
			out.println(parser.getHeader());
			for(ScoredPair p : parser.getPairs()){
				if(p.getMatched3pNum() > num3pPaired && p.getMatched5pNum() > num5pPaired){
					out.println(p);
					sum ++;
					if(p.getOverHang() == 2) o2 ++;
				}
			}
			System.out.println(sum + " " + (double)o2/sum * 100);
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	// clip bed file
	// db fasta
	// parameter file
	
	static public void main(String[] args) throws IOException{
		String sample = "h19x24";
		String parentFolder = "/media/kyowon/Data1/Dropbox/fCLIP/new24/";		
		String bedFileName = "/media/kyowon/Data1/fCLIP/samples/sample3/bed/" + sample + ".sorted.bed";
		String dbFasta = "/media/kyowon/Data1/RPF_Project/genomes/hg19.fa";
		String parameterFileName = "/media/kyowon/Data1/fCLIP/samples/sample3/bed/" + sample + ".sorted.param";
		
		String trainOutFileName = parentFolder + sample + ".train.csv";
		String arffTrainOutFileName = parentFolder + sample + ".train.arff";
		
		String cisOutFileName = parentFolder + sample + ".csv";
		String transOutFileNameM = parentFolder + sample + ".pair.M.csv";
		String transOutFileNameU = parentFolder + sample + ".pair.U.csv";
		String transControlOutFileNameM = parentFolder + sample + ".pair.M.AntiSense.csv";
		String transControlOutFileNameU = parentFolder + sample + ".pair.U.AntiSense.csv";
		String rmskBed = "/media/kyowon/Data1/fCLIP/Data/cat.rmsk.bed";
		String siControlBed = "/media/kyowon/Data1/fCLIP/RNAseq/siControl_R1_Aligned_Sorted.bed";
		String siKDDroshaBed = "/media/kyowon/Data1/fCLIP/RNAseq/siDrosha_R1_Aligned_Sorted.bed";
		String siKDDicerBed = "/media/kyowon/Data1/fCLIP/RNAseq/siDicer_R1_Aligned_Sorted.bed";
		
		String cisBed5pFileName = cisOutFileName + ".5p.bed";
		String cisBed3pFileName = cisOutFileName + ".3p.bed";
		
		double unpairedScoreThreshold = 0.25; 
		double pairedScoreThreshold = unpairedScoreThreshold/5;
	//	double motifScoreThreshold = 0.15;
		
		int num3pPaired = 10;
		int num5pPaired = 20;
		
		String filteredTransOutFileName = transOutFileNameM + ".filtered.csv";
		
		int blatHitThreshold = 100000000;
		int transPairSeqLength = FCLIP_Scorer.getFlankingNTNumber() *2 + 80;
		
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser(dbFasta);
		MirGff3FileParser mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/genomes/hsa_hg19.gff3");
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg19.refFlat.txt");
		FCLIP_ScorerTrainer.train(parameterFileName, bedFileName, mirParser, annotationParser);
		
		train(trainOutFileName, arffTrainOutFileName, bedFileName, fastaParser, annotationParser, mirParser, parameterFileName, unpairedScoreThreshold, pairedScoreThreshold);
		runCis(cisOutFileName, arffTrainOutFileName, bedFileName, fastaParser, annotationParser, mirParser, parameterFileName, unpairedScoreThreshold, pairedScoreThreshold);

		ScoredPositionOutputParser.generateBedFromCsv(cisOutFileName, cisBed5pFileName, cisBed3pFileName, false); // for fold change..
	//	ScoredPositionOutputParser.generateFastaForMotif(cisOutFileName, cisOutFileName + ".5p.motif.fa", cisOutFileName + ".3p.motif.fa",cisOutFileName + ".motif.m", "M");
		
		runTrans(transOutFileNameM, transOutFileNameU, cisOutFileName, arffTrainOutFileName, blatHitThreshold, transPairSeqLength, false);
		CheckRepeat.generate(cisOutFileName, transOutFileNameM, cisBed5pFileName, cisBed3pFileName, rmskBed);
		CheckRepeat.generate(cisOutFileName, transOutFileNameU, cisBed5pFileName, cisBed3pFileName, rmskBed);
		ScoredPairOutputParser.generateFastaForMotif(transOutFileNameM, transOutFileNameM + ".5p.motif.fa", transOutFileNameM + ".3p.motif.fa");
		//GenerateCircosLinkFiles.run(transOutFileNameM + ".rmsk.csv", transOutFileNameM + ".link.txt");
	//	CheckCoverage.generate(cisOutFileName, transOutFileNameM, cisBed5pFileName, cisBed3pFileName, siKDDroshaBed, siControlBed);
	///	CheckCoverage.generate(cisOutFileName, transOutFileNameM, cisBed5pFileName, cisBed3pFileName, siKDDicerBed, siControlBed);
		
		runTrans(transControlOutFileNameM, transControlOutFileNameU, cisOutFileName, arffTrainOutFileName, blatHitThreshold, transPairSeqLength, true);		
		CheckRepeat.generate(cisOutFileName, transControlOutFileNameM, cisBed5pFileName, cisBed3pFileName, rmskBed);
		CheckRepeat.generate(cisOutFileName, transControlOutFileNameU, cisBed5pFileName, cisBed3pFileName, rmskBed);
		ScoredPairOutputParser.generateFastaForMotif(transControlOutFileNameM, transControlOutFileNameM + ".5p.motif.fa", transControlOutFileNameM + ".3p.motif.fa");
	//	GenerateCircosLinkFiles.run(transControlOutFileNameM + ".rmsk.csv", transControlOutFileNameM + ".link.txt");
		//CheckCoverage.generate(cisOutFileName, transControlOutFileNameM, cisBed5pFileName, cisBed3pFileName, siKDDroshaBed, siControlBed);
	//	CheckCoverage.generate(cisOutFileName, transControlOutFileNameM, cisBed5pFileName, cisBed3pFileName, siKDDicerBed, siControlBed);
		
		filterTransPairs(transOutFileNameM, filteredTransOutFileName, num3pPaired, num5pPaired);
		CheckRepeat.generate(cisOutFileName, filteredTransOutFileName, cisBed5pFileName, cisBed3pFileName, rmskBed);
		ScoredPairOutputParser.generateFastaForMotif(filteredTransOutFileName, filteredTransOutFileName + ".5p.motif.fa", filteredTransOutFileName + ".3p.motif.fa");
		//GenerateCircosLinkFiles.run(filteredTransOutFileName + ".rmsk.csv", filteredTransOutFileName + ".link.txt");
	///	CheckCoverage.generate(cisOutFileName, filteredTransOutFileName, cisBed5pFileName, cisBed3pFileName, siKDDroshaBed, siControlBed);
	//	CheckCoverage.generate(cisOutFileName, filteredTransOutFileName, cisBed5pFileName, cisBed3pFileName, siKDDicerBed, siControlBed);
		
		
		// generate background for drosha and dicer..
		
		
		// motif different conditioin?? 
		
				
		//GenerateDepthsForEncodeDataSets.generate(outFileName, outFileName + ".encode.csv", annotationParser);
		
		
	}

}
