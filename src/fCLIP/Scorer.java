package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import fCLIP.MirGff3FileParser.MiRNA;
import fCLIP.parser.ScoringOutputParser;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.MafParser;
import parser.ZeroBasedFastaParser;
import util.GenerateMaffileWithSpecificSpecies;
public class Scorer {
	
	static private int leftWindowSize = 40;
	static private int rightWindowSize = 40;
	static private int maxDepthThreshold = 0;
	static public int flankingNTNumber = 13; // inclusive.. 
	private Bed12Parser bedParser;
	private double[] filter5p;
	private double[] filter3p;
	private double filter5pNorm;
	private double filter3pNorm;
	private ZeroBasedFastaParser fastaParser;
	private MirGff3FileParser mirParser;
	static private int minReadDiff = 50;
	static private int maxReadDiff = 100;
	
	
	private void read(String filename){
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(filename);
			String s;
			int mode = 0;
			int i = 0;
			while((s=in.readLine())!=null){
				String[] token = null;
				if(s.startsWith("#")) token = s.split("\t");
				
				if(s.startsWith("#LEFT")){
					leftWindowSize = Integer.parseInt(token[1]);
				}else if(s.startsWith("#RIGHT")){
					rightWindowSize = Integer.parseInt(token[1]);
				}else if(s.startsWith("#MAXDEPTHTHRESHOLD")){
					maxDepthThreshold = Integer.parseInt(token[1]);
				}else if(s.startsWith("#FILTER5P")){
					filter5p = new double[Integer.parseInt(token[1])];
					mode = 1;
					i = 0;
				}else if(s.startsWith("#FILTER3P")){
					filter3p = new double[Integer.parseInt(token[1])];
					mode = 2;
					i = 0;
				}else if(s.startsWith("#SIGNAL5P")){
					mode = 3;
					i = 0;
				}else if(s.startsWith("#NOISE5P")){
					mode = 4;
					i = 0;
				}else if(s.startsWith("#SIGNAL3P")){
					mode = 5;
					i = 0;
				}else if(s.startsWith("#NOISE3P")){
					mode = 6;
					i = 0;
				}else{
					if(mode == 1){
						filter5p[i++] = Double.parseDouble(s);
					}else if(mode == 2){
						filter3p[i++] = Double.parseDouble(s);
					}
				}
			}
			
			filter3pNorm = rpf.Scorer.getNorm(filter3p);
			filter5pNorm = rpf.Scorer.getNorm(filter5p);
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static double getRawScore(double[] filter, double[] cov, double filterNorm){
		if(rpf.Scorer.sum(cov) < maxDepthThreshold) return 0;
		double norm = rpf.Scorer.getNorm(cov);
		//System.out.println(norm);
		if(norm <= 0) return 0;
				
		double ip = rpf.Scorer.getInnerProduct(filter, cov);
		
	//	System.out.println("hh " + ip + " " + filterNorm + " " + norm);
		
		return ip / filterNorm / norm;
	}
	
	
	private HashMap<Integer, ArrayList<ScoredPosition>> getScoredPositionMap(AnnotationFileParser parser, boolean isPlusStrand, boolean is5p, double scoreThreshold){
		Iterator<Integer> iterator = is5p ? bedParser.getNonZero5pPositionIterator(isPlusStrand) : bedParser.getNonZero3pPositionIterator(isPlusStrand);
		HashMap<Integer, ArrayList<ScoredPosition>> scoredPositionMap = new HashMap<Integer, ArrayList<ScoredPosition>>();
		ArrayList<ScoredPosition> scoredPositions = new ArrayList<ScoringOutputParser.ScoredPosition>();
		HashSet<Integer> considered = new HashSet<Integer>();
		
		while(iterator.hasNext()){
			int position = iterator.next();
			//for(int offset=0;offset<1;offset++){
				int currentPosition = isPlusStrand == is5p? position - 1 : position + 1;
				if(considered.contains(currentPosition)) continue;
				
				considered.add(currentPosition);
				
				for(ArrayList<Integer> co : bedParser.getCoordinates(currentPosition + (isPlusStrand == is5p ? 1 :  -1), rightWindowSize, isPlusStrand, is5p)){ // should be left inclusive..
					if(!Bed12Parser.getSplices(co).isEmpty()) continue;
					ArrayList<Integer> coordinate = new ArrayList<Integer>(co);					
					
					for(int i=0;i<leftWindowSize;i++){
						if(is5p) coordinate.add(0, currentPosition + (isPlusStrand? -i : i));
						else coordinate.add(currentPosition + (isPlusStrand? i : -i));
					}
					
					double[] cov = is5p? bedParser.get5pCoverages(isPlusStrand, coordinate) : bedParser.get3pCoverages(isPlusStrand, coordinate);
					double[] cov2 = is5p? bedParser.get3pCoverages(isPlusStrand, coordinate) : bedParser.get5pCoverages(isPlusStrand, coordinate);
					
					double score = 0;
					if(rpf.Scorer.max(cov) < maxDepthThreshold && rpf.Scorer.max(cov2) < maxDepthThreshold) score = 0;
					else{
						double[] depth = is5p ? bedParser.get5pSignalForfCLIP(isPlusStrand, coordinate) :
							bedParser.get3pSignalForfCLIP(isPlusStrand, coordinate);
						score = getRawScore(is5p? filter5p : filter3p, depth, is5p? filter5pNorm : filter3pNorm);
					}
				 
					if(score > scoreThreshold){
						ScoredPosition scoredPosition = is5p? new ScoredPosition(bedParser.getContig(), isPlusStrand, currentPosition, coordinate, 0, null, parser).setFivePScore(score):
							new ScoredPosition(bedParser.getContig(), isPlusStrand, 0, null, currentPosition, coordinate, parser).setThreePScore(score);
						scoredPositions.add(scoredPosition);
					//	if(!scoredPositionMap.containsKey(currentPosition)) scoredPositionMap.put(currentPosition, new ArrayList<ScoredPosition>());
					//	scoredPositionMap.get(currentPosition).add(scoredPosition);
						
						/*if(!is5p){
							System.out.println(isPlusStrand + " " + is5p + " " + currentPosition);
							System.out.println(coordinate);
							for(double c : cov) System.out.print(c + " ");
							System.out.println();
						}*/
					}					
				}
			}			
		//}
		Collections.sort(scoredPositions);
		scoredPositions = windowFilter(scoredPositions, leftWindowSize, 1, is5p);
		
		for(ScoredPosition s: scoredPositions){
			int position = is5p? s.getFivePposition() : s.getThreePposition();
			if(!scoredPositionMap.containsKey(position)) scoredPositionMap.put(position, new ArrayList<ScoredPosition>());
			scoredPositionMap.get(position).add(s);
		}
		
		return scoredPositionMap;
	}
	
	public ArrayList<ScoredPosition> windowFilter(ArrayList<ScoredPosition> s, int window, int top, boolean is5p) {
	    
	    // select each peak if it is top n within window (-window,+window) around it
		ArrayList<ScoredPosition> retPositions = new ArrayList<ScoringOutputParser.ScoredPosition>();
	   
	    for(int positionIndex = 0; positionIndex < s.size(); positionIndex++) {
	      int rank = 1;
	      
	      ScoredPosition thisPeak = s.get(positionIndex);
	      int thisPosition = is5p? thisPeak.getFivePposition() : thisPeak.getThreePposition() ;
	      double thisScore = is5p? thisPeak.getFivePScore() : thisPeak.getThreePScore();
	      
	      // move left
	      int prevIndex = positionIndex-1;
	      while(prevIndex >= 0) {
	    	ScoredPosition prevPeak = s.get(prevIndex);
	    	int prevPosition = is5p? prevPeak.getFivePposition() : prevPeak.getThreePposition();
	    	double prevScore = is5p? prevPeak.getFivePScore() : prevPeak.getThreePScore();
	        if(thisPosition - prevPosition > window)    break;
	        if(prevScore > thisScore)            rank++;
	        prevIndex--;
	      }

	      // move right
	      int nextIndex = positionIndex+1;
	      while(nextIndex < s.size()) {
	    	ScoredPosition nextPeak = s.get(nextIndex);
	    	int nextPosition = is5p? nextPeak.getFivePposition() : nextPeak.getThreePposition();
	    	double nextScore = is5p? nextPeak.getFivePScore() : nextPeak.getThreePScore();
	    	 if(nextPosition - thisPosition > window)    break;
		     if(nextScore > thisScore)            rank++;   
		     nextIndex++;
	      }
	    
	      if(rank <= top)   retPositions.add(thisPeak);
	    }
	    return retPositions;
	  }
	
	private void updatePaired(HashSet<ScoredPosition> positionSet, int p3p, int p5p, HashMap<Integer, ArrayList<ScoredPosition>> sp3p, HashMap<Integer, ArrayList<ScoredPosition>> sp5p,
			double energyThreshold, double depthThreshold,  boolean isPlusStrand){
		String seq = isPlusStrand? fastaParser.getSequence(bedParser.getContig(), p3p - flankingNTNumber + 1, p5p + flankingNTNumber)
		    : ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(bedParser.getContig(), p5p-flankingNTNumber + 1, p3p + flankingNTNumber), true);
		
		RNAfoldLauncher fold = new RNAfoldLauncher(seq, flankingNTNumber);
		if(fold.getEnergy() > energyThreshold) return;
		if(fold.getDepth() < depthThreshold) return;
		if(fold.getOverHang() < 0 || fold.getOverHang() > 5) return;
		for(ScoredPosition s3 : sp3p.get(p3p)){
			for(ScoredPosition s5 : sp5p.get(p5p)){
				ScoredPosition sp = new ScoredPosition(bedParser.getContig(), isPlusStrand, p5p, s5.getFivePcoordinate(), p3p, s3.getThreePcoordinate(), null);
				sp.set(fold).setFivePScore(s5.getFivePScore()).setThreePScore(s3.getThreePScore()).setSeq(seq);
				sp.addGenes(s5);sp.addGenes(s3);
				positionSet.add(sp);
			}
		}
	}
	
	private void updateUnPaired(HashSet<ScoredPosition> positionSet, int p, HashMap<Integer, ArrayList<ScoredPosition>> sp, Classifier classifier, double unpairedScoreThreshold, double energyThreshold, double depthThreshold, boolean isPlusStrand, boolean is3p){
		boolean isLeft = isPlusStrand == is3p;
		ArrayList<ScoredPosition> sps = new ArrayList<ScoredPosition>();
		for(ScoredPosition s: sp.get(p)){
			if(is3p && s.getThreePScore() > unpairedScoreThreshold) sps.add(s);
			else if(!is3p && s.getFivePScore() > unpairedScoreThreshold) sps.add(s);
		}		
		if(sps.isEmpty()) return;
		
		String seq = isLeft? fastaParser.getSequence(bedParser.getContig(), p - flankingNTNumber + 1, p + maxReadDiff + flankingNTNumber)
		 :fastaParser.getSequence(bedParser.getContig(), p - maxReadDiff - flankingNTNumber + 1, p + flankingNTNumber);
		if(!isPlusStrand) seq = ZeroBasedFastaParser.getComplementarySequence(seq, true);
		
	//	System.out.println(is3p + " " + isPlusStrand + " " + (2 * flankingNTNumber + maxReadDiff - 1) + "  " + seq.length()); 
		if(seq.length() < 2 * flankingNTNumber + maxReadDiff - 1) return;
		double minEnergy = 1000;
		String maxSeq = null;
		RNAfoldLauncher maxFold = null;
		for(int k=2*flankingNTNumber+minReadDiff; k<2*flankingNTNumber+maxReadDiff; k++){
			String subSeq = is3p? seq.substring(0, k) : seq.substring(seq.length() - k);
			RNAfoldLauncher fold = new RNAfoldLauncher(subSeq, flankingNTNumber);
			if(fold.getOverHang() > 5 || fold.getOverHang() < 2) continue; //TODO
			if(fold.getEnergy() > energyThreshold) continue;
			if(fold.getDepth() < depthThreshold) continue;
			if(minEnergy > fold.getEnergy()){
				minEnergy = fold.getEnergy();
				maxFold = fold;
				maxSeq = subSeq;
			}
		}
		if(maxFold != null){
			double prediction = classifier.classify(maxFold.toInstance(classifier.getDataset()));
			if(prediction>0){
				int len = maxSeq.length() - 2*flankingNTNumber;				
				for(ScoredPosition s : sps){
					s.set(maxFold).setSeq(maxSeq);
					if(is3p)
						s.setFivePposition(s.getThreePposition() + (isPlusStrand? len + 1: - len - 1));
					else
						s.setThreePposition(s.getFivePposition() + (isPlusStrand? -len - 1: len + 1));
					positionSet.add(s);
				}
			}else{
				RNAfoldLauncher fold = new RNAfoldLauncher(seq, flankingNTNumber);
				for(ScoredPosition s : sps){
					s.set(fold).setSeq(seq);
					positionSet.add(s);
				}
			}
		}else{
			RNAfoldLauncher fold = new RNAfoldLauncher(seq, flankingNTNumber);
			for(ScoredPosition s : sps){
				s.set(fold).setSeq(seq);
				positionSet.add(s);
			}
		}
	}
	
	private void updateScoredPositions(HashSet<ScoredPosition> positionSet, HashMap<Integer, ArrayList<ScoredPosition>> sp3p, HashMap<Integer, ArrayList<ScoredPosition>> sp5p, 
			Classifier classifier, double unpairedScoreThreshold, double energyThreshold, double depthThreshold, boolean isPlusStrand){
		ArrayList<Integer> positions5p = new ArrayList<Integer>(sp5p.keySet());
		ArrayList<Integer> positions3p = new ArrayList<Integer>(sp3p.keySet());
		Collections.sort(positions5p);
		Collections.sort(positions3p);
		
		HashSet<ScoredPosition> tPositionSet = new HashSet<ScoringOutputParser.ScoredPosition>();
		HashSet<Integer> paired3ps = new HashSet<Integer>();
		HashSet<Integer> paired5ps = new HashSet<Integer>();
		
		for(int p3p : positions3p){
			int i = Collections.binarySearch(positions5p, p3p + (isPlusStrand? minReadDiff : -maxReadDiff));
			i = i<0? -i-1 : i;
			boolean pairExists = false;
			while(i < positions5p.size()){
				int p5p = positions5p.get(i);
				if(isPlusStrand?  p5p > p3p + maxReadDiff : p3p < p5p + minReadDiff) break;
				pairExists = true;
				updatePaired(tPositionSet, p3p, p5p, sp3p, sp5p, energyThreshold, depthThreshold, isPlusStrand);
				i++;
			}
			if(classifier!= null && !pairExists) updateUnPaired(tPositionSet, p3p, sp3p, classifier, unpairedScoreThreshold, energyThreshold, depthThreshold, isPlusStrand, true);
		}
		for(int p5p : positions5p){
			int i = Collections.binarySearch(positions3p, p5p + (isPlusStrand? - maxReadDiff : minReadDiff));
			i = i<0? -i-1 : i;
			boolean pairExists = false;
			while(i < positions3p.size()){
				int p3p = positions3p.get(i);
				if(isPlusStrand? p5p < p3p + minReadDiff : p3p > p5p + maxReadDiff) break;
				pairExists = true;
				updatePaired(tPositionSet, p3p, p5p, sp3p, sp5p, energyThreshold, depthThreshold, isPlusStrand);
				i++;							
			}
			if(classifier!= null && !pairExists) updateUnPaired(tPositionSet, p5p, sp5p, classifier, unpairedScoreThreshold, energyThreshold, depthThreshold, isPlusStrand, false);
		}
		
		for(ScoredPosition p : tPositionSet){
			if(!p.isPaired()) continue;
			paired3ps.add(p.getThreePposition());// 91
			paired5ps.add(p.getFivePposition());
			positionSet.add(p);
		}
		
		for(ScoredPosition p : tPositionSet){
			if(p.isPaired()) continue;
			int off = p.isPlusStrand()? 1 : -1;
			if(p.is3pScored()){
				if(paired5ps.contains(p.getThreePposition() - off )) continue;
			}else{
				if(paired3ps.contains(p.getFivePposition() + off )) continue;
			}
			positionSet.add(p);
		}
		
		if(classifier != null){
			for(ScoredPosition p : positionSet){
				classifier.setClassification(p);
			}
		}
		
	}
	
	private ArrayList<ScoredPosition> getScoredPositions(AnnotationFileParser parser, Classifier classifier, double unpairedScoreThreshold, double pairedScoreThreshold, double depthThreshold, double energyThreshold, boolean isPlusStrand){
		HashSet<ScoredPosition> positionSet = new HashSet<ScoredPosition>();
		System.out.println("Scoring for " + bedParser.getContig() + " " + (isPlusStrand? "+" : "-") + " strand");
		
		HashMap<Integer, ArrayList<ScoredPosition>> sp5p = getScoredPositionMap(parser, isPlusStrand, true, pairedScoreThreshold);
		HashMap<Integer, ArrayList<ScoredPosition>> sp3p = getScoredPositionMap(parser, isPlusStrand, false, pairedScoreThreshold);
		
		updateScoredPositions(positionSet, sp3p, sp5p, classifier, unpairedScoreThreshold, energyThreshold, depthThreshold, isPlusStrand);
		
		System.out.println("\t"+positionSet.size()+" points to be scored");
		for(ScoredPosition sp : positionSet){
			ArrayList<MiRNA> miRNAs = null;
			if(sp.getFivePcoordinate() != null){
				miRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
				ArrayList<MiRNA> matched = mirParser.getMatchingMiRNAs(bedParser.getContig(), isPlusStrand, sp.getFivePposition(), sp.getFivePcoordinate());
				if(matched != null) miRNAs.addAll(matched);
			}
			if(sp.getThreePcoordinate() != null){
				if(miRNAs == null){
					miRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
					ArrayList<MiRNA> matched = mirParser.getMatchingMiRNAs(bedParser.getContig(), isPlusStrand, sp.getThreePposition(), sp.getThreePcoordinate());
					if(matched != null) miRNAs.addAll(matched);
				}else{
					ArrayList<MiRNA> intersectedMiRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
					ArrayList<MiRNA> matched = mirParser.getMatchingMiRNAs(bedParser.getContig(), isPlusStrand, sp.getThreePposition(), sp.getThreePcoordinate());
					if(matched != null){
						for(MiRNA miRNA : matched){
							if(miRNAs.contains(miRNA)){
								intersectedMiRNAs.add(miRNA);
							}
						}
					}
					miRNAs = intersectedMiRNAs;
				}
			}
			if(miRNAs != null && !miRNAs.isEmpty()){
				sp.setMiRNAs(miRNAs);
			}
		}
		ArrayList<ScoredPosition> positions = new ArrayList<ScoredPosition>(positionSet);
		Collections.sort(positions);
		return positions;
	}
	
	public Scorer(Bed12Parser bedParser, ZeroBasedFastaParser fastaParser, MirGff3FileParser mirParser, String parameterFileName){
		this.bedParser = bedParser;
		this.fastaParser = fastaParser;
		this.mirParser = mirParser;
		read(parameterFileName);
	}
	
	static public void main(String[] args) throws IOException{
		GenerateMaffileWithSpecificSpecies.main(null);
		
		String sample = "x4";
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		MirGff3FileParser mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/genomes/hsa.gff3");
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg19.refFlat.txt");
		String bedFileName = "/media/kyowon/Data1/fCLIP/samples/sample2/bed/" + sample + ".sorted.bed";
		String parameterFileName = "/media/kyowon/Data1/fCLIP/samples/sample2/bed/" + sample + ".sorted.param";
		String outFileName = "/media/kyowon/Data1/Dropbox/" + sample + ".sorted.out.csv";
		String arffTrainOutFileName = "/media/kyowon/Data1/Dropbox/" + sample + ".sorted.out.train.arff";
		String arffOutFileName = "/media/kyowon/Data1/Dropbox/" + sample + ".sorted.out.arff";
		MafParser mafParser = new MafParser("/media/kyowon/Data1/fCLIP/genomes/maf/");
		double unpairedScoreThreshold = 0.095;
		double pairedScoreThreshold = 0.095;
		double energyThreshold = 1000;
		double depthThreshold = -1;
		mafParser.generateIndexFile();
		mafParser.readIndexFile();
		
		String[] args1 = new String[1];
		args1[0] = sample;
		ScorerTrainer.main(args1);
				
		if(!new File(arffTrainOutFileName).exists()){
			System.out.println("Training Classifier for " + sample);
			HashSet<ScoredPosition> trainingPositions = new HashSet<ScoredPosition>();
			int annotated = 0;
			for(String contig : fastaParser.getContigs()){
				if(contig.length() > 5 || contig.equals("chrM")) continue;
				Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
				Scorer scorer = new Scorer(bedParser, fastaParser, mirParser, parameterFileName);
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, null, unpairedScoreThreshold, 0.075, depthThreshold, energyThreshold, true)){
					trainingPositions.add(sp);
				}
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, null, unpairedScoreThreshold, 0.075, depthThreshold, energyThreshold, false)){
					trainingPositions.add(sp);
				}		
			}	
			PrintStream outTrainArff = new PrintStream(arffTrainOutFileName);
			outTrainArff.println(ScoringOutputParser.ScoredPosition.getArffHeader());
			for(ScoredPosition sp : trainingPositions){
				if(sp.getMiRNAs() != null) 	annotated++;
			}
			
			System.out.println("Training size : " + trainingPositions.size() + " Annotated : " + annotated);
			int numUnannotated = 100;
			for(ScoredPosition sp : trainingPositions){
				if(sp.getMiRNAs() == null) numUnannotated--;
				if(numUnannotated < 0 && sp.getMiRNAs() == null) continue;
				
				sp.setRNAzScores(mafParser);
				outTrainArff.println(sp.toArffString());
			}
			outTrainArff.close(); // training done			
		}		
		
		if(!new File(outFileName).exists()){
			Classifier classifier = new Classifier(arffTrainOutFileName);
			PrintStream outArff = new PrintStream(arffOutFileName);
			
			PrintStream out = new PrintStream(outFileName);
			
			outArff.println(ScoringOutputParser.ScoredPosition.getArffHeader());
			out.println(ScoredPosition.getHeader());
			
			for(String contig : fastaParser.getContigs()){
				if(contig.length() > 5 || contig.equals("chrM")) continue;
				
				Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
				Scorer scorer = new Scorer(bedParser, fastaParser, mirParser, parameterFileName);
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, classifier, unpairedScoreThreshold, pairedScoreThreshold, depthThreshold, energyThreshold, true)){
					if(sp.getClassification().equals("M")){ 
						sp.setRNAzScores(mafParser);
						sp.setDnDs(mafParser);
						classifier.setClassification(sp);
					}
					outArff.println(sp.toArffString());
					out.println(sp);
				}
				for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, classifier, unpairedScoreThreshold, pairedScoreThreshold, depthThreshold, energyThreshold, false)){
					if(sp.getClassification().equals("M")){ 
						sp.setRNAzScores(mafParser);
						sp.setDnDs(mafParser);
						classifier.setClassification(sp);
					}
					outArff.println(sp.toArffString());
					out.println(sp);
				}
			}			
			outArff.close(); // training done
			out.close();
		}
		ScoringOutputParser sparser = new ScoringOutputParser(outFileName);
		ScorePairs.run(sample, sparser, mafParser, arffTrainOutFileName,  false);
		ScorePairs.run(sample, sparser, mafParser, arffTrainOutFileName, true);
	}	
}
