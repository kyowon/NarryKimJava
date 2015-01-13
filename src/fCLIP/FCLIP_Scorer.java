package fCLIP;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import launcher.RNAfoldLauncher;

import fCLIP.MirGff3FileParser.MiRNA;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.ZeroBasedFastaParser;
public class FCLIP_Scorer {
	
	static private int leftWindowSize = 40;
	static private int rightWindowSize = 40;
	static private int maxDepthThreshold = 0;
	static public int flankingNTNumber = 15; // inclusive.. 
	private Bed12Parser bedParser;
	private double[] filter5p;
	private double[] filter3p;
	private double filter5pNorm;
	private double filter3pNorm;
	private ZeroBasedFastaParser fastaParser;
	private MirGff3FileParser mirParser;
	final static private int minReadDiff = 40;
	final static private int maxReadDiff = 150;
	
	
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
	
	
	private double getScore(ArrayList<Integer> co, int currentPosition, boolean isPlusStrand, boolean is5p){
		ArrayList<Integer> coordinate = new ArrayList<Integer>(co);
		for(int i=0;i<leftWindowSize;i++){
			if(is5p) coordinate.add(0, currentPosition + (isPlusStrand? -i : i));
			else coordinate.add(currentPosition + (isPlusStrand? i : -i));
		}
		
		double[] depth = is5p ? bedParser.get5pSignalForfCLIP(isPlusStrand, coordinate) :
			bedParser.get3pSignalForfCLIP(isPlusStrand, coordinate);
		return getRawScore(is5p? filter5p : filter3p, depth, is5p? filter5pNorm : filter3pNorm);
		
	}
	
	
	private HashMap<Integer, ArrayList<ScoredPosition>> getScoredPositionMap(AnnotationFileParser parser, boolean isPlusStrand, boolean is5p, double scoreThreshold){
		HashMap<Integer, ArrayList<ScoredPosition>> scoredPositionMap = new HashMap<Integer, ArrayList<ScoredPosition>>();
		ArrayList<ScoredPosition> scoredPositions = new ArrayList<ScoredPositionOutputParser.ScoredPosition>();
		HashSet<Integer> positionsToConsider = new HashSet<Integer>();
		
		Iterator<Integer> iterator = is5p ? bedParser.getNonZero5pPositionIterator(isPlusStrand) : bedParser.getNonZero3pPositionIterator(isPlusStrand);
		while(iterator.hasNext()){
			int position = iterator.next() + (isPlusStrand == is5p? -1 : 1);
			positionsToConsider.add(position);
		}
		Iterator<Integer> iterator2 = is5p ? bedParser.getNonZero3pPositionIterator(isPlusStrand) : bedParser.getNonZero5pPositionIterator(isPlusStrand);
		while(iterator2.hasNext()){
			int position = iterator2.next();
			positionsToConsider.add(position);
		}
		ArrayList<Integer> positionListToConsider = new ArrayList<Integer>(positionsToConsider);
		Collections.sort(positionListToConsider);
		

		//windowFilter(positionListToConsider, 3, 1, isPlusStrand, is5p);
		
		for(int currentPosition : positionListToConsider){			
			double cov = is5p? bedParser.get5pDepth(isPlusStrand, currentPosition + (isPlusStrand == is5p ? 1 :  -1)) : 
				bedParser.get3pDepth(isPlusStrand, currentPosition + (isPlusStrand == is5p ? 1 :  -1));//     get5pCoverages(isPlusStrand, coordinate) : bedParser.get3pCoverages(isPlusStrand, coordinate);
			double cov2 = is5p? bedParser.get3pDepth(isPlusStrand, currentPosition) : bedParser.get5pDepth(isPlusStrand, currentPosition);
			if(cov + cov2 < maxDepthThreshold) continue;
			
			for(ArrayList<Integer> co : bedParser.getCoordinates(currentPosition + (isPlusStrand == is5p ? 1 :  -1), rightWindowSize, isPlusStrand, is5p)){ // should be left inclusive..
				if(!Bed12Parser.getSplices(co).isEmpty()) continue;
				double score = getScore(co, currentPosition, isPlusStrand, is5p);
				if(score > scoreThreshold){
					ScoredPosition scoredPosition = is5p? new ScoredPosition(bedParser.getContig(), isPlusStrand, currentPosition, -1, score, -1, parser):
						new ScoredPosition(bedParser.getContig(), isPlusStrand, -1, currentPosition, -1, score, parser);
					scoredPositions.add(scoredPosition);
				}					
			}
		}			
		
		Collections.sort(scoredPositions);
		scoredPositions = windowFilter(scoredPositions, 5, 1, is5p);
		
		for(ScoredPosition s: scoredPositions){
			int position = is5p? s.getFivePposition() : s.getThreePposition();
			if(!scoredPositionMap.containsKey(position)) scoredPositionMap.put(position, new ArrayList<ScoredPosition>());
			scoredPositionMap.get(position).add(s);
		}
		
		return scoredPositionMap;
	}
	
	private ArrayList<ScoredPosition> windowFilter(ArrayList<ScoredPosition> s, int window, int top, boolean is5p) {
	    
	    // select each peak if it is top n within window (-window,+window) around it
		ArrayList<ScoredPosition> retPositions = new ArrayList<ScoredPositionOutputParser.ScoredPosition>();
	   
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
	
	private ArrayList<Integer> windowFilter(ArrayList<Integer> s, int window, int top, boolean isPlusStrand, boolean is5p) {
	    
	    // select each peak if it is top n within window (-window,+window) around it
		ArrayList<Integer> retPositions = new ArrayList<Integer>();
	   
	    for(int positionIndex = 0; positionIndex < s.size(); positionIndex++) {
	      int rank = 1;
	      
	      int thisPosition = s.get(positionIndex);
	      double thisScore = is5p? bedParser.get5pDepth(isPlusStrand, thisPosition) : bedParser.get3pDepth(isPlusStrand, thisPosition);
	      // move left
	      int prevIndex = positionIndex-1;
	      while(prevIndex >= 0) {
	    	int prevPosition =  s.get(prevIndex);
	    	double prevScore = is5p? bedParser.get5pDepth(isPlusStrand, prevPosition) : bedParser.get3pDepth(isPlusStrand, prevPosition);
	        if(thisPosition - prevPosition > window)    break;
	        if(prevScore > thisScore)            rank++;
	        prevIndex--;
	      }

	      // move right
	      int nextIndex = positionIndex+1;
	      while(nextIndex < s.size()) {
	    	int nextPosition = s.get(nextIndex);
	    	double nextScore = is5p? bedParser.get5pDepth(isPlusStrand, nextPosition) : bedParser.get3pDepth(isPlusStrand, nextPosition);
	    	 if(nextPosition - thisPosition > window)    break;
		     if(nextScore > thisScore)            rank++;   
		     nextIndex++;
	      }
	    
	      if(rank <= top)   retPositions.add(thisPosition);
	    }
	    return retPositions;
	  }

//CGCCGCCGCCGC CGCCCCGG CACCCGCCTCCCGGCGCTGACGGTCTCGTACGAAGCCGGCGAGGGGGAG CCAGCAGC  GGCGGTCGCCGG
	private void updatePaired(HashSet<ScoredPosition> positionSet, int p3p, int p5p, HashMap<Integer, ArrayList<ScoredPosition>> sp3p, HashMap<Integer, ArrayList<ScoredPosition>> sp5p,
			double unpairedScoreThreshold, boolean isPlusStrand, AnnotationFileParser parser){
		String seq = isPlusStrand? fastaParser.getSequence(bedParser.getContig(), p3p - flankingNTNumber, p5p + flankingNTNumber + 1)
		    : ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(bedParser.getContig(), p5p-flankingNTNumber, p3p + flankingNTNumber + 1), true);
		
		RNAfoldLauncher fold = new RNAfoldLauncher(seq, flankingNTNumber);
	//	if(fold.getEnergyPerNT() > energyThreshold) return;
	//	if(fold.getDepth() < depthThreshold) return;
		//if(checkFoldFeasibility && !fold.isFeasibleFold()) return;
		for(ScoredPosition s3 : sp3p.get(p3p)){
			for(ScoredPosition s5 : sp5p.get(p5p)){
				if(s5.getFivePScore() < unpairedScoreThreshold && s3.getThreePScore() < unpairedScoreThreshold) continue;
				ScoredPosition sp = new ScoredPosition(bedParser.getContig(), isPlusStrand, p5p, p3p, s5.getFivePScore(), s3.getThreePScore(), parser);
				sp.set(fold).setSeq(seq);
				sp.addGenes(s5);
				sp.addGenes(s3);
				//sp.se
				sp.setGeneInfo(parser);
				positionSet.add(sp);
			}
		}
	}
	
	private void updateUnPaired(HashSet<ScoredPosition> positionSet, int p, HashMap<Integer, ArrayList<ScoredPosition>> sp, double unpairedScoreThreshold,
			boolean isPlusStrand, boolean is3p, AnnotationFileParser annotationParser){
		boolean isLeft = isPlusStrand == is3p;
		ArrayList<ScoredPosition> sps = new ArrayList<ScoredPosition>();
		for(ScoredPosition s: sp.get(p)){
			if(is3p && s.getThreePScore() > unpairedScoreThreshold) sps.add(s);
			else if(!is3p && s.getFivePScore() > unpairedScoreThreshold) sps.add(s);
		}		
		if(sps.isEmpty()) return;
		
		String seq = isLeft? fastaParser.getSequence(bedParser.getContig(), p - flankingNTNumber, p + maxReadDiff + flankingNTNumber + 1)
		 :fastaParser.getSequence(bedParser.getContig(), p - maxReadDiff - flankingNTNumber, p + flankingNTNumber + 1);
		if(!isPlusStrand) seq = ZeroBasedFastaParser.getComplementarySequence(seq, true);
		
	//	System.out.println(is3p + " " + isPlusStrand + " " + (2 * flankingNTNumber + maxReadDiff - 1) + "  " + seq.length()); 
		if(seq.length() < 2 * flankingNTNumber + maxReadDiff - 1) return;
		double minEnergy = 1000;
		String maxSeq = null;
		RNAfoldLauncher maxFold = null;
		//for(int k=2*flankingNTNumber+minReadDiff; k<2*flankingNTNumber+maxReadDiff; k++) // TODO
		{
			String subSeq = seq;//is3p? seq.substring(0, k) : seq.substring(2*flankingNTNumber+maxReadDiff - k - 1);
			RNAfoldLauncher fold = new RNAfoldLauncher(subSeq, flankingNTNumber);
			//if(!fold.isFeasibleFold()) continue;
			//if(fold.getEnergy() > energyThreshold) continue;
			//if(fold.getDepth() < depthThreshold) continue;
			if(minEnergy > fold.getEnergy()){
				minEnergy = fold.getEnergy();
				maxFold = fold;
				maxSeq = subSeq;
			}
		}
		if(maxFold != null){
			int len = maxSeq.length() - 2*flankingNTNumber;				
			for(ScoredPosition s : sps){
				s.set(maxFold).setSeq(maxSeq);
				if(is3p)
					s.setFivePposition(s.getThreePposition() + (isPlusStrand? len - 1: - len + 1), annotationParser);
				else
					s.setThreePposition(s.getFivePposition() + (isPlusStrand? - len + 1: len - 1), annotationParser);
				positionSet.add(s);
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
			Classifier classifier, double unpairedScoreThreshold, boolean isPlusStrand, AnnotationFileParser annotationParser){
		ArrayList<Integer> positions5p = new ArrayList<Integer>(sp5p.keySet());
		ArrayList<Integer> positions3p = new ArrayList<Integer>(sp3p.keySet());
		Collections.sort(positions5p);
		Collections.sort(positions3p);
		
		HashSet<ScoredPosition> tPositionSet = new HashSet<ScoredPositionOutputParser.ScoredPosition>();
		
		for(int p3p : positions3p){
			int i = Collections.binarySearch(positions5p, p3p + (isPlusStrand? minReadDiff : -maxReadDiff));
			i = i<0? -i-1 : i;
			ArrayList<Integer> paired5pList = new ArrayList<Integer>();
			
			while(i < positions5p.size()){
				int p5p = positions5p.get(i);
				if(isPlusStrand?  p5p > p3p + maxReadDiff : p3p < p5p + minReadDiff) break;
				paired5pList.add(p5p);
				i++;
			}
			
			if(!paired5pList.isEmpty())
				updatePaired(tPositionSet, p3p, isPlusStrand? paired5pList.get(0) : paired5pList.get(paired5pList.size()-1), sp3p, sp5p, unpairedScoreThreshold, isPlusStrand, annotationParser);
			else if(classifier!= null) updateUnPaired(tPositionSet, p3p, sp3p, unpairedScoreThreshold, isPlusStrand, true, annotationParser);
		}
				
		for(int p5p : positions5p){
			int i = Collections.binarySearch(positions3p, p5p + (isPlusStrand? - maxReadDiff : minReadDiff));
			i = i<0? -i-1 : i;
			ArrayList<Integer> paired3pList = new ArrayList<Integer>();
			
			while(i < positions3p.size()){
				int p3p = positions3p.get(i);
				if(isPlusStrand? p5p < p3p + minReadDiff : p3p > p5p + maxReadDiff) break;
				paired3pList.add(p3p);
				i++;							
			}
			
			if(!paired3pList.isEmpty())
				updatePaired(tPositionSet, isPlusStrand? paired3pList.get(paired3pList.size()-1) : paired3pList.get(0), p5p, sp3p, sp5p, unpairedScoreThreshold, isPlusStrand, annotationParser);
			else if(classifier!= null) updateUnPaired(tPositionSet, p5p, sp5p, unpairedScoreThreshold, isPlusStrand, false, annotationParser);
		}
		
		for(ScoredPosition p : tPositionSet){ // set read length variance
		
			int offset = isPlusStrand ? 1 : -1;
						
			int ft = bedParser.get3pDepth(isPlusStrand, p.getFivePposition());
			int ff = bedParser.get5pDepth(isPlusStrand, p.getFivePposition() + offset);
			int tf = bedParser.get5pDepth(isPlusStrand, p.getThreePposition());
			int tt = bedParser.get3pDepth(isPlusStrand, p.getThreePposition() - offset);
		
			p.set5pReads(ft, ff);
			p.set3pReads(tt, tf);
			
			double sum = 0;
			
			if(isPlusStrand){
				for(int i=p.getThreePposition() + 1;i<=p.getFivePposition() - 1;i++){
					sum += bedParser.get5pDepth(isPlusStrand, i);
					sum += bedParser.get3pDepth(isPlusStrand, i);
				}				
			}else{
				for(int i=p.getFivePposition() + 1;i<=p.getThreePposition() - 1;i++){
					sum += bedParser.get5pDepth(isPlusStrand, i);
					sum += bedParser.get3pDepth(isPlusStrand, i);
				}
			}
			
			p.setHetero(Math.log10(sum / (ft + ff + tf + tt) + .01));
			
		}
		
		HashMap<Integer, Double> fivePScores = new HashMap<Integer, Double>();
		HashMap<Integer, Double> threePScores = new HashMap<Integer, Double>();
		HashSet<Integer> paired3ps = new HashSet<Integer>();
		HashSet<Integer> paired5ps = new HashSet<Integer>();
		
		
		for(ScoredPosition p : tPositionSet){
			if(!p.isPaired()){
				if(p.is5pScored()) fivePScores.put(p.getFivePposition(), p.getFivePScore());
				else if(p.is3pScored()) threePScores.put(p.getThreePposition(), p.getThreePScore());
				continue;
			}
			paired3ps.add(p.getThreePposition());// 91
			paired5ps.add(p.getFivePposition());
			positionSet.add(p);
		}
		
		HashSet<Integer> mClassified3ps = new HashSet<Integer>();
		HashSet<Integer> mClassified5ps = new HashSet<Integer>();
		
		
		//String mafString = mafParser.getSeqsInMafFormat(contig, position.getCoordinate(), position.getPosition(), isPlusStrand, positionQuantityChangeLength);
		//this.phyloP = new PhyloPLauncher(mafString).getPvalConservation();
		
		
		if(classifier != null){
			for(ScoredPosition p : tPositionSet){
				//classifier.setClassification(p); // 1st classification
				//if(p.getClassification().toUpperCase().equals("M")){ // 2nd classification
			//		p.setRNAzScoresNSeedConservation(mafParser);
			//		p.setDnDs(mafParser);
					classifier.setClassification(p);
				//}
				
				if(p.isPaired()) continue;
				if(p.getClassification().equals("M")){
					if(p.is3pScored()) mClassified3ps.add(p.getThreePposition());
					else mClassified5ps.add(p.getFivePposition());
				}
			}
		}
		
		mClassified5ps.remove(0);
		mClassified3ps.remove(0);
		
		for(ScoredPosition p : tPositionSet){
			if(p.isPaired()) continue;
			int off = p.isPlusStrand()? 1 : -1;
			if(p.is3pScored()){
				if(paired5ps.contains(p.getThreePposition() - off )) continue;
				Double fs = fivePScores.get(p.getThreePposition() - off);
				if(fs!=null && fs > p.getThreePScore()) continue;
				if(!mClassified3ps.contains(p.getThreePposition()) && mClassified5ps.contains(p.getThreePposition() - off)) continue;
			}else{
				if(paired3ps.contains(p.getFivePposition() + off )) continue;
				Double ts = threePScores.get(p.getFivePposition() + off); 
				if(ts!=null && ts > p.getFivePScore()) continue;
				if(!mClassified5ps.contains(p.getFivePposition()) && mClassified3ps.contains(p.getFivePposition() + off)) continue;
			}
			positionSet.add(p);
		}
	}
	
	ArrayList<ScoredPosition> getScoredPositions(AnnotationFileParser parser, Classifier classifier, double unpairedScoreThreshold, double pairedScoreThreshold, boolean isPlusStrand){
		HashSet<ScoredPosition> positionSet = new HashSet<ScoredPosition>();
		System.out.println("Scoring for " + bedParser.getContig() + " " + (isPlusStrand? "+" : "-") + " strand");
		
		HashMap<Integer, ArrayList<ScoredPosition>> sp5p = getScoredPositionMap(parser, isPlusStrand, true, pairedScoreThreshold);
		HashMap<Integer, ArrayList<ScoredPosition>> sp3p = getScoredPositionMap(parser, isPlusStrand, false, pairedScoreThreshold);
		
		updateScoredPositions(positionSet, sp3p, sp5p, classifier, unpairedScoreThreshold, isPlusStrand, parser);
		
		System.out.println("\t"+positionSet.size()+" points to be scored");
		for(ScoredPosition sp : positionSet){
			ArrayList<MiRNA> miRNAs = null;
			if(sp.is5pScored()){
				miRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
				ArrayList<MiRNA> matched = mirParser.getContainingMiRNAs(bedParser.getContig(), isPlusStrand, sp.getFivePposition());
				if(matched != null) miRNAs.addAll(matched);
			}
			if(sp.is3pScored()){
				if(miRNAs == null){
					miRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
					ArrayList<MiRNA> matched = mirParser.getContainingMiRNAs(bedParser.getContig(), isPlusStrand, sp.getThreePposition());
					if(matched != null) miRNAs.addAll(matched);
				}else{
					ArrayList<MiRNA> intersectedMiRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
					ArrayList<MiRNA> matched = mirParser.getContainingMiRNAs(bedParser.getContig(), isPlusStrand, sp.getThreePposition());
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
	
	public FCLIP_Scorer(Bed12Parser bedParser, ZeroBasedFastaParser fastaParser, MirGff3FileParser mirParser, String parameterFileName){
		this.bedParser = bedParser;
		this.fastaParser = fastaParser;
		this.mirParser = mirParser;
		read(parameterFileName);
	}
	
	public static void main(String[] args){
		Bed12Parser bedParser = new Bed12Parser("/media/kyowon/Data1/fCLIP/samples/sample4/bed/merged.se.bed", 
				new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt"),
				"chr19", true);
		//ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/fCLIP/genomes/hg38.fa");
		//MirGff3FileParser mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3");
		String parameterFileName = "/media/kyowon/Data1/fCLIP/samples/sample4/results/training/merged.param";
		FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, null, null, parameterFileName);
		
		
		for(int o = - 5; o<5; o++){
			int currentPosition = 58349275 + o;
			boolean isPlusStrand = true;
			boolean is5p = true;
			ArrayList<ArrayList<Integer>> cos = bedParser.getCoordinates(currentPosition + (isPlusStrand == is5p ? 1 :  -1), FCLIP_Scorer.rightWindowSize, isPlusStrand, is5p);			
			if(cos == null) System.out.println("null");
			else System.out.println(currentPosition + " " + scorer.getScore(cos.get(0), currentPosition, isPlusStrand, is5p));
		}
		System.out.println();
		for(int o = - 5; o<5; o++){
			int currentPosition = 58349275 + o;
			boolean isPlusStrand = true;
			boolean is5p = false;
			ArrayList<ArrayList<Integer>> cos = bedParser.getCoordinates(currentPosition + (isPlusStrand == is5p ? 1 :  -1), FCLIP_Scorer.rightWindowSize, isPlusStrand, is5p);			
			if(cos == null) System.out.println("null");
			else System.out.println(currentPosition + " " + scorer.getScore(cos.get(0), currentPosition, isPlusStrand, is5p));
		}
		System.out.println();
		for(int o = - 5; o<5; o++){
			int currentPosition = 58349236 + o;
			boolean isPlusStrand = true;
			boolean is5p = true;
			ArrayList<ArrayList<Integer>> cos = bedParser.getCoordinates(currentPosition + (isPlusStrand == is5p ? 1 :  -1), FCLIP_Scorer.rightWindowSize, isPlusStrand, is5p);			
			if(cos == null) System.out.println("null");
			else System.out.println(currentPosition + " " + scorer.getScore(cos.get(0), currentPosition, isPlusStrand, is5p));
		}
		System.out.println();
		for(int o = - 5; o<5; o++){
			int currentPosition = 58349236 + o;
			boolean isPlusStrand = true;
			boolean is5p = false;
			ArrayList<ArrayList<Integer>> cos = bedParser.getCoordinates(currentPosition + (isPlusStrand == is5p ? 1 :  -1), FCLIP_Scorer.rightWindowSize, isPlusStrand, is5p);			
			if(cos == null) System.out.println("null");
			else System.out.println(currentPosition + " " + scorer.getScore(cos.get(0), currentPosition, isPlusStrand, is5p));
		}
		System.out.println();
		
	}
	
}
