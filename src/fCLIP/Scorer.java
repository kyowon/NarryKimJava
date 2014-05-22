package fCLIP;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import fCLIP.MirGff3FileParser.MiRNA;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.ZeroBasedFastaParser;
public class Scorer {
	
	static private int leftWindowSize = 40;
	static private int rightWindowSize = 40;
	static private int readCountThreshold = 0;
 
	private Bed12Parser bedParser;
	private double[] filter5p;
	private double[] filter3p;
	private double filter5pNorm;
	private double filter3pNorm;
	private ZeroBasedFastaParser fastaParser;
	private MirGff3FileParser mirParser;
	
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
				}else if(s.startsWith("#READCOUNTTHRESHOLD")){
					readCountThreshold = Integer.parseInt(token[1]);
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private static double getRawScore(double[] filter, double[] cov, double filterNorm){
		if(rpf.Scorer.sum(cov) < readCountThreshold) return 0;
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
		HashSet<Integer> considered = new HashSet<Integer>();
		
		while(iterator.hasNext()){
			int position = iterator.next();
			for(int offset=0;offset<30;offset++){
				int currentPosition = isPlusStrand == is5p? position + offset : position - offset;
				if(considered.contains(currentPosition)) continue;
				
				considered.add(currentPosition);
				
				for(ArrayList<Integer> co : bedParser.getCoordinates(currentPosition, rightWindowSize, isPlusStrand, is5p)){
					ArrayList<Integer> coordinate = new ArrayList<Integer>(co);					
			
					for(int i=1;i<=leftWindowSize;i++){
						if(is5p) coordinate.add(0, currentPosition + (isPlusStrand? -i : i));
						else coordinate.add(currentPosition + (isPlusStrand? i : -i));
					}
					
					double[] cov = is5p? bedParser.get5pCoverages(isPlusStrand, coordinate) : bedParser.get3pCoverages(isPlusStrand, coordinate);
					double score = 0;
					if(rpf.Scorer.sum(cov) < readCountThreshold) score = 0;
					else{
						double[] depth = bedParser.getDepths(isPlusStrand, coordinate);
						score = getRawScore(is5p? filter5p : filter3p, depth, is5p? filter5pNorm : filter3pNorm);
					}
				 
					if(score > scoreThreshold){
						ScoredPosition scoredPosition = is5p? new ScoredPosition(bedParser.getContig(), isPlusStrand, currentPosition, coordinate, 0, null, parser).setFivePScore(score):
							new ScoredPosition(bedParser.getContig(), isPlusStrand, 0, null, currentPosition, coordinate, parser).setThreePScore(score);						
						if(!scoredPositionMap.containsKey(currentPosition)) scoredPositionMap.put(currentPosition, new ArrayList<ScoredPosition>());
						scoredPositionMap.get(currentPosition).add(scoredPosition);
						
						/*if(!is5p){
							System.out.println(isPlusStrand + " " + is5p + " " + currentPosition);
							System.out.println(coordinate);
							for(double c : cov) System.out.print(c + " ");
							System.out.println();
						}*/
					}					
				}
			}			
		}
		return scoredPositionMap;
	}
	
	
	//TODO 
	
	int min = 20;
	int seqlength = 60;
	int max = 150;
	
	private ArrayList<ScoredPosition> getScoredPositions(AnnotationFileParser parser, double scoreThreshold, double depthThreshold, double energyThreshold, boolean isPlusStrand){
		HashSet<ScoredPosition> positionSet = new HashSet<ScoredPosition>();
		System.out.println("Scoring for " + bedParser.getContig() + " " + (isPlusStrand? "+" : "-") + " strand");
		
		HashMap<Integer, ArrayList<ScoredPosition>> sp5p = getScoredPositionMap(parser, isPlusStrand, true, scoreThreshold);
		HashMap<Integer, ArrayList<ScoredPosition>> sp3p = getScoredPositionMap(parser, isPlusStrand, false, scoreThreshold);
		
		ArrayList<Integer> positions5p = new ArrayList<Integer>(sp5p.keySet());
		ArrayList<Integer> positions3p = new ArrayList<Integer>(sp3p.keySet());
		Collections.sort(positions5p);
		Collections.sort(positions3p);
		
		//20 - 50?
		if(isPlusStrand){ // TODO opt later
			for(int p3p : positions3p){
				int i = Collections.binarySearch(positions5p, p3p + min);
				i = i<0? -i-1 : i;
				boolean pairExists = false;
				while(i < positions5p.size()){
					int p5p = positions5p.get(i);
					if(p5p > p3p + max) break;
					
					String seq = fastaParser.getSequence(bedParser.getContig(), (p3p+p5p+1-seqlength)/2, (p3p+p5p+1+seqlength)/2);
					ArrayList<Double> scores = RunRNAfold.run(seq);
					pairExists = true;
					i++;
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s3 : sp3p.get(p3p)){
						for(ScoredPosition s5 : sp5p.get(p5p)){
							ScoredPosition sp = new ScoredPosition(bedParser.getContig(), isPlusStrand, p5p, s5.getFivePcoordinate(), p3p, s3.getThreePcoordinate(), null);
							sp.setFivePScore(s5.getFivePScore()).setThreePScore(s3.getThreePScore()).setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
							sp.addGenes(s5);sp.addGenes(s3);
							positionSet.add(sp);
						}						
					}					
				}
				if(!pairExists){
					String seq = fastaParser.getSequence(bedParser.getContig(), p3p, p3p+seqlength);
					ArrayList<Double> scores = RunRNAfold.run(seq);	
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s3 : sp3p.get(p3p)){
						s3.setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
						positionSet.add(s3);
					}
				}
			}
			
			
			for(int p5p : positions5p){
				int i = Collections.binarySearch(positions3p, p5p - max);
				i = i<0? -i-1 : i;
				boolean pairExists = false;
				while(i < positions3p.size()){
					int p3p = positions3p.get(i);
					if(p5p < p3p + min) break;
					
					String seq = fastaParser.getSequence(bedParser.getContig(), (p3p+p5p+1-seqlength)/2, (p3p+p5p+1+seqlength)/2);
					ArrayList<Double> scores = RunRNAfold.run(seq);
					pairExists = true;
					i++;
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s3 : sp3p.get(p3p)){
						for(ScoredPosition s5 : sp5p.get(p5p)){
							ScoredPosition sp = new ScoredPosition(bedParser.getContig(), isPlusStrand, p5p, s5.getFivePcoordinate(), p3p, s3.getThreePcoordinate(), null);
							sp.setFivePScore(s5.getFivePScore()).setThreePScore(s3.getThreePScore()).setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
							sp.addGenes(s5);sp.addGenes(s3);
							positionSet.add(sp);
						}						
					}					
				}
				if(!pairExists){
					String seq = fastaParser.getSequence(bedParser.getContig(), p5p-seqlength, p5p);
					ArrayList<Double> scores = RunRNAfold.run(seq);	
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s5 : sp5p.get(p5p)){
						s5.setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
						positionSet.add(s5);
					}
				}
			}
		}else{
			for(int p5p : positions5p){
				int i = Collections.binarySearch(positions3p, p5p + min);
				i = i<0? -i-1 : i;
				boolean pairExists = false;
				while(i < positions3p.size()){
					int p3p = positions3p.get(i);
					if(p3p > p5p + max) break;
					String seq = ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(bedParser.getContig(), (p3p+p5p+1-seqlength)/2, (p3p+p5p+1+seqlength)/2), true);
					ArrayList<Double> scores = RunRNAfold.run(seq);
					pairExists = true;
					i++;
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s3 : sp3p.get(p3p)){
						for(ScoredPosition s5 : sp5p.get(p5p)){
							ScoredPosition sp = new ScoredPosition(bedParser.getContig(), isPlusStrand, p5p, s5.getFivePcoordinate(), p3p, s3.getThreePcoordinate(), null);
							sp.setFivePScore(s5.getFivePScore()).setThreePScore(s3.getThreePScore()).setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
							sp.addGenes(s5);sp.addGenes(s3);
							positionSet.add(sp);
						}						
					}
					
				}
				if(!pairExists){
					String seq = ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(bedParser.getContig(), p5p, p5p+seqlength), true);
					ArrayList<Double> scores = RunRNAfold.run(seq);	
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s5 : sp5p.get(p5p)){
						s5.setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
						positionSet.add(s5);
					}
				}
			}
			
			for(int p3p : positions3p){
				int i = Collections.binarySearch(positions5p, p3p - max);
				i = i<0? -i-1 : i;
				boolean pairExists = false;
				while(i < positions5p.size()){
					int p5p = positions5p.get(i);
					if(p3p < p5p + min) break;
					String seq = ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(bedParser.getContig(), (p3p+p5p+1-seqlength)/2, (p3p+p5p+1+seqlength)/2), true);
					ArrayList<Double> scores = RunRNAfold.run(seq);
					pairExists = true;
					i++;
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s3 : sp3p.get(p3p)){
						for(ScoredPosition s5 : sp5p.get(p5p)){
							ScoredPosition sp = new ScoredPosition(bedParser.getContig(), isPlusStrand, p5p, s5.getFivePcoordinate(), p3p, s3.getThreePcoordinate(), null);
							sp.setFivePScore(s5.getFivePScore()).setThreePScore(s3.getThreePScore()).setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
							sp.addGenes(s5);sp.addGenes(s3);
							positionSet.add(sp);
						}						
					}
					
				}
				if(!pairExists){
					String seq = ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(bedParser.getContig(), p3p-seqlength, p3p), true);
					ArrayList<Double> scores = RunRNAfold.run(seq);	
					if(scores.get(0) > energyThreshold) continue;
					if(scores.get(1) < depthThreshold) continue;
					for(ScoredPosition s3 : sp3p.get(p3p)){
						s3.setDepth(scores.get(1)).setEnergy(scores.get(0)).setSeq(seq);
						positionSet.add(s3);
					}
				}
			}
		}
		
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
	//	System.out.println("Haha");
		
		String sample = "Drosha2";
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		MirGff3FileParser mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/genomes/hsa.gff3");
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg19.refFlat.txt");
		String bedFileName = "/media/kyowon/Data1/fCLIP/samples/sample1/bed/" + sample + ".sorted.bed";
		String parameterFileName = "/media/kyowon/Data1/fCLIP/samples/sample1/bed/Drosha2.sorted.param";
		String outFileName = "/media/kyowon/Data1/fCLIP/samples/sample1/bed/" + sample + ".sorted.out.csv";
		PrintStream out = new PrintStream(outFileName);
		double scoreThreshold = 0.1;
		double energyThreshold = 1000;
		double depthThreshold = 0;
		out.println(ScoredPosition.getHeader());
		for(String contig : fastaParser.getContigs()){
			//System.out.println("gg1");
			Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
			//System.out.println("gg");
			Scorer scorer = new Scorer(bedParser, fastaParser, mirParser, parameterFileName);
			for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, scoreThreshold, depthThreshold, energyThreshold, true)){
				out.println(sp);
			}
			for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, scoreThreshold, depthThreshold, energyThreshold, false)){
				out.println(sp);
			}			
		}		
		out.close();
	}
	
}
