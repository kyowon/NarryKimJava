package util;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import fCLIP.parser.ScoringOutputParser;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.Bed12Parser.Bed12Read;

public class FoldChangerCalculator {
	public static class FoldChangeFilter{
		private int t1, t2;
		boolean shouldSatisfyBothThresholds;
		
		private FoldChangeFilter(int t1, int t2, boolean shouldSatisfyBothThresholds){
			this.t1 = t1;
			this.t2 = t2;
			this.shouldSatisfyBothThresholds = shouldSatisfyBothThresholds;
		}
		
		public final static FoldChangeFilter nullFilter = new FoldChangeFilter(-1,-1, false);
		public final static FoldChangeFilter unMappedGeneFilter = new FoldChangeFilter(1, 1, true);
		public final static FoldChangeFilter singleSampleMappedGeneFilter = new FoldChangeFilter(1, 1, false);
		public static FoldChangeFilter getFilterWithReadThresholds(int t1, int t2, boolean shouldSatisfyBothThresholds) { return new FoldChangeFilter(t1, t2, shouldSatisfyBothThresholds); }
		
		private boolean pass(int r1, int r2){
			if(shouldSatisfyBothThresholds) return r1 >= t1 && r2 >= t2;			
			return r1 >= t1 || r2 >= t2;
		}
		
	}
	
	
	//private Bed12Parser parser1, parser2;
	private AnnotationFileParser annotationFileParser;
	private String bed1, bed2;
	private int sudoCount = 1;
	private double median = 0;
	private boolean invertStrand = false;
	private static FoldChangeFilter filter = FoldChangeFilter.unMappedGeneFilter;
	//private HashMap<String, Bed12Parser> bedParserMap1;
	//private HashMap<String, Bed12Parser> bedParserMap2;
	private HashMap<AnnotatedGene, Double> cMap;
	private HashMap<AnnotatedGene, Double> tMap;
	private Bed12Parser bedParser1 = null, bedParser2 = null;
	
	public double getMedian(){
		return median;
	}
	
	public HashMap<AnnotatedGene, double[]> getReadCounts(){
		HashMap<AnnotatedGene, double[]> ret = new HashMap<AnnotatedGene, double[]>();
		for(AnnotatedGene gene : cMap.keySet()){
			double[] v = { cMap.get(gene), tMap.get(gene)};
			ret.put(gene, v);
		}
		return ret;
	}
	
	public FoldChangerCalculator(String bed1, String bed2, AnnotationFileParser annotationFileParser, boolean invertStrand){
		//this.annotationFileParser = annotationFileParser;
		//bedParserMap1 = new HashMap<String, Bed12Parser>();
		//bedParserMap2 = new HashMap<String, Bed12Parser>();
		System.out.println("Initiating fold change calculation...");
		this.invertStrand = invertStrand;
		this.bed1 = bed1;
		this.bed2 = bed2;
		this.annotationFileParser = annotationFileParser;
		
		if(false){
			cMap = new HashMap<AnnotationFileParser.AnnotatedGene, Double>();
			tMap = new HashMap<AnnotationFileParser.AnnotatedGene, Double>();
			
			ArrayList<Double> fcs = new ArrayList<Double>();
			for(String contig : annotationFileParser.getContigs()){
				if(contig.length() > 5 || contig.equals("chrM")) continue;
				
				System.out.println(contig);
				Bed12Parser parser1 = new Bed12Parser(bed1, annotationFileParser, contig, false, invertStrand);
				Bed12Parser parser2 = new Bed12Parser(bed2, annotationFileParser, contig, false, invertStrand);
				System.out.println("Bed reading done");
				
			//	bedParserMap1.put(contig, parser1);
			//	bedParserMap2.put(contig, parser2);
							
				Iterator<AnnotatedGene> iterator = annotationFileParser.getAnnotatedGeneIterator(contig);
				while(iterator.hasNext()){
					AnnotatedGene gene = iterator.next();
					//System.out.println(gene);
					double r1 = 0;
					for(double r : parser1.getCoverages(gene)) r1 += r;
					
					double r2 = 0;
					for(double r : parser2.getCoverages(gene)) r2 += r;
				
					if(!filter.pass((int)r1, (int)r2)) continue;
					double fc = Math.log10((r1+sudoCount)/(r2 + sudoCount));
					
					fcs.add(fc);
					cMap.put(gene, r2);
					tMap.put(gene, r1);
					//cMap.put(gene, fc);
				}			
			}
			
			Collections.sort(fcs);
			median = fcs.get(fcs.size()/2);
			System.out.println("Fold change median : " + median); // 0.04717907769546925
		}
	}
	
	//public double getFoldChange(String contig, int p1, int p2, boolean isPlusStrand){
		
	//}
	
	public double getControlReadCount(AnnotatedGene gene){
		Double fc = cMap.get(gene);
		if(fc == null) return 0.0;
		return cMap.get(gene);
	}
	
	
	public double getTargetReadCount(AnnotatedGene gene){
		Double fc = tMap.get(gene);
		if(fc == null) return 0.0;
		return tMap.get(gene);
	}
	
	public double getTargetReadCount(String contig, int start, int end, boolean isPlusStrand){ // inclusive
		if(bedParser1 == null || !bedParser1.getContig().equals(contig)){
			System.out.print("Target Reading " + contig);
			bedParser1 =  new Bed12Parser(bed1, annotationFileParser, contig, true, invertStrand);
			System.out.println(" Done");
		}
		
		double c1 = 0;
		
		if(isPlusStrand){
			for(int i=start;i<=end;i++){
				ArrayList<Bed12Read> reads1 = bedParser1.get5pReads(isPlusStrand, i);
				if(reads1 != null){
					for(Bed12Read read : reads1){
						if(read.get3p()<=end) c1++;
					}
					//System.out.println(i + " " + reads1.size() + " " + c1);
				}				
			}
		}else{
			for(int i=start;i>=end;i--){
				ArrayList<Bed12Read> reads1 = bedParser1.get5pReads(isPlusStrand, i);
				if(reads1 != null){
					for(Bed12Read read : reads1){
						if(read.get3p()>=start) c1++;
					}
				}
			}
		}
		//System.out.println(contig + " " + bedParser1.getContig() + " " +  start + " " + end + " " + isPlusStrand + " " + c1 + " " + c2);
		return c1;
	}
	
	public double getControlReadCount(String contig, int start, int end, boolean isPlusStrand){ // inclusive
		if(bedParser2 == null || !bedParser2.getContig().equals(contig)){
			System.out.print("Control Reading " + contig);
			bedParser2 =  new Bed12Parser(bed2, annotationFileParser, contig, true, invertStrand);
			System.out.println(" Done");
		}
		
		double c1 = 0;
		
		if(isPlusStrand){
			for(int i=start;i<=end;i++){
				ArrayList<Bed12Read> reads1 = bedParser2.get5pReads(isPlusStrand, i);
				if(reads1 != null){
					for(Bed12Read read : reads1){
						if(read.get3p()<=end) c1++;
					}
					//System.out.println(i + " " + reads1.size() + " " + c1);
				}				
			}
		}else{
			for(int i=start;i>=end;i--){
				ArrayList<Bed12Read> reads1 = bedParser2.get5pReads(isPlusStrand, i);
				if(reads1 != null){
					for(Bed12Read read : reads1){
						if(read.get3p()>=start) c1++;
					}
				}
			}
		}
		//System.out.println(contig + " " + bedParser1.getContig() + " " +  start + " " + end + " " + isPlusStrand + " " + c1 + " " + c2);
		return c1;
	}
	
	public static void setFilter(FoldChangeFilter filter){
		FoldChangerCalculator.filter = filter;
	}
	
	

}
