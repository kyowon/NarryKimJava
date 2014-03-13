package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import parser.AnnotationFileParser.AnnotatedGene;
import net.sf.samtools.util.BufferedLineReader;

public class BedCovFileParser {
	
	private HashMap<String, HashMap<Integer, Integer>> coverageMap;
	private AnnotationFileParser annotationParser;
	private int totalReadCount = 0;
	
	public BedCovFileParser(String bedCovFile, AnnotationFileParser annotationParser){
		read(bedCovFile);
		this.annotationParser = annotationParser;
	}
	
	public int getTotalReadCount() { return totalReadCount; }
	
	public ArrayList<String> getContigs(){
		ArrayList<String> contigs = new ArrayList<String>(coverageMap.keySet());
		Collections.sort(contigs);
		return contigs;
	}
	
	public Iterator<Integer> getNonZeroCoveragePositionIterator(String contig){
		ArrayList<Integer> positions = new ArrayList<Integer>();
		HashMap<Integer, Integer> subCoverageMap = coverageMap.get(contig);
		for(int key : subCoverageMap.keySet()){
			positions.add(key);
		}
		Collections.sort(positions);
		return positions.iterator();
	}
	
	private Iterator<Integer> getBedCoverageIterator(String contig, int start, int length, boolean isPlusStrand){// start inclusive
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		ArrayList<Integer> coveragesToreturn = new ArrayList<Integer>();
		
		if(coverages == null) return coveragesToreturn.iterator();		
		ArrayList<ArrayList<Integer>> positions = annotationParser.getLiftOverPositions(contig, isPlusStrand, start, length);			
		if(positions == null || positions.isEmpty()) return coveragesToreturn.iterator();
		
		for(ArrayList<Integer> subPositions : positions){
			for(int position : subPositions){
				Integer coverage = coverages.get(position);
				coveragesToreturn.add(coverage == null ? 0 : coverage);
			}
			//System.out.println(position + " " + coverage);
		}	
		//System.out.println(positions);
		return coveragesToreturn.iterator();
	}
	
	
	private Iterator<Integer> getBedCoverageIterator(String contig, int start, boolean isPlusStrand, int maxLength, ZeroBasedFastaParser fastaParser){// start inclusive
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		ArrayList<Integer> coveragesToreturn = new ArrayList<Integer>();
		
		if(coverages == null) return coveragesToreturn.iterator();		
		
		ArrayList<ArrayList<Integer>> positions = annotationParser.getLiftOverPositionsTillNextStopCodon(contig, isPlusStrand, start, maxLength, fastaParser);
		if(positions == null || positions.isEmpty()) return coveragesToreturn.iterator();
		
		for(ArrayList<Integer> subPositions : positions){
			for(int position : subPositions){
				Integer coverage = coverages.get(position);
				coveragesToreturn.add(coverage == null ? 0 : coverage);
			}
			//System.out.println(position + " " + coverage);
		}	
		//System.out.println(positions);
		return coveragesToreturn.iterator();
	}
	
	
	public double[] getSqrtCoverages(String contig, int position, int leftWindowSize, int rightWindowSize, boolean isPlusStrand){
		double[] cov = getCoverages(contig, position, leftWindowSize, rightWindowSize, isPlusStrand);
		double[] sqrtCov = new double[cov.length];
		for(int i=0; i<cov.length; i++){
			if(cov[i]>0) sqrtCov[i] = Math.sqrt(cov[i]);
		}
		return sqrtCov;
	}
	
	public double[] getCoverages(String contig, int position, int leftWindowSize, int rightWindowSize, boolean isPlusStrand){
		double[] coverages = new double[leftWindowSize + rightWindowSize];
		int start;
		if(isPlusStrand){
			start = position - leftWindowSize;
		}else{
			start = position + leftWindowSize;
		}
		Iterator<Integer> iterator = getBedCoverageIterator(contig, start, coverages.length, isPlusStrand);
		
		int index = 0;
		while(iterator.hasNext()){
			int coverage = iterator.next();
			coverages[index++] = coverage;
		}		
		return coverages;
	}
	
	public double[] getSqrtCoveragesTillnextStopCodon(String contig, int position, int leftWindowSize, boolean isPlusStrand, int maxLength, ZeroBasedFastaParser fastaParser){
		double[] cov = getCoveragesTillnextStopCodon(contig, position, leftWindowSize, isPlusStrand, maxLength, fastaParser);
		double[] sqrtCov = new double[cov.length];
		for(int i=0; i<cov.length; i++){
			if(cov[i]>0) sqrtCov[i] = Math.sqrt(cov[i]);
		}
		return sqrtCov;
	}
	
	public double[] getCoveragesTillnextStopCodon(String contig, int position, int leftWindowSize, boolean isPlusStrand, int maxLength, ZeroBasedFastaParser fastaParser){
		ArrayList<Integer> coverageList = new ArrayList<Integer>();
		int start;
		if(isPlusStrand){
			start = position - leftWindowSize;
		}else{
			start = position + leftWindowSize;
		}
		Iterator<Integer> iterator = getBedCoverageIterator(contig, start, isPlusStrand, maxLength, fastaParser);
		
		while(iterator.hasNext()){
			int coverage = iterator.next();
			coverageList.add(coverage);
		}	
		double[] coverages = new double[coverageList.size()];
		for(int i=0;i<coverages.length;i++){
			coverages[i] = coverageList.get(i);
		}		
		return coverages;
	}
	
	
	public double getTotalCDSCoverage(AnnotatedGene gene, boolean normalizeByLength){
		HashMap<Integer, Integer> coverages = coverageMap.get(gene.getContig());
		if(coverages == null) return 0;
		double totalCoverages = 0;
		int tl = 0;
		ArrayList<ArrayList<Integer>> positions = annotationParser.getLiftOverCDSPositions(gene);
	
		for(ArrayList<Integer> subPositions : positions){
			for(int position : subPositions){
				Integer coverage = coverages.get(position);
				tl++;
				if(coverage == null) continue;
				totalCoverages += coverage;
			}	
		}
		if(normalizeByLength) totalCoverages/=tl;
		return totalCoverages;
	}
	
	public double getTotalCoverageTillnextStopCodon(String contig, boolean isPlusStrand, int position, int maxLength, ZeroBasedFastaParser fastaParser, boolean normalizeByLength){ // use getCoveragesTillnextStopCodon - fix later TODO
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		if(coverages == null) return 0;
		double totalCoverages = 0;
		int tl = 0;
		ArrayList<ArrayList<Integer>> positions = annotationParser.getLiftOverPositionsTillNextStopCodon(contig, isPlusStrand, position, maxLength, fastaParser);
		//System.out.println(positions.size());
		for(ArrayList<Integer> subPositions : positions){
			for(int p : subPositions){
				Integer coverage = coverages.get(p);
				tl++;
				if(coverage == null) continue;
				totalCoverages += coverage;
			}	
		}
		if(normalizeByLength) totalCoverages/=tl;
		return totalCoverages;
	}
	
	public ArrayList<Double> getCoverageBeforeNAfternextStopCodonNStopCodonPosition(String contig, boolean isPlusStrand, int position, int length, int maxLength, ZeroBasedFastaParser fastaParser, boolean normalizeByLength){
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		ArrayList<Double> ret = new ArrayList<Double>();
		if(coverages == null){
			return null;
		}
		ArrayList<ArrayList<Integer>> positions = annotationParser.getLiftOverPositionsTillNextStopCodon(contig, isPlusStrand, position, maxLength, fastaParser);
		double before = 0;
		int n = 0;
		int stopPosition = -1;
		for(int i=positions.size()-1;i>=0;i--){
			ArrayList<Integer> subPositions = positions.get(i);
			for(int j = subPositions.size()-1;j>=0;j--){
				if(stopPosition < 0) stopPosition = subPositions.get(j);
				Integer coverage = coverages.get(subPositions.get(j));
				if(n++>=length) break;
				if(coverage == null) continue;
				before += coverage;
			}
			if(n>=length) break;
		}
		if(normalizeByLength) before /= length;
		ret.add(before);
		double after = 0;
		
		positions = annotationParser.getLiftOverPositions(contig, isPlusStrand, stopPosition + (isPlusStrand? 1 : -1), length);
		n = 0;
		for(ArrayList<Integer> subPositions : positions){
			for(int p : subPositions){
				Integer coverage = coverages.get(p);
				if(n++>=length) break;
				if(coverage == null) continue;
				after += coverage;
			}
			if(n>=length) break;
		}
		if(normalizeByLength) after /= length;
		ret.add(after);
		ret.add((double)stopPosition);
		return ret;
	}
	
	
	public double getTotalCoverage(String contig, boolean isPlusStrand, int position, int length, boolean normalizeByLength){
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		if(coverages == null) return 0;
		double totalCoverages = 0;
		int tl = 0;
		ArrayList<ArrayList<Integer>> positions = annotationParser.getLiftOverPositions(contig, isPlusStrand, position, length);
		for(ArrayList<Integer> subPositions : positions){
			for(int p : subPositions){
				Integer coverage = coverages.get(p);
				tl++;
				if(coverage == null) continue;
				totalCoverages += coverage;
			}
		}	
		if(normalizeByLength) totalCoverages/=tl;
		return totalCoverages;
	}
	
	private void read(String bedCovFile){
		coverageMap = new HashMap<String, HashMap<Integer, Integer>>();
		totalReadCount = 0;
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(bedCovFile));
			String s;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				String contig = token[0];
				int position = Integer.parseInt(token[1]);
				int coverage = (int)Double.parseDouble(token[2]);
				if(!coverageMap.containsKey(contig)) coverageMap.put(contig, new HashMap<Integer, Integer>());
				coverageMap.get(contig).put(position, coverage);
				totalReadCount += coverage;
			}
			in.close();					
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		BedCovFileParser test = new BedCovFileParser("/media/kyowon/Data1/RPF_Project/samples/sample3/coverages/RPF-30_1-uncollapsed.plus.cov", 
				new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.refFlat.txt"));
	//chr10	127852386
		
		int i=0;
		for(double cov : test.getCoverages("chr5", 3803296, 30,100, true))
			System.out.println(i++ + " " + cov);
	}
	
}
