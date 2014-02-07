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
	
	private Iterator<Integer> getBedCoverageIterator(String contig, int start, int length, boolean isPlusStrand){// start inclusive TODO
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		ArrayList<Integer> coveragesToreturn = new ArrayList<Integer>();
		
		if(coverages == null) return coveragesToreturn.iterator();		
		ArrayList<Integer> positions = annotationParser.getLiftOverPositions(contig, isPlusStrand, start, length);			
		if(positions == null) return coveragesToreturn.iterator();
		
		for(int position : positions){
			Integer coverage = coverages.get(position);
			coveragesToreturn.add(coverage == null ? 0 : coverage);
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
	
	public double getTotalCDSCoverage(AnnotatedGene gene, boolean normalizeByLength){
		HashMap<Integer, Integer> coverages = coverageMap.get(gene.getContig());
		if(coverages == null) return 0;
		double totalCoverages = 0;
		ArrayList<Integer> positions = annotationParser.getLiftOverCDSPositions(gene);
	
		for(int position : positions){
			Integer coverage = coverages.get(position);
			if(coverage == null) continue;
			totalCoverages += coverage;
		}	
		if(normalizeByLength) totalCoverages/=positions.size();
		return totalCoverages;
	}
	
	public double getTotalCoverageTillnextStopCodon(String contig, boolean isPlusStrand, int position, int maxLength, ZeroBasedFastaParser fastaParser, boolean normalizeByLength){
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		if(coverages == null) return 0;
		double totalCoverages = 0;
		ArrayList<Integer> positions = annotationParser.getLiftOverPositionsTillNextStopCodon(contig, isPlusStrand, position, maxLength, fastaParser);
		//System.out.println(positions.size());
		for(int p : positions){
			Integer coverage = coverages.get(p);
			if(coverage == null) continue;
			totalCoverages += coverage;
		}	
		if(normalizeByLength) totalCoverages/=positions.size();
		return totalCoverages;
	}
	
	public double getTotalCoverage(String contig, boolean isPlusStrand, int position, int length, boolean normalizeByLength){
		HashMap<Integer, Integer> coverages = coverageMap.get(contig);
		if(coverages == null) return 0;
		double totalCoverages = 0;
		ArrayList<Integer> positions = annotationParser.getLiftOverPositions(contig, isPlusStrand, position, length);
		for(int p : positions){
			Integer coverage = coverages.get(p);
			if(coverage == null) continue;
			totalCoverages += coverage;
		}	
		if(normalizeByLength) totalCoverages/=positions.size();
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
		BedCovFileParser test = new BedCovFileParser("/media/kyowon/Data1/RPF_Project/samples/sample3/coverages/Harr_C-uncollapsed.plus.cov", 
				new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.refFlat.txt"));
	//chr10	127852386

		int i=0;
		for(double cov : test.getCoverages("chr10", 127852386, 30,10, true))
			System.out.println(i++ + " " + cov);
	}
	
}
