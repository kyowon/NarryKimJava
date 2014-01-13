package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import net.sf.samtools.util.BufferedLineReader;

public class BedCovFileParser {
	
	private HashMap<String, HashMap<Integer, Integer>> coverageMap;
	private AnnotationFileParser annotationParser;
	
	public BedCovFileParser(String bedCovFile, String annotationFile){
		read(bedCovFile);
		annotationParser = new AnnotationFileParser(annotationFile);
	}
	
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
		}		
		return coveragesToreturn.iterator();
	}
	
	public double[] getSqrtCoverages(String contig, int position, int leftWindowSize, int rightWindowSize, boolean isPlusStrand){
		double[] cov = getCoverages(contig, position, leftWindowSize, rightWindowSize, isPlusStrand);
		double[] sqrtCov = new double[cov.length];
		for(int i=0; i<cov.length; i++){
			sqrtCov[i] = Math.sqrt(cov[i]);
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
	
	
	private void read(String bedCovFile){
		coverageMap = new HashMap<String, HashMap<Integer, Integer>>();
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
			}
			in.close();					
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		BedCovFileParser test = new BedCovFileParser("/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/Noco_Harr_10mHsum-uncollapsed.plus.cov", "/media/kyowon/Data1/RPF_Project/genomes/refFlatHuman.txt");
	
		for(double cov : test.getCoverages("chr1", 566744, 5,10, false))
			System.out.println(cov);
	}
	
}
