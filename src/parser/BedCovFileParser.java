package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import net.sf.samtools.util.BufferedLineReader;

public class BedCovFileParser {
	
	private class BedCoverage implements Comparable<BedCoverage>{ // contigs are not considered for memory usage
		private int position;
		private int coverage;
		private double sqrtCoverage;
		
		private BedCoverage(int position, int coverage){
			this.position = position;
			this.coverage = coverage;
			this.sqrtCoverage = coverage == 0? 0.0 : Math.sqrt(coverage);
		}
		
		@Override
		public boolean equals(Object o){
			if(o instanceof BedCoverage){
				BedCoverage other = (BedCoverage)o;
				return this.position == other.position;
			}
			return false;
		}
		
		@Override
		public int hashCode(){
			return new Long(this.position).hashCode();			
		}
		
		public double getSqrtCoverage(){
			return sqrtCoverage;
		}
		
		public int getCoverage(){
			return coverage;
		}

		public int getPosition(){
			return position;
		}
		
		public int compareTo(BedCoverage o) {
			return new Long(this.position).compareTo(new Long(o.position)); 
		}
		
		@Override
		public String toString(){
			return position + " " + coverage;
		}		
		
	}
	
	private HashMap<String, ArrayList<BedCoverage>> bedCoverageMap;
	
	
	public BedCovFileParser(String bedCovFile){
		read(bedCovFile);
	}
	
	public ArrayList<String> getContigs(){
		ArrayList<String> contigs = new ArrayList<String>(bedCoverageMap.keySet());
		Collections.sort(contigs);
		return contigs;
	}
	
	public Iterator<Integer> getNonZeroCoveragePositionIterator(String contig){
		ArrayList<Integer> positions = new ArrayList<Integer>();
		Iterator<BedCoverage> iterator = getBedCoverageIterator(contig, 0, Integer.MAX_VALUE);
		while(iterator.hasNext()){
			positions.add(iterator.next().getPosition());
		}
		return positions.iterator();
	}
	
	private Iterator<BedCoverage> getBedCoverageIterator(String contig, int start, int end){// start inclusive, end exclusive
		ArrayList<BedCoverage> bedCoverages = bedCoverageMap.get(contig);
		ArrayList<BedCoverage> bedCoveragesToreturn = new ArrayList<BedCoverage>();
		
		if(bedCoverages == null) return bedCoveragesToreturn.iterator();
		int startIndex = Collections.binarySearch(bedCoverages, new BedCoverage(start, 0));
		startIndex = startIndex < 0? - startIndex - 1 : startIndex;
		int endIndex = Collections.binarySearch(bedCoverages, new BedCoverage(end, 0));
		endIndex = endIndex < 0? - endIndex - 1 : endIndex;
			for(int i=startIndex;i<Math.min(endIndex, bedCoverages.size());i++)
			bedCoveragesToreturn.add(bedCoverages.get(i));
	
		return bedCoveragesToreturn.iterator();
	}
	
	public double[] getSqrtCoverages(String contig, int position, int leftWindowSize, int rightWindowSize, boolean isPlusStrand){
		double[] coverages = new double[leftWindowSize + rightWindowSize];
		int start, end;
		if(isPlusStrand){
			start = position - leftWindowSize;
			end = position + rightWindowSize;
		}else{
			end = position + leftWindowSize + 1;
			start = position - rightWindowSize + 1;
		}
		Iterator<BedCoverage> iterator = getBedCoverageIterator(contig, start, end);
		
		while(iterator.hasNext()){
			BedCoverage bedCoverage = iterator.next();
			int index = bedCoverage.getPosition() - start;
			index = isPlusStrand? index : coverages.length - index - 1;
				coverages[index] = bedCoverage.getSqrtCoverage();
		}		
		return coverages;
	}
	
	private void read(String bedCovFile){
		bedCoverageMap = new HashMap<String, ArrayList<BedCoverage>>();
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(bedCovFile));
			String s;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				String contig = token[0];
				int position = Integer.parseInt(token[1]);
				int coverage = (int)Double.parseDouble(token[2]);
				if(!bedCoverageMap.containsKey(contig)) bedCoverageMap.put(contig, new ArrayList<BedCoverage>());
				bedCoverageMap.get(contig).add(new BedCoverage(position, coverage));
			}
			in.close();
			
			for(String contig : bedCoverageMap.keySet()){
				Collections.sort(bedCoverageMap.get(contig));
			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		BedCovFileParser test = new BedCovFileParser("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.plus.cov");
	
		for(double cov : test.getSqrtCoverages("chr1", 566743, 5,10, true))
			System.out.println(cov);
	}
	
}
