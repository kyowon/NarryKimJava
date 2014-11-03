package fCLIP;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

import net.sf.picard.util.IntervalTree;
import net.sf.picard.util.IntervalTree.Node;
import parser.Bed12Parser;
import parser.BufferedLineReader;

public class MirGff3FileParser {
	
	public class fivePComparator implements Comparator<MiRNA>{
		public int compare(MiRNA o1, MiRNA o2) {
			int i= new Integer(o1.fivepPosition).compareTo(new Integer(o2.fivepPosition));
			if(i!=0) return i;
			return new Integer(o1.threepPosition).compareTo(new Integer(o2.threepPosition));
		}
	}
	public class threePComparator implements Comparator<MiRNA>{
		public int compare(MiRNA o1, MiRNA o2) {
			int i= new Integer(o1.threepPosition).compareTo(new Integer(o2.threepPosition));
			if(i!=0) return i;
			return new Integer(o1.fivepPosition).compareTo(new Integer(o2.fivepPosition));
		}
	}
	
	public static class MiRNA{
		private String contig;
		private String name;
		private int fivepPosition; // inclusive
		private int threepPosition;
		private boolean isPlusStrand;
		//private int[][] pres;
		public String getContig() {
			return contig;
		}

		public String getName() {
			return name;
		}

		public int get5p() {
			return fivepPosition;
		}

		public int get3p() {
			return threepPosition;
		}
		
		
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((contig == null) ? 0 : contig.hashCode());
			result = prime * result + fivepPosition;
			result = prime * result + (isPlusStrand ? 1231 : 1237);
			result = prime * result + ((name == null) ? 0 : name.hashCode());
			result = prime * result + threepPosition;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof MiRNA))
				return false;
			MiRNA other = (MiRNA) obj;
			if (contig == null) {
				if (other.contig != null)
					return false;
			} else if (!contig.equals(other.contig))
				return false;
			if (fivepPosition != other.fivepPosition)
				return false;
			if (isPlusStrand != other.isPlusStrand)
				return false;
			if (name == null) {
				if (other.name != null)
					return false;
			} else if (!name.equals(other.name))
				return false;
			if (threepPosition != other.threepPosition)
				return false;
			return true;
		}

		public ArrayList<Integer> get5pCoordinate(int leftWindowSize, int rightWindowSize){
			return get5pCoordinate(leftWindowSize, rightWindowSize, 0);
		}
		
		// left window inclusive
		public ArrayList<Integer> get5pCoordinate(int leftWindowSize, int rightWindowSize, int offset){
			ArrayList<Integer> coordinate = new ArrayList<Integer>();
			if(isPlusStrand){
				int start = fivepPosition - leftWindowSize;
				int end = fivepPosition + rightWindowSize;
				for(int i=start;i<end;i++){
					coordinate.add(i + offset);
				}
			}else{
				int start = fivepPosition + leftWindowSize;
				int end = fivepPosition - rightWindowSize;
				for(int i=start;i>end;i--){
					coordinate.add(i + offset);
				}
			}
			return coordinate;
		}
		
		public ArrayList<Integer> get3pCoordinate(int leftWindowSize, int rightWindowSize){
			return get3pCoordinate(leftWindowSize, rightWindowSize, 0);
		}

		public ArrayList<Integer> get3pCoordinate(int leftWindowSize, int rightWindowSize, int offset){
			ArrayList<Integer> coordinate = new ArrayList<Integer>();
			if(isPlusStrand){
				int start = threepPosition - leftWindowSize;
				int end = threepPosition + rightWindowSize;
				for(int i=start;i<end;i++){
					coordinate.add(i + offset);
				}
			}else{
				int start = threepPosition + leftWindowSize;
				int end = threepPosition - rightWindowSize;
				for(int i=start;i>end;i--){
					coordinate.add(i + offset);
				}
			}
			return coordinate;
		}

		public boolean isPlusStrand() {
			return isPlusStrand;
		}
		
		MiRNA(String s){
			String[] lines = s.split("\n");
			String[] token = lines[0].split("\t");
			this.contig = token[0];
			this.name = token[token.length-1].substring(token[token.length-1].lastIndexOf("Name=")+5);
			this.isPlusStrand = token[6].equals("+");
			if(this.isPlusStrand){
				this.fivepPosition = Integer.parseInt(token[3]);
				this.threepPosition = Integer.parseInt(token[4]);
			}else{
				this.fivepPosition = Integer.parseInt(token[4]);
				this.threepPosition = Integer.parseInt(token[3]);
			}
			//if(lines.length > 1){
				//TODO if necessary..
			//}
		}
		
		private MiRNA(int fivepPosition, int threepPosition){ // constructor for comparison
			this.fivepPosition = fivepPosition;
			this.threepPosition = threepPosition;
		}
		
		@Override
		public String toString(){
			return contig + ";" + name + ";" + fivepPosition + ";" + threepPosition + ";" + (isPlusStrand? '+' : '-'); 
		}
			
		
	}
	
	private HashMap<String, ArrayList<MiRNA>> miRNAMap;
	private HashMap<String, IntervalTree<ArrayList<Integer>>> miRNAIntervalMap;
	
	private void read(String filename){
		miRNAMap = new HashMap<String, ArrayList<MiRNA>>();
		try {
			BufferedLineReader in = new BufferedLineReader(filename);
			String s;
			String c = null;
			String contig = null;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				contig = token[0];
				if(token[2].equals("miRNA_primary_transcript")){
					if(c != null){
						if(!miRNAMap.containsKey(contig)) miRNAMap.put(contig, new ArrayList<MiRNA>());
						MiRNA mi = new MiRNA(c);
						//System.out.println(mi);
						if(mi !=null) miRNAMap.get(contig).add(mi);
					}
					c = s;
					continue;
				}
				c += "\n" + s;
			}
			if(c != null){
				if(!miRNAMap.containsKey(contig)) miRNAMap.put(contig, new ArrayList<MiRNA>());
				MiRNA mi = new MiRNA(c);
				//System.out.println(mi);
				if(mi !=null) miRNAMap.get(contig).add(mi);
			}
			
			miRNAIntervalMap = new HashMap<String, IntervalTree<ArrayList<Integer>>>();
			for(String con : miRNAMap.keySet()){
				miRNAIntervalMap.put(con, new IntervalTree<ArrayList<Integer>>());
				IntervalTree<ArrayList<Integer>> tree = miRNAIntervalMap.get(con);
				
				for(int i=0;i<miRNAMap.get(con).size();i++){
					MiRNA miRNA = miRNAMap.get(con).get(i);
					int start = miRNA.get5p();
					int end = miRNA.get3p();
					if(!miRNA.isPlusStrand){
						int tmp = start;
						start = end;
						end = tmp;
					}
					Node<ArrayList<Integer>> node = tree.find(start, end);
					ArrayList<Integer> v;
					
					if(node == null)	
						v = new ArrayList<Integer>();
					else 
						v = node.getValue();
					
					v.add(i);
					tree.put(start, end, v);
				}
			}
			
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public MirGff3FileParser(String filename){
		read(filename);
	}
	
	public Iterator<MiRNA> getMiRNAIterator(){
		ArrayList<MiRNA> allmiRNAs = new ArrayList<MiRNA>();
		for(String contig : miRNAMap.keySet())
			allmiRNAs.addAll(miRNAMap.get(contig));
		return allmiRNAs.iterator();
	}
	
	public Iterator<MiRNA> getMiRNAIterator(String contig){
		ArrayList<MiRNA> allmiRNAs = miRNAMap.get(contig);
		if(allmiRNAs != null)
			return allmiRNAs.iterator();
		return null;
	}
	
	public static void main(String[] args){
		new MirGff3FileParser("/media/kyowon/Data1/fCLIP/Genome/hsa.gff3");
	}
	
	/**
	    * Get the containing miRNA
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the tx start position
	    * @return the containing miRNA
	    */
	public ArrayList<MiRNA> getContainingMiRNAs(String contig, boolean isPlusStrand, int position){
		ArrayList<MiRNA> ret = new ArrayList<MiRNA>();
		ArrayList<MiRNA> miRNAs = miRNAMap.get(contig);
		if(miRNAs != null)// ret = null;
		{
			Iterator<Node<ArrayList<Integer>>> it = miRNAIntervalMap.get(contig).overlappers(position, position);
			while(it.hasNext()){
				Node<ArrayList<Integer>> is = it.next();
				for(int i : is.getValue()){
					MiRNA miRNA = miRNAs.get(i);
					if(miRNA.isPlusStrand() == isPlusStrand)
						ret.add(miRNA);
				}
			}
		}
		return ret.isEmpty()? null : ret;
	}

	
	
}
