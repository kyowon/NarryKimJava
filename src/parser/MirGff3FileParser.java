package parser;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

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
		
		private int prifp;
		private int pritp;
		private Integer fivepPosition = null; // inclusive, pre position
		private Integer threepPosition = null;// inclusive, pre position
		private boolean isPlusStrand;
		//private int[][] pres;
		public String getContig() {
			return contig;
		}

		public String getName() {
			return name;
		}

		public Integer getPri5p(){
			return prifp;
		}
		
		public Integer getPri3p(){
			return pritp;
		}
				
		public Integer get5p() {
			return fivepPosition;
		}

		public Integer get3p() {
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
		//	int prifivepPosition = ;
			this.pritp = Integer.parseInt(token[isPlusStrand? 4:3]);
			this.prifp = Integer.parseInt(token[isPlusStrand? 3:4]);
			
			
			if(lines.length > 2){
				String[] stoken1 = lines[1].split("\t");
				String[] stoken2 = lines[2].split("\t");
				
				if(stoken1[stoken1.length-1].split(";")[2].endsWith("-5p")){
					this.fivepPosition = Integer.parseInt(stoken1[isPlusStrand? 3:4]) - 1;
					this.threepPosition = Integer.parseInt(stoken2[isPlusStrand? 4:3]) - 1;
				}else{
					this.fivepPosition = Integer.parseInt(stoken2[isPlusStrand? 3:4]) - 1;
					this.threepPosition = Integer.parseInt(stoken1[isPlusStrand? 4:3]) - 1;
				}				
			}else{
				String[] stoken = lines[1].split("\t");
				String name = stoken[stoken.length-1].split(";")[2];
				if(name.endsWith("-5p")){
					this.fivepPosition = Integer.parseInt(stoken[isPlusStrand? 3:4]) - 1;
				}else if(name.endsWith("-3p")){
					this.threepPosition = Integer.parseInt(stoken[isPlusStrand? 4:3]) - 1;
				}else{
					int n1 = Integer.parseInt(stoken[3]) - 1;
					int n2 = Integer.parseInt(stoken[4]) - 1;
					int p1 = Integer.parseInt(token[3]) - 1;
					int p2 = Integer.parseInt(token[4]) - 1;
					
					if(this.isPlusStrand == n1 - p1 > p2 - n2){
						this.threepPosition = isPlusStrand? n2 : n1;
					}else{
						this.fivepPosition = isPlusStrand? n1 : n2;
					}
				}
			}
			
			/*if(isPlusStrand){
				if(this.fivepPosition != null) this.fivepPosition--;
			}else{
				if(this.threepPosition != null) this.threepPosition--;
			}*/
			
				
			//if(lines.length > 1){
				//TODO if necessary..
			//}
			//System.out.println(s);
		//	System.out.println(this);
		//	System.exit(0);
			
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
				contig = token[0];				
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
					//if(miRNA.get3p() == null || miRNA.get5p() == null) continue;
					int start = miRNA.prifp;
					int end = miRNA.pritp;
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
	
	
	public Iterator<MiRNA> getMiRNAIterator(String contig){
		ArrayList<MiRNA> allmiRNAs = miRNAMap.get(contig);
		if(allmiRNAs != null)
			return allmiRNAs.iterator();
		return null;
	}
	
	public MiRNA getMiRNA(String contig, String name){ // very slow..
		ArrayList<MiRNA> allmiRNAs = miRNAMap.get(contig);
		if(allmiRNAs != null){
			for(MiRNA mr : allmiRNAs){
				if(mr.name.equals(name)) return mr;
			}
		}
		return null;
	}
	
	
	public Iterator<MiRNA> getMiRNAIterator(){
		ArrayList<MiRNA> allmiRNAs = new ArrayList<MiRNA>();
		
		for(String contig : miRNAMap.keySet()){
			for(MiRNA mi : miRNAMap.get(contig)){
				allmiRNAs.add(mi);
			}
		}
		return allmiRNAs.iterator();
	}
	
	

	/**
	    * Get the containing miRNA
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the tx start position
	    * @return the containing miRNA
	    */
	public ArrayList<MiRNA> getContainingMiRNAs(String contig, boolean isPlusStrand, int position, int margin){
		ArrayList<MiRNA> ret = new ArrayList<MiRNA>();
		ArrayList<MiRNA> miRNAs = miRNAMap.get(contig);
		if(miRNAs != null)// ret = null;
		{
			Iterator<Node<ArrayList<Integer>>> it = miRNAIntervalMap.get(contig).overlappers(position-margin, position+margin);
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

	public void correctAnnotationUsignMirGeneDB(String file){
		HashMap<String, HashMap<String, String>> miRNAsInMDB = new HashMap<String, HashMap<String,String>>();
		//HashSet<String> miRNANamesInMDB = new HashSet<String>();
		try {
			BufferedLineReader in = new BufferedLineReader(file);
			String s;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				String name = token[1].trim();
				//miRNANamesInMDB.add(name);
				String contig = token[5].trim();
				int l = Integer.parseInt(token[6].trim());
				int r = Integer.parseInt(token[7].trim());
				boolean isPlusStrand = token[8].trim().equals("+");
				int fp = isPlusStrand? l : r;
				int tp = isPlusStrand? r : l;
				if(!miRNAsInMDB.containsKey(contig)) miRNAsInMDB.put(contig, new HashMap<String, String>());
				miRNAsInMDB.get(contig).put(name,fp + "\t" + tp);						
			}
			in.close();
			
			for(String contig : miRNAsInMDB.keySet()){
				ArrayList<MiRNA> miRNAs = miRNAMap.get(contig);
				ArrayList<MiRNA> filteredmiRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
				HashMap<String, String> mirRNAsInMDBsub = miRNAsInMDB.get(contig);
				for(MiRNA m : miRNAs){
					if(mirRNAsInMDBsub.containsKey(m.getName())){
						String mm = mirRNAsInMDBsub.get(m.getName());
						String[] tt = mm.split("\t");
						m.fivepPosition = Integer.parseInt(tt[0])-1;
						m.threepPosition = Integer.parseInt(tt[1])-1;
						filteredmiRNAs.add(m);
					}
				}
				miRNAMap.put(contig, filteredmiRNAs);				
			}			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException{
		MirGff3FileParser mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3");
		
		Iterator<MiRNA> mirit = mirParser.getMiRNAIterator();
		int i = 0;
		int j = 0;
		while(mirit.hasNext()){
			MiRNA mirna = mirit.next();
			i++;
			if(mirna.get3p() == null || mirna.get5p() == null) j++;
		}
		
		System.out.println(i + " " + j);
		/*String csv = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.csv";
		PrintStream outCsv = new PrintStream(csv+"modified.csv");
		ScoredPositionOutputParser sparser = new ScoredPositionOutputParser(csv);
		int mirBaseSearchMargin = 5;
		outCsv.println(sparser.getHeader());
		for(ScoredPosition sp : sparser.getPositions()){
			ArrayList<MiRNA> miRNAs = null;
			boolean isPlusStrand = sp.isPlusStrand();
			if(sp.is5pScored()){
				miRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
				ArrayList<MiRNA> matched = mirParser.getContainingMiRNAs(sp.getContig(), isPlusStrand, sp.getFivePPosition(),mirBaseSearchMargin);
				if(matched != null) miRNAs.addAll(matched);
			}
			if(sp.is3pScored()){
				if(miRNAs == null){
					miRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
					ArrayList<MiRNA> matched = mirParser.getContainingMiRNAs(sp.getContig(), isPlusStrand, sp.getThreePPosition(),mirBaseSearchMargin);
					if(matched != null) miRNAs.addAll(matched);
				}else{
					ArrayList<MiRNA> intersectedMiRNAs = new ArrayList<MirGff3FileParser.MiRNA>();
					ArrayList<MiRNA> matched = mirParser.getContainingMiRNAs(sp.getContig(), isPlusStrand, sp.getThreePPosition(),mirBaseSearchMargin);
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
			
			outCsv.println(sp);
			
		}
		
		
		outCsv.close();
		*/
		
	}
	
	
}
