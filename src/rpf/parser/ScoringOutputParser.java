package rpf.parser;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BufferedLineReader;
import rpf.Quantifier;

public class ScoringOutputParser {
	//private HashMap<String, HashMap<Integer, >>
	
	public class ScoredPosition implements Comparable<ScoredPosition>{
	

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((contig == null) ? 0 : contig.hashCode());
			result = prime * result + (isPlusStrand ? 1231 : 1237);
			result = prime
					* result
					+ ((mappedRegionStartsEnds == null) ? 0
							: mappedRegionStartsEnds.hashCode());
			result = prime * result + position;
			return result;
		}
		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			ScoredPosition other = (ScoredPosition) obj;		
			if (contig == null) {
				if (other.contig != null)
					return false;
			} else if (!contig.equals(other.contig))
				return false;
			if (isPlusStrand != other.isPlusStrand)
				return false;
			if (mappedRegionStartsEnds == null) {
				if (other.mappedRegionStartsEnds != null)
					return false;
			} else if (!mappedRegionStartsEnds
					.equals(other.mappedRegionStartsEnds))
				return false;
			if (position != other.position)
				return false;
			return true;
		}
		private String contig;		
		private boolean isPlusStrand;
		private int position;
		private double startScore;
		//private ArrayList<Integer> coordinate; // save memory!!
		//private double stopScore;
		public ArrayList<ArrayList<Integer>> mappedRegionStartsEnds;
		private String seq;
		private boolean isAnnotated = false;
		private AnnotatedGene gene = null;
		private String genomicRegion;
		private String frameShift;
		private HashSet<Integer> rpfrnaIndices = new HashSet<Integer>();
		private HashSet<Integer> harrIndices = new HashSet<Integer>();
	
		public HashSet<Integer> getRPFRNAIndices() {
			return rpfrnaIndices;
		}
		public void setRPFRNAIndices(HashSet<Integer> conditionIndices) {
			this.rpfrnaIndices = conditionIndices;
		}
		
		public String getContig() {
			return contig;
		}
		public boolean isPlusStrand() {
			return isPlusStrand;
		}
		public int getPosition() {
			return position;
		}
		public double getScore() {
			return startScore;
		}
		public String getSequence() {
			return seq;
		}
		public boolean isAnnotated(){
			return isAnnotated;
		}
		public AnnotatedGene getGene(){
			return gene;
		}
		
		public String getGenomicRegion(){
			return genomicRegion;
		}
		public String getFrameShift(){
			return frameShift;
		}
		
		/**
		 * @return the coordinate
		 */
		public ArrayList<Integer> getCoordinate() {
			return Bed12Parser.getCoordinate(mappedRegionStartsEnds, isPlusStrand);
		}
		private ScoredPosition(String s){ 
			int i=0;
			String[] token = s.split("\t");
			contig = token[i++];
			position = Integer.parseInt(token[i++]);
			String[] mrStarts = token[i++].split(",");
			String[] mrEnds = token[i++].split(",");
			mappedRegionStartsEnds = new ArrayList<ArrayList<Integer>>();
			for(int j=0;j<mrStarts.length;j++){
				ArrayList<Integer> exom = new ArrayList<Integer>();
				exom.add(Integer.parseInt(mrStarts[j]));
				exom.add(Integer.parseInt(mrEnds[j]));
				mappedRegionStartsEnds.add(exom);
			}
			isPlusStrand = token[i++].equals("+");
		
			
			startScore = Double.parseDouble(token[i++]);
			//stopScore = Double.parseDouble(token[4]);
			seq = token[i++];
			genomicRegion = token[i++];
			frameShift = token[i++];
			if(!token[i].equals("_")){				
				isAnnotated = token[i++].equals("T");
				StringBuffer gsb = new StringBuffer();
				for(;i<token.length;i++){
					gsb.append(token[i]);
					gsb.append('\t');
				}
				gene = new AnnotatedGene(gsb.toString());
			}			
		}
		
		
		
		
		public ScoredPosition(String contig, int position, ArrayList<Integer> coordinate, boolean isPlusStrand, double startScore, String seq, AnnotatedGene gene, boolean isAnnotated, String genomicRegion, String frameShift){
			this.contig = contig;
			this.position = position;
			this.isPlusStrand = isPlusStrand;
			this.startScore = startScore;
			this.mappedRegionStartsEnds = Bed12Parser.getCoordinateStartsEnds(coordinate, isPlusStrand);
			//this.stopScore = stopScore;
			this.seq = seq;
			this.genomicRegion = genomicRegion;
			this.frameShift = frameShift;
			if(gene != null){
				this.isAnnotated = isAnnotated;
				this.gene = gene;			
			}
			
		}
			
		public int getLength(){
			return getCoordinate().size();
		}
		
		public int compareTo(ScoredPosition o) {
			return new Integer(this.getPosition()).compareTo(new Integer(o.getPosition()));
		}
		

		public String toMappedRegionString(StringBuffer sb){
			//StringBuilder sb = new StringBuilder();
			//ArrayList<ArrayList<Integer>> mappedRegionStartsEnds = Bed12Parser.getCoordinateStartsEnds(coordinate, isPlusStrand);
			for(int j=0;j<mappedRegionStartsEnds.size()-1;j++){
				sb.append(mappedRegionStartsEnds.get(j).get(0));
				sb.append(',');
			}
			sb.append(mappedRegionStartsEnds.get(mappedRegionStartsEnds.size()-1).get(0));
			sb.append('\t');
			for(int j=0;j<mappedRegionStartsEnds.size()-1;j++){
				sb.append(mappedRegionStartsEnds.get(j).get(1));
				sb.append(',');
			}
			sb.append(mappedRegionStartsEnds.get(mappedRegionStartsEnds.size()-1).get(1));
			//sb.append('\t');
			return sb.toString();
		}
	
		@Override
		public String toString(){
			StringBuffer sb = new StringBuffer();
			sb.append(contig); sb.append('\t');
			sb.append(position); sb.append('\t');
			toMappedRegionString(sb);
			sb.append('\t');
			sb.append((isPlusStrand? '+': '-')); sb.append('\t');
			sb.append(startScore); sb.append('\t');
			//sb.append(stopScore); sb.append('\t');
			sb.append(seq);sb.append('\t');
			sb.append(genomicRegion);sb.append('\t');
			sb.append(frameShift);sb.append('\t');
			if(gene !=null){	
				sb.append(isAnnotated()? 'T' : 'F');
				sb.append('\t');
				sb.append(gene);
			}else{ 
				sb.append('_');
				sb.append('\t');
				sb.append(AnnotatedGene.getEmptyString());
			}
			return sb.toString(); 
		}		
	
		
	}
	
	private HashMap<String, ArrayList<ScoredPosition>> positionMap;
	
	public ScoringOutputParser(String inFile){
		read(inFile);
	}
	
	public ScoringOutputParser(String inFile, String contig) {
		read(inFile, contig);
	}
	
	public ScoringOutputParser(){
		
	}

	private void read(String outFile){
		read(outFile, null);
	}
	
	private void read(String outFile, String contig){
		positionMap = new HashMap<String, ArrayList<ScoredPosition>>();
		try {
			BufferedLineReader in  = new BufferedLineReader((outFile));
			String s;
			while((s=in.readLine())!=null){
				ScoredPosition position = new ScoredPosition(s);	
				if(contig != null && !position.getContig().equals(contig)) continue;
				if(!positionMap.containsKey(position.getContig()))
					positionMap.put(position.getContig(), new ArrayList<ScoredPosition>());
				positionMap.get(position.getContig()).add(position);
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public Set<String> getContigs(){
		return positionMap.keySet();
	}
	
	public ArrayList<ScoredPosition> getPositions(String contig){
		if(positionMap.containsKey(contig)) return positionMap.get(contig);
		return new ArrayList<ScoredPosition>();
	}
	
	public ArrayList<ScoredPosition> getPositions(){
		ArrayList<ScoredPosition> positions = new ArrayList<ScoredPosition>();
		for(String key : positionMap.keySet())
			positions.addAll(positionMap.get(key));
		return positions;
	}
	
	/*public static ArrayList<Integer> getIntersectedPositions(ScoringOutputParser a, ScoringOutputParser b, String contig){
		HashSet<Integer> positions = new HashSet<Integer>();
		for(ScoredPosition position : a.positionMap.get(contig)){
			positions.add(position.getPosition());
		}
		for(ScoredPosition position : b.positionMap.get(contig)){
			positions.remove(position.getPosition());			
		}
		ArrayList<Integer> positionList = new ArrayList<Integer>(positions);
		Collections.sort(positionList);
		return positionList;
	}
	*/
	
	public static ScoredPosition getScoredPosition(String contig, int position, ArrayList<Integer> coordinate, boolean isPlusStrand, double startScore, 
			String seq, AnnotatedGene gene, boolean isAnnotated, String genomicRegion, String frameShift){
		return  new ScoringOutputParser().new ScoredPosition(contig, position, coordinate, isPlusStrand, startScore, seq, gene, isAnnotated, genomicRegion, frameShift);			
	}
	
	public static ArrayList<ScoredPosition> getUnionPositions(ScoringOutputParser[] parsers,
			AnnotationFileParser annotationParser, String contig, double scoreThreshold){
		//HashMap<Integer, ArrayList<ScoredPosition>> positions = new HashMap<Integer, ArrayList<ScoredPosition>>();
		
		ArrayList<ScoredPosition> positionList = new ArrayList<ScoredPosition>();
		HashMap<ScoredPosition, HashSet<Integer>> positions = new HashMap<ScoredPosition, HashSet<Integer>>();
		
		for(int i=0; i<parsers.length;i++){
			for(ScoredPosition position : parsers[i].getPositions(contig)){
				if(position.getScore() < scoreThreshold) continue;
				if( position.getCoordinate().isEmpty()){
					System.out.println("Empty : " + position );
					continue;
				}
				if(!positions.containsKey(position)) positions.put(position, new HashSet<Integer>());
				positions.get(position).add(i); // rpf index
			}
		}
		
		for(ScoredPosition position : positions.keySet()){
			position.setRPFRNAIndices(positions.get(position)); // set rpf index		
		}	
		positionList.addAll(positions.keySet());
	
		Collections.sort(positionList);
		System.out.println("Union positions done.. " + positionList.size() + " positions");
		return positionList;
	}	
	
	static public ArrayList<ScoredPosition> getWindowFilteredPositionis(ArrayList<ScoredPosition> positions, int window){
		Collections.sort(positions);
		ArrayList<ScoredPosition> ret = new ArrayList<ScoringOutputParser.ScoredPosition>();
		for(int i = 0; i < positions.size(); i++) {
	      int rank = 1;
	      
	      ScoredPosition thisPosition = positions.get(i);
	      
	      // move left
	      int prevIndex = i-1;
	      while(prevIndex >= 0) {
	    	ScoredPosition prevPosition = positions.get(prevIndex);
	        if(thisPosition.getPosition() - prevPosition.getPosition() > window)    break;
	        if(prevPosition.getScore() > thisPosition.getScore()) rank++;
	        prevIndex--;
	      }

	      // move right
	      int nextIndex = i+1;
	      while(nextIndex < positions.size()) {
	    	ScoredPosition nextPosition = positions.get(nextIndex);
	        if(nextPosition.getPosition() - thisPosition.getPosition() > window)    break;
	        if(nextPosition.getScore() > thisPosition.getScore()) rank++;
	        nextIndex++;
	      }
	    
	      if(rank <= 1) ret.add(thisPosition);
		}
		return ret;
	}
	
	static public ArrayList<ScoredPosition> getQuantityFilteredPositions(ArrayList<ScoredPosition> positions, Quantifier quantifier, double RPKM){
		Quantifier[] quantifiers = new Quantifier[1];
		quantifiers[0] = quantifier;
		return getQuantityFilteredPositions(positions, quantifiers, RPKM);
	}
	
	static public ArrayList<ScoredPosition> getQuantityFilteredPositions(ArrayList<ScoredPosition> positions, Quantifier[] quantifiers, double RPKM){
		ArrayList<ScoredPosition> ret = new ArrayList<ScoredPosition>();
		//int i = 0;
		for(ScoredPosition position : positions){
			//System.out.print(i++);
			boolean isAbundant = false;
			for(Quantifier q : quantifiers){
				if(q.isAbundant(position, RPKM)){
					isAbundant = true;
					break;
				}
			}
			//System.out.println(" " + isAbundant);
			if(isAbundant)				
				ret.add(position);
		}
		return ret;
	}
	
	public static void main(String[] args) throws FileNotFoundException{
		ScoringOutputParser test = new ScoringOutputParser("/media/kyowon/Data1/RPF_Project/samples/sample5/bed/RPF-C_1-uncollapsed.bed.score0.3.tsv");
		//PrintStream out = new PrintStream("/media/kyowon/Data1/RPF_Project/samples/sample5/results/out_0.3_smORFcopy.csv");
		HashSet<ScoredPosition> set = new HashSet<ScoringOutputParser.ScoredPosition>(test.getPositions());
		System.out.println(set.size());
		System.out.println(test.getPositions().get(0).equals(test.getPositions().get(1)));
		//out.close();
	}
}
