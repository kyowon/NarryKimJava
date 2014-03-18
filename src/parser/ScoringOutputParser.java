package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import parser.AnnotationFileParser.AnnotatedGene;
import rpf.Quantifier;
import net.sf.samtools.util.BufferedLineReader;

public class ScoringOutputParser {
	//private HashMap<String, HashMap<Integer, >>
	
	public class ScoredPosition implements Comparable<ScoredPosition>{
		private String contig;		
		private boolean isPlusStrand;
		private int position;
		private double startScore;
		//private double stopScore;
		private String codon;
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
		public HashSet<Integer> getHarrIndices() {
			return harrIndices;
		}
		public void setHarrIndices(HashSet<Integer> conditionIndices) {
			this.harrIndices = conditionIndices;
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
		public String getCodon() {
			return codon;
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
		private ScoredPosition(String s){ 
			String[] token = s.split("\t");
			contig = token[0];
			position = Integer.parseInt(token[1]);
			isPlusStrand = token[2].equals("+");
			startScore = Double.parseDouble(token[3]);
			//stopScore = Double.parseDouble(token[4]);
			codon = token[4];
			genomicRegion = token[5];
			frameShift = token[6];
			if(!token[7].equals("_")){				
				isAnnotated = token[7].equals("T");
				StringBuffer gsb = new StringBuffer();
				for(int i=8;i<token.length;i++){
					gsb.append(token[i]);
					gsb.append('\t');
				}
				gene = new AnnotatedGene(gsb.toString());
			}
			
		}
		
		private ScoredPosition(String contig, int position, boolean isPlusStrand, double startScore, String codon, AnnotatedGene gene, boolean isAnnotated, String genomicRegion, String frameShift){
			this.contig = contig;
			this.position = position;
			this.isPlusStrand = isPlusStrand;
			this.startScore = startScore;
			//this.stopScore = stopScore;
			this.codon = codon;
			this.genomicRegion = genomicRegion;
			this.frameShift = frameShift;
			if(gene != null){
				this.isAnnotated = isAnnotated;
				this.gene = gene;			
			}
			
		}
				
		public int compareTo(ScoredPosition o) {
			return new Integer(this.getPosition()).compareTo(new Integer(o.getPosition()));
		}
		
		@Override
		public boolean equals(Object o){
			if(o instanceof ScoredPosition){
				ScoredPosition other = (ScoredPosition)o;
				return this.position == other.position && this.isPlusStrand == other.isPlusStrand && this.contig.equals(other.contig);
			}
			return false;
		}
		
		@Override
		public int hashCode(){
			return new Integer(position).hashCode() * new Boolean(isPlusStrand).hashCode() * contig.hashCode();
		}
		
		@Override
		public String toString(){
			StringBuffer sb = new StringBuffer();
			sb.append(contig); sb.append('\t');
			sb.append(position); sb.append('\t');
			sb.append((isPlusStrand? '+': '-')); sb.append('\t');
			sb.append(startScore); sb.append('\t');
			//sb.append(stopScore); sb.append('\t');
			sb.append(codon);sb.append('\t');
			sb.append(genomicRegion);sb.append('\t');
			sb.append(frameShift);sb.append('\t');
			if(gene !=null){	
				sb.append(isAnnotated()? 'T' : 'F');
				sb.append('\t');
				sb.append(gene);
			}else{ 
				sb.append('_');
				sb.append('\t');
				sb.append(AnnotatedGene.getEmptyGeneString());
			}
			return sb.toString(); 
		}		
	}
	
	private HashMap<String, ArrayList<ScoredPosition>> positionMap;
	
	public ScoringOutputParser(String outFile){
		read(outFile);
	}
	
	public ScoringOutputParser() {
	}

	private void read(String outFile){
		positionMap = new HashMap<String, ArrayList<ScoredPosition>>();
		try {
			BufferedLineReader in  = new BufferedLineReader(new FileInputStream(outFile));
			String s;
			while((s=in.readLine())!=null){
				ScoredPosition position = new ScoredPosition(s);				
				if(!positionMap.containsKey(position.getContig()))
					positionMap.put(position.getContig(), new ArrayList<ScoredPosition>());
				positionMap.get(position.getContig()).add(position);
			}
			in.close();
		} catch (FileNotFoundException e) {
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
	
	public static ScoredPosition getScoredPosition(String contig, int position, boolean isPlusStrand, double startScore, String codon, AnnotatedGene gene, boolean isAnnotated, String genomicRegion, String frameShift){
		return  new ScoringOutputParser().new ScoredPosition(contig, position, isPlusStrand, startScore, codon, gene, isAnnotated, genomicRegion, frameShift);			
	}
	
	public static ArrayList<ScoredPosition> getUnionPositions(ScoringOutputParser[] parsers, Quantifier[] rpfQuantifiers, Quantifier[] rnaQuantifiers, AnnotationFileParser annotationParser, String contig, double scoreThreshold, double rpfRPKMThreshold, double rnaRPKMThreshold){
		//HashMap<Integer, ArrayList<ScoredPosition>> positions = new HashMap<Integer, ArrayList<ScoredPosition>>();
		
		ArrayList<ScoredPosition> positionList = new ArrayList<ScoredPosition>();
		HashMap<ScoredPosition, HashSet<Integer>> positions = new HashMap<ScoredPosition, HashSet<Integer>>();
		
		for(int i=0; i<parsers.length;i++){
			for(ScoredPosition position : parsers[i].getPositions(contig)){
				if(position.getScore() < scoreThreshold) continue;
				if(!positions.containsKey(position)) positions.put(position, new HashSet<Integer>());
				positions.get(position).add(i); // harr index
			}
		}
		
		for(ScoredPosition position : positions.keySet()){
			position.setHarrIndices(positions.get(position)); // set harr index
			
			HashSet<Integer> rpfAbundantIndices = null;
			AnnotatedGene gene = null;//(annotationParser != null? annotationParser.getContainingGene(contig, position.isPlusStrand(), position.getPosition()) : null);
			
			if(rpfQuantifiers!=null){
			//	boolean abundant = false;	
				rpfAbundantIndices = new HashSet<Integer>();
				for(int j=0;j<rpfQuantifiers.length; j++){
					double rpfRPKM = (gene == null ? rpfQuantifiers[j].getPositionRPKM(contig, position.getPosition(), position.isPlusStrand(), 150) : rpfQuantifiers[j].getCDSRPKM(gene));
					if(rpfRPKM > rpfRPKMThreshold){
					//	abundant = true;
						rpfAbundantIndices.add(j);
						//break;
					}
				}
				//if(!abundant) continue;					
			}
			
			HashSet<Integer> rnaAbundantIndices = null;
			if(rnaQuantifiers!=null){
			//	boolean abundant = false;
				rnaAbundantIndices = new HashSet<Integer>();
				for(int j=0;j<rnaQuantifiers.length; j++){
					double rnaRPKM = (gene == null ? rnaQuantifiers[j].getPositionRPKM(contig, position.getPosition(), position.isPlusStrand(), 150) : rnaQuantifiers[j].getCDSRPKM(gene));
					if(rnaRPKM > rnaRPKMThreshold){
					//	abundant = true;
						rnaAbundantIndices.add(j);
						//break;
					}
				}
				//if(!abundant) continue;					
			}
			
			HashSet<Integer> indices = new HashSet<Integer>();
			
			if(rpfAbundantIndices == null && rnaAbundantIndices == null){
				for(int j=0; j<parsers.length;j++) indices.add(j);
			}else if(rpfAbundantIndices != null && rnaAbundantIndices == null){
				indices = rpfAbundantIndices;
			}else if(rpfAbundantIndices == null && rnaAbundantIndices != null){
				indices = rnaAbundantIndices;
			}else{
				for(int j : rpfAbundantIndices){
					if(rnaAbundantIndices.contains(j)) indices.add(j);
				}
			}
			
			position.setRPFRNAIndices(indices);
			//positions.add(position);
			
			//if(!positions.containsKey(position.getPosition())) positions.put(position.getPosition(), new ArrayList<ScoredPosition>());
		//	positions.get(position.getPosition()).add(position);
				
		}	
		positionList.addAll(positions.keySet());
	//	ArrayList<ScoredPosition> positionList = new ArrayList<ScoredPosition>();
		
	//	for(int position : positions.keySet()){
	//		positionList.addAll(positions.get(position));
	//		ArrayList<ScoredPosition> ps = positions.get(position);
			//ScoredPosition tp = null;
		//	for(ScoredPosition p : ps){
				//if(scoreThreshold>0){
		//			if(p.getScore() >= scoreThreshold){
			//			positionList.add(p);
		//				break;
		//			}
				//}
				/*else{
					tp = p;
					if(p.getScore() >= -scoreThreshold){
						tp = null;
						break;
					}
				}*/
			//	if(tp!=null) positionList.add(tp);
		//	}			
		//}		
		Collections.sort(positionList);
		System.out.println("Union positions done.. " + positionList.size() + " positions");
		return positionList;
	}	
}
