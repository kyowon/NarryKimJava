package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import parser.AnnotationFileParser.AnnotatedGene;
import net.sf.samtools.util.BufferedLineReader;

public class ScoringOutputParser {
	//private HashMap<String, HashMap<Integer, >>
	
	public class ScoredPosition implements Comparable<ScoredPosition>{
		private String contig;		
		private boolean isPlusStrand;
		private int position;
		private double score;
		private String codon;
		private boolean isAnnotated = false;
		private AnnotatedGene gene = null;
		private String genomicRegion;
		private String frameShift;
		
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
			return score;
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
			score = Double.parseDouble(token[3]);
			codon = token[4];
			genomicRegion = token[5];
			frameShift = token[6];
			if(!token[8].equals("_")){				
				isAnnotated = token[7].equals("T");
				StringBuffer gsb = new StringBuffer();
				for(int i=8;i<token.length;i++){
					gsb.append(token[i]);
					gsb.append('\t');
				}
				gene = new AnnotatedGene(gsb.toString());
			}
			
		}
		
		private ScoredPosition(String contig, int position, boolean isPlusStrand, double score, String codon, AnnotatedGene gene, boolean isAnnotated, String genomicRegion, String frameShift){
			this.contig = contig;
			this.position = position;
			this.isPlusStrand = isPlusStrand;
			this.score = score;
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
			sb.append(score); sb.append('\t');
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
	
	public static ScoredPosition getScoredPosition(String contig, int position, boolean isPlusStrand, double score, String codon, AnnotatedGene gene, boolean isAnnotated, String genomicRegion, String frameShift){
		return  new ScoringOutputParser().new ScoredPosition(contig, position, isPlusStrand, score, codon, gene, isAnnotated, genomicRegion, frameShift);			
	}
	
	public static ArrayList<ScoredPosition> getUnionPositions(ScoringOutputParser[] parsers, String contig, double scoreThreshold){
		HashSet<ScoredPosition> positions = new HashSet<ScoredPosition>();
		for(int i=0; i<parsers.length;i++){
			positions.addAll(parsers[i].getPositions(contig));
		}	
		ArrayList<ScoredPosition> positionList = new ArrayList<ScoredPosition>();
		for(ScoredPosition p : positions){
			if(scoreThreshold > 0 && p.getScore() > scoreThreshold) positionList.add(p);
			else if(scoreThreshold <= 0 && p.getScore() < -scoreThreshold) positionList.add(p);
		}		
		Collections.sort(positionList);
		return positionList;
	}
	
}
