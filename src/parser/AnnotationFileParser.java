package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import net.sf.picard.util.IntervalTree;
import net.sf.picard.util.IntervalTree.Node;
import net.sf.samtools.util.BufferedLineReader;

public class AnnotationFileParser {
	public class txStartComparator implements Comparator<AnnotatedGene>{
		public int compare(AnnotatedGene o1, AnnotatedGene o2) {
			int i= new Integer(o1.txStart).compareTo(new Integer(o2.txStart));
			if(i!=0) return i;
			return new Integer(o1.txEnd).compareTo(new Integer(o2.txEnd));
		}
	}
	public class txEndComparator implements Comparator<AnnotatedGene>{
		public int compare(AnnotatedGene o1, AnnotatedGene o2) {
			int i= new Integer(o1.txEnd).compareTo(new Integer(o2.txEnd));
			if(i!=0) return i;
			return new Integer(o1.txStart).compareTo(new Integer(o2.txStart));
		}
	}
	
	public static class AnnotatedGene{
		private String geneName;		
		private String accession;
		private String contig;
		private boolean isPlusStrand;
		private int txStart;
		private int txEnd;
		private int cdsStart;
		private int cdsEnd;
		private int exonCount;		
		private int[] exonStarts;
		private int[] exonEnds;
		private int[][] introns; // intron inclusive, increasing
		private HashMap<Integer, Integer> fivepSplices = null;
		private HashMap<Integer, Integer> threepSplices = null;
		
		public AnnotatedGene(String s){
			String[] token = s.split("\t");
			geneName = token[0];
			accession = token[1];
			contig = token[2];
			isPlusStrand = token[3].equals("+");
			txStart = Integer.parseInt(token[4]);
			txEnd = Integer.parseInt(token[5]);
			cdsStart = Integer.parseInt(token[6]);
			cdsEnd = Integer.parseInt(token[7]);
			exonCount = Integer.parseInt(token[8]);
			exonStarts = new int[exonCount];
			int i=0;
			for(String t : token[9].split(",")){
				if(t.isEmpty()) continue;
				exonStarts[i++] = Integer.parseInt(t);
			}
			i=0;
			exonEnds = new int[exonCount];
			for(String t : token[10].split(",")){
				if(t.isEmpty()) continue;
				exonEnds[i++] = Integer.parseInt(t);
			}
			introns = new int[getExonCount()-1][2];
			if(introns.length > 0){
				fivepSplices = new HashMap<Integer, Integer>();
				threepSplices = new HashMap<Integer, Integer>();
			}
			for(int j=1;j<exonCount;j++){
				introns[j-1][0] = getExonEnds()[j-1];
				introns[j-1][1] = getExonStarts()[j]-1;
				
				if(isPlusStrand){
					threepSplices.put(introns[j-1][0]-1, introns[j-1][1] + 1);
					fivepSplices.put(introns[j-1][1] + 1, introns[j-1][0]-1);
				}else{
					fivepSplices.put(introns[j-1][0]-1, introns[j-1][1] + 1);
					threepSplices.put(introns[j-1][1] + 1, introns[j-1][0]-1);
				}
				
			}
		}	
		
		
		/*// construct sudo annotated gene
		private AnnotatedGene(ScoredPosition position){ 
			genomeBrowserGeneName = position.getGBGeneName() != null? position.getGBGeneName() : "Putative";
			geneName = position.getGeneName() != null? position.getGeneName() : ((Integer)position.getPosition()).toString();
			contig = position.getContig();
			isPlusStrand = position.isPlusStrand();
			txStart = position.getGeneName() != null ? position.getTxStart() : position.getPosition() - 100;
			txEnd = position.getGeneName() != null? position.getTxEnd() : position.getPosition() + 100;
			cdsStart = isPlusStrand? position.getPosition() : txStart;
			cdsEnd = isPlusStrand? txEnd : position.getPosition() + 1;
			exonCount = 1;
			exonStarts = new int[exonCount];
			exonEnds = new int[exonCount];			
		}*/
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((accession == null) ? 0 : accession.hashCode());
			return result;
		}


		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			AnnotatedGene other = (AnnotatedGene) obj;
			if (accession == null) {
				if (other.accession != null)
					return false;
			} else if (!accession.equals(other.accession))
				return false;
			return true;
		}

		public Integer getDistFromClosestExon(int position){ // 0 if within an exon strand matters 
			if(exonCount==0) return null;		
			for(int j=0; j<exonCount;j++){
				if(position < exonStarts[j]){
					int d1 = exonStarts[j] - position; 
					int d2 = j>0? position - exonEnds[j-1] + 1: d1 + 1;
					if(isPlusStrand) d1 = -d1;
					else d2 = -d2;
					return Math.abs(d1) < Math.abs(d2)? d1 : d2;
				}
				
				if(position < exonEnds[j]){
					return 0;
				}
			}
			int d2 = position - exonEnds[exonCount-1] + 1;
			if(!isPlusStrand) d2 = -d2;
			//else d1 = -d1;
			return d2;
			//int d1 = position - exonStarts[exonCount-1];
			//int d2 = exonEnds[exonCount-1] - position - 1;
			//if(isPlusStrand) d2 = -d2;
			//else d1 = -d1;
			//return Math.abs(d1) < Math.abs(d2)? d1 : d2;
		}

		private AnnotatedGene(int txStart, int txEnd){ // constructor for comparison
			this.txStart = txStart;
			this.txEnd = txEnd;
		}
		
		@Override
		public String toString(){
			StringBuffer sb = new StringBuffer();			
			sb.append(geneName); sb.append('\t');
			sb.append(accession); sb.append('\t');
			sb.append(contig); sb.append('\t');
			sb.append(isPlusStrand? '+' : '-'); sb.append('\t');
			sb.append(txStart); sb.append('\t');
			sb.append(txEnd); sb.append('\t');
			sb.append(cdsStart); sb.append('\t');
			sb.append(cdsEnd); sb.append('\t');
			sb.append(exonCount); sb.append('\t');
			for(int t : exonStarts){
				sb.append(t);sb.append(',');
			}
			sb.append('\t');
			for(int t : exonEnds){
				sb.append(t);sb.append(',');
			}
			return sb.toString();
		}
		
		public String toSimpleString(){
			StringBuffer sb = new StringBuffer();			
			sb.append(geneName); sb.append('\t');
			sb.append(accession); 
			
			return sb.toString();
		}
		
		public String toDetailedInfoString(){
			StringBuffer sb = new StringBuffer();			
			//sb.append(genomeBrowserGeneName); sb.append('\t');
			//sb.append(geneName); sb.append('\t');
			//sb.append(contig); sb.append('\t');
			sb.append(isPlusStrand? '+' : '-'); sb.append('\t');
			sb.append(txStart); sb.append('\t');
			sb.append(txEnd); sb.append('\t');
			sb.append(cdsStart); sb.append('\t');
			sb.append(cdsEnd); sb.append('\t');
			sb.append(exonCount); sb.append('\t');
			for(int t : exonStarts){
				sb.append(t);sb.append(',');
			}
			sb.append('\t');
			for(int t : exonEnds){
				sb.append(t);sb.append(',');
			}
			return sb.toString();
		}
		
		public static String getEmptyString(){
			StringBuffer sb = new StringBuffer();			
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); 
			return sb.toString();
		}
		
		public static String getEmptySimpleString(){
			StringBuffer sb = new StringBuffer();			
			sb.append('_'); sb.append('\t');
			//sb.append('_'); sb.append('\t');
			sb.append('_'); //sb.append('\t');
			//sb.append('_'); 
			return sb.toString();
		}
		
		public static String getEmptyInfoString(){
			StringBuffer sb = new StringBuffer();			
		//	sb.append('_'); sb.append('\t');
		//	sb.append('_'); sb.append('\t');
		//	sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); sb.append('\t');
			sb.append('_'); 
			return sb.toString();
		}
		
		public static String getHeader(){
			StringBuffer sb = new StringBuffer();
			sb.append("GeneName");sb.append('\t');
			sb.append("Accession");sb.append('\t');
			sb.append("Contig");sb.append('\t');
			sb.append("Strand");sb.append('\t');
			sb.append("txStart");sb.append('\t');
			sb.append("txEnd");sb.append('\t');
			sb.append("cdsStart");sb.append('\t');
			sb.append("cdsEnd");sb.append('\t');
			sb.append("ExonCount");sb.append('\t');
			sb.append("ExonStarts");sb.append('\t');
			sb.append("ExonEnds");		
			return sb.toString();
		}
		
		public static String getSimpleHeader(){
			StringBuffer sb = new StringBuffer();
			sb.append("GeneName");sb.append('\t');
			sb.append("Accession");
			return sb.toString();
		}
		
		public static String getDetailedInfoHeader(){
			StringBuffer sb = new StringBuffer();
			sb.append("Strand");sb.append('\t');
			sb.append("txStart");sb.append('\t');
			sb.append("txEnd");sb.append('\t');
			sb.append("cdsStart");sb.append('\t');
			sb.append("cdsEnd");sb.append('\t');
			sb.append("ExonCount");sb.append('\t');
			sb.append("ExonStarts");sb.append('\t');
			sb.append("ExonEnds");		
			return sb.toString();
		}
		
		public String getGeneName() {
			return geneName;
		}
		public String getAccession() {
			return accession;
		}
		public String getContig() {
			return contig;
		}
		public boolean isPlusStrand() {
			return isPlusStrand;
		}
		public int getTxStart() {
			return txStart;
		}
		public int getTxEnd() {
			return txEnd;
		}
		public int getCdsStart() {
			return cdsStart;
		}
		public int getCdsEnd() {
			return cdsEnd;
		}
		public int getExonCount() {
			return exonCount;
		}
		public int[] getExonStarts() {
			return exonStarts;
		}
		public int[] getExonEnds() {
			return exonEnds;
		}
		public int[][] getIntrons(){
			return introns;
		}
		public boolean isAnnotated(int position){
			return isPlusStrand? position == this.cdsStart : position == this.cdsEnd - 1;
		}
		
		
		private boolean matching(boolean isPlusStrand, ArrayList<Integer> coordinates, ArrayList<ArrayList<Integer>> splices){ // TODO
			if(isPlusStrand != this.isPlusStrand) return false;
			if(introns == null || introns.length == 0){
				if(splices == null || splices.isEmpty()) return true;
			}
			
			//System.out.println(splices); 
			ArrayList<ArrayList<Integer>> intronSet = new ArrayList<ArrayList<Integer>>();
			if(isPlusStrand){
				for(int i=0;i<introns.length;i++){
					int[] intron = introns[i]; // check order.. TODO check inclusive / exclusive
					if(coordinates.get(0) > intron[1] || coordinates.get(coordinates.size()-1) < intron[0]) continue;  
					
					ArrayList<Integer> it = new ArrayList<Integer>();
					it.add(intron[0]);
					it.add(intron[1]);
					intronSet.add(it);
					//if(!splices.contains(il)) return false;
				}
			}else{
				for(int i=introns.length-1;i>=0;i--){
					int[] intron = introns[i]; // check order.. TODO check inclusive / exclusive
					if(coordinates.get(coordinates.size()-1) > intron[1] || coordinates.get(0) < intron[0]) continue;

					ArrayList<Integer> it = new ArrayList<Integer>();
					it.add(intron[0]);
					it.add(intron[1]);
					intronSet.add(it);
					//if(!splices.contains(il)) return false;
				}
			}
			
			ArrayList<ArrayList<Integer>> splicesInGene = new ArrayList<ArrayList<Integer>>();
			for(ArrayList<Integer> splice : splices){
				if(splice.get(0) >= this.txEnd || splice.get(1)<this.txStart) continue;
				splicesInGene.add(splice);	
			}
			
			return intronSet.equals(splicesInGene);
		}
		
		
		
		
		// start : inclusive
		public ArrayList<Integer> getLiftOverPositions(int position, int leftWindowSize, int rightWindowSize, boolean stopAtCDSEnd){
			ArrayList<Integer> ret = new ArrayList<Integer>();		
				// 5p direction
			int cp = position;
			int next = cp;
			while(leftWindowSize-- > 0){
				if(fivepSplices != null && fivepSplices.containsKey(cp)){
					next = fivepSplices.get(cp);
				}else next += isPlusStrand? -1 : 1;
				cp = next;
				ret.add(0, cp);	
			}
				// 3p direction
			cp = position;
			next = cp;
			while(rightWindowSize-- > 0){
				if(threepSplices != null && threepSplices.containsKey(cp)){
					next = threepSplices.get(cp);
				}else next += isPlusStrand? 1 : -1;
				ret.add(cp);
				cp = next;
				if(stopAtCDSEnd){
					if(isPlusStrand){
						if(cp >= this.cdsEnd -1) break;
					}else{
						if(cp <= this.cdsStart) break;
					}
				}
			}	
					
			return ret;
		}
		
		
		public ArrayList<Integer> getLiftOverPositionsTillNextStopCodon(int start, int minLength, int maxLength, ZeroBasedFastaParser fastaParser){
			ArrayList<Integer> positions = getLiftOverPositions(start, 0, maxLength, false);
			ArrayList<Integer> ret = new ArrayList<Integer>();
		//	boolean stop = false;
		//	int len = 0;
			
			for(int j=0;j<positions.size()-2;j+=3){				
				String tc = fastaParser.getSequence(contig, positions.subList(j, j+3));
				if(!isPlusStrand) tc = ZeroBasedFastaParser.getComplementarySequence(tc, false);
				tc = tc.toUpperCase();
				ret.addAll(positions.subList(j, j+3));
				if(ret.size() < minLength) continue;
				if(ret.size() > maxLength) break;
				if(tc.equals("TAG") || tc.equals("TAA") || tc.equals("TGA")){ 					
					break;
				}
			}	
			
			return ret;
		}
		
	}
	
	public Set<String> getContigs(){
		return annotatedGeneSetMap.keySet();
	}
	
	private HashMap<String, ArrayList<AnnotatedGene>> annotatedGeneSetMap;
	private HashMap<String, AnnotatedGene> accessionMap;
	private HashMap<String, IntervalTree<ArrayList<Integer>>> txIntervalMap;
	private HashMap<String, HashMap<Boolean, HashMap<Integer, ArrayList<AnnotatedGene>>>> annotatedGenePositionMap;
	
//	public static AnnotatedGene getSudoAnnotatedGene(ScoredPosition position){
//		return new AnnotationFileParser().new AnnotatedGene(position);
//	}
	
	
	public AnnotationFileParser(String annotationFile){
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(annotationFile));
			annotatedGeneSetMap = new HashMap<String, ArrayList<AnnotatedGene>>();
			accessionMap = new HashMap<String, AnnotationFileParser.AnnotatedGene>();
			annotatedGenePositionMap = new HashMap<String, HashMap<Boolean,HashMap<Integer,ArrayList<AnnotatedGene>>>>();
			String s;
			while((s=in.readLine())!=null){
				AnnotatedGene gene = new AnnotatedGene(s);
				accessionMap.put(gene.accession, gene);
				String contig = gene.getContig();
				if(!annotatedGeneSetMap.containsKey(contig)){
					annotatedGeneSetMap.put(contig, new ArrayList<AnnotatedGene>());
					annotatedGenePositionMap.put(contig, new HashMap<Boolean, HashMap<Integer, ArrayList<AnnotatedGene>>>());
					annotatedGenePositionMap.get(contig).put(true, new HashMap<Integer, ArrayList<AnnotatedGene>>());
					annotatedGenePositionMap.get(contig).put(false, new HashMap<Integer, ArrayList<AnnotatedGene>>());
				}
				annotatedGeneSetMap.get(contig).add(gene);	
				int position = gene.isPlusStrand()? gene.getCdsStart() : gene.getCdsEnd();
				if(!annotatedGenePositionMap.get(contig).get(gene.isPlusStrand()).containsKey(position))
					annotatedGenePositionMap.get(contig).get(gene.isPlusStrand()).put(position, new ArrayList<AnnotatedGene>());
				annotatedGenePositionMap.get(contig).get(gene.isPlusStrand()).get(position).add(gene);
			}
			for(String contig : annotatedGeneSetMap.keySet()){
				Collections.sort(annotatedGeneSetMap.get(contig), new txStartComparator());
			}
			in.close();
			
			txIntervalMap = new HashMap<String, IntervalTree<ArrayList<Integer>>>();
			for(String contig : annotatedGeneSetMap.keySet()){
				txIntervalMap.put(contig, new IntervalTree<ArrayList<Integer>>());
				IntervalTree<ArrayList<Integer>> tree = txIntervalMap.get(contig);
				
				for(int i=0;i<annotatedGeneSetMap.get(contig).size();i++){
					AnnotatedGene gene = annotatedGeneSetMap.get(contig).get(i);
					int start = gene.txStart;
					int end = gene.txEnd-1;
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
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	/**
	    * Get the annotated gene
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the coding start position ( see below for the examples)
	    * @return the annotated genes with the position as the CDS start or as the CDS end.
	    * For example, 
	    * 
	    * Gene1 : chr1   +    cds start 0 cds end 100    will be returned for getAnnotatedGene(chr1, true, 0)
	    * Gene2:  chr3   -    cds start 10  cds end 123  will be returned for getAnnotatedGene(chr3, false, 122), not for getAnnotatedGene(chr3, false, 123) 
	    * 
	    * This is because cds end is exclusive in the ref flat file
	    */
	/*public ArrayList<AnnotatedGene> getAnnotatedGenes(String contig, boolean isPlusStrand, int position){
		HashMap<Boolean, HashMap<Integer, ArrayList<AnnotatedGene>>> sub = annotatedGenePositionMap.get(contig);
		if(sub == null) return null;
		HashMap<Integer, ArrayList<AnnotatedGene>> ssub = sub.get(isPlusStrand);
		if(ssub == null) return null;
		return ssub.get(position + (isPlusStrand? 0 : 1));
	}*/
	
	/**
	    * Get the annotated gene
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the coding start position ( see below for the examples)
	    * @param coordianteSet the coordinates (considering read splices)
	    * @return the annotated genes with the position as the CDS start or as the CDS end.
	    * For example, 
	    * 
	    * Gene1 : chr1   +    cds start 0 cds end 100    will be returned for getAnnotatedGene(chr1, true, 0)
	    * Gene2:  chr3   -    cds start 10  cds end 123  will be returned for getAnnotatedGene(chr3, false, 122), not for getAnnotatedGene(chr3, false, 123) 
	    * 
	    * This is because cds end is exclusive in the ref flat file
	    */
	/*public ArrayList<AnnotatedGene> getAnnotatedGenes(String contig, boolean isPlusStrand, int position, HashSet<Integer> coordinateSet){
		ArrayList<AnnotatedGene> genes = getAnnotatedGenes(contig, isPlusStrand, position);
		if(genes == null) return null;
		ArrayList<AnnotatedGene> ret = new ArrayList<AnnotationFileParser.AnnotatedGene>();
		for(AnnotatedGene gene : genes){
			if(gene.containing(coordinateSet)) ret.add(gene);
		}
		return ret.isEmpty()? null : ret;
	}*/
	
	
	/**
	    * Get the matching genes
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the genomic position
	    * @param coordinate the coordiate considering read splices
	    * @return the matching genes
	    */
	public ArrayList<AnnotatedGene> getMatchingGenes(String contig, boolean isPlusStrand, int position, ArrayList<Integer> coordinate){
		ArrayList<AnnotatedGene> genes = getContainingGenes(contig, isPlusStrand, position);
		return getMatchingGenes(genes, isPlusStrand, coordinate);
	}
	
	public ArrayList<AnnotatedGene> getMatchingGenes(ArrayList<AnnotatedGene> containingGenes, boolean isPlusStrand, ArrayList<Integer> coordinate){
		if(containingGenes == null) return null;
		ArrayList<AnnotatedGene> ret = new ArrayList<AnnotationFileParser.AnnotatedGene>();
		ArrayList<ArrayList<Integer>> splices = Bed12Parser.getSplices(coordinate);
		
		for(AnnotatedGene gene: containingGenes){
			if(gene.matching(isPlusStrand, coordinate, splices))
				ret.add(gene);
		}
		
		return ret.isEmpty()? null : ret;
	}
	
	/**
	    * Get the containing genes
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the genomic position
	    * @return the containing genes
	    */
	public ArrayList<AnnotatedGene> getContainingGenes(String contig, boolean isPlusStrand, int position){
		ArrayList<AnnotatedGene> ret = new ArrayList<AnnotationFileParser.AnnotatedGene>();
		ArrayList<AnnotatedGene> genes = annotatedGeneSetMap.get(contig);
		if(genes != null)// ret = null;
		{
			Iterator<Node<ArrayList<Integer>>> it = txIntervalMap.get(contig).overlappers(position, position+1);
			while(it.hasNext()){
				Node<ArrayList<Integer>> is = it.next();
				for(int i : is.getValue()){
					AnnotatedGene gene = genes.get(i);
					if(gene.isPlusStrand() == isPlusStrand)
						ret.add(gene);
				}
			}
		}
		return ret.isEmpty()? null : ret;
	}
	
	/**
	    * Get the genomic region and frameshift information
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the position
	     * @param gene the gene not only containing position but matching the coordinate 
	    * @return ArrayList of String : 0 genomic region, 1 frameshift (0/1/2)
	    */
	public ArrayList<String> getGenomicRegionNameAndFrameShift(String contig, boolean isPlusStrand, int position, AnnotatedGene gene, boolean isMatchingGene){ 
		ArrayList<String> ret = new ArrayList<String>();
		String name;
		String frameShift = "_";
		
		if(gene == null){
			name = "InterGenic";
		}else{
			if(gene.getGeneName().startsWith("LINC")) name = "LINC";
			else if(gene.getAccession().startsWith("NR")) name = "NR";
			else{
				name = "NM";
				if(position >= gene.getCdsStart() && position < gene.getCdsEnd()){
					name += "_ORF";			
				}else{
					if(isPlusStrand && position < gene.getCdsStart()) name += "_5_UTR";
					else if(!isPlusStrand && position >= gene.getCdsEnd()) name += "_5_UTR";
					else name += "_3_UTR";
				}				
			}
			int[][] introns = gene.getIntrons();
			if(introns != null){
				for(int i=0;i<introns.length;i++){//106129273-5
					if(position >=introns[i][0] && position <= introns[i][1]){
						name += "_Intron";
						break;
					}
				}	
			}
			if(!isMatchingGene){				
				if(!name.endsWith("_Intron")) name += "_ISOFORM";
			}
			
			if(name.endsWith("_ORF") || name.endsWith("_UTR")){
				Integer j = 0;
				if(isPlusStrand){
					int i = 0;
					for(;i<gene.getExonCount()-1;i++){
						if(position < gene.getExonEnds()[i]) break;
						j+=gene.getExonEnds()[i] - gene.getExonStarts()[i];
					}
					assert(position >= gene.getExonStarts()[i]);
					j += position - gene.getExonStarts()[i];
					
					i = 0;
					int startPosition = gene.getCdsStart();
					for(;i<gene.getExonCount()-1;i++){
						if(startPosition < gene.getExonEnds()[i]) break;					
						j-=gene.getExonEnds()[i] - gene.getExonStarts()[i];
					}
					assert(startPosition >= gene.getExonStarts()[i]);
					j -= startPosition - gene.getExonStarts()[i];
					j = j%3;
					if(j<0) j+=3;
				}else{				
					int i = gene.getExonCount()-1;
					for(;i>0;i--){
						if(position >= gene.getExonStarts()[i]) break;
						j+=gene.getExonEnds()[i] - gene.getExonStarts()[i];
					}
					assert(position < gene.getExonEnds()[i]);
					j += gene.getExonEnds()[i] - position;
					
					i = gene.getExonCount()-1;
					int startPosition = gene.getCdsEnd() - 1;
					for(;i>0;i--){
						if(startPosition >= gene.getExonStarts()[i]) break;
						j-=gene.getExonEnds()[i] - gene.getExonStarts()[i];
					}
					assert(startPosition < gene.getExonEnds()[i]);
					j -= gene.getExonEnds()[i] - startPosition;
					j = j%3;
					if(j<0) j+=3;
				}	
				frameShift = j.toString();
			}
		}
		
		ret.add(name); ret.add(frameShift);
		return ret;
	}
	
	public AnnotatedGene getGeneByAccession(String accession){
		return accessionMap.get(accession);
	}
	
/*	public ArrayList<ArrayList<Integer>> getLiftOverCDSPositions(AnnotatedGene gene){	
		return getLiftOverPositions(gene, gene.isPlusStrand(), gene.isPlusStrand()? gene.cdsStart : gene.cdsEnd - 1, gene == null? 0 : Integer.MAX_VALUE, true);
	}
	*/
	//private int pstart = -1;
	//private int pminLength = 0;
	private String pcontig = "";
	//private boolean pIsPlusStrand = false;
	//private ArrayList<ArrayList<Integer>> pPositions = null;
	private HashMap<Integer, HashMap<Integer, HashMap<Boolean, ArrayList<ArrayList<Integer>>>>> pmap = null;
	
	// super stupid .. fix later
	/**/
	
	/*
	public ArrayList<ArrayList<Integer>> getLiftOverPositions(String contig, boolean isPlusStrand, int start, int length){
		AnnotatedGene gene = getContainingGene(contig, isPlusStrand, start);
		//if(gene == null){
		//	if(isPlusStrand) gene = getContainingGene(contig, isPlusStrand, start + length - 1);
		//	else gene = getContainingGene(contig, isPlusStrand, start - length + 1);			
		//}
		return getLiftOverPositions(gene, isPlusStrand, start, length, false);
	}*/
	

	public Iterator<AnnotatedGene> getAnnotatedGeneIterator(){
		ArrayList<AnnotatedGene> allGenes = new ArrayList<AnnotatedGene>();
		for(String contig : annotatedGeneSetMap.keySet())
			allGenes.addAll(annotatedGeneSetMap.get(contig));
		return allGenes.iterator();
	}

	public Iterator<AnnotatedGene> getAnnotatedGeneIterator(String contig){
		if(!annotatedGeneSetMap.containsKey(contig)) return new ArrayList<AnnotatedGene>().iterator();
		return annotatedGeneSetMap.get(contig).iterator();
	}
	
	public ArrayList<AnnotatedGene> getAnnotatedGenes(String contig){
		if(!annotatedGeneSetMap.containsKey(contig)) return new ArrayList<AnnotatedGene>();
		return annotatedGeneSetMap.get(contig);
	}
	
	public void toBedFile(String bed){
		try {
			PrintStream out = new PrintStream(bed);
			Iterator<AnnotatedGene> iterator = getAnnotatedGeneIterator();
			while(iterator.hasNext()){
				AnnotatedGene gene = iterator.next();
				StringBuilder sb = new StringBuilder();
				sb.append(gene.getContig());sb.append('\t');
				sb.append(gene.getTxStart());sb.append('\t');
				sb.append(gene.getTxEnd());sb.append('\t');
				sb.append(gene.getAccession());sb.append('\t');
				sb.append('0');sb.append('\t');
				sb.append(gene.isPlusStrand()? '+' : '-');
				out.println(sb.toString());
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	

	public static void main(String[] args){
		new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt").toBedFile("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.bed");
	}
}
