package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import parser.ScoringOutputParser.ScoredPosition;
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
		private String genomeBrowserGeneName;		
		private String geneName;
		private String contig;
		private boolean isPlusStrand;
		private int txStart;
		private int txEnd;
		private int cdsStart;
		private int cdsEnd;
		private int exonCount;		
		private int[] exonStarts;
		private int[] exonEnds;
		private int[][] introns;
		
		public AnnotatedGene(String s){
			String[] token = s.split("\t");
			genomeBrowserGeneName = token[0];
			geneName = token[1];
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
		
		private AnnotatedGene(int txStart, int txEnd){ // constructor for comparison
			this.txStart = txStart;
			this.txEnd = txEnd;
		}
		
		public String toString(){
			StringBuffer sb = new StringBuffer();			
			sb.append(genomeBrowserGeneName); sb.append('\t');
			sb.append(geneName); sb.append('\t');
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
		
		public static String getEmptyGeneString(){
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
		
		public static String getHeader(){
			StringBuffer sb = new StringBuffer();
			sb.append("GenomeBrowserGeneName");sb.append('\t');
			sb.append("GeneName");sb.append('\t');
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
		
		public String getGenomeBrowserGeneName() {
			return genomeBrowserGeneName;
		}
		public String getGeneName() {
			return geneName;
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
			if(introns == null){
				introns = new int[getExonCount()-1][2];
				for(int i=1;i<getExonCount();i++){
					introns[i-1][0] = getExonEnds()[i-1];
					introns[i-1][1] = getExonStarts()[i]-1;
				}
			}
			return introns;
		}
		
	}
	
	private HashMap<String, ArrayList<AnnotatedGene>> annotatedGeneSetMap;
	private HashMap<String, HashMap<Boolean, HashMap<Integer, AnnotatedGene>>> annotatedGenePositionMap;
	
//	public static AnnotatedGene getSudoAnnotatedGene(ScoredPosition position){
//		return new AnnotationFileParser().new AnnotatedGene(position);
//	}
	
	
	public AnnotationFileParser(String annotationFile){
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(annotationFile));
			annotatedGeneSetMap = new HashMap<String, ArrayList<AnnotatedGene>>();
			annotatedGenePositionMap = new HashMap<String, HashMap<Boolean, HashMap<Integer, AnnotatedGene>>>();
			String s;
			while((s=in.readLine())!=null){
				AnnotatedGene gene = new AnnotatedGene(s);
				String contig = gene.getContig();
				if(!annotatedGeneSetMap.containsKey(contig)){
					annotatedGeneSetMap.put(contig, new ArrayList<AnnotatedGene>());
					annotatedGenePositionMap.put(contig, new HashMap<Boolean, HashMap<Integer, AnnotatedGene>>());
					annotatedGenePositionMap.get(contig).put(true, new HashMap<Integer, AnnotatedGene>());
					annotatedGenePositionMap.get(contig).put(false, new HashMap<Integer, AnnotatedGene>());
				}
				annotatedGeneSetMap.get(contig).add(gene);	
				annotatedGenePositionMap.get(contig).get(gene.isPlusStrand()).put((gene.isPlusStrand()? gene.getCdsStart() : gene.getCdsEnd()), gene);
			}
			for(String contig : annotatedGeneSetMap.keySet()){
				Collections.sort(annotatedGeneSetMap.get(contig), new txStartComparator());
			}
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	/**
	    * Get the annotated gene
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the coding start position ( see below for the examples)
	    * @return the annotated gene with the position as the CDS start or as the CDS end.
	    * For example, 
	    * 
	    * Gene1 : chr1   +    cds start 0 cds end 100    will be returned for getAnnotatedGene(chr1, true, 0)
	    * Gene2:  chr3   -    cds start 10  cds end 123  will be returned for getAnnotatedGene(chr3, false, 122), not for getAnnotatedGene(chr3, false, 123) 
	    * 
	    * This is because cds end is exclusive in the ref flat file
	    */
	public AnnotatedGene getAnnotatedGene(String contig, boolean isPlusStrand, int position){
		HashMap<Boolean, HashMap<Integer, AnnotatedGene>> sub = annotatedGenePositionMap.get(contig);
		if(sub == null) return null;
		HashMap<Integer, AnnotatedGene> ssub = sub.get(isPlusStrand);
		if(ssub == null) return null;
		return ssub.get(position + (isPlusStrand? 0 : 1));
	}
	
	/**
	    * Get the containing gene
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the tx start position
	    * @return the containing gene
	     * For example, 
	    * 
	    * Gene1 : chr1   +    tx start 0 tx end 100    will be returned for getContainingGene(chr1, true, 50)
	    * Gene2:  chr3   -    tx start 10  tx end 123  will be returned for getContainingGene(chr3, false, 120), but not for getContainingGene(chr3, false, 123) 
	    * 
	    * This is because tx end is exclusive in the ref flat file
	    */
	public AnnotatedGene getContainingGene(String contig, boolean isPlusStrand, int position){
		AnnotatedGene ret = getAnnotatedGene(contig, isPlusStrand, position);
		if(ret != null) return ret;
		ArrayList<AnnotatedGene> genes = annotatedGeneSetMap.get(contig);
		if(genes == null) ret = null;
		else{
			int index = Collections.binarySearch(genes, new AnnotatedGene(position, Integer.MAX_VALUE), new txStartComparator());
			if(index < 0) index = -index - 1;
			if(index - 1 >= genes.size() || index - 1 < 0) ret = null;
			else{
				AnnotatedGene gene = genes.get(index - 1);
				if(gene.getTxEnd() > position && gene.isPlusStrand() == isPlusStrand) ret = gene;
			}
			if(ret == null){
				index = Collections.binarySearch(genes, new AnnotatedGene(0, position + 1), new txEndComparator());
				if(index < 0) index = -index - 1;
				if(index >= genes.size() || index < 0) ret = null;
				else{
					AnnotatedGene gene = genes.get(index);
					if(gene.getTxStart() <= position && gene.isPlusStrand() == isPlusStrand) ret = gene;
				}
			}
		}
		return ret;
	}
	
	/**
	    * Get the genomic region and frameshift information
	    * @param contig the contig
	    * @param isPlusStrand is plus strand?
	    * @param position the position
	    * @return ArrayList of String : 0 genomic region, 1 frameshift (0/1/2)
	    */
	public ArrayList<String> getGenomicRegionNameAndFrameShift(String contig, boolean isPlusStrand, int position){
		ArrayList<String> ret = new ArrayList<String>();
		AnnotatedGene gene = getAnnotatedGene(contig, isPlusStrand, position);
		if(gene != null){
			//ret.add("NM_Start"); ret.add("0"); 
			//return ret;
		}
		else gene = getContainingGene(contig, isPlusStrand, position);								
		
		String name;
		String frameShift;
		
		if(gene == null){
			name = "InterGenic";
			frameShift = "_";
		}else{
			if(gene.getGenomeBrowserGeneName().startsWith("LINC")) name = "LINC";
			else if(gene.getGeneName().startsWith("NR")) name = "NR";
			else name = "NM";
			
			if(position >= gene.getCdsStart() && position < gene.getCdsEnd()){
				name += "_ORF";
				int[][] introns = gene.getIntrons();
				for(int i=0;i<introns.length;i++){
					if(position >=introns[i][0] && position <= introns[i][1]){
						name += "_Intron";
						break;
					}
				}				
			}else{
				if(isPlusStrand && position < gene.getCdsStart()) name += "_5_UTR";
				else if(!isPlusStrand && position >= gene.getCdsEnd()) name += "_5_UTR";
				else name += "_3_UTR";
			}	
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
		
		ret.add(name); ret.add(frameShift);
		return ret;
	}
	
	
	public ArrayList<ArrayList<Integer>> getLiftOverCDSPositions(AnnotatedGene gene){	
		return getLiftOverPositions(gene, gene.isPlusStrand(), gene.isPlusStrand()? gene.cdsStart : gene.cdsEnd - 1, gene == null? 0 : Integer.MAX_VALUE, true);
	}
	
	public ArrayList<ArrayList<Integer>> getLiftOverPositionsTillNextStopCodon(String contig, boolean isPlusStrand, int start, int maxLength, ZeroBasedFastaParser fastaParser){
		AnnotatedGene gene = getContainingGene(contig, isPlusStrand, start);
		//AnnotatedGene tgene = getAnnotatedGene(contig, isPlusStrand, start);
		
		ArrayList<ArrayList<Integer>> positions = getLiftOverPositions(gene, isPlusStrand, start, maxLength, true);
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();
		HashSet<String> stopCodons = new HashSet<String>();
		stopCodons.add("TAG");
		stopCodons.add("TAA");
		stopCodons.add("TGA");
		boolean stop = false;
		ArrayList<Integer> tp = new ArrayList<Integer>();
		for(int j=0;j<positions.size();j++){
			ArrayList<Integer> subPositions = positions.get(j);
			ArrayList<Integer> subRet = new ArrayList<Integer>();
			for(int i=0;i<subPositions.size();i++){				
				subRet.add(subPositions.get(i));
				tp.add(subPositions.get(i));
				if(tp.size()%3 ==0 && tp.size()>=3){
					String tc = fastaParser.getSequence(contig, tp.subList(tp.size()-3, tp.size()));
					if(!isPlusStrand) tc = ZeroBasedFastaParser.getComplementarySequence(tc, false);
					if(stopCodons.contains(tc.toUpperCase())){ 					
						ret.add(subRet);
						stop = true;
						break;
					}
				}
			}
			if(stop) break;
			//Collections.sort(subRet);
			ret.add(subRet);
		}
		return ret;
	}
	
	
	public ArrayList<ArrayList<Integer>> getLiftOverPositions(String contig, boolean isPlusStrand, int start, int length){
		AnnotatedGene gene = getContainingGene(contig, isPlusStrand, start);
		//if(gene == null){
		//	if(isPlusStrand) gene = getContainingGene(contig, isPlusStrand, start + length - 1);
		//	else gene = getContainingGene(contig, isPlusStrand, start - length + 1);			
		//}
		return getLiftOverPositions(gene, isPlusStrand, start, length, false);
	}
	
	// start : inclusive
	private ArrayList<ArrayList<Integer>> getLiftOverPositions(AnnotatedGene gene, boolean isPlusStrand, int start, int length, boolean stopAtCDSEnd){
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();
		int[][] introns = gene==null? new int[0][0] : gene.getIntrons();
		if(isPlusStrand){
			int intronIndex = 0;
			int c = start;
			for(;intronIndex<introns.length;intronIndex++){
				if(introns[intronIndex][1] >= c){
					//if(introns[intronIndex][0] <= c) return null;
					break;
				}
			}	
			//System.out.println(gene);
			//System.out.println(start + " " + gene.getCdsEnd());
			boolean startLift = false;
			ArrayList<Integer> subPositions = new ArrayList<Integer>();
			int tlength = 0;
			while(true){
				subPositions.add(c++);
				tlength++;
				if(stopAtCDSEnd && gene != null && (start >= gene.getCdsEnd() - 3? c >=gene.getTxEnd() : c >= gene.getCdsEnd())) break;
				else if(tlength >= length) break;
				
				if(intronIndex >= 0 && intronIndex < introns.length && c >= introns[intronIndex][0]){ // if c is within intron
					if(!startLift) continue;
					c = introns[intronIndex][1] + 1;
					intronIndex++;
					if(!subPositions.isEmpty())ret.add(subPositions);
					subPositions = new ArrayList<Integer>();
				}else startLift = true;
			}
			if(!subPositions.isEmpty())ret.add(subPositions);
		//	for(ArrayList<Integer> sss : ret)
		//		System.out.println(sss);
		}else{
			int intronIndex = introns.length-1;
			int c = start;
			for(;intronIndex>=0;intronIndex--){
				if(introns[intronIndex][0] <= c){
					//if(introns[intronIndex][1] >= c) return null;
					break;
				}
			}		
			boolean startLift = false; 
			ArrayList<Integer> subPositions = new ArrayList<Integer>();
			int tlength = 0;
			while(true){
				subPositions.add(c--);
				tlength++;
				if(stopAtCDSEnd && gene != null && (start < gene.getCdsStart() + 3 ? c < gene.getTxStart() : c < gene.getCdsStart())) break;
				else if(tlength >= length || c < 0) break;
							
				//if(stopAtCDSEnd && gene!=null && c < gene.getCdsStart()) break;
				if(intronIndex >= 0 && intronIndex < introns.length && c <= introns[intronIndex][1]){					
					if(!startLift) continue;
					c = introns[intronIndex][0] - 1;
					intronIndex--;
					if(!subPositions.isEmpty())ret.add(subPositions);
					subPositions = new ArrayList<Integer>();
				}else startLift = true;
			}	
			if(!subPositions.isEmpty())ret.add(subPositions);
		}		
		return ret;
	}
	
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
	

	public static void main(String[] args){
		AnnotationFileParser test = new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.refFlat.txt");
		//System.out.println();
		
		
		ZeroBasedFastaParser fasta = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.fa");
	//	System.out.println(test.getLiftOverPositions("chr17", false, 28695253+12, 5));
	//
		//getLiftOverPositionsTillNextStopCodon(//)//test.getAnnotatedGene("chr5", false, 102410577))
		
		for(ArrayList<Integer> positions : test.getLiftOverPositionsTillNextStopCodon("chr5", false, 102410577, 150 , fasta)){
			System.out.println(positions.get(0) + " " + positions.get(positions.size()-1));
			//System.out.println();
			System.out.println(ZeroBasedFastaParser.getComplementarySequence(fasta.getSequence("chr5", positions),false));
			//break;
		}
		//SREK1	NM_001270493	chr5	+	65440045	65453568	65440204	65451900	3	65440045,65449290,65451892,	65440365,65449424,65453568,

		//GATGAACATGGTGAAGAGGATCATGGGGCGGCCTCGGCAGGAGGAGTGCAGCCCGCAAGACAACGCCTTAGGCCTGATGCACCTCCGCCGGCTCTTCACCGAGCTGTGCCACCCTCCGAGGCACATGACCCAGAAGGAGCAGGAGGAGAA
		//System.out.println(fasta.getSequence("chr5", 102410577, 102410570));
	//	System.out.println(fasta.getSequence("chr5", test.getLiftOverPositions("chr5", false, 3648025, 3)));
	//	System.out.println(fasta.getSequence("chr5", test.getLiftOverPositions("chr5", false, 4086564, 3)));
	}//4086564
//GYPE	NM_198682	chr4	-	144792018	144826716	144797907	144826660	4	144792018,144797898,144801563,144826623,	144793405,144798008,144801662,144826716,
	//PTH	NM_000315	chr11	-	13513600	13517567	13513951	13514402	3	13513600,13514316,13517457,	13514213,13514407,13517567,

}
