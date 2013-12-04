package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

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
	
	public class AnnotatedGene{
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
		
		public AnnotatedGene(){
			
		}
		
		AnnotatedGene(int txStart, int txEnd){ // constructor for comparison
			this.txStart = txStart;
			this.txEnd = txEnd;
		}
		
		public String toString(){
			StringBuffer sb = new StringBuffer();
			
			sb.append(genomeBrowserGeneName); sb.append('\t');
			sb.append(geneName); sb.append('\t');
			sb.append(contig); sb.append('\t');
			sb.append(isPlusStrand); sb.append('\t');
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
	}
	
	private HashMap<String, ArrayList<AnnotatedGene>> annotatedGeneSetMap;
	private HashMap<String, HashMap<Boolean, HashMap<Integer, AnnotatedGene>>> annotatedGenePositionMap;
	
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
				annotatedGenePositionMap.get(contig).get(gene.isPlusStrand()).put((gene.isPlusStrand()? gene.getCdsStart() : gene.getCdsEnd() - 1), gene);
			}
			for(String contig : annotatedGeneSetMap.keySet()){
				Collections.sort(annotatedGeneSetMap.get(contig), new txStartComparator());
			}
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	
	public AnnotatedGene getAnnotatedGene(String contig, boolean isPlusStrand, int position){
		HashMap<Boolean, HashMap<Integer, AnnotatedGene>> sub = annotatedGenePositionMap.get(contig);
		if(sub == null) return null;
		HashMap<Integer, AnnotatedGene> ssub = sub.get(isPlusStrand);
		if(ssub == null) return null;
		return ssub.get(position);
	}
	
	public AnnotatedGene getContainingGene(String contig, boolean isPlusStrand, int position){
		ArrayList<AnnotatedGene> genes = annotatedGeneSetMap.get(contig);
		AnnotatedGene ret = null;
		if(genes == null) ret = null;
		else{
			int index = Collections.binarySearch(genes, new AnnotatedGene(position, Integer.MAX_VALUE), new txStartComparator());
			if(index < 0) index = -index - 1;
			if(index - 1 >= genes.size() || index - 1 < 0) ret = null;
			else{
				AnnotatedGene gene = genes.get(index - 1);
				if(gene.getTxEnd() >= position && gene.isPlusStrand() == isPlusStrand) ret = gene;
			}
			if(ret == null){
				index = Collections.binarySearch(genes, new AnnotatedGene(0, position), new txEndComparator());
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
		AnnotationFileParser test = new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/data/refFlatHuman.txt");
		//ZeroBasedFastaParser fasta = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/data/hg19.fa");
		Iterator<AnnotatedGene> iterator = test.getAnnotatedGeneIterator("chr1");
		while(iterator.hasNext()){
			AnnotatedGene gene = iterator.next();
			//if(gene.txStart == 235272657)System.out.println(test.annotatedGeneSetMap.get("chr1").indexOf(gene) + " " +  gene.txStart + " " + gene.txEnd + " " + gene.geneName + " " + gene.getGenomeBrowserGeneName());
			//if(gene.isPlusStrand())
				//System.out.println(fasta.getSequence(gene.getContig(), gene.getCdsStart(), gene.getCdsStart()+3));
			//else 
				//System.out.println(fasta.getSequence(gene.getContig(), gene.getCdsEnd()-3, gene.getCdsEnd()));
		}
		//System.out.println(test.getContainingGene("chr5", 36152364).getGeneName());
	}
}
