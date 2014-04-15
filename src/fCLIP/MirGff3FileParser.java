package fCLIP;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

import parser.BufferedLineReader;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.AnnotationFileParser.txEndComparator;
import parser.AnnotationFileParser.txStartComparator;

public class MirGff3FileParser {
	
	public class startComparator implements Comparator<MiRNA>{
		public int compare(MiRNA o1, MiRNA o2) {
			int i= new Integer(o1.start).compareTo(new Integer(o2.start));
			if(i!=0) return i;
			return new Integer(o1.end).compareTo(new Integer(o2.end));
		}
	}
	public class endComparator implements Comparator<MiRNA>{
		public int compare(MiRNA o1, MiRNA o2) {
			int i= new Integer(o1.end).compareTo(new Integer(o2.end));
			if(i!=0) return i;
			return new Integer(o1.start).compareTo(new Integer(o2.start));
		}
	}
	
	public static class MiRNA{
		private String contig;
		private String name;
		private int start;
		private int end;
		private boolean isPlusStrand;
		//private int[][] pres;
		public String getContig() {
			return contig;
		}

		public String getName() {
			return name;
		}

		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}

		public boolean isPlusStrand() {
			return isPlusStrand;
		}
		
		MiRNA(String s){
			String[] lines = s.split("\n");
			String[] token = lines[0].split("\t");
			this.contig = token[0];
			this.name = token[token.length-1].substring(token[token.length-1].lastIndexOf("Name=")+5);
			this.start = Integer.parseInt(token[3]);
			this.end = Integer.parseInt(token[4]);
			this.isPlusStrand = token[6].equals("+");
			//if(lines.length > 1){
				//TODO if necessary..
			//}
		}
		
		private MiRNA(int start, int end){ // constructor for comparison
			this.start = start;
			this.end = end;
		}
		
		public String toString(){
			return contig + ";" + name + ";" + start + ";" + end + ";" + (isPlusStrand? '+' : '-'); 
		}
	}
	
	private HashMap<String, ArrayList<MiRNA>> miRNAMap;
	
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
	public MiRNA getContainingMiRNA(String contig, boolean isPlusStrand, int position){
		MiRNA ret = null;
		ArrayList<MiRNA> miRNAs = miRNAMap.get(contig);
		if(miRNAs == null) ret = null;
		else{
			int index = Collections.binarySearch(miRNAs, new MiRNA(position, Integer.MAX_VALUE), new startComparator());
			if(index < 0) index = -index - 1;
			if(index - 1 >= miRNAs.size() || index - 1 < 0) ret = null;
			else{
				MiRNA miRNA = miRNAs.get(index - 1);
				if(miRNA.getEnd() > position && miRNA.isPlusStrand() == isPlusStrand) ret = miRNA;
			}
			if(ret == null){
				index = Collections.binarySearch(miRNAs, new MiRNA(0, position + 1), new endComparator());
				if(index < 0) index = -index - 1;
				if(index >= miRNAs.size() || index < 0) ret = null;
				else{
					MiRNA miRNA = miRNAs.get(index);
					if(miRNA.getStart() <= position && miRNA.isPlusStrand() == isPlusStrand) ret = miRNA;
				}
			}
		}
		return ret;
	}
	
}
