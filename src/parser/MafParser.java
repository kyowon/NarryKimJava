package parser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import util.DnDsCalculator;
import util.Nucleotide;
import net.sf.samtools.util.BufferedLineReader;

public class MafParser {
	private String folderName;
	private String indexFilename;
	private HashMap<String, HashMap<Integer, Long>> indexMap;
	private HashMap<Integer, ArrayList<String>> map;
	private ArrayList<Integer> sps;
	private ArrayList<Integer> eps;
	private HashMap<String, ArrayList<Integer>> indexSPositions;
	private HashMap<String, String> fileNameMap;
//	private Bed12Parser parser;
	public static int minNumSpecies = 0;

	
	
	
	public MafParser(String folderName){
		this.folderName = folderName;
		this.indexFilename = folderName + "mafIndex.txt";		
		//this.parser = parser;
	}	
	
	public void readIndexFile(){
		indexMap = new HashMap<String, HashMap<Integer, Long>>();
		indexSPositions = new HashMap<String, ArrayList<Integer>>();
		fileNameMap = new HashMap<String, String>();
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(indexFilename));
			String s;
			HashMap<Integer, Long> indexSubMap = null;
			ArrayList<Integer> indexSubSPositions = null;
			
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				if(s.startsWith("#CONTIG")){
					String contig = token[1];
					fileNameMap.put(contig, token[2]);
					indexMap.put(contig, new HashMap<Integer, Long>());
					indexSPositions.put(contig, new ArrayList<Integer>());
					indexSubMap = indexMap.get(contig);
					
					
					
					
					
					indexSubSPositions = indexSPositions.get(contig);
					continue;
				}
				if(indexSubMap!=null){
					indexSubMap.put(Integer.parseInt(token[0]), Long.parseLong(token[1]));
					indexSubSPositions.add(Integer.parseInt(token[0]));
				}
			}
			for(ArrayList<Integer> i : indexSPositions.values())
				Collections.sort(i);
			in.close();
		}catch (IOException e) {			
			e.printStackTrace();
		}
		
	}
	
	public void generateIndexFile(){
		try {	
			if(new File(indexFilename).exists()) return;
			PrintStream out = new PrintStream(indexFilename);
			
			for(File file : new File(folderName).listFiles()){
				if(!file.getName().endsWith(".maf")) continue;
				if(!file.getName().startsWith("chr")) continue;
				String filename = file.getAbsolutePath();
				System.out.println("Generating maf index file for " + file.getName());
				String s;
				RandomAccessFile accessFile = new RandomAccessFile(filename,    "r");
				String contig = null;
				
				for(long i=0; i<accessFile.length(); i+=100000){
					accessFile.seek(i);
					boolean start = false;
					boolean write = false;
					int sp = -1;
					while((s=accessFile.readLine())!=null){
						if(s.startsWith("#")) continue;
						if(s.startsWith("s ") && start && !write){
							String[] token = s.replaceAll(" +", "\t").split("\t");
							sp = Integer.parseInt(token[2]);
							if(contig == null){
								contig = token[1].substring(token[1].indexOf('.')+1);
								out.println("#CONTIG\t"+contig+"\t"+file.getAbsolutePath());
							}
							write = true;						
						}
						if(s.startsWith("a score")){ 							
							start = true;
							if(write){
								out.println(sp + "\t" + accessFile.getFilePointer());
								break;
							}
						}						
					}					
				}				
				accessFile.close();
			}	
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	   
	}
	
	private void getSubSeqs(int position, int length, int prevIndex, ArrayList<StringBuffer> ret){
		if(length <= 0) return;
		int index = Collections.binarySearch(sps, position);
		//
		if(index < 0) index = - index - 2;
		//
		index = Math.max(prevIndex, index);
		int cposition = 0;
		ArrayList<String> seqs = map.get(sps.get(index < 0? 0 : index));
		//System.out.println(seqs);
		if(ret.isEmpty()) for(int i=0;i<seqs.size();i++) ret.add(new StringBuffer());
		
		if(index >= 0){
			cposition = sps.get(index);			
			for(int i=0;i<seqs.get(0).length();i++){
				boolean isInserted = false;
				if(seqs.get(0).charAt(i) == '-'){
					isInserted = true;
					//continue; // TODO 
				}
				
				if(!isInserted) cposition++;
				if(cposition <= position) continue;
				for(int j=0;j< ret.size();j++){
					ret.get(j).append(seqs.size() <= j? '-' : seqs.get(j).charAt(i));
				}
				if(isInserted) continue;
				if(--length<=0) break;
			}	
	    }
			
		//System.out.println(ret.get(0).length() + " " + length + " " +sps.size() + " " + (index + 1));
		if(length>0 && sps.size() > index + 1){			
			int bposition = cposition;
			cposition = sps.get(index + 1);
			if(bposition < cposition){
				
				length -= cposition - position;
				int append = cposition - position;
				for(int j=0;j< ret.size();j++){
					for(int i=0;i<(length > 0 ? append : append + length);i++)
						ret.get(j).append('-');
				}	
			}
			getSubSeqs(cposition, length, index + 1, ret);
		}
		return;
	}
		
	
	public String getSeqsInMafFormat(String contig, int position, boolean isPlusStrand, int length){
		ArrayList<StringBuffer> seqs = getSubSeqs(contig, position, isPlusStrand, length);
		StringBuffer out = new StringBuffer();
		out.append("a score=1\n");
		for(int i=0;i<seqs.size();i++){
			if(seqs.get(i).length() == 0) continue;
			out.append("s s " + i + " 1 "); 
			out.append(isPlusStrand? "+ 1 " : "- 1 ");
			out.append(seqs.get(i));
			out.append('\n');
		}
		return out.toString();
	}
	
	public String[] getSeqs(String contig, int position, boolean isPlusStrand, int length){
		ArrayList<StringBuffer> seqs = getSubSeqs(contig, position, isPlusStrand, length);
		String[] ret = new String[seqs.size()];
		for(int i=0;i<seqs.size();i++){
			ret[i] = seqs.get(i).toString();
		}
		return ret;
	}
	
	public String getSeqsInMafFormat(String[] contigs, int[] positions, boolean[] strands, int length){
		ArrayList<ArrayList<StringBuffer>> subSeqs = new ArrayList<ArrayList<StringBuffer>>();
		for(int i=0;i<contigs.length;i++){
			ArrayList<StringBuffer> subSeq = getSubSeqs(contigs[i], positions[i], strands[i], length);
			subSeqs.add(subSeq);
		}
		ArrayList<StringBuffer> seqs = mergeSubSeqs(subSeqs);
		StringBuffer out = new StringBuffer();
		out.append("a score=1\n");
		for(int i=0;i<seqs.size();i++){
			out.append("s s " + i + " 1 "); 
			out.append(strands[0]? "+ 1 " : "- 1 ");
			out.append(seqs.get(i));
			out.append('\n');
		}
		return out.toString();
	}
	
	public String[] getSeqs(String[] contigs, int[] positions, boolean[] strands, int length){
		ArrayList<ArrayList<StringBuffer>> subSeqs = new ArrayList<ArrayList<StringBuffer>>();
		for(int i=0;i<contigs.length;i++){
			ArrayList<StringBuffer> subSeq = getSubSeqs(contigs[i], positions[i], strands[i], length);
			subSeqs.add(subSeq);
		}
		ArrayList<StringBuffer> seqs = mergeSubSeqs(subSeqs);
		String[] ret = new String[seqs.size()];
		for(int i=0;i<seqs.size();i++){
			ret[i] = seqs.get(i).toString();
		}
		return ret;
	}
	
	private ArrayList<StringBuffer> getSubSeqs(String contig, int position, boolean isPlusStrand, int length){
		ArrayList<StringBuffer> tmp = new ArrayList<StringBuffer>();
		//System.out.println(position + " " + length);
		if(!isPlusStrand) position -= length-1;
		if(!indexSPositions.containsKey(contig)) return tmp; 
		int index = Collections.binarySearch(indexSPositions.get(contig), position);
		index = index<0? -index-2: index ; 
		index--;
		try {
			RandomAccessFile accessFile = new RandomAccessFile(fileNameMap.get(contig),    "r");
			//System.out.println(indexMap.get(contig).isEmpty());
			if(index >= 0) accessFile.seek(indexMap.get(contig).get(indexSPositions.get(contig).get(index)));
			readChunk(accessFile, position + length);				
			accessFile.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
		
		getSubSeqs(position, length, -1, tmp);
		ArrayList<StringBuffer> ret = new ArrayList<StringBuffer>();
		for(StringBuffer t : tmp) if(t.length()>0) ret.add(t);
		return ret;
	}

	public String[] getSeqs(String contig, ArrayList<Integer> coordinate, int position, boolean isPlusStrand, int length){
		ArrayList<ArrayList<StringBuffer>> subSeqs = new ArrayList<ArrayList<StringBuffer>>();
		//int maxNum = 0;
		//boolean startFlag = false;
		for(ArrayList<Integer> subPositions : Bed12Parser.getCoordinateStartsEnds(coordinate, isPlusStrand)){
			int start = subPositions.get(0);
			int end = subPositions.get(1);
			//System.out.println("* " + start + " " + end);
			if(isPlusStrand){
				if(position >= end) continue;
				start = Math.max(position, start);
			}else{
				if(position < start) continue;
				end = Math.min(position + 1, end);
			}
			
			length -= (end - start);
			ArrayList<StringBuffer> subSeq = getSubSeqs(contig, isPlusStrand? start : end - 1, isPlusStrand, end - start + (length<0? length : 0));
			
			//System.out.println(subSeq.size());
			subSeqs.add(subSeq);
			//maxNum = Math.max(maxNum, subSeq.size());
			if(length <=0) break;
		}
		ArrayList<StringBuffer> seqs = mergeSubSeqs(subSeqs);
		String[] ret = new String[seqs.size()];
		for(int i=0;i<ret.length;i++){
			ret[i] = seqs.get(i).toString();
			if(!isPlusStrand) ret[i] = Nucleotide.getComplementarySeq(ret[i]);
		}
		return ret;
	}
	
	
	private ArrayList<StringBuffer> mergeSubSeqs(ArrayList<ArrayList<StringBuffer>> subSeqs){
		ArrayList<StringBuffer> seqs = new ArrayList<StringBuffer>();
		int max = 0;
		for(ArrayList<StringBuffer> subSeq : subSeqs){
			max = Math.max(max, subSeq.size());
		}
		for(int j=0;j<max;j++){
			seqs.add(new StringBuffer());
		}
		
		for(int i=0;i<subSeqs.size();i++){
			ArrayList<StringBuffer> subSeq = subSeqs.get(i);
			for(int j=0;j<max;j++){
				if(j<subSeq.size()){
					seqs.get(j).append(subSeq.get(j));
				}else if(!subSeq.isEmpty()){
					for(int k=0; k<subSeq.get(0).length();k++)
						seqs.get(j).append('-');
				}
			}
		}
		return seqs;
	}
	
	private void readChunk(RandomAccessFile accessFile, int lastPosition){
		map = new HashMap<Integer, ArrayList<String>>();
		sps = new ArrayList<Integer>();
		eps = new ArrayList<Integer>();
		
		String s;
		ArrayList<String> seqs = new ArrayList<String>();
		int sposition = 0;
		int eposition = 0;
		long spointer = 0;
		try {
			spointer = accessFile.getFilePointer();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		boolean isPlusStrand = true;
		try {
			while((s=accessFile.readLine())!=null){
				//System.out.println(s);
				if(s.startsWith("#")) continue;
				if(s.startsWith("s")){
					String[] token = s.replaceAll(" +", "\t").split("\t");								
					if(seqs.isEmpty()){
						int i=1;								
						i++;
						
						sposition = Integer.parseInt(token[i++]);
						eposition = sposition + Integer.parseInt(token[i++]);
						
						if(!token[i].equals("+")){
							isPlusStrand = false;
						}else isPlusStrand = true;
					}
					seqs.add(token[token.length-1]);
				}				
				if(s.isEmpty() && isPlusStrand){
					if(seqs.size() >= minNumSpecies){					
						map.put(sposition, seqs);
						sps.add(sposition);
						eps.add(eposition);
					}
					seqs = new ArrayList<String>();
				}	
				
				if(s.isEmpty() && (eposition > lastPosition || accessFile.getFilePointer() - spointer > 3000000)) break;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}			
	}
	
	
	
	public static void main(String[] args) throws IOException{
		
		/*RandomAccessFile acc = new RandomAccessFile("/media/kyowon/Data1/RPF_Project/genomes/mm9/maf/chr17.maf", "r");
		acc.seek(1543);
		String s;
		while(true){
			System.out.println(s=acc.readLine());
			if(s.isEmpty()) break;
		}
		*/
		
		String file = "/media/kyowon/Data1/RPF_Project/genomes/hg9/maf/";
		//MafParser.minNumSpecies
		MafParser test = new MafParser(file);
		test.generateIndexFile();
		test.readIndexFile();
		boolean isPlusStrand = false;
		ArrayList<Integer> coordinate = new ArrayList<Integer>();
		coordinate.add(25853545);
		coordinate.add(25853544);
		coordinate.add(25853543);
		coordinate.add(25853542);
		coordinate.add(25853541);	
		coordinate.add(25853540);
		coordinate.add(25853540-1);
		coordinate.add(25853540-2);
		coordinate.add(25853540-3);
		coordinate.add(25853540-4);
		//coordinate.add(25853540-1);
		//coordinate.add(25853540-2);
		//coordinate.add(25853540-3);
		String[] seqs =  test.getSeqs("chr7", coordinate, 25853540+5, isPlusStrand, 10); // why not 150???? TODO
		for(String seq : seqs){
			System.out.println(seq + " " + seq.length());
		}
		//chrX	131254689 chr13	23673470
	
/*chr9	120864455	+
*/

//GATGAACATGGTGAAGAGGATCATGGGGCGGCCTCGGCAGGAGGAGTGCAGCCCGCAAGACAACGCCTTAGGCCTGATGCACCTCCGCCGGCTCTTCACCGAGCTGTGCCACCCTCCGAGGCACATGACCCAGAAGGAGCAGGAGGAGAA
		System.out.println(DnDsCalculator.calculate(seqs));
		//ATGAACATGGTGAAGAGGATCATGGGGCGGCCTCGGCAGGAGGAGTGCAGCCCGCAAGACAACGCCTTAGGCCTGATGCACCTCCGCCGGCTCTTCACCGAGCTGTGCCACCCTCCGAGGCACATGACCCAGAAGGAGCAGGAGGAGAAG
		
		//ATTGCACGCTGTGCCGGCCCTGGAGAAATGGCAGATAAATTATTACTCACTACTCCCTCCAAAAAATTTACATGTCAAGGTCCCGTGGATATCACTATTCAAGCCAAGTGTAATCCCTGCTTATCAAATCCATGTAAAAATGATGGCACC 150
		//ATTGCACGCTGTGCCGGCCCTGGAGAAATGGCAGATAAATTATTACTCACTACTCCCTCCAAAAAATTTACATGTCAAGGTCCCGTGGATATCACTATTCAAGCCAAGTGTAATCCCTGCTTATCAAATCCATGTAAAAATGATGGCACC
	}
}
