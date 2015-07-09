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

public class MafParser {
	private String folderName;
	private String indexFilename;
	private HashMap<String, HashMap<Integer, Long>> indexMap;
	private HashMap<Integer, HashMap<String, String>> map;
	private String mainSpecies = null;
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
			BufferedLineReader in = new BufferedLineReader((indexFilename));
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
	
	private void getSubSeqs(int position, int length, int prevIndex, HashMap<String, StringBuffer> ret){
		if(length <= 0) return;
		int index = Collections.binarySearch(sps, position);
		//
		if(index < 0) index = - index - 2;
		//
		index = Math.max(prevIndex, index);
		int cposition = 0;
		HashMap<String, String> seqs = map.get(sps.get(index < 0? 0 : index));
	
		for(String species : seqs.keySet()){
			if(!ret.containsKey(species)){ 
				ret.put(species, new StringBuffer());
				if(ret.containsKey(mainSpecies)){
					for(int i=0;i<ret.get(mainSpecies).length();i++)
						ret.get(species).append('-');
				}
			}
		}
	
		if(index >= 0){
			cposition = sps.get(index);			
			for(int i=0;i<seqs.get(mainSpecies).length();i++){
				boolean isInserted = false;
				if(seqs.get(mainSpecies).charAt(i) == '-'){
					isInserted = true;
					//continue; // TODO 
				}
				
				if(!isInserted) cposition++;
				if(cposition <= position) continue;
				
				
				for(String species : ret.keySet()){
					ret.get(species).append(!seqs.containsKey(species)? '-' : seqs.get(species).charAt(i));
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
				for(String species : ret.keySet()){
					for(int i=0;i<(length > 0 ? append : append + length);i++)
						ret.get(species).append('-');
				}	
			}
			getSubSeqs(cposition, length, index + 1, ret);
		}
		return;
	}
		
	
	public String getSeqsInMafFormat(String contig, int position, boolean isPlusStrand, int length){
		HashMap<String, StringBuffer> seqs = getSubSeqs(contig, position, isPlusStrand, length);
		StringBuffer out = new StringBuffer();
		out.append("##maf version=1 scoring=autoMZ.v1\na score=1\n");
		
		out.append("s ");out.append(mainSpecies);
		out.append(".");out.append(contig);
		out.append(" ");out.append(position);
		out.append(" ");out.append(length);
		out.append(" + 1 ");
		out.append(seqs.get(mainSpecies).toString());
		out.append('\n');
		
		for(String species : seqs.keySet()){
			if(species.equals(mainSpecies)) continue;
			out.append("s ");out.append(species);
			out.append(".");out.append(contig);
			out.append(" ");out.append(position);
			out.append(" ");out.append(length);
			out.append(" + 1 ");
			out.append(seqs.get(species).toString());
			out.append('\n');
		}
		return out.toString();
	}
	
	public String getSeqsInMafFormat(String contig, ArrayList<Integer> coordinate, int position, boolean isPlusStrand, int length){
		ArrayList<HashMap<String, StringBuffer>> subSeqs = new ArrayList<HashMap<String, StringBuffer>>();
		//int maxNum = 0;
		int tlength = length;
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
			
			tlength -= (end - start);
			HashMap<String, StringBuffer> subSeq = getSubSeqs(contig, isPlusStrand? start : end - 1, isPlusStrand, end - start + (tlength<0? tlength : 0));
			
			//System.out.println(subSeq.size());
			subSeqs.add(subSeq);
			//maxNum = Math.max(maxNum, subSeq.size());
			if(tlength <=0) break;
		}
		HashMap<String, StringBuffer> seqs = mergeSubSeqs(subSeqs);
		StringBuffer out = new StringBuffer();
		out.append("##maf version=1 scoring=autoMZ.v1\na score=1\n");
		
		out.append("s ");out.append(mainSpecies);
		out.append(".");out.append(contig);
		out.append(" ");out.append(position);
		out.append(" ");out.append(length);
		out.append(" + 1 ");
		out.append(seqs.get(mainSpecies).toString());
		out.append('\n');
		
		for(String species : seqs.keySet()){
			if(species.equals(mainSpecies)) continue;
			out.append("s ");out.append(species);
			out.append(".");out.append(contig);
			out.append(" ");out.append(position);
			out.append(" ");out.append(length);
			out.append(" + 1 ");
			out.append(seqs.get(species).toString());
			out.append('\n');
		}
		return out.toString();
	}
	
	public String getSeqsInMafFormat(String[] contigs, int[] positions, boolean[] strands, int length){
		ArrayList<HashMap<String, StringBuffer>> subSeqs = new ArrayList<HashMap<String, StringBuffer>>();
		for(int i=0;i<contigs.length;i++){
			HashMap<String, StringBuffer> subSeq = getSubSeqs(contigs[i], positions[i], strands[i], length); //TODO length should take effect!!!
			subSeqs.add(subSeq);
		}
		HashMap<String, StringBuffer> seqs = mergeSubSeqs(subSeqs);
		StringBuffer out = new StringBuffer();
		out.append("##maf version=1 scoring=autoMZ.v1\na score=1\n");
		
		out.append("s ");out.append(mainSpecies);
		out.append(".");out.append(contigs[0]);
		out.append(" ");out.append(positions[0]);
		out.append(" ");out.append(length);
		out.append(" + 1 ");
		out.append(seqs.get(mainSpecies).toString());
		out.append('\n');
		
		for(String species : seqs.keySet()){
			if(species.equals(mainSpecies)) continue;
			out.append("s ");out.append(species);
			out.append(".");out.append(contigs[0]);
			out.append(" ");out.append(positions[0]);
			out.append(" ");out.append(length);
			out.append(" + 1 ");
			out.append(seqs.get(species).toString());
			out.append('\n');
		}
		return out.toString();
	}
	
	public String[] getSeqs(String contig, int position, boolean isPlusStrand, int length){
		HashMap<String, StringBuffer> seqs = getSubSeqs(contig, position, isPlusStrand, length);
		String[] ret = new String[seqs.size()];
		int i=0;
		ret[i++] = seqs.get(mainSpecies).toString();
		for(String species : seqs.keySet()){
			if(species.equals(mainSpecies)) continue;
			ret[i++] = seqs.get(species).toString();
		}
		return ret;
	}
	
	public String[] getSeqs(String[] contigs, int[] positions, boolean[] strands, int length){
		ArrayList<HashMap<String, StringBuffer>> subSeqs = new ArrayList<HashMap<String, StringBuffer>>();
		for(int i=0;i<contigs.length;i++){
			HashMap<String, StringBuffer> subSeq = getSubSeqs(contigs[i], positions[i], strands[i], length);
			subSeqs.add(subSeq);
		}
		HashMap<String, StringBuffer> seqs = mergeSubSeqs(subSeqs);
		String[] ret = new String[seqs.size()];
		int i=0;
		ret[i++] = seqs.get(mainSpecies).toString();
		for(String species : seqs.keySet()){
			if(species.equals(mainSpecies)) continue;
			ret[i++] = seqs.get(species).toString();
		}
		return ret;
	}
	
	private HashMap<String, StringBuffer> getSubSeqs(String contig, int position, boolean isPlusStrand, int length){
		HashMap<String, StringBuffer> tmp = new HashMap<String, StringBuffer>();
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
		HashMap<String, StringBuffer> ret = new HashMap<String, StringBuffer>();
		for(String species : tmp.keySet()){
			StringBuffer t = tmp.get(species);
			if(!isPlusStrand) t = Nucleotide.getComplementarySeq(t);
			if(t.length()>0) 
				ret.put(species, t);
			//else
			//	ret.put(key, value);
		}
		removeAllIntertedPositions(ret);
		return ret;
	}
	
	private void removeAllIntertedPositions(HashMap<String, StringBuffer> map){
		ArrayList<Integer> indexToRemove = new ArrayList<Integer>();
		StringBuffer s = map.get(mainSpecies);
		for(int i=0;i<s.length();i++){
			boolean toRemove = false;
			for(StringBuffer t : map.values()){
				if(t.charAt(i) != '-'){
					toRemove = false;
					break;
				}
				toRemove = true;
			}
			if(toRemove) indexToRemove.add(i);
		}
		int off = 0;
		for(int i : indexToRemove){
			for(StringBuffer t : map.values()){
				t.deleteCharAt(i - off);
			}
			off++;
		}
		
	}
	

	public String[] getSeqs(String contig, ArrayList<Integer> coordinate, int position, boolean isPlusStrand, int length){
		ArrayList<HashMap<String, StringBuffer>> subSeqs = new ArrayList<HashMap<String, StringBuffer>>();
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
			HashMap<String, StringBuffer> subSeq = getSubSeqs(contig, isPlusStrand? start : end - 1, isPlusStrand, end - start + (length<0? length : 0));
			
			//System.out.println(subSeq.size());
			subSeqs.add(subSeq);
			//maxNum = Math.max(maxNum, subSeq.size());
			if(length <=0) break;
		}
		HashMap<String, StringBuffer> seqs = mergeSubSeqs(subSeqs);
		String[] ret = new String[seqs.size()];
		
		int i=0;
		ret[i++] = seqs.get(mainSpecies).toString();

		for(String species : seqs.keySet()){
			if(species.equals(mainSpecies)) continue;
			ret[i++] = seqs.get(species).toString();
		}
		return ret;
	}
	
	
	
	
	
	private HashMap<String, StringBuffer> mergeSubSeqs(ArrayList<HashMap<String, StringBuffer>> subSeqs){
		HashMap<String, StringBuffer> tseqs = new HashMap<String, StringBuffer>();
		
		for(HashMap<String, StringBuffer> subSeq : subSeqs){
			for(String species : subSeq.keySet())
				tseqs.put(species, new StringBuffer());
		}
		
		for(int i=0;i<subSeqs.size();i++){
			HashMap<String, StringBuffer> subSeq = subSeqs.get(i);
			for(String species : tseqs.keySet()){
				if(subSeq.containsKey(species)){
					tseqs.get(species).append(subSeq.get(species));
				}else if(!subSeq.isEmpty()){
					for(int k=0; k<subSeq.get(mainSpecies).length();k++)
						tseqs.get(species).append('-');
				}
			}
		}
		
		HashMap<String, StringBuffer> seqs = new HashMap<String, StringBuffer>();
		
		int len = 0;
		for(String k : tseqs.keySet()){
			len = tseqs.get(k).length();
			seqs.put(k, new StringBuffer());
		}
		
		for(int i=0;i<len;i++){
			boolean toWrite = false;
			for(String k : tseqs.keySet()){
				if(tseqs.get(k).charAt(i)!='-'){
					toWrite = true; 
					break;
				}
			}
			if(toWrite){
				for(String k : tseqs.keySet()){
					seqs.get(k).append(tseqs.get(k).charAt(i));
				}
			}
			
		}
		
		return seqs;
	}
	
	private void readChunk(RandomAccessFile accessFile, int lastPosition){
		map = new HashMap<Integer, HashMap<String, String>>();
		sps = new ArrayList<Integer>();
		eps = new ArrayList<Integer>();
		
		String s;
		HashMap<String, String> seqs = new HashMap<String, String> ();
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
						mainSpecies = token[1].substring(0, token[1].indexOf('.'));
						isPlusStrand = token[i].equals("+");
					}
					seqs.put(token[1].substring(0, token[1].indexOf('.')), token[token.length-1]);
				}				
				if(s.isEmpty() && isPlusStrand){
					if(seqs.size() >= minNumSpecies){							
						map.put(sposition, seqs);
						sps.add(sposition);
						eps.add(eposition);
					}
					seqs = new HashMap<String, String>();
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
		
		String file = "/media/kyowon/Data1/RPF_Project/genomes/mm9/maf/";
		//MafParser.minNumSpecies
		MafParser test = new MafParser(file);
		test.generateIndexFile();
		test.readIndexFile();
		boolean isPlusStrand = true;
		ArrayList<Integer> coordinate = new ArrayList<Integer>();
		
		for(int i=0;i<150;i++)
			coordinate.add(150158242+i);
		
	
		
		System.out.println(test.getSeqsInMafFormat("chrX", 150158242, isPlusStrand, 120));
		
		//System.out.println(test.getSeqsInMafFormat("chr2", 136711414, !isPlusStrand, 150));

		
		//DnDsCalculator.numSpecies = 12;
		
		
		String[] seqs =  test.getSeqs("chrX", coordinate, 150158242, isPlusStrand, 120); 
		for(String seq : seqs){
			System.out.println(seq + " " + seq.length());
		}
	
		System.out.println(DnDsCalculator.calculate(seqs));
	}
}
