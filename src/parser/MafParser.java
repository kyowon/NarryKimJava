package parser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import util.Nucleotide;
import net.sf.samtools.util.BufferedLineReader;

public class MafParser {
	private HashMap<String, HashMap<Integer, ArrayList<String>>> map;
	private HashMap<String, ArrayList<Integer>> sPositionMap;
	private HashMap<String, ArrayList<Integer>> ePositionMap;
	
	public MafParser(String name){
		map = new HashMap<String, HashMap<Integer, ArrayList<String>>>();
		sPositionMap = new HashMap<String, ArrayList<Integer>>();
		ePositionMap = new HashMap<String, ArrayList<Integer>>();
		
		try {
			BufferedLineReader in = new BufferedLineReader(new GZIPInputStream(new FileInputStream(name)));
			String s;
			ArrayList<String> seqs = new ArrayList<String>();
			int sposition = 0;
			int eposition = 0;
			String contig = null;
			boolean isPlusStrand = true;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				if(s.startsWith("s")){
					String[] token = s.split(" ");								
					if(seqs.isEmpty()){
						int i=1;		
						contig = token[i].substring(token[i].indexOf('.')+1);
						i++;
						for(;i<token.length;i++) if(!token[i].isEmpty()) break;
						
						sposition = Integer.parseInt(token[i++]);
						for(;i<token.length;i++) if(!token[i].isEmpty()) break;
						eposition = sposition + Integer.parseInt(token[i++]);
						for(;i<token.length;i++) if(!token[i].isEmpty()) break;
						
						if(!token[i].equals("+")){
							//System.out.println("Damn, -" + name);
							//System.out.println(contig + " " + s);
							//System.out.println(token[i]);
							//System.exit(1);
							isPlusStrand = false;
						}else isPlusStrand = true;
					}
					seqs.add(token[token.length-1]);
				}				
				if(s.isEmpty() && isPlusStrand){
					if(seqs.size() >= 10){
						if(!map.containsKey(contig)){
							map.put(contig, new HashMap<Integer, ArrayList<String>>());
							sPositionMap.put(contig, new ArrayList<Integer>());
							ePositionMap.put(contig, new ArrayList<Integer>());
						}
						map.get(contig).put(sposition, seqs);
						sPositionMap.get(contig).add(sposition);
						ePositionMap.get(contig).add(eposition);
					}
					seqs = new ArrayList<String>();
				}				
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		for(String contig : sPositionMap.keySet()){
			Collections.sort(sPositionMap.get(contig));
			Collections.sort(ePositionMap.get(contig));
		}
		
	}
	
	private void getSeqs(String contig, int position, int length, int prevIndex, ArrayList<StringBuffer> ret){
	//	System.out.println(position);
		if(length <= 0) return;
		ArrayList<Integer> positions = sPositionMap.get(contig);
		if(positions == null) return;
		int index = Collections.binarySearch(positions, position);
		//
		if(index < 0) index = - index - 2;
		//
		index = Math.max(prevIndex, index);
		int cposition = 0;
		ArrayList<String> seqs = map.get(contig).get(positions.get(index < 0? 0 : index));
		if(ret.isEmpty()) for(int i=0;i<seqs.size();i++) ret.add(new StringBuffer());
		
		if(index >= 0){
			cposition = positions.get(index);
			
			//System.out.println(position);
			for(int i=0;i<seqs.get(0).length();i++){
				if(seqs.get(0).charAt(i) == '-'){
					continue;
				}
				cposition++;
				//System.out.println(cposition + " " + position);
				if(cposition <= position) continue;
				for(int j=0;j< ret.size();j++){
					ret.get(j).append(seqs.size() <= j? '-' : seqs.get(j).charAt(i));
				}
				if(--length<=0) break;
			}	
	    }
		if(length>0 && positions.size() > index + 1){		
			int bposition = cposition;
			cposition = positions.get(index + 1);
			//System.out.println(cposition + " " + position + " " + length);
			if(bposition < cposition){
				length -= cposition - position;
				int append = cposition - position;
			//	System.out.println(append + " " + length);
				for(int j=0;j< ret.size();j++){
					for(int i=0;i<(length > 0 ? append : append + length);i++)
						ret.get(j).append('-');
				}	
			}
			//}		
		//	System.out.println(cposition + " " + position);
			getSeqs(contig, cposition, length, index + 1, ret);
		}
		return;
	}
	
	public String[] getSeqs(String contig, int position, boolean isPlusStrand, int length){
		ArrayList<StringBuffer> ret = new ArrayList<StringBuffer>();
		if(!isPlusStrand) position -= length;
		getSeqs(contig, position, length, -1, ret);
		String[] seqs = new String[ret.size()];
		for(int i=0;i<seqs.length;i++){
			// seqs[i] = ret.get(i).reverse().toString();
			seqs[i] = ret.get(i).toString();
			if(!isPlusStrand) seqs[i] = Nucleotide.getComplementarySeq(seqs[i]);
		}
		return seqs;
	}
	
	public Set<String> getContigs(){
		return sPositionMap.keySet();
	}
	
	public ArrayList<Integer> getStartPositions(String contig){
		return sPositionMap.get(contig);
	}
	
	public ArrayList<Integer> getEndPositions(String contig){
		return ePositionMap.get(contig);
	}
	
	
	static public void main(String[] args){
		String mafFileDir = "/media/kyowon/Data1/RPF_Project/genomes/mm9/maf";
		for(File mafFile : new File(mafFileDir).listFiles()){
			if(!mafFile.getName().endsWith("18.maf.gz"))continue;
			MafParser mp = new MafParser(mafFile.getAbsolutePath());
			System.out.println(mp.getSeqs("chr18", 3115284+1, true, 30)[0]);
			System.out.println(mafFile.getName());
		}		
	}
	//3115282
}
