package parser;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import util.Nucleotide;
import net.sf.samtools.util.BufferedLineReader;

public class MafParser {
	private HashMap<String, HashMap<Integer, ArrayList<String>>> map;
	private HashMap<String, ArrayList<Integer>> positionMap;
	
	public MafParser(String name){
		map = new HashMap<String, HashMap<Integer, ArrayList<String>>>();
		positionMap = new HashMap<String, ArrayList<Integer>>();
		
		try {
			BufferedLineReader in = new BufferedLineReader(new GZIPInputStream(new FileInputStream(name)));
			String s;
			ArrayList<String> seqs = new ArrayList<String>();
			int position = 0;
			String contig = null;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				if(s.startsWith("s")){
					String[] token = s.split(" ");								
					if(seqs.isEmpty()){
						int i=1;		
						contig = token[i].substring(token[i].indexOf('.')+1);
						i++;
						for(;i<token.length;i++) if(!token[i].isEmpty()) break;
						
						position = Integer.parseInt(token[i++]);
						for(;i<token.length;i++) if(!token[i].isEmpty()) break;
						i++;
						for(;i<token.length;i++) if(!token[i].isEmpty()) break;
						
						if(!token[i].equals("+")){
							System.out.println("Damn, -");
							System.out.println(token[i]);
						}
					}
					seqs.add(token[token.length-1]);
				}				
				if(s.isEmpty()){
					if(!map.containsKey(contig)){
						map.put(contig, new HashMap<Integer, ArrayList<String>>());
						positionMap.put(contig, new ArrayList<Integer>());
					}
					map.get(contig).put(position, seqs);
					positionMap.get(contig).add(position);
					seqs = new ArrayList<String>();
				}				
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		for(String contig : positionMap.keySet())
			Collections.sort(positionMap.get(contig));
	}
	
	private void getSeqs(String contig, int position, int length, int prevIndex, ArrayList<StringBuffer> ret){
		if(length <= 0) return;
		ArrayList<Integer> positions = positionMap.get(contig);
		if(positions == null) return;
		int index = Collections.binarySearch(positions, position);
		//
		if(index < 0) index = - index - 2;
		//
		index = Math.max(prevIndex, index);
		int cposition = 0;
		ArrayList<String> seqs = map.get(contig).get(positions.get(0));
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
		if(length>0){		
			
			int bposition = cposition;
			cposition = positions.get(index + 1);
			//System.out.println(cposition + " " + position + " " + length);
			if(bposition < cposition){
				length -= cposition - position;
				int append = cposition - position;
				System.out.println(append + " " + length);
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
	
	
	static public void main(String[] args){
		String maf = "/media/kyowon/Data1/RPF_Project/genomes/mm9/maf/chr13_random.maf.gz";
		MafParser test = new MafParser(maf);
		
	//	ArrayList<StringBuffer> ret = new ArrayList<StringBuffer>();
		//589 -586
		//2 1
		//test.getSeqs("chr13_random", 1578, 4, 0, ret);
		for(String seq : test.getSeqs("chr13_random", 6, false, 3)){
			System.out.println(seq.toString());
			//System.out.println(seq.length());
		}
	}
	
}
