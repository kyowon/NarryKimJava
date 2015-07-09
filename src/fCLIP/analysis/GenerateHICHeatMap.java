package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;

public class GenerateHICHeatMap {
	private static HashMap<String, HashMap<String, HashMap<Integer, HashMap<Integer, Integer>>>> map;	
	public static void main(String[] args) {
		try {
			BufferedLineReader sam1 = new BufferedLineReader("/media/kyowon/Data1/fCLIP/hic/SRR027962_1.split.sam");
			BufferedLineReader sam2 = new BufferedLineReader("/media/kyowon/Data1/fCLIP/hic/SRR027962_2.split.sam");
			String outFolder = "/media/kyowon/Data1/fCLIP/hic/result2";
			String s1, s2;
			int binsize = 1000000;
			
			while((s1=sam1.readLine()).startsWith("@"));
			while((s2=sam2.readLine()).startsWith("@"));
			
			HashMap<Integer, ArrayList<String>> read1Map = new HashMap<Integer, ArrayList<String>>(); // read number, chr, position strings
			HashMap<Integer, ArrayList<String>> read2Map = new HashMap<Integer, ArrayList<String>>(); // read number, chr, position strings
			
			while((s1=sam1.readLine())!=null){
				String[] token1 = s1.split("\t");
				int readNumber1 = Integer.parseInt(token1[0].substring(token1[0].indexOf('.') + 1));
				int binnumber1 = Integer.parseInt(token1[3])/binsize;
				String chr1 = token1[2];
				
				if(!read1Map.containsKey(readNumber1)) read1Map.put(readNumber1, new ArrayList<String>());
				read1Map.get(readNumber1).add(chr1 + ";" + binnumber1);
			}
			
			sam1.close();
			
			while((s2=sam2.readLine())!=null){
				String[] token2 = s2.split("\t");
				int readNumber2 = Integer.parseInt(token2[0].substring(token2[0].indexOf('.') + 1));
				int binnumber2 = Integer.parseInt(token2[3])/binsize;
				String chr2 = token2[2];
				
				if(!read2Map.containsKey(readNumber2)) read2Map.put(readNumber2, new ArrayList<String>());
				read2Map.get(readNumber2).add(chr2 + ";" + binnumber2);
			}
			sam2.close();
			
			map = new HashMap<String, HashMap<String,HashMap<Integer,HashMap<Integer,Integer>>>>(); // chr1 chr2 pos1 pos2
			
			for(int readNumber1 : read1Map.keySet()){
				if(!read2Map.containsKey(readNumber1)) continue;
				ArrayList<String> pos1 = read1Map.get(readNumber1);
				ArrayList<String> pos2 = read2Map.get(readNumber1);
				
				for(String p1 : pos1){
					String[] t1 = p1.split(";");
					String key1 = t1[0];
					int binnumber1 = Integer.parseInt(t1[1]);
					for(String p2 : pos2){
						String[] t2 = p2.split(";");
						String key2 = t2[0];						
						int binnumber2 = Integer.parseInt(t2[1]); // binnumber 1 2 exchangable -_-
						
						if(!map.containsKey(key1)) map.put(key1, new HashMap<String, HashMap<Integer,HashMap<Integer,Integer>>>());
						HashMap<String, HashMap<Integer,HashMap<Integer,Integer>>> smap = map.get(key1);
						if(!smap.containsKey(key2)) smap.put(key2, new HashMap<Integer, HashMap<Integer,Integer>>());
						HashMap<Integer, HashMap<Integer,Integer>> ssmap = smap.get(key2);
						if(!ssmap.containsKey(binnumber1)) ssmap.put(binnumber1, new HashMap<Integer,Integer>());
						HashMap<Integer,Integer> sssmap = ssmap.get(binnumber1);
						if(!sssmap.containsKey(binnumber2)) sssmap.put(binnumber2, 0);
						sssmap.put(binnumber2, sssmap.get(binnumber2) + 1);
					}
				}
			}
			
			HashMap<String, Integer> binNumber = new HashMap<String, Integer>();
			BufferedLineReader genomeSizeReader = new BufferedLineReader("/media/kyowon/Data1/fCLIP/genomes/genome.size");
			String s;
			while((s=genomeSizeReader.readLine())!=null){
				String[] token = s.split("\t");
				binNumber.put(token[0], 1 + Integer.parseInt(token[1]) / binsize);
			}
			genomeSizeReader.close();
			
			for(String key1 : map.keySet()){
				if(key1.equals("chrM") || key1.length()>5) continue;
				//System.out.println(key1);
				Integer binNumber1 = binNumber.get(key1);
				if(binNumber1 == null) continue;
				for(String key2 : map.get(key1).keySet()){
					if(key2.equals("chrM") || key2.length()>5) continue;
					Integer binNumber2 = binNumber.get(key2);
					if(binNumber2 == null) continue;
					String fileKey = (key1.equals("chrX")? "chr23" : (key1.equals("chrY") ? "chr24" : key1)) + "_" + (key2.equals("chrX")? "chr23" : (key2.equals("chrY") ? "chr24" : key2));
					StringBuffer fileValue = new StringBuffer();
					 HashMap<Integer, HashMap<Integer, Integer>> smap = map.get(key1).get(key2);
					for(int bn1 = 0; bn1<binNumber1;bn1++){
						if(!smap.containsKey(bn1)) smap.put(bn1, new HashMap<Integer, Integer>());
						HashMap<Integer, Integer> ssmap = smap.get(bn1);
						for(int bn2 = 0; bn2<binNumber2;bn2++){
							if(!ssmap.containsKey(bn2)) ssmap.put(bn2, 0);	
							if(bn1 < bn2) fileValue.append(0); //TODO
							else fileValue.append(ssmap.get(bn2));
							if(bn2<binNumber2-1) fileValue.append("\t");
						}
						fileValue.append("\n");
					}				
					
					PrintStream out = new PrintStream(outFolder + "/" + fileKey + ".txt");
					out.print(fileValue.toString());
					out.close();					
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

}
