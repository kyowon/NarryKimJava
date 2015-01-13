package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class GetRegionDistribution {

	public static void getScoredPositionRegion(String file, String outFile){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(file);
		HashMap<String, HashMap<String, Integer>> map5 = new HashMap<String, HashMap<String, Integer>>(); 
		HashMap<String, HashMap<String, Integer>> map3 = new HashMap<String, HashMap<String, Integer>>(); 
		
		for(ScoredPosition p : parser.getPositions()){
			String key = p.hasMatchingMiRNA() ? "miRNA" :  p.getClassification();
			if(key.equals("M")) key += p.isPaired()? "T" : "F";
			
			String r5 = p.getGenomicRegions3p() == null? "_" : p.getGenomicRegions3p().get(0);
			String v5 = "";
			if(r5.equals("_")) v5 = "InterGenic";
			else if(r5.endsWith("Intron")) v5 = "Intron";
			else if(r5.equals("NR")) v5 = "NR";
			else if(r5.endsWith("UTR")) v5 = "UTR";
			else if(r5.equals("LINC")) v5 = "LINC";
			else if(r5.endsWith("ORF")) v5 = "ORF";
			else v5 = r5;
			
			if(!map5.containsKey(key)) map5.put(key, new HashMap<String, Integer>());
			HashMap<String, Integer> smap5 = map5.get(key);
			if(!smap5.containsKey(v5)) smap5.put(v5, 0);
			smap5.put(v5, smap5.get(v5) + 1);
			
			String r3 = p.getGenomicRegions5p() == null? "_" : p.getGenomicRegions5p().get(0);
			String v3 = "";
			if(r3.equals("_")) v3 = "InterGenic";
			else if(r3.endsWith("Intron")) v3 = "Intron";
			else if(r3.equals("NR")) v3 = "NR";
			else if(r3.endsWith("UTR")) v3 = "UTR";
			else if(r3.equals("LINC")) v3 = "LINC";
			else if(r3.endsWith("ORF")) v3 = "ORF";
			else v3 = r3;
			
			if(!map3.containsKey(key)) map3.put(key, new HashMap<String, Integer>());
			HashMap<String, Integer> smap3 = map3.get(key);
			if(!smap3.containsKey(v3)) smap3.put(v3, 0);
			smap3.put(v3, smap3.get(v3) + 1);
		}
		
		try {
			PrintStream out = new PrintStream(outFile);
			
			out.println("5p");
			for(String key : map5.keySet()){
				int sum = 0;
				for(String skey : map5.get(key).keySet()){
					out.println(key+"\t"+skey+"\t"+map5.get(key).get(skey));
					sum += map5.get(key).get(skey);
				}
				out.println("\t\t\t\t" + sum);
			}
			out.println("3p");
			for(String key : map3.keySet()){
				for(String skey : map3.get(key).keySet()){
					out.println(key+"\t"+skey+"\t"+map3.get(key).get(skey));
				}
				out.println();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void getScoredPairRegion(String file, String outFile){
		ScoredPairOutputParser parser = new ScoredPairOutputParser(file);
		HashMap<String, HashMap<String, Integer>> map5 = new HashMap<String, HashMap<String, Integer>>(); 
		HashMap<String, HashMap<String, Integer>> map3 = new HashMap<String, HashMap<String, Integer>>(); 
		
		for(ScoredPair p : parser.getPairs()){
			String key = p.getClassification();
			
			String r5 = p.getGenomicRegions3p() == null? "_" : p.getGenomicRegions3p().get(0);
			String v5 = "";
			if(r5.equals("_")) v5 = "InterGenic";
			else if(r5.endsWith("Intron")) v5 = "Intron";
			else if(r5.equals("NR")) v5 = "NR";
			else if(r5.endsWith("UTR")) v5 = "UTR";
			else if(r5.equals("LINC")) v5 = "LINC";
			else if(r5.endsWith("ORF")) v5 = "ORF";
			else v5 = r5;
			
			if(!map5.containsKey(key)) map5.put(key, new HashMap<String, Integer>());
			HashMap<String, Integer> smap5 = map5.get(key);
			if(!smap5.containsKey(v5)) smap5.put(v5, 0);
			smap5.put(v5, smap5.get(v5) + 1);
			
			String r3 = p.getGenomicRegions5p() == null? "_" : p.getGenomicRegions5p().get(0);
			String v3 = "";
			if(r3.equals("_")) v3 = "InterGenic";
			else if(r3.endsWith("Intron")) v3 = "Intron";
			else if(r3.equals("NR")) v3 = "NR";
			else if(r3.endsWith("UTR")) v3 = "UTR";
			else if(r3.equals("LINC")) v3 = "LINC";
			else if(r3.endsWith("ORF")) v3 = "ORF";
			else v3 = r3;
			
			if(!map3.containsKey(key)) map3.put(key, new HashMap<String, Integer>());
			HashMap<String, Integer> smap3 = map3.get(key);
			if(!smap3.containsKey(v3)) smap3.put(v3, 0);
			smap3.put(v3, smap3.get(v3) + 1);
		}
		
		try {
			PrintStream out = new PrintStream(outFile);
			
			out.println("5p");
			for(String key : map5.keySet()){
				int sum = 0;
				for(String skey : map5.get(key).keySet()){
					out.println(key+"\t"+skey+"\t"+map5.get(key).get(skey));
					sum += map5.get(key).get(skey);
				}
				out.println("\t\t\t\t" + sum);
			}
			out.println("3p");
			for(String key : map3.keySet()){
				for(String skey : map3.get(key).keySet()){
					out.println(key+"\t"+skey+"\t"+map3.get(key).get(skey));
				}
				out.println();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		String file = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.csv";
		String outFile = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.region.csv";
		String pairFile = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.pair.multihit.U.200.csv";
		String pairOutFile = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.pair.multihit.U.200.region.csv";
		
		GetRegionDistribution.getScoredPositionRegion(file, outFile);		
		GetRegionDistribution.getScoredPairRegion(pairFile, pairOutFile);
	}

}
