package fCLIP.analysis;


import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class DivideTransDuplexResults {

	public static void main(String[] args) {
		String inFile = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/trans/merged.trans.complete.M.csv";
		String outFile = inFile.replace('.', '_') + "_alu.m";
		
		ScoredPairOutputParser parser = new ScoredPairOutputParser(inFile); 
		HashMap<String, HashMap<ScoredPosition, int[]>> map3p = new HashMap<String, HashMap<ScoredPosition, int[]>>();
		HashMap<String, HashMap<ScoredPosition, int[]>> map3pinverted = new HashMap<String, HashMap<ScoredPosition, int[]>>();
		HashMap<String, HashMap<ScoredPosition, int[]>> map5p = new HashMap<String, HashMap<ScoredPosition, int[]>>();
		HashMap<String, HashMap<ScoredPosition, int[]>> map5pinverted = new HashMap<String, HashMap<ScoredPosition, int[]>>();
		
		String[] htoken = parser.getHeader().split("\t");
		int index = 0;
		for(int i=0;i<htoken.length;i++){
			if(htoken[i].equals("5p rmsk")){
				index = i;
				break;
			}
		}
	//	System.out.println("1\t2\t\t\t".split("\t").length);
	//	System.out.println(index);
		index -= ScoredPair.getHeader().split("\t").length;
	
		for(ScoredPair pair : parser.getPairs()){
			if(pair.getOverHang() !=2) continue;
			String[] token = pair.getMisc().split("\t");
			String rmsk5p = token[index];
			//System.out.println(pair.getMisc());
			
			String rmsk3p = token[index+1];
			int[] dist5p = new int[2];
			int[] dist3p = new int[2];
			dist5p[0] = token[index+2].isEmpty() || token[index+2].equals(" ")? 0 : Integer.parseInt(token[index+2]);
			dist5p[1] = token[index+3].isEmpty() || token[index+3].equals(" ")? 0 : Integer.parseInt(token[index+3]);
			dist3p[0] = token[index+4].isEmpty() || token[index+4].equals(" ")? 0 : Integer.parseInt(token[index+4]);
			dist3p[1] = token[index+5].isEmpty() || token[index+5].equals(" ")? 0 : Integer.parseInt(token[index+5]);
			
			boolean ts = pair.getFivePStrand(); 
			boolean fs = pair.getThreePStrand(); // read three p thus five p in duplex
			
			boolean alufs = dist5p[1] > dist5p[0];
			boolean aluts = dist3p[1] > dist3p[0];
			
			if(alufs){
				dist5p[0] = -dist5p[0];
				dist5p[1] = -dist5p[1];
			}
			if(aluts){
				dist3p[0] = -dist3p[0];
				dist3p[1] = -dist3p[1];
			}
			boolean isInverted5p = alufs == fs;
			boolean isInverted3p = aluts == ts;
			
			if(rmsk5p.startsWith("_")){
				rmsk5p = "NoRepeat";
			}else{
				//System.out.println(rmsk5p);
				String[] rt = rmsk5p.split("\\|");
				rmsk5p = rt[0]+"_"+rt[1].substring(0, Math.min(rt[1].length(), 4));
				rmsk5p = rmsk5p.replace('(', '_');
				rmsk5p = rmsk5p.replace(')', '_');
				rmsk5p = rmsk5p.replace('-', '_');
				//System.out.println(rmsk5p);
			}
			
			if(rmsk3p.startsWith("_")){
				rmsk3p = "NoRepeat";
			}else{
				String[] rt = rmsk3p.split("\\|");
				rmsk3p = rt[0]+"_"+rt[1].substring(0, Math.min(rt[1].length(), 4));
				rmsk3p = rmsk3p.replace('(', '_');
				rmsk3p = rmsk3p.replace(')', '_');
				rmsk3p = rmsk3p.replace('-', '_');
				
			}
			
			ScoredPosition s5 = pair.getPairedScoredPositions()[0];
			ScoredPosition s3 = pair.getPairedScoredPositions()[1];
			
			if(isInverted5p){
				String prefix = "fivePInv_";
				if(!map5pinverted.containsKey(prefix + rmsk5p)) map5pinverted.put(prefix + rmsk5p, new HashMap<ScoredPosition, int[]>());
				HashMap<ScoredPosition, int[]> sm = map5pinverted.get(prefix + rmsk5p);
				if(!sm.containsKey(s5)) sm.put(s5, dist5p);
			}else{
				String prefix = "fiveP_";
				if(!map5p.containsKey(prefix + rmsk5p)) map5p.put(prefix + rmsk5p, new HashMap<ScoredPosition, int[]>());
				HashMap<ScoredPosition, int[]> sm = map5p.get(prefix + rmsk5p);
				if(!sm.containsKey(s5)) sm.put(s5, dist5p);
				
			}
			
			if(isInverted3p){
				String prefix = "threePInv_";
				if(!map3pinverted.containsKey(prefix + rmsk3p)) map3pinverted.put(prefix + rmsk3p, new HashMap<ScoredPosition, int[]>());
				HashMap<ScoredPosition, int[]> sm = map3pinverted .get(prefix + rmsk3p);
				if(!sm.containsKey(s3)) sm.put(s3, dist3p);
			}else{
				String prefix = "threeP_";				
				if(!map3p.containsKey(prefix + rmsk3p)) map3p.put(prefix + rmsk3p, new HashMap<ScoredPosition, int[]>());
				HashMap<ScoredPosition, int[]> sm = map3p.get(prefix + rmsk3p);
				if(!sm.containsKey(s3)) sm.put(s3, dist3p);
			}
		}
		
		try {
			PrintStream out = new PrintStream(outFile);
			
			for(String k : map5p.keySet()){
				out.println(k+"_from_origin=[");
				for(int[] d : map5p.get(k).values()){
					out.print(d[0] + " ");
				}
				out.println("];");
			}
			
			for(String k : map5p.keySet()){
				out.println(k+"_from_end=[");
				for(int[] d : map5p.get(k).values()){
					out.print(d[1] + " ");
				}
				out.println("];");
			}
			
			
			for(String k : map3p.keySet()){
				out.println(k+"_from_origin=[");
				for(int[] d : map3p.get(k).values()){
					out.print(d[0] + " ");
				}
				out.println("];");
			}
			
			for(String k : map3p.keySet()){
				out.println(k+"_from_end=[");
				for(int[] d : map3p.get(k).values()){
					out.print(d[1] + " ");
				}
				out.println("];");
			}
			
			for(String k : map5pinverted.keySet()){
				out.println(k+"_from_origin=[");
				for(int[] d : map5pinverted.get(k).values()){
					out.print(d[0] + " ");
				}
				out.println("];");
			}
			
			for(String k : map5pinverted.keySet()){
				out.println(k+"_from_end=[");
				for(int[] d : map5pinverted.get(k).values()){
					out.print(d[1] + " ");
				}
				out.println("];");
			}
			
			for(String k : map3pinverted.keySet()){
				out.println(k+"_from_origin=[");
				for(int[] d : map3pinverted.get(k).values()){
					out.print(d[0] + " ");
				}
				out.println("];");
			}
			
			for(String k : map3pinverted.keySet()){
				out.println(k+"_from_end=[");
				for(int[] d : map3pinverted.get(k).values()){
					out.print(d[1] + " ");
				}
				out.println("];");
			}
			
			
			
			
			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

}










