package fCLIP.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import launcher.ShellLauncher;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;
import parser.BufferedLineReader;

public class CheckRepeat {

	//bedtools intersect -S -split -wb -a /media/kyowon/Data1/Dropbox/fCLIP/new/h19x2.csv.5p.sorted.bed -b /media/kyowon/Data1/fCLIP/Data/cat.rmsk.bed -sorted > /media/kyowon/Data1/Dropbox/fCLIP/new/h19x2.rmsk.5p.bed
	//bedtools intersect -S -split -wb -a /media/kyowon/Data1/Dropbox/fCLIP/new/h19x2.csv.5p.bed -b /media/kyowon/Data1/Dropbox/fCLIP/new/h19x2.rmsk.5p.bed -sorted > /media/kyowon/Data1/Dropbox/fCLIP/new/h19x2.csv.5p.bed.rmsk.bed

	private static String getGenomicRegion(ScoredPosition p, boolean is5p, String suffix){
		ArrayList<String> rs = is5p? p.getGenomicRegions3p() : p.getGenomicRegions5p();
		String r = rs == null? "_" : rs.get(0);
		String k = "";
		if(r.equals("_")) k = "InterGenic";
		else if(r.endsWith("Intron")) k = "Intron";
		else if(r.equals("NR")) k = "NR";
		else if(r.endsWith("5_UTR")) k = "5_UTR";
		else if(r.endsWith("3_UTR")) k = "3_UTR";
		else if(r.equals("LINC")) k = "LINC";
		else if(r.endsWith("ORF")) k = "ORF";
		else k = r;
		
		if(!suffix.isEmpty()) suffix = "_" + suffix;
		return k  + suffix;
	}
	
	private static String getGenomicRegion(ScoredPair p, boolean is5p, String suffix){
		ArrayList<String> rs = is5p? p.getGenomicRegions3p() : p.getGenomicRegions5p();
		String r = rs == null? "_" : rs.get(0);
		String k = "";
		if(r.equals("_")) k = "InterGenic";
		else if(r.endsWith("Intron")) k = "Intron";
		else if(r.equals("NR")) k = "NR";
		else if(r.endsWith("5_UTR")) k = "5_UTR";
		else if(r.endsWith("3_UTR")) k = "3_UTR";
		else if(r.equals("LINC")) k = "LINC";
		else if(r.endsWith("ORF")) k = "ORF";
		else k = r;
		
		if(!suffix.isEmpty()) suffix = "_" + suffix;
		return k  + suffix;
	}
	
	
	public static void generateForTrans(String pairOut, String statOut, String pairCsv, String rmskBed5p, String rmskBed3p){
		HashMap<String, String> repeat5pMap = new HashMap<String, String>();
		HashMap<String, String> repeat3pMap = new HashMap<String, String>();
		HashMap<String, int[]> originEnd5pMap = new HashMap<String, int[]>();
		HashMap<String, int[]> originEnd3pMap = new HashMap<String, int[]>();
		try {
			BufferedLineReader in5 = new BufferedLineReader(rmskBed5p); // bed
			BufferedLineReader in3 = new BufferedLineReader(rmskBed3p);
			String s;
			
			while((s=in5.readLine())!=null){
				String[] token = s.split("\t");
				String key = token[0] + "\t" + token[3];
				repeat5pMap.put(key, token[10]);
				int p1 = Integer.parseInt(token[8]);
				int p2 = Integer.parseInt(token[9]);
				
				int[] v = new int[2];
				if(token[12].equals("+")){
					v[0] = p1; v[1] = p2;
				}else{
					v[1] = p1; v[0] = p2;
				}
				originEnd5pMap.put(key, v);
			}
		
			while((s=in3.readLine())!=null){
				String[] token = s.split("\t");
				String key = token[0] + "\t" + token[3];
				repeat3pMap.put(key, token[10]);
				int p1 = Integer.parseInt(token[8]);
				int p2 = Integer.parseInt(token[9]);
				
				int[] v = new int[2];
				if(token[12].equals("+")){
					v[0] = p1; v[1] = p2;
				}else{
					v[1] = p1; v[0] = p2;
				}
				originEnd3pMap.put(key, v);
			}
			
			in5.close();
			in3.close();
			
			HashMap<String, HashMap<String,Integer>> stat5p = new HashMap<String, HashMap<String,Integer>>();
			HashMap<String, HashMap<String,Integer>> stat3p = new HashMap<String, HashMap<String,Integer>>();
				
			HashMap<String, HashMap<String,Integer>> stat5pp = new HashMap<String, HashMap<String,Integer>>();
			HashMap<String, HashMap<String,Integer>> stat3pp = new HashMap<String, HashMap<String,Integer>>();
				
			HashSet<ScoredPosition> positions5p = new HashSet<ScoredPosition>();
			HashSet<ScoredPosition> positions3p = new HashSet<ScoredPosition>();
			
			ScoredPairOutputParser pparser = new ScoredPairOutputParser(pairCsv);
			PrintStream p2 = new PrintStream(pairOut);
			p2.println(pparser.getHeader() + "\t5p rmsk\t3p rmsk\t5p dist from origin\t5p dist from end\t3p dist from origin\t3p dist from end");
			for(ScoredPair p : pparser.getPairs()){
				//if(p.getMatched3pNum() < 20 && p.getMatched5pNum() < 20) continue;
				String key5 = p.getThreePContig() + "\t" + p.getThreePPosition();
				String key3 = p.getFivePContig() + "\t" + p.getFivePPosition();
				String v5 = repeat5pMap.get(key5);
				String v3 = repeat3pMap.get(key3);
				int[] p5 = originEnd5pMap.get(key5);
				int[] p3 = originEnd3pMap.get(key3);
				
				p2.println(p + "\t" + (v5 == null? "_" : v5) +  "\t" + (v3 == null? "_" : v3) 
						+ "\t" + (p5 == null? "" : p5[0] - p.getThreePPosition()) + "\t" + (p5 == null? "" : p5[1] - p.getThreePPosition())
						+ "\t" + (p3 == null? "" : p3[0] - p.getFivePPosition()) + "\t" + (p3 == null? "" : p3[1] - p.getFivePPosition()));
				
			
				String mKey = p.getClassification();
				
				if(!stat5p.containsKey(mKey)) stat5p.put(mKey, new HashMap<String, Integer>());
				if(!stat3p.containsKey(mKey)) stat3p.put(mKey, new HashMap<String, Integer>());	
				
				HashMap<String, Integer> mv5p = stat5p.get(mKey);
				String suffix5 = "";
				if(v5 != null){
					if(v5.startsWith("SINE")) suffix5 = "SINE";
					else if(v5.startsWith("LINE")) suffix5 = "LINE";
					else suffix5 = "R";
				}
				String k5 = getGenomicRegion(p, true, suffix5);
				if(!mv5p.containsKey(k5)) mv5p.put(k5, 0);
				mv5p.put(k5, mv5p.get(k5)+1);
				
				if(!positions5p.contains(p.getPairedScoredPositions()[0])){
					positions5p.add(p.getPairedScoredPositions()[0]);
					if(!stat5pp.containsKey(mKey)) stat5pp.put(mKey, new HashMap<String, Integer>());
					HashMap<String, Integer> mv5pp = stat5pp.get(mKey);				
					
					if(!mv5pp.containsKey(k5)) mv5pp.put(k5, 0);
					mv5pp.put(k5, mv5pp.get(k5)+1);
				}
				HashMap<String, Integer> mv3p = stat3p.get(mKey);
				
				String suffix3 = "";
				if(v3 != null){
					if(v3.startsWith("SINE")) suffix3 = "SINE";
					else if(v3.startsWith("LINE")) suffix3 = "LINE";
					else suffix3 = "R";
				}
				String k3 = getGenomicRegion(p, false, suffix3);
				if(!mv3p.containsKey(k3)) mv3p.put(k3, 0);
				mv3p.put(k3, mv3p.get(k3)+1);	
				
				if(!positions3p.contains(p.getPairedScoredPositions()[1])){
					positions3p.add(p.getPairedScoredPositions()[1]);
					if(!stat3pp.containsKey(mKey)) stat3pp.put(mKey, new HashMap<String, Integer>());
					HashMap<String, Integer> mv3pp = stat3pp.get(mKey);				
					
					if(!mv3pp.containsKey(k3)) mv3pp.put(k3, 0);
					mv3pp.put(k3, mv3pp.get(k3)+1);
				}
			}
				
			p2.close();
			PrintStream statp = new PrintStream(statOut);
			
			statp.println("Trans duplexes\n");
			for(String mKey : stat5p.keySet()){
				HashMap<String, Integer> mv5p = stat5p.get(mKey);
				statp.println(mKey + " 5p");
				int sum = 0;
				for(String k5 : mv5p.keySet()){
					sum += mv5p.get(k5);
				}
				
				for(String k5 : mv5p.keySet()){
					statp.println(k5 + "\t" + mv5p.get(k5) + "\t" + ((double)mv5p.get(k5) * 100 / sum));
				}
				//	
				statp.println("total\t" + sum + "\n\n");
			}
			
			for(String mKey : stat3p.keySet()){
				HashMap<String, Integer> mv3p = stat3p.get(mKey);
				statp.println(mKey + " 3p");
				int sum = 0;
				for(String k3 : mv3p.keySet()){
					sum += mv3p.get(k3);
				}
				
				for(String k3 : mv3p.keySet()){
					statp.println(k3 + "\t" + mv3p.get(k3) + "\t" + ((double)mv3p.get(k3) * 100 / sum));
				}
				//	
				statp.println("total\t" + sum + "\n\n");
			}
			
			statp.println("Paired Positions\n");
			for(String mKey : stat5pp.keySet()){
				HashMap<String, Integer> mv5p = stat5pp.get(mKey);
				statp.println(mKey + " 5p");
				int sum = 0;
				for(String k5 : mv5p.keySet()){
					sum += mv5p.get(k5);
				}
				
				for(String k5 : mv5p.keySet()){
					statp.println(k5 + "\t" + mv5p.get(k5) + "\t" + ((double)mv5p.get(k5) * 100 / sum));
				}
				//	
				statp.println("total\t" + sum + "\n\n");
			}
			
			for(String mKey : stat3pp.keySet()){
				HashMap<String, Integer> mv3p = stat3pp.get(mKey);
				statp.println(mKey + " 3p");
				int sum = 0;
				for(String k3 : mv3p.keySet()){
					sum += mv3p.get(k3);
				}
				
				for(String k3 : mv3p.keySet()){
					statp.println(k3 + "\t" + mv3p.get(k3) + "\t" + ((double)mv3p.get(k3) * 100 / sum));
				}
				//	
				statp.println("total\t" + sum + "\n\n");
			}
			
			statp.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void generateForCis(String cisOut, String statOut, String cisCsv, String rmskBed5p, String rmskBed3p){
		HashMap<String, String> repeat5pMap = new HashMap<String, String>();
		HashMap<String, String> repeat3pMap = new HashMap<String, String>();
		HashMap<String, int[]> originEnd5pMap = new HashMap<String, int[]>();
		HashMap<String, int[]> originEnd3pMap = new HashMap<String, int[]>();
		try {
			BufferedLineReader in5 = new BufferedLineReader(rmskBed5p); // bed
			BufferedLineReader in3 = new BufferedLineReader(rmskBed3p);
			String s;
			
			while((s=in5.readLine())!=null){
				String[] token = s.split("\t");
				String key = token[0] + "\t" + token[3];
				repeat5pMap.put(key, token[10]);
				int p1 = Integer.parseInt(token[8]);
				int p2 = Integer.parseInt(token[9]);
				
				int[] v = new int[2];
				if(token[12].equals("+")){
					v[0] = p1; v[1] = p2;
				}else{
					v[1] = p1; v[0] = p2;
				}
				originEnd5pMap.put(key, v);
			}
		
			while((s=in3.readLine())!=null){
				String[] token = s.split("\t");
				String key = token[0] + "\t" + token[3];
				repeat3pMap.put(key, token[10]);
				int p1 = Integer.parseInt(token[8]);
				int p2 = Integer.parseInt(token[9]);
				
				int[] v = new int[2];
				if(token[12].equals("+")){
					v[0] = p1; v[1] = p2;
				}else{
					v[1] = p1; v[0] = p2;
				}
				originEnd3pMap.put(key, v);
			}
			
			in5.close();
			in3.close();
			
			ScoredPositionOutputParser parser = new ScoredPositionOutputParser(cisCsv);
			PrintStream p1 = new PrintStream(cisOut);
			
			p1.println(parser.getHeader() + "\t5p rmsk\t3p rmsk\t5p dist from origin\t5p dist from end\t3p dist from origin\t3p dist from end");
			HashMap<String, HashMap<String, Integer>> stat5p = new HashMap<String, HashMap<String, Integer>>();
			HashMap<String, HashMap<String, Integer>> stat3p = new HashMap<String, HashMap<String, Integer>>();
			
			for(ScoredPosition p : parser.getPositions()){
				String key5 = p.getContig() + "\t" + p.getThreePposition();
				String key3 = p.getContig() + "\t" + p.getFivePposition();
				String v5 = repeat5pMap.get(key5);
				String v3 = repeat3pMap.get(key3);
				int[] p5 = originEnd5pMap.get(key5);
				int[] p3 = originEnd3pMap.get(key3);
				
				p1.println(p + "\t" + (v5 == null? "_" : v5) +  "\t" + (v3 == null? "_" : v3) 
						+ "\t" + (p5 == null? "" : p5[0] - p.getThreePposition()) + "\t" + (p5 == null? "" : p5[1] - p.getThreePposition())
						+ "\t" + (p3 == null? "" : p3[0] - p.getFivePposition()) + "\t" + (p3 == null? "" : p3[1] - p.getFivePposition()));
				
				String mKey = "";
				
				if(p.hasMatchingMiRNA()){
					mKey += "MiRNA";
				}else{
					mKey += p.getClassification() + (p.isPaired() ? "T" : "F");
				}
				if(!stat5p.containsKey(mKey)) stat5p.put(mKey, new HashMap<String, Integer>());
				if(!stat3p.containsKey(mKey)) stat3p.put(mKey, new HashMap<String, Integer>());
				
				HashMap<String, Integer> mv5p = stat5p.get(mKey);
				
				String suffix5 = "";
				if(v5 != null){
					if(v5.startsWith("SINE")) suffix5 = "SINE";
					else if(v5.startsWith("LINE")) suffix5 = "LINE";
					else suffix5 = "R";
				}
				String k5 = getGenomicRegion(p, true, suffix5);
				if(!mv5p.containsKey(k5)) mv5p.put(k5, 0);
				mv5p.put(k5, mv5p.get(k5)+1);
				
				
				HashMap<String, Integer> mv3p = stat3p.get(mKey);
				
				String suffix3 = "";
				if(v3 != null){
					if(v3.startsWith("SINE")) suffix3 = "SINE";
					else if(v3.startsWith("LINE")) suffix3 = "LINE";
					else suffix3 = "R";
				}
				String k3 = getGenomicRegion(p, false, suffix3);
				if(!mv3p.containsKey(k3)) mv3p.put(k3, 0);
				mv3p.put(k3, mv3p.get(k3)+1);						
			}
			p1.close();
			PrintStream statp = new PrintStream(statOut);
			
			statp.println("CIS\n");
			for(String mKey : stat5p.keySet()){
				HashMap<String, Integer> mv5p = stat5p.get(mKey);
				statp.println(mKey + " 5p");
				int sum = 0;
				for(String k5 : mv5p.keySet()){
					sum += mv5p.get(k5);
				}
				
				for(String k5 : mv5p.keySet()){
					statp.println(k5 + "\t" + mv5p.get(k5) + "\t" + ((double)mv5p.get(k5) * 100 / sum));
				}
				//	
				statp.println("total\t" + sum + "\n\n");
			}
			
			for(String mKey : stat3p.keySet()){
				HashMap<String, Integer> mv3p = stat3p.get(mKey);
				statp.println(mKey + " 3p");
				int sum = 0;
				for(String k3 : mv3p.keySet()){
					sum += mv3p.get(k3);
				}
				
				for(String k3 : mv3p.keySet()){
					statp.println(k3 + "\t" + mv3p.get(k3) + "\t" + ((double)mv3p.get(k3) * 100 / sum));
				}
				//	
				statp.println("total\t" + sum + "\n\n");
			}
			statp.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	public static void generate(String csv, String pairCsv, String bed5p, String bed3p, String rmskBed) throws IOException {
		HashMap<String, String> repeat5pMap = new HashMap<String, String>();
		HashMap<String, String> repeat3pMap = new HashMap<String, String>();
		HashMap<String, int[]> originEnd5pMap = new HashMap<String, int[]>();
		HashMap<String, int[]> originEnd3pMap = new HashMap<String, int[]>();
		
		if(!new File(bed5p + ".rmsk.bed").exists() || !new File(bed3p + ".rmsk.bed").exists()){
			ShellLauncher.run("bedtools intersect -split -wb -a " + bed5p+" -b " + rmskBed + " -sorted > " + bed5p + ".rmsk.bed");
			ShellLauncher.run("bedtools intersect -split -wb -a " + bed3p+" -b " + rmskBed + " -sorted > " + bed3p + ".rmsk.bed");
		}
		
		String sfileOut = csv + ".rmsk.csv";	
		String pairFileOut = pairCsv + ".rmsk.csv";
		
		String statOut = pairCsv + ".regionStat.csv";
		
		BufferedLineReader in5 = new BufferedLineReader(bed5p + ".rmsk.bed"); // bed
		BufferedLineReader in3 = new BufferedLineReader(bed3p + ".rmsk.bed");
		String s;
		
		while((s=in5.readLine())!=null){
			String[] token = s.split("\t");
			String key = token[0] + "\t" + token[3];
			repeat5pMap.put(key, token[9]);
			int p1 = Integer.parseInt(token[7]);
			int p2 = Integer.parseInt(token[8]);
			
			int[] v = new int[2];
			if(token[11].equals("+")){
				v[0] = p1; v[1] = p2;
			}else{
				v[1] = p1; v[0] = p2;
			}
			originEnd5pMap.put(key, v);
		}
	
		while((s=in3.readLine())!=null){
			String[] token = s.split("\t");
			String key = token[0] + "\t" + token[3];
			repeat3pMap.put(key, token[9]);
			int p1 = Integer.parseInt(token[7]);
			int p2 = Integer.parseInt(token[8]);
			
			int[] v = new int[2];
			if(token[11].equals("+")){
				v[0] = p1; v[1] = p2;
			}else{
				v[1] = p1; v[0] = p2;
			}
			originEnd3pMap.put(key, v);
		}
		
		in5.close();
		in3.close();
		
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		PrintStream p1 = new PrintStream(sfileOut);
		
		p1.println(ScoredPosition.getHeader() + "\t5p rmsk\t3p rmsk\t5p dist from origin\t5p dist from end\t3p dist from origin\t3p dist from end");
		HashMap<String, HashMap<String, Integer>> stat5p = new HashMap<String, HashMap<String, Integer>>();
		HashMap<String, HashMap<String, Integer>> stat3p = new HashMap<String, HashMap<String, Integer>>();
		
		for(ScoredPosition p : parser.getPositions()){
			String key5 = p.getContig() + "\t" + p.getThreePposition();
			String key3 = p.getContig() + "\t" + p.getFivePposition();
			String v5 = repeat5pMap.get(key5);
			String v3 = repeat3pMap.get(key3);
			int[] p5 = originEnd5pMap.get(key5);
			int[] p3 = originEnd3pMap.get(key3);
			
			p1.println(p + "\t" + (v5 == null? "_" : v5) +  "\t" + (v3 == null? "_" : v3) 
					+ "\t" + (p5 == null? "" : p5[0] - p.getThreePposition()) + "\t" + (p5 == null? "" : p5[1] - p.getThreePposition())
					+ "\t" + (p3 == null? "" : p3[0] - p.getFivePposition()) + "\t" + (p3 == null? "" : p3[1] - p.getFivePposition()));
			
			String mKey = "";
			
			if(p.hasMatchingMiRNA()){
				mKey += "MiRNA";
			}else{
				mKey += p.getClassification() + (p.isPaired() ? "T" : "F");
			}
			if(!stat5p.containsKey(mKey)) stat5p.put(mKey, new HashMap<String, Integer>());
			if(!stat3p.containsKey(mKey)) stat3p.put(mKey, new HashMap<String, Integer>());
			
			HashMap<String, Integer> mv5p = stat5p.get(mKey);
			
			String suffix5 = "";
			if(v5 != null){
				if(v5.startsWith("SINE")) suffix5 = "SINE";
				else if(v5.startsWith("LINE")) suffix5 = "LINE";
				else suffix5 = "R";
			}
			String k5 = getGenomicRegion(p, true, suffix5);
			if(!mv5p.containsKey(k5)) mv5p.put(k5, 0);
			mv5p.put(k5, mv5p.get(k5)+1);
			
			
			HashMap<String, Integer> mv3p = stat3p.get(mKey);
			
			String suffix3 = "";
			if(v3 != null){
				if(v3.startsWith("SINE")) suffix3 = "SINE";
				else if(v3.startsWith("LINE")) suffix3 = "LINE";
				else suffix3 = "R";
			}
			String k3 = getGenomicRegion(p, false, suffix3);
			if(!mv3p.containsKey(k3)) mv3p.put(k3, 0);
			mv3p.put(k3, mv3p.get(k3)+1);						
		}
		p1.close();
		PrintStream statp = new PrintStream(statOut);
		
		statp.println("CIS\n");
		for(String mKey : stat5p.keySet()){
			HashMap<String, Integer> mv5p = stat5p.get(mKey);
			statp.println(mKey + " 5p");
			int sum = 0;
			for(String k5 : mv5p.keySet()){
				sum += mv5p.get(k5);
			}
			
			for(String k5 : mv5p.keySet()){
				statp.println(k5 + "\t" + mv5p.get(k5) + "\t" + ((double)mv5p.get(k5) * 100 / sum));
			}
			//	
			statp.println("total\t" + sum + "\n\n");
		}
		
		for(String mKey : stat3p.keySet()){
			HashMap<String, Integer> mv3p = stat3p.get(mKey);
			statp.println(mKey + " 3p");
			int sum = 0;
			for(String k3 : mv3p.keySet()){
				sum += mv3p.get(k3);
			}
			
			for(String k3 : mv3p.keySet()){
				statp.println(k3 + "\t" + mv3p.get(k3) + "\t" + ((double)mv3p.get(k3) * 100 / sum));
			}
			//	
			statp.println("total\t" + sum + "\n\n");
		}
		
		
		stat5p = new HashMap<String, HashMap<String,Integer>>();
		stat3p = new HashMap<String, HashMap<String,Integer>>();
				
		ScoredPairOutputParser pparser = new ScoredPairOutputParser(pairCsv);
		PrintStream p2 = new PrintStream(pairFileOut);
		p2.println(ScoredPair.getHeader() + "\t5p rmsk\t3p rmsk\t5p dist from origin\t5p dist from end\t3p dist from origin\t3p dist from end");
		for(ScoredPair p : pparser.getPairs()){
			//if(p.getMatched3pNum() < 20 && p.getMatched5pNum() < 20) continue;
			String key5 = p.getThreePContig() + "\t" + p.getThreePPosition();
			String key3 = p.getFivePContig() + "\t" + p.getFivePPosition();
			String v5 = repeat5pMap.get(key5);
			String v3 = repeat3pMap.get(key3);
			int[] p5 = originEnd5pMap.get(key5);
			int[] p3 = originEnd3pMap.get(key3);
			
			p2.println(p + "\t" + (v5 == null? "_" : v5) +  "\t" + (v3 == null? "_" : v3) 
					+ "\t" + (p5 == null? "" : p5[0] - p.getThreePPosition()) + "\t" + (p5 == null? "" : p5[1] - p.getThreePPosition())
					+ "\t" + (p3 == null? "" : p3[0] - p.getFivePPosition()) + "\t" + (p3 == null? "" : p3[1] - p.getFivePPosition()));
			
			
		//	p2.println(p + "\t" + (v5 == null? "_" : v5) +  "\t" + (v3 == null? "_" : v3));
			
			String mKey = p.getClassification();
			
			if(!stat5p.containsKey(mKey)) stat5p.put(mKey, new HashMap<String, Integer>());
			if(!stat3p.containsKey(mKey)) stat3p.put(mKey, new HashMap<String, Integer>());
			
			HashMap<String, Integer> mv5p = stat5p.get(mKey);
			
			String suffix5 = "";
			if(v5 != null){
				if(v5.startsWith("SINE")) suffix5 = "SINE";
				else if(v5.startsWith("LINE")) suffix5 = "LINE";
				else suffix5 = "R";
			}
			String k5 = getGenomicRegion(p, true, suffix5);
			if(!mv5p.containsKey(k5)) mv5p.put(k5, 0);
			mv5p.put(k5, mv5p.get(k5)+1);
			
			
			HashMap<String, Integer> mv3p = stat3p.get(mKey);
			
			String suffix3 = "";
			if(v3 != null){
				if(v3.startsWith("SINE")) suffix3 = "SINE";
				else if(v3.startsWith("LINE")) suffix3 = "LINE";
				else suffix3 = "R";
			}
			String k3 = getGenomicRegion(p, false, suffix3);
			if(!mv3p.containsKey(k3)) mv3p.put(k3, 0);
			mv3p.put(k3, mv3p.get(k3)+1);	
		}
			
		p2.close();
		
		statp.println("Trans\n");
		for(String mKey : stat5p.keySet()){
			HashMap<String, Integer> mv5p = stat5p.get(mKey);
			statp.println(mKey + " 5p");
			int sum = 0;
			for(String k5 : mv5p.keySet()){
				sum += mv5p.get(k5);
			}
			
			for(String k5 : mv5p.keySet()){
				statp.println(k5 + "\t" + mv5p.get(k5) + "\t" + ((double)mv5p.get(k5) * 100 / sum));
			}
			//	
			statp.println("total\t" + sum + "\n\n");
		}
		
		for(String mKey : stat3p.keySet()){
			HashMap<String, Integer> mv3p = stat3p.get(mKey);
			statp.println(mKey + " 3p");
			int sum = 0;
			for(String k3 : mv3p.keySet()){
				sum += mv3p.get(k3);
			}
			
			for(String k3 : mv3p.keySet()){
				statp.println(k3 + "\t" + mv3p.get(k3) + "\t" + ((double)mv3p.get(k3) * 100 / sum));
			}
			//	
			statp.println("total\t" + sum + "\n\n");
		}
		
		statp.close();
		
	}
	
	public static void main(String[] args) {	
		String out = args[0];
		String stat = args[1];
		String csv = args[2];
		String rmskBed5p = args[3];
		String rmskBed3p = args[4];
		if(args[5].equals("cis")){			
			CheckRepeat.generateForCis(out, stat, csv, rmskBed5p, rmskBed3p);
		}
		if(args[5].equals("trans")){
			CheckRepeat.generateForTrans(out, stat, csv, rmskBed5p, rmskBed3p);
		}
	}
}
