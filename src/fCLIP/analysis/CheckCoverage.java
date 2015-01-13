package fCLIP.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.HashSet;
import launcher.ShellLauncher;
import parser.BufferedLineReader;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class CheckCoverage {
	public static class MapRunner implements Runnable {
		private String[] beds; 
		private ArrayList<Hashtable<String, Double>>  readMaps; 
		private Hashtable<String, Integer>  regionMap; 
		
		private int maxSpan;
		
		MapRunner(ArrayList<Hashtable<String, Double>> readMaps, Hashtable<String, Integer> regionMap, String[] beds, int maxSpan){
			this.readMaps = readMaps;
			this.beds = beds;
			this.maxSpan = maxSpan;		
			this.regionMap = regionMap;
		}
		
		MapRunner(ArrayList<Hashtable<String, Double>> readMaps, Hashtable<String, Integer> regionMap,  ArrayList<String> beds, int maxSpan){
			this.readMaps = readMaps;
			this.beds = new String[beds.size()];
			for(int i=0;i<beds.size();i++) this.beds[i] = beds.get(i);
			this.maxSpan = maxSpan;
			this.regionMap = regionMap;
		}
		
		public void run() {
			for(int i=0;i< beds.length;i++){
				String s;				
				Hashtable<String, Double> readMap = new Hashtable<String, Double>(); 
				BufferedLineReader tin;
				try {
					tin = new BufferedLineReader(beds[i]);
				
				while((s=tin.readLine())!=null){
					String[] token = s.split("\t");
					if(Integer.parseInt(token[4]) - Integer.parseInt(token[2]) > maxSpan) continue;
					if(Integer.parseInt(token[1]) == Integer.parseInt(token[12]) - 1) continue;
					if(Integer.parseInt(token[5]) - 1 == Integer.parseInt(token[11])) continue;
					StringBuilder sb = new StringBuilder();
					sb.append(token[10]); sb.append('\t');
					sb.append(token[13]); sb.append('\t');
					sb.append(token[15].equals("+"));
					
					String k = sb.toString();
					
					Double n = readMap.get(k);
					if(n == null)
						n =  0.0;
					readMap.put(k, n + 1.0);						
					regionMap.put(k, Integer.parseInt(token[16]));					
				}
				readMaps.add(readMap);
				tin.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}			
		}		
	}
	
	public static void generateForBackground(String mOut, String[] intersectedBeds, String[] keys, int cutThreshold, int maxSpan, int numThreads){
		try {	
			ArrayList<Hashtable<String, Double>> readCntMaps = new ArrayList<Hashtable<String, Double>>();
			ArrayList<ArrayList<Hashtable<String, Double>>> treadCntMaps = new ArrayList<ArrayList<Hashtable<String, Double>>>();
			
			Hashtable<String, Integer> regionMap = new Hashtable<String, Integer>();
			//ArrayList<ArrayList<Hashtable<String, Integer>>> tregionMaps = new ArrayList<ArrayList<Hashtable<String, Integer>>>();
					
			
			ArrayList<Thread> threads = new ArrayList<Thread>();
			
			for(int i=0;i<numThreads;i++){
				treadCntMaps.add(new ArrayList<Hashtable<String, Double>>());
			//	tregionMaps.add(new ArrayList<Hashtable<String,Integer>>());
				
				ArrayList<String> beds = new ArrayList<String>();
				
				for(int j=0;j<intersectedBeds.length;j++){
					if(j%numThreads == i) beds.add(intersectedBeds[j]);
				}
				MapRunner runner = new MapRunner(treadCntMaps.get(treadCntMaps.size() - 1), regionMap, beds, maxSpan);
				Thread thread = new Thread(runner, "" + i);      
				thread.start();				
				threads.add(thread);
			}
			
			try {
				for(Thread thread : threads){
					thread.join();
				}
			} catch (InterruptedException e) {						
				e.printStackTrace();
			}
			
			HashSet<String> accs = new HashSet<String>();
			
			for(int i=0;i<=intersectedBeds.length/numThreads;i++){
				for(int j=0;j<treadCntMaps.size();j++){
					readCntMaps.add(treadCntMaps.get(j).get(i));
					//regionMap.add(tregionMaps.get(j).get(i));
					if(readCntMaps.size() == intersectedBeds.length) break;					
				}
				if(readCntMaps.size() == intersectedBeds.length) break;
			}
			
			for(Hashtable<String, Double> cr : readCntMaps){
				accs.addAll(cr.keySet());
			}
			
			PrintStream outm = new PrintStream(mOut);
			Hashtable<String, StringBuilder> mStringMap = new Hashtable<String, StringBuilder>();
			Hashtable<String, ArrayList<Double>> ctm = new Hashtable<String, ArrayList<Double>>();
			
			for(String acc : accs){
				for(int i=0; i<keys.length;i++){
					Double c5 = readCntMaps.get(i).get(acc);
					Integer regionID = regionMap.get(acc);
					if(regionID == null) regionID = 0;
					String mKey = "BG_" + new String(keys[i])+"_"+regionID;
					
					if(c5 == null) c5 = .0;
					if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
					if(!ctm.containsKey(keys[i]+"_"+regionID)) ctm.put(keys[i]+"_"+regionID, new ArrayList<Double>());
					
					mStringMap.get(mKey).append((c5) + " ");
					ctm.get(keys[i]+"_"+regionID).add(c5);
				}
			}
		
			
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("];");
				//outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=20 & " + mkey + "(:,2)>= 20),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			
			for(String key1 : ctm.keySet()){
				for(String key2 : ctm.keySet()){
					if(key1.equals(key2)) continue;
					if(!key1.substring(key1.length()-2).equals(key2.substring(key2.length()-2))) continue;
					outm.print("%"+key1 + ";"+key2+";");
					outm.println(getMedianFC(ctm.get(key1), ctm.get(key2), cutThreshold));
				}
			}
			
			outm.close();
			
			
		}catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	private static double getMedianFC(ArrayList<Double> c, ArrayList<Double> t, int th){
		ArrayList<Double> fcs = new ArrayList<Double>();
		for(int i=0;i<c.size();i++){
			double n = c.get(i);
			double m = t.get(i);
			if(n + m <= th) continue;
			
			fcs.add(Math.log((m+1)/(n+1))/Math.log(2));
		}
		Collections.sort(fcs);
		return fcs.get(fcs.size()/2);
	}
	
	public static void generateForCis(String cisOut, String mOut, String cisCsv, 
			String[] cBeds5p, String[] cBeds3p, String[] tBeds5p, String[] tBeds3p, 
			String[] cKeys, String[] tKeys, int maxSpan, int fcThreshold, Hashtable<String, Double> medians){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(cisCsv);
		ArrayList<Hashtable<String, Double>> creadCnt3pMaps = new ArrayList<Hashtable<String, Double>>();
		ArrayList<Hashtable<String, Double>> creadCnt5pMaps = new ArrayList<Hashtable<String, Double>>();
		ArrayList<Hashtable<String, Double>> treadCnt3pMaps = new ArrayList<Hashtable<String, Double>>();
		ArrayList<Hashtable<String, Double>> treadCnt5pMaps = new ArrayList<Hashtable<String, Double>>();
		
		Hashtable<String, Integer> region3pMap = new Hashtable<String, Integer>();
		Hashtable<String, Integer> region5pMap = new Hashtable<String, Integer>();
			
		try {			
			//System.out.println(cKeys.length +  " " + tKeys.length);
			PrintStream out = new PrintStream(cisOut);
			out.print(parser.getHeader());
			for(int i=0; i<cKeys.length; i++){
				String ckey = cKeys[i];
				out.print("\t" + ckey);
				for(int j=0; i + j*cKeys.length<tKeys.length; j++){
					//System.out.println(i + j*cKeys.length);
					String tkey = tKeys[i + j*cKeys.length];
					out.print("\t" + tkey + "\t" + tkey + "_FC");
				}
			}
			//System.out.println(out);
			out.println();
			MapRunner runner1 = new MapRunner(creadCnt3pMaps, region3pMap, cBeds3p, maxSpan);
			MapRunner runner2 = new MapRunner(creadCnt5pMaps, region5pMap, cBeds5p, maxSpan);
			MapRunner runner3 = new MapRunner(treadCnt3pMaps, region3pMap, tBeds3p, maxSpan);
			MapRunner runner4 = new MapRunner(treadCnt5pMaps, region5pMap, tBeds5p, maxSpan);
			
			Thread thread1 = new Thread(runner1, "" + 1);       
			Thread thread2 = new Thread(runner2, "" + 2);
			Thread thread3 = new Thread(runner3, "" + 3);
			Thread thread4 = new Thread(runner4, "" + 4);
			
			thread1.start();
			thread2.start();
			thread3.start();
			thread4.start();
			
			try {
				thread1.join();
				thread2.join();
				thread3.join();
				thread4.join();
			} catch (InterruptedException e) {						
				e.printStackTrace();
			}
					
			PrintStream outm = new PrintStream(mOut);
			Hashtable<String, StringBuilder> mStringMap = new Hashtable<String, StringBuilder>();
			
			for(ScoredPosition sp : parser.getPositions()){
				if(sp.getOverHang() !=2) continue;
			//	if(sp.getGenomicRegions3p() == null || !sp.getGenomicRegions3p().get(0).endsWith("Intron")) continue;
			//	if(sp.getGenomicRegions5p() == null || !sp.getGenomicRegions5p().get(0).endsWith("Intron")) continue;
				
				out.print(sp);
				String ksuf = "";
				if(sp.getClassification().equals("M") && sp.hasMatchingMiRNA()){
					 
					ksuf += "MiRNA";
				}else{
					ksuf += sp.getClassification() + (sp.isPaired() ? "T" : "F");
				}
				String k3 = sp.getContig() + "\t" + sp.getThreePposition() + "\t" + sp.isPlusStrand();
				String k5 = sp.getContig() + "\t" + sp.getFivePposition() + "\t" + sp.isPlusStrand();
				Integer regionID5 = region5pMap.get(k3);
				Integer regionID3 = region3pMap.get(k5);
				
				if(regionID5 == null) regionID5 = 0;
				if(regionID3 == null) regionID3 = 0;
				int regionID = Math.max(regionID5, regionID3);
				
				for(int i=0; i<cKeys.length;i++){
					
					String cKey = ksuf + "_" +  cKeys[i] + "_" + regionID;				
					
					Double c5 = creadCnt5pMaps.get(i).get(k3); 
					Double c3 = creadCnt3pMaps.get(i).get(k5);
									
					if(c5 == null) c5 = .0;
					if(c3 == null) c3 = .0;
					
					out.print("\t" + (c5 + c3));
					
					if(!mStringMap.containsKey(cKey)) mStringMap.put(cKey, new StringBuilder());
					mStringMap.get(cKey).append((c5 + c3) + " ");
					
					for(int j=0; i + j*cKeys.length<tKeys.length; j++){
						String tKey = ksuf + "_" + tKeys[i + j*cKeys.length] + "_" + regionID;						
						Double t5 = treadCnt5pMaps.get(i + j*cKeys.length).get(k3); 
						Double t3 = treadCnt3pMaps.get(i + j*cKeys.length).get(k5);
						double median = medians.get(cKeys[i]+"_" + regionID+";"+tKeys[i + j*cKeys.length]+"_" + regionID);		
						if(t5 == null) t5 = .0;
						if(t3 == null) t3 = .0;
						
						out.print("\t" + (t5 + t3) +  "\t" + ((t5+t3 + c5+c3 > fcThreshold)?  Math.log((t5+t3)/(c5+c3))/Math.log(2) - median : " ") );
						
						if(!mStringMap.containsKey(tKey)) mStringMap.put(tKey, new StringBuilder());
						mStringMap.get(tKey).append((t5 + t3) + " ");
					}
				}
				out.println();				
			}
		
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("];");
				//outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=20 & " + mkey + "(:,2)>= 20),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			outm.close();
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	public static void generateForTrans(String transOut, String mOut, String transCsv, 
			String[] cBeds5p, String[] cBeds3p, String[] tBeds5p, String[] tBeds3p, 
			String[] cKeys, String[] tKeys,  int maxSpan, boolean filtered, int fcThreshold, Hashtable<String, Double> medians){
		ScoredPairOutputParser parser = new ScoredPairOutputParser(transCsv);
		ArrayList<Hashtable<String, Double>> creadCnt3pMaps = new ArrayList<Hashtable<String, Double>>();
		ArrayList<Hashtable<String, Double>> creadCnt5pMaps = new ArrayList<Hashtable<String, Double>>();
		ArrayList<Hashtable<String, Double>> treadCnt3pMaps = new ArrayList<Hashtable<String, Double>>();
		ArrayList<Hashtable<String, Double>> treadCnt5pMaps = new ArrayList<Hashtable<String, Double>>();
		
		Hashtable<String, Integer> region3pMap = new Hashtable<String, Integer>();
		Hashtable<String, Integer> region5pMap = new Hashtable<String, Integer>();
			
		try {	
			
			PrintStream out = new PrintStream(transOut);
			out.print(parser.getHeader());
			//out.print(parser.getHeader());
			for(int i=0; i<cKeys.length; i++){
				String ckey = cKeys[i];
				out.print("\t" + ckey);
				for(int j=0; i + j*cKeys.length<tKeys.length; j++){
					String tkey = tKeys[i + j*cKeys.length];
					out.print("\t" + tkey + "\t" + tkey + "_FC");
				}
			}
			out.println();
		
			MapRunner runner1 = new MapRunner(creadCnt3pMaps, region3pMap, cBeds3p, maxSpan);
			MapRunner runner2 = new MapRunner(creadCnt5pMaps, region5pMap, cBeds5p, maxSpan);
			MapRunner runner3 = new MapRunner(treadCnt3pMaps, region3pMap, tBeds3p, maxSpan);
			MapRunner runner4 = new MapRunner(treadCnt5pMaps, region5pMap, tBeds5p, maxSpan);
			
			Thread thread1 = new Thread(runner1, "" + 1);       
			Thread thread2 = new Thread(runner2, "" + 2);
			Thread thread3 = new Thread(runner3, "" + 3);
			Thread thread4 = new Thread(runner4, "" + 4);
			
			thread1.start();
			thread2.start();
			thread3.start();
			thread4.start();
			
			try {
				thread1.join();
				thread2.join();
				thread3.join();
				thread4.join();
			} catch (InterruptedException e) {						
				e.printStackTrace();
			}
			
			PrintStream outm = new PrintStream(mOut);
			Hashtable<String, StringBuilder> mStringMap = new Hashtable<String, StringBuilder>();
			
			for(ScoredPair pair : parser.getPairs()){
				if(pair.getOverHang() !=2) continue;
			//	if(pair.getGenomicRegions3p() == null || !pair.getGenomicRegions3p().get(0).endsWith("Intron")) continue;
			//	if(pair.getGenomicRegions5p() == null || !pair.getGenomicRegions5p().get(0).endsWith("Intron")) continue;
				
				out.print(pair);
				String ksuf = pair.getClassification();
				String contig1 = pair.getThreePContig();
				String contig2 = pair.getFivePContig();
				boolean isPlusStrand1 = pair.getThreePStrand();
				boolean isPlusStrand2 = pair.getFivePStrand();
				
				String k3 = contig1  + "\t" + pair.getThreePPosition() + "\t" + isPlusStrand1; 
				String k5 = contig2 + "\t" + pair.getFivePPosition() + "\t" + isPlusStrand2;
				Integer regionID5 = region5pMap.get(k3);
				Integer regionID3 = region3pMap.get(k5);
				
				if(regionID5 == null) regionID5 = 0;
				if(regionID3 == null) regionID3 = 0;
				int regionID = Math.max(regionID5, regionID3);
				
				for(int i=0; i<cKeys.length;i++){
					
					String cKey = (filtered? "TransFiltered" : "Trans") + ksuf + "_" + cKeys[i] + "_" + regionID;				
					
					Double c5 = creadCnt5pMaps.get(i).get(k3); 
					Double c3 = creadCnt3pMaps.get(i).get(k5);
									
					if(c5 == null) c5 = .0;
					if(c3 == null) c3 = .0;
					
					out.print("\t" + (c5 + c3));
					
					if(!mStringMap.containsKey(cKey)) mStringMap.put(cKey, new StringBuilder());
					mStringMap.get(cKey).append((c5 + c3) + " ");
					
					for(int j=0; i + j*cKeys.length<tKeys.length; j++){
						String tKey = (filtered? "TransFiltered" : "Trans") + ksuf + "_" + tKeys[i + j*cKeys.length]  + "_" + regionID;						
						Double t5 = treadCnt5pMaps.get(i + j*cKeys.length).get(k3); 
						Double t3 = treadCnt3pMaps.get(i + j*cKeys.length).get(k5);
						double median = medians.get(cKeys[i]+"_" + regionID+";"+tKeys[i + j*cKeys.length]+"_" + regionID);		
									
						if(t5 == null) t5 = .0;
						if(t3 == null) t3 = .0;
						
						out.print("\t" + (t5 + t3) +  "\t" + ((t5+t3 + c5+c3 > fcThreshold)?  Math.log((t5+t3)/(c5+c3))/Math.log(2) - median: " ") );
						
						if(!mStringMap.containsKey(tKey)) mStringMap.put(tKey, new StringBuilder());
						mStringMap.get(tKey).append((t5 + t3) + " ");
					}
				}
				
				out.println();				
			}
		
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("];"); 
				//outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=20 & " + mkey + "(:,2)>= 20),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			outm.close();
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}	
	}	
	//bedtools intersect -f 0.6 -c -S -split -a h19x2.sorted.out.5p.bed -b /media/kyowon/Data1/fCLIP/RNAseq/siControl_R1_Aligned_Sorted.bed -sorted > h19x2.sorted.out.Control.5p.csv
	
	public static void generate(String csv, String pairCsv, String bed5p, String bed3p, String siKDBed, String siControlBed) {		
		int si = siKDBed.lastIndexOf(System.getProperty("file.separator"));
		String key = siKDBed.substring(si + 1, si + 6);
		
		String threeBed = csv + "." + key + ".3p.bed";
		String fiveBed = csv + "." + key + ".5p.bed";
		
		String threeCBed = csv + ".siControl.3p.bed";
		String fiveCBed = csv + ".siControl.5p.bed";
		
		if(!new File(threeCBed).exists()){
			ShellLauncher.run("bedtools intersect -f 0.6 -c -S -split -a " + bed5p+" -b " + siControlBed + " -sorted > " + fiveCBed);
			ShellLauncher.run("bedtools intersect -f 0.6 -c -S -split -a " + bed3p+" -b " + siControlBed + " -sorted > " + threeCBed);
		}
		if(!new File(threeBed).exists()){
			ShellLauncher.run("bedtools intersect -f 0.6 -c -S -split -a " + bed5p+" -b " + siKDBed + " -sorted > " + fiveBed);
			ShellLauncher.run("bedtools intersect -f 0.6 -c -S -split -a " + bed3p+" -b " + siKDBed + " -sorted > " + threeBed);
		}
		
		String outCsv = csv + "." + key + ".csv";
		String outPairCsv = pairCsv + "." + key + ".csv";
		
		String outMfile = new File(outPairCsv).getParent() + System.getProperty("file.separator") + new File(outPairCsv).getName().replace('.', '_') + ".m";  
	//	System.out.println(outMfile);
		
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		Hashtable<String, Double> readCCnt3pMap = new Hashtable<String, Double>();
		Hashtable<String, Double> readTCnt3pMap = new Hashtable<String, Double>();
		
		Hashtable<String, Double> readCCnt5pMap = new Hashtable<String, Double>();
		Hashtable<String, Double> readTCnt5pMap = new Hashtable<String, Double>();
		
		try {				
			PrintStream outm = new PrintStream(outMfile);
			Hashtable<String, StringBuilder> mStringMap = new Hashtable<String, StringBuilder>();
			String mKeyPrefix = key;
			String s;
			BufferedLineReader tin = new BufferedLineReader(threeBed);
			while((s=tin.readLine())!=null){
				String[] token = s.split("\t");
				String k = token[0] + "\t" + token[3] + "\t" + (token[5].equals("+"));
				readTCnt3pMap.put(k, Double.parseDouble(token[6]));
			}
			tin.close();
			
			BufferedLineReader fin = new BufferedLineReader(fiveBed);
			while((s=fin.readLine())!=null){
				String[] token = s.split("\t");
				String k = token[0] + "\t" + token[3] + "\t" + (token[5].equals("+"));
				readTCnt5pMap.put(k, Double.parseDouble(token[6]));
				
			}
			fin.close();
			
			
			
			BufferedLineReader tcin = new BufferedLineReader(threeCBed);
			while((s=tcin.readLine())!=null){
				String[] token = s.split("\t");
				String k = token[0] + "\t" + token[3] + "\t" + (token[5].equals("+"));
				readCCnt3pMap.put(k,  Double.parseDouble(token[6]));				
			}
			tcin.close();
			
			BufferedLineReader fcin = new BufferedLineReader(fiveCBed);
			while((s=fcin.readLine())!=null){
				String[] token = s.split("\t");
				String k = token[0] + "\t" + token[3] + "\t" + (token[5].equals("+"));
				readCCnt5pMap.put(k,  Double.parseDouble(token[6]));
			
			}
			fcin.close();			
			
			PrintStream out = new PrintStream(outCsv);
			out.println(ScoredPosition.getHeader() + "\tControl5p\tControl3p\tTarget5p\tTarget3p");
			for(ScoredPosition sp : parser.getPositions()){
				String mKey = mKeyPrefix;
				
				if(sp.hasMatchingMiRNA()){
					mKey += "MiRNA";
				}else{
					mKey += sp.getClassification() + (sp.isPaired() ? "T" : "F");
				}
				out.print(sp);
				String k3 = sp.getContig() + "\t" + sp.getThreePposition() + "\t" + sp.isPlusStrand();
				String k5 = sp.getContig() + "\t" + sp.getFivePposition() + "\t" + sp.isPlusStrand();
				
				Double c5 = readCCnt5pMap.get(k3); 
				Double c3 = readCCnt3pMap.get(k5);
				
				Double t5 = readTCnt5pMap.get(k3); 
				Double t3 = readTCnt3pMap.get(k5);
				
				
				
				if(c5 == null) c5 = .0;
				if(c3 == null) c3 = .0;
				if(t5 == null) t5 = .0;
				if(t3 == null) t3 = .0;
				
				out.println("\t" + c5 + "\t" + c3 + "\t" + t5 + "\t" + t3);
				
				if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
				mStringMap.get(mKey).append((c5 + c3) + " " + (t5 + t3) + ";");
			}
			out.close();
			
			ScoredPairOutputParser pairParser = new ScoredPairOutputParser(pairCsv);
			PrintStream pairOut = new PrintStream(outPairCsv);
			HashSet<String> k3s = new HashSet<String>();
			HashSet<String> k5s = new HashSet<String>();
			
			pairOut.println(ScoredPair.getHeader() + "\tControl5p\tControl3p\tTarget5p\tTarget3p");
			
			for(ScoredPair pair : pairParser.getPairs()){
				String mKey = mKeyPrefix + "Pair" + pair.getClassification();
				
				String contig1 = pair.getThreePContig();
				String contig2 = pair.getFivePContig();
				boolean isPlusStrand1 = pair.getThreePStrand();
				boolean isPlusStrand2 = pair.getFivePStrand();
				
				String k3 = contig1  + "\t" + pair.getThreePPosition() + "\t" + isPlusStrand1; 
				String k5 = contig2 + "\t" + pair.getFivePPosition() + "\t" + isPlusStrand2;
			//	if(k3s.contains(k3)) continue;
			//	if(k5s.contains(k5)) continue;
				
				k3s.add(k3);
				k5s.add(k5);
				pairOut.print(pair);
				
				Double c5 = readCCnt5pMap.get(k3); 
				Double c3 = readCCnt3pMap.get(k5);
				
				Double t5 = readTCnt5pMap.get(k3); 
				Double t3 = readTCnt3pMap.get(k5);
				
				if(c5 == null) c5 = .0;
				if(c3 == null) c3 = .0;
				if(t5 == null) t5 = .0;
				if(t3 == null) t3 = .0;
				//c5 = t5 = .0;
				pairOut.print("\t" + c5 + "\t" + c3	+ "\t" + t5 + "\t" + t3);
				
				if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
				mStringMap.get(mKey).append((c5 + c3) + " " + (t5 + t3) + ";");
				pairOut.println();
			}
			
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("];");
				outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=20 & " + mkey + "(:,2)>= 20),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			outm.close();
			
			pairOut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	private static Hashtable<String, Double> getMedians(String file){
		Hashtable<String, Double> medians = new Hashtable<String, Double>();
		try {
			BufferedLineReader in = new BufferedLineReader(file);
			String s;
			
			while((s=in.readLine())!=null){
				if(!s.startsWith("%")) continue;
				String[] token = s.split(";");
				medians.put(token[0].substring(1)+";"+token[1], Double.parseDouble(token[2]));
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return medians;
	}
	
	public static void main(String[] args) {
		if(args[args.length - 1].equals("CIS")){ // cis
			String cisOut = args[0];
			//int regionID = Integer.parseInt(args[11]);
			
			String cisMOut = args[1] + new File(cisOut).getName().replace('.', '_');
			cisMOut = cisMOut.substring(0, cisMOut.length() - 4) + ".m";			
			String cisCsv = args[2];			
			String[] cBeds5p = args[3].split(",");
			String[] cBeds3p = args[4].split(",");
			String[] tBeds5p = args[5].split(",");
			String[] tBeds3p = args[6].split(",");
			String[] cKeys = args[7].split(",");
			String[] tKeys = args[8].split(",");
			int maxSpan = Integer.parseInt(args[9]);	
			int fcThreshold = Integer.parseInt(args[10]);
			generateForCis(cisOut, cisMOut, cisCsv, cBeds5p, cBeds3p, tBeds5p, tBeds3p, cKeys, tKeys, maxSpan, fcThreshold, getMedians(args[11]));
		}
		if(args[args.length - 1].equals("TRANS")){
			String transOut = args[0];
			String transMOut =args[1] + new File(transOut).getName().replace('.', '_');
		//	int regionID = Integer.parseInt(args[12]);
			
			transMOut = transMOut.substring(0, transMOut.length() - 4)  + ".m";
			String transCsv = args[2];
			String[] cBeds5p = args[3].split(",");
			String[] cBeds3p = args[4].split(",");
			String[] tBeds5p = args[5].split(",");
			String[] tBeds3p = args[6].split(",");
			String[] cKeys = args[7].split(",");
			String[] tKeys = args[8].split(",");
			int maxSpan = Integer.parseInt(args[9]);		
			boolean filtered = args[10].equals("True");
			int fcThreshold = Integer.parseInt(args[11]);
			generateForTrans(transOut, transMOut, transCsv, cBeds5p, cBeds3p, tBeds5p, tBeds3p, cKeys, tKeys, maxSpan, filtered, fcThreshold, getMedians(args[12]));
		}
		
		if(args[args.length - 1].equals("BG")){
			String mOut = args[0];
			String[] intersectedBeds = args[1].split(",");
			String[] keys = args[2].split(",");	
			int maxSpan = Integer.parseInt(args[3]);		
			int cutThreshold = Integer.parseInt(args[4]);
			int numThreads = Integer.parseInt(args[5]);		
		//	int regionID = Integer.parseInt(args[6]);
			generateForBackground(mOut, intersectedBeds, keys, cutThreshold, maxSpan, numThreads);
		}
		
	}
}
