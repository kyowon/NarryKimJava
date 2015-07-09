package fCLIP.analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;
import parser.AnnotationFileParser;
import parser.BufferedLineReader;
import parser.ZeroBasedFastaParser;

public class NewCheckGeneFoldChange {
	// trans
	
	public static void generateForBackground(String mOut, String[] beds, String[] keys){
		
		try {	
			PrintStream out = new PrintStream(mOut);		
			Hashtable<String, Vector<Double>> ctm = new Hashtable<String, Vector<Double>>(); // key , coverages
			ArrayList<String> accessions = new ArrayList<String>(getGeneCoverage(beds[0]).keySet());			
			
			for(int i=0;i<keys.length;i++){
				Vector<Double> covs = new Vector<Double>();
				HashMap<String, Integer> cov = getGeneCoverage(beds[i]);
				out.println("GeneBG_" + keys[i].replace('-', '_') + "=[");
				for(String accession : accessions){
					covs.add((double)cov.get(accession));
					out.println(cov.get(accession));
				}
				ctm.put(keys[i], covs);
				out.println("]';");
			}
					
			for(String key1 : ctm.keySet()){
				for(String key2 : ctm.keySet()){
					if(key1.equals(key2)) continue;					
					out.print("%"+key1 + ";"+key2+";");
					out.println(NewCheckCoverage.getMedianFC(ctm.get(key1), ctm.get(key2), 1));
				}
			}
			
			
			out.close();
		}catch(IOException e){
			
		}
	}
	
	public static void generateForTrans(String transOut, String mOut, String transCsv,//
			String cisCsv, String[] cKeys, String[] tKeys, boolean filtered, boolean control, Hashtable<String, Double> medians){
		ScoredPairOutputParser transParser = new ScoredPairOutputParser(transCsv);
		ScoredPositionOutputParser cisParser = new ScoredPositionOutputParser(cisCsv);
		Hashtable<String, Integer> headerIndex = new Hashtable<String, Integer>();
		Hashtable<String, Hashtable<String, Integer>> cntrTable = new Hashtable<String, Hashtable<String,Integer>>();
		
		String[] cisHeader = cisParser.getHeader().split("\t");
		int oi = ScoredPosition.getHeader().split("\t").length;
		for(int i=oi ;i<cisHeader.length;i++){
			headerIndex.put(cisHeader[i], i-oi);
		}
		
		for(ScoredPosition sp : cisParser.getPositions()){
			if(sp.getContainingGeneAccessions() == null || sp.getContainingGeneAccessions().isEmpty()) continue;
			String[] misc = sp.getMiscInfo().split("\t");
			
			for(int n=0;n<sp.getContainingGeneAccessions().size();n++){
				String accession = sp.getContainingGeneAccessions().get(n);
				if(!cntrTable.containsKey(accession)) cntrTable.put(accession, new Hashtable<String, Integer>());
				for(String ck : cKeys){
					int i = headerIndex.get("Gene_"+ck);
					int r = Integer.parseInt(misc[i].split(";")[n]);				
					cntrTable.get(accession).put(ck, r);
				}				
				for(String tk : tKeys){
					int i = headerIndex.get("Gene_"+tk);
					int r = Integer.parseInt(misc[i].split(";")[n]);
					cntrTable.get(accession).put(tk, r);
				}
			}				
		}
		try {	
			PrintStream out = new PrintStream(transOut);
			out.print(transParser.getHeader());
			for(int i=0; i<cKeys.length; i++){
				String ckey = "Gene_" + cKeys[i];
				out.print("\t" + ckey + "_5p\t" + ckey + "_3p");
				for(int j=0; i + j*cKeys.length<tKeys.length; j++){
					String tkey = "Gene_" + tKeys[i + j*cKeys.length];
					out.print("\t" + tkey + "_5p\t" + tkey + "_3p\t" + tkey + "_FC");
				}
			}
			out.println();
			
			PrintStream outm = new PrintStream(mOut);
			Hashtable<String, StringBuffer> mStringMap = new Hashtable<String, StringBuffer>();
			
			for(ScoredPair pair : transParser.getPairs()){
				int region1 = 0;
				if(pair.getGenomicRegions3p() != null && !pair.getGenomicRegions3p().isEmpty()){
					region1 = 1;
					for(String r : pair.getGenomicRegions3p()){
						if(r.endsWith("ORF") || r.endsWith("UTR")){
							region1 = 2;
							break;
						}
					}
				}
					
				int region2 = 0;
				if(pair.getGenomicRegions5p() != null && !pair.getGenomicRegions5p().isEmpty()){
					region2 = 1;
					for(String r : pair.getGenomicRegions5p()){
						if(r.endsWith("ORF") || r.endsWith("UTR")){
							region2 = 2;
							break;
						}
					}
				}
			
				
				out.print(pair);
				String ksuf = (control? "Control" : "") + pair.getClassification();
				
				for(int i=0; i<cKeys.length;i++){
					
					String cKey = (filtered? "TransFiltered" : "Trans") + ksuf + "_" + cKeys[i];				
					if(!mStringMap.containsKey(cKey)) mStringMap.put(cKey, new StringBuffer());
					StringBuffer sbc = mStringMap.get(cKey);
					
					String c5str = "";
					String c3str = "";
					for(String accession5 : pair.getPairedScoredPositions()[0].getContainingGeneAccessions()){
						for(String accession3 : pair.getPairedScoredPositions()[1].getContainingGeneAccessions()){
							int c5 = cntrTable.get(accession5).get(cKeys[i]);
							int c3 = cntrTable.get(accession3).get(cKeys[i]);			
							c5str += c5 + ";";
							c3str += c3 +";" ;
							
							sbc.append(c5+c3);sbc.append(",");
							//sbc.append(c3);sbc.append(",");
							sbc.append(pair.getOverHang());sbc.append(",");
							sbc.append(pair.getEnergy());sbc.append(",");
							sbc.append(region1);sbc.append(",");
							sbc.append(region2);sbc.append(",");
							sbc.append(pair.isRepeat5p()? 1 :0);sbc.append(",");
							sbc.append(pair.isRepeat3p()? 1 :0);sbc.append(",");						
							sbc.append(0);sbc.append(",");
							sbc.append(pair.getPredictionScore());sbc.append(",");
									
							sbc.append("%");
							sbc.append(pair.getThreePContig());sbc.append(" ");
							sbc.append(pair.getThreePPosition());sbc.append(" ");
							sbc.append(pair.getFivePContig());sbc.append(" ");
							sbc.append(pair.getFivePPosition());sbc.append(" ");
							sbc.append("\n");
						}
					}
					out.print("\t" + (c5str.isEmpty()? "_" : c5str) +  "\t" + (c3str.isEmpty()? "_" : c3str)); 
					
					for(int j=0; i + j*cKeys.length<tKeys.length; j++){
						String tKey = (filtered? "TransFiltered" : "Trans") + ksuf + "_" + tKeys[i + j*cKeys.length];	
						if(!mStringMap.containsKey(tKey)) mStringMap.put(tKey, new StringBuffer());
						StringBuffer sbt = mStringMap.get(tKey);
						double median = medians.get(cKeys[i]+";"+tKeys[i + j*cKeys.length]);								
						
						String t5str = "";
						String t3str = "";						
						String fcstr = "";						
						
						for(String accession5 : pair.getPairedScoredPositions()[0].getContainingGeneAccessions()){
							for(String accession3 : pair.getPairedScoredPositions()[1].getContainingGeneAccessions()){
								int c5 = cntrTable.get(accession5).get(cKeys[i]);
								int c3 = cntrTable.get(accession3).get(cKeys[i]);					
								
								int t5 = cntrTable.get(accession5).get(tKeys[i]);
								int t3 = cntrTable.get(accession3).get(tKeys[i]);		
								t5str += t5 + ";";
								t3str += t3 + ";";
								fcstr += (Math.log(((double)t5+t3)/((double)c5+c3))/Math.log(2) - median) + ";";
								sbt.append(t5+t3);sbt.append(",");
								//sbt.append(t3);sbt.append(",");
								sbt.append(pair.getOverHang());sbt.append(",");
								sbt.append(pair.getEnergy());sbt.append(",");
								sbt.append(region1);sbt.append(",");
								sbt.append(region2);sbt.append(",");
								sbt.append(pair.isRepeat5p()? 1 :0);sbt.append(",");
								sbt.append(pair.isRepeat3p()? 1 :0);sbt.append(",");						
								sbt.append(0);sbt.append(",");
								sbt.append(pair.getPredictionScore());sbt.append(",");
										
								sbt.append("%");
								sbt.append(pair.getThreePContig());sbt.append(" ");
								sbt.append(pair.getThreePPosition());sbt.append(" ");
								sbt.append(pair.getFivePContig());sbt.append(" ");
								sbt.append(pair.getFivePPosition());sbt.append(" ");
								sbt.append("\n");
							}
						}	
						out.print("\t" + (t5str.isEmpty()? "_" : t5str) + "\t" + (t3str.isEmpty()? "_" : t3str) +  "\t" + (fcstr.isEmpty()? "_" : fcstr) );						
					}
				}				
				out.println();				
			}
		
			for(String mkey : mStringMap.keySet()){
				outm.println("Gene" + mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("]';"); 
			}
			outm.close();
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	
	public static void generateForCis(String cisOut, String mOut, String cisCsv, 
			String[] cBeds, String[] tBeds,	String[] cKeys, String[] tKeys, Hashtable<String, Double> medians){ 
		
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(cisCsv);		
		try {	
			PrintStream out = new PrintStream(cisOut);
			
			out.print(parser.getHeader());
			for(int i=0; i<cKeys.length; i++){
				String ckey = "Gene_" + cKeys[i];
				out.print("\t" + ckey);
				for(int j=0; i + j*cKeys.length<tKeys.length; j++){
					String tkey = "Gene_" +  tKeys[i + j*cKeys.length];
					out.print("\t" + tkey + "\t" + tkey + "_FC");
				}
			}
			out.println(); // print header
			
			HashMap<String, HashMap<String, Integer>> covMap = new HashMap<String, HashMap<String,Integer>>();
			HashMap<String, StringBuffer> mStringMTMap = new HashMap<String, StringBuffer>();
			HashMap<String, StringBuffer> mStringMiRNAMap = new HashMap<String, StringBuffer>();
			HashMap<String, StringBuffer> mStringMap = null;
			
			for(int i=0;i<cKeys.length;i++){
				covMap.put(cKeys[i], getGeneCoverage(cBeds[i]));
			}
			
			for(int i=0;i<tKeys.length;i++){
				covMap.put(tKeys[i], getGeneCoverage(tBeds[i]));
			}
			
			for(ScoredPosition sp : parser.getPositions()){
				
				if(sp.hasMatchingMiRNA()) mStringMap = mStringMiRNAMap;
				else mStringMap = mStringMTMap;
				
				int region1 = 0;
				if(sp.getGenomicRegions3p() != null && !sp.getGenomicRegions3p().isEmpty()){
					region1 = 1;
					for(String r : sp.getGenomicRegions3p()){
						if(r.endsWith("ORF") || r.endsWith("UTR")){
							region1 = 2;
							break;
						}
					}
				}
					
				int region2 = 0;
				if(sp.getGenomicRegions5p() != null && !sp.getGenomicRegions5p().isEmpty()){
					region2 = 1;
					for(String r : sp.getGenomicRegions5p()){
						if(r.endsWith("ORF") || r.endsWith("UTR")){
							region2 = 2;
							break;
						}
					}
				}
				
				out.print(sp); out.print('\t');
				int[] ccovs = null;
				int[] tcovs = null;
				for(int i=0; i<cKeys.length; i++){
					String ckey = cKeys[i];		
					if(!mStringMap.containsKey(ckey)) mStringMap.put(ckey, new StringBuffer());
					if(sp.getContainingGeneAccessions() == null || sp.getContainingGeneAccessions().isEmpty()){
						out.print("_\t");
					}else{
						
						StringBuffer sbc = mStringMap.get(ckey);
						HashMap<String, Integer> cmap = covMap.get(ckey);
						ccovs = new int[sp.getContainingGeneAccessions().size()];
						for(int k=0; k<sp.getContainingGeneAccessions().size(); k++){
							String accession = sp.getContainingGeneAccessions().get(k);
							ccovs[k] = cmap.get(accession);
							out.print(ccovs[k]+";");
							if(sp.getClassification().equals("M")){	
								sbc.append(ccovs[k]);sbc.append(",");
								sbc.append(sp.getOverHang());sbc.append(",");
								sbc.append(sp.getEnergy());sbc.append(",");
								sbc.append(region1);sbc.append(",");
								sbc.append(region2);sbc.append(",");
								sbc.append(sp.isRepeat5p()? 1 :0);sbc.append(",");
								sbc.append(sp.isRepeat3p()? 1 :0);sbc.append(",");	
								sbc.append(sp.getPreLength());sbc.append(",");
								sbc.append(sp.getPredictionScore());sbc.append(",");
							
								sbc.append("%");
								sbc.append(sp.getContig());sbc.append(" ");
								sbc.append(sp.getThreePPosition());sbc.append(" ");
								sbc.append(sp.getFivePPosition());sbc.append(" ");
								sbc.append("\n");
							}
						}
						out.print("\t");
					}
					for(int j=0; i + j*cKeys.length<tKeys.length; j++){
						String tkey = tKeys[i + j*cKeys.length];
						if(!mStringMap.containsKey(tkey)) mStringMap.put(tkey, new StringBuffer());
						
						double median = medians.get(ckey+";"+tkey);								
						if(sp.getContainingGeneAccessions() == null || sp.getContainingGeneAccessions().isEmpty()){
							out.print("_\t_\t");
						}else{
							StringBuffer sbt = mStringMap.get(tkey);							
							HashMap<String, Integer> tmap = covMap.get(tkey);
							tcovs = new int[sp.getContainingGeneAccessions().size()];
							for(int k=0; k<sp.getContainingGeneAccessions().size(); k++){
								String accession = sp.getContainingGeneAccessions().get(k);
								tcovs[k] = tmap.get(accession);
								out.print(tcovs[k]+";");
								
								if(sp.getClassification().equals("M")){
									sbt.append(tcovs[k]);sbt.append(",");
									sbt.append(sp.getOverHang());sbt.append(",");
									sbt.append(sp.getEnergy());sbt.append(",");
									sbt.append(region1);sbt.append(",");
									sbt.append(region2);sbt.append(",");
									sbt.append(sp.isRepeat5p()? 1 :0);sbt.append(",");
									sbt.append(sp.isRepeat3p()? 1 :0);sbt.append(",");	
									sbt.append(sp.getPreLength());sbt.append(",");
									sbt.append(sp.getPredictionScore());sbt.append(",");
																	
									sbt.append("%");
									sbt.append(sp.getContig());sbt.append(" ");
									sbt.append(sp.getThreePPosition());sbt.append(" ");
									sbt.append(sp.getFivePPosition());sbt.append(" ");
									sbt.append("\n");
								}
							}
							out.print("\t");
							for(int k=0;k<tcovs.length;k++){
								double c = ccovs[k];
								double t = tcovs[k];
								double fc = Math.log(t/c)/Math.log(2) - median;
								out.print(fc+";");
							}
							out.print("\t");
						}						
					}
				}
				out.println();
			}			
			out.close();
			PrintStream outm = new PrintStream(mOut);
			for(String mkey : mStringMTMap.keySet()){
				outm.println("GeneMT_" + mkey.replace('-', '_') + "= [");
				outm.println(mStringMTMap.get(mkey));
				outm.println("]';");
			}
			
			for(String mkey : mStringMiRNAMap.keySet()){
				outm.println("GeneMiRNA_" + mkey.replace('-', '_') + "= [");
				outm.println(mStringMiRNAMap.get(mkey));
				outm.println("]';");
			}
			outm.close();			
		
		}catch (IOException e) {
			
		}
	}
	
	
	private static HashMap<String, Integer> getGeneCoverage(String covBed){ // key : accession value : coverage
		BufferedLineReader in;	
		HashMap<String, Integer> ret = new HashMap<String, Integer>();
		try {
			in = new BufferedLineReader(covBed);
			String s;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				//3 vs 12
				ret.put(token[3], Integer.parseInt(token[12]));
			}
			in.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
		return ret;
	}
	
	public static void main(String[] args){		
		if(args[args.length - 1].equals("CIS")){ // cis
			String cisOut = args[0];
			String cisMOut = args[1] + new File(cisOut).getName().replace('.', '_');
			cisMOut = cisMOut.substring(0, cisMOut.length() - 4) + "_gene.m";			
			String cisCsv = args[2];			
			String[] cBeds = args[3].split(",");
			String[] tBeds = args[4].split(",");
			String[] cKeys = args[5].split(",");
			String[] tKeys = args[6].split(",");
			String bc = args[7];		    
		    generateForCis(cisOut, cisMOut, cisCsv, cBeds, tBeds, cKeys, tKeys, NewCheckCoverage.getMedians(bc));	
		}
		if(args[args.length - 1].equals("TRANS")){
			String transOut = args[0];
			String transMOut =args[1] + new File(transOut).getName().replace('.', '_');
			transMOut = transMOut.substring(0, transMOut.length() - 4)  + "_gene.m";
			String transCsv = args[2];
			String cisCsv = args[3];
			String[] cKeys = args[4].split(",");
			String[] tKeys = args[5].split(",");
		//	int maxSpan = Integer.parseInt(args[6]);		
			boolean filtered = args[6].equals("True");
			boolean control = args[7].equals("True");
			String bc = args[8];
			generateForTrans(transOut, transMOut, transCsv, cisCsv, cKeys, tKeys, filtered, control, NewCheckCoverage.getMedians(bc));
		}
		
		if(args[args.length - 1].equals("BG")){
			String mOut = args[0];
			String[] beds = args[1].split(",");
			String[] keys = args[2].split(",");	
			generateForBackground(mOut, beds, keys);
		}
	}
}
