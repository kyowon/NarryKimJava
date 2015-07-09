package fCLIP.analysis;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Random;
import java.util.Vector;

import launcher.RNAfoldLauncher;
import parser.AnnotationFileParser;
import parser.BufferedLineReader;
import parser.ZeroBasedFastaParser;
import parser.AnnotationFileParser.AnnotatedGene;
import fCLIP.FCLIP_Scorer;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;


public class NewCheckCoverage {
	
	public static class RunnerForCis implements Runnable{
		private String outfile;
		private String[] cKeys, tKeys;
		private SamReader[] creaders1, treaders1;
		private Hashtable<String, StringBuffer> mStringMap;
		private Hashtable<String, Double> medians;
		private double fcThreshold;
		private Vector<ScoredPosition> positions;
		private int maxSpan;
		
		public RunnerForCis(String outfile, ArrayList<ScoredPosition> sps, String[] cKeys, String[] tKeys, 
				SamReader[] creaders1, SamReader[] treaders1, Hashtable<String, StringBuffer> mStringMap, Hashtable<String, Double> medians,
				double fcThreshold, int maxSpan, int nt, int tn)
		{
			this.outfile = outfile;
			this.positions = new Vector<ScoredPositionOutputParser.ScoredPosition>();
			for(int i=0;i<sps.size();i++){
				if(i%nt == tn) this.positions.add(sps.get(i));
			}
			
			this.cKeys = cKeys;
			this.tKeys = tKeys;
			this.creaders1 = creaders1;
			this.treaders1 = treaders1;
			this.mStringMap = mStringMap;
			this.medians = medians;
			this.fcThreshold = fcThreshold;
			this.maxSpan = maxSpan;
		}
		
		public void run(){
			PrintStream out;
			try {
				out = new PrintStream(outfile);		 
				int total = positions.size();
				int cntr = 0;
				int prev = 1;
				
				for(ScoredPosition sp : positions){
					//if(sp.getOverHang() !=2) continue;
					out.print(sp);
					String ksuf = "";
					if(sp.getClassification().equals("M") && sp.hasMatchingMiRNA()){						 
						ksuf += "MiRNA";
					}else{
						ksuf += sp.getClassification() + (sp.isPaired() ? "T" : "F");
					}
								
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
					
					for(int i=0; i<cKeys.length;i++){
						String cKey = ksuf + "_" +  cKeys[i];				
						
						int[] c = getPassingThroughReadCount(sp, creaders1[i]);
						int c5 = c[0];
						int c3 = c[1];
						
						out.print("\t" + c5 + "\t" + c3);
						
						if(!mStringMap.containsKey(cKey)) mStringMap.put(cKey, new StringBuffer());
						//mStringMap.get(cKey).append((c5 + c3) + ";");
						
						StringBuffer sbc = new StringBuffer();
						sbc.append((c5));sbc.append(",");
						sbc.append((c3));sbc.append(",");
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
						mStringMap.get(cKey).append(sbc);
						
						for(int j=0; i + j*cKeys.length<tKeys.length; j++){
							String tKey = ksuf + "_" + tKeys[i + j*cKeys.length];						
							
							int[] t =  getPassingThroughReadCount(sp, treaders1[i + j*cKeys.length]);
							int t5 = t[0];
							int t3 = t[1];
							
							double median = medians.get(cKeys[i]+";"+tKeys[i + j*cKeys.length]);		
							
							out.print("\t" + t5 +  "\t" + t3 +  "\t" + ((t5+t3 + c5+c3 > fcThreshold)?  Math.log(((double)t5+t3)/((double)c5+c3))/Math.log(2) - median : " ") );
							
							if(!mStringMap.containsKey(tKey)) mStringMap.put(tKey, new StringBuffer());
							StringBuffer sbt =  new StringBuffer();
							sbt.append((t5));sbt.append(",");
							sbt.append((t3));sbt.append(",");
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
							mStringMap.get(tKey).append(sbt);
						}
					}
					if(100 * cntr++ / total > prev){
						System.out.println("Thread " +Thread.currentThread().getName() + ": CIS " + (prev) + " % done..");
						prev += 1;
					}	
					out.println();	
				}
					
			} catch (IOException e) {
			e.printStackTrace();
			}
		}
	}
	
	
	public static class RunnerForBackground implements Runnable{
		private Vector<String> positions = null;
		private String[] keys;
		private SamReader[] readers1;
		private Hashtable<String, Vector<Double>> ctm;	 
		private ZeroBasedFastaParser fastaParser;
	    private Hashtable<String, StringBuffer> mStringMap;
	  //  private int maxSpan;
	    
		public RunnerForBackground(ArrayList<String> positions,  String[] keys, SamReader[] readers1, 
				Hashtable<String, StringBuffer> mStringMap, 
				Hashtable<String, Vector<Double>> ctm,
				int maxSpan, int nt, int tn , ZeroBasedFastaParser fastaParser){
			this.positions = new Vector<String>();
			for(int i=0;i<positions.size();i++){
				if(i%nt == tn) this.positions.add(positions.get(i));
			}
			this.keys = keys;
			this.readers1 = readers1;
			this.ctm = ctm;
			this.mStringMap = mStringMap;
			//this.maxSpan = maxSpan;	
			this.fastaParser = fastaParser;
		}
		
		public void run() {						
			int total = positions.size();
			//mStringMap = new Hashtable<String, StringBuffer>();
			
			int cntr = 0;
			int prev = 1;
			//System.out.println("Thread " +Thread.currentThread().getName() + ": BG start");
			for(String position : positions){
				String[] token = position.split("\t");
				int p1 = Integer.parseInt(token[1]);
				int p2 = Integer.parseInt(token[2]); 
				int fn = FCLIP_Scorer.getFlankingNTNumber();
				
				String seq = fastaParser.getSequence(token[0], p1 - fn, p2 + fn + 1);
				RNAfoldLauncher fold = new RNAfoldLauncher(seq, fn);			
				
				for(int i=0; i<keys.length;i++){
					int[] c = getPassingThroughReadCount(token[0], p1, p2, readers1[i]);					
					String mKey = "BG_" + new String(keys[i]);
					if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuffer());
					if(!ctm.containsKey(keys[i])) ctm.put(keys[i], new Vector<Double>());					
					StringBuffer sb = new StringBuffer();
					sb.append(c[0]);sb.append(",");
					sb.append(c[1]);sb.append(",");
					sb.append(fold.getOverHang());sb.append(",");
					sb.append(fold.getEnergyPerNT());sb.append(",");
					sb.append(token[3]);sb.append(",");
					sb.append(token[4]);sb.append(",");					
					sb.append(" %");
					sb.append(token[0]);
					sb.append(";");
					sb.append(p1);
					sb.append(";");
					sb.append(p2);					
					sb.append("\n");
					mStringMap.get(mKey).append(sb);
					ctm.get(keys[i]).add((double)c[0] + c[1]); // TODO change later..
				//	System.out.print(Thread.currentThread().getName() + " " + sb);
				}
				//System.out.println("Thread " +Thread.currentThread().getName() + ": " + cntr + " out of " + positions.size());
				if(100 * cntr++ / total > prev){
					System.out.println("Thread " +Thread.currentThread().getName() + ": BG " + (prev) + " % done..");
					prev += 1;
				}
			}				
			
		
		}
	}
	
	//p1 < p2 and 0 based
	public static boolean isPassingThrough(SAMRecord rec, int p1){
		int p2 = p1 + 1;
		if(p1 <= 0 || p2 <= 0) return false;
		for(AlignmentBlock b : rec.getAlignmentBlocks()){
			int s = b.getReferenceStart() - 1; // zero based
			int e = s + b.getLength() - 1;
			if(e < p1) continue;
			if(s > p2) return false;
			if(s <= p1 && e >= p2) return true;
			else return false;			
		}
		return false;
	}
	
	// rec should be the read 1 and paired
	/*private static boolean isPairPassingThrough(SAMRecord mrec, SAMRecord rec, int p1, int maxSpan){ // p1 to p1+1
		if(p1 <= 0) return false;
		boolean r1passing = isPassingThrough(rec, p1);
		if(r1passing) return true;
		
		//SAMRecord mrec = reader.queryMate(rec);
		if(mrec == null) return r1passing;
		
		int le, rs;
		if(rec.getAlignmentStart() < mrec.getAlignmentStart()){
			le = rec.getAlignmentEnd();
			rs = mrec.getAlignmentStart();
		}else{
			le = mrec.getAlignmentEnd();
			rs = rec.getAlignmentStart();
		}
		le--; rs--; // zero based
		if(rs>=le & rs-le<=maxSpan & le<= p1 & p1+1<=rs) return true; 
		
		return isPassingThrough(mrec, p1);		
	}
	
  
	public static int[] getPassingThroughReadPairCount(String contig, int p1, int p2, SamReader reader1, int maxSpan){
		SAMRecordIterator iterator = reader1.query(contig, p1+1-maxSpan, p2+1+maxSpan, false);	
		int[] c = new int[2];
		ArrayList<SAMRecord> recs1 = new ArrayList<SAMRecord>();
		Hashtable<String, SAMRecord> recs2 = new Hashtable<String, SAMRecord>();
		
		int cntr = 0;
		
		while(true){
			boolean hasNext = iterator.hasNext();
			if(hasNext){
				SAMRecord rec = iterator.next();
		        //if(!rec.getReadPairedFlag()) continue;
				if(rec.getFirstOfPairFlag()) 			
					recs1.add(rec);
				else recs2.put(rec.getReadName(), rec);
				
				cntr++;
			} 
			if(!hasNext || cntr > 1000000){
				cntr = 0;
				//System.out.print(contig + " " + p1 + " " + p2 + " : read number : " + recs1.size());
				for(SAMRecord tr : recs1){
					SAMRecord mtr = recs2.get(tr.getReadName());
					if(isPairPassingThrough(mtr, tr, p1-1, maxSpan)) c[0]++;
					if(isPairPassingThrough(mtr, tr, p2, maxSpan)) c[1]++;
					recs2.remove(tr.getReadName());
				}
				for(SAMRecord tr : recs2.values()){
					if(isPairPassingThrough(null, tr, p1-1, maxSpan)) c[0]++;
					if(isPairPassingThrough(null, tr, p2, maxSpan)) c[1]++;
				}				
				recs1.clear();	
				recs2.clear();
				//System.out.println(" - counting done");
			}
			if(!hasNext) break;
		}						
		iterator.close();
		return c;
	}*/
	
	// is3p = true : left side of cleavage if strand == '+'
	private static int[] getPassingThroughReadCount(String contig, int p1, int p2, SamReader reader){
		
		SAMRecordIterator iterator = reader.query(contig, p1+1, p2+1, false);
		try{
			int[] c = new int[2];
			while(iterator.hasNext()){
				SAMRecord rec = iterator.next(); 
				//if(!rec.getReadPairedFlag()) continue;
				//if(countRead1Only && !rec.getFirstOfPairFlag()) continue;				
				if(isPassingThrough(rec, p1-1)) c[0]++;
				if(isPassingThrough(rec, p2)) c[1]++;
			}
			return c;
		}finally{
			iterator.close();
		}
	}
	
	/*private static int[] getPassingThroughReadPairCount(ScoredPosition position, SamReader reader, int maxSpan){
		boolean strand = position.isPlusStrand();
		int p1 = strand? position.getThreePPosition() : position.getFivePPosition();
		int p2 = strand? position.getFivePPosition() : position.getThreePPosition();	
		return getPassingThroughReadPairCount(position.getContig(),  (p1<0?p2 : p1), (p2<0?p1 : p2), reader, maxSpan);
		
	}*/	
		
	// output[0] = 3p output[1] = 5p
	private static int[] getPassingThroughReadCount(ScoredPosition position, SamReader reader){
		boolean strand = position.isPlusStrand();
		int p1 = strand? position.getThreePPosition() : position.getFivePPosition();
		int p2 = strand? position.getFivePPosition() : position.getThreePPosition();	
	    		
		SAMRecordIterator iterator = reader.query(position.getContig(), (p1<0?p2 : p1)+1, (p2<0?p1 : p2)+1, false);
	
		try{
			int[] c = new int[2];
			while(iterator.hasNext()){
				SAMRecord rec = iterator.next(); 
				//if(!rec.getReadPairedFlag()) continue;
				//if(countRead1Only && !rec.getFirstOfPairFlag()) continue;				
				
				if(isPassingThrough(rec, p1-1)) c[strand? 0 : 1]++;
				if(isPassingThrough(rec, p2)) c[strand? 1 : 0]++;
			}
		//	System.out.println(position + " " + (p1<0?p2 : p1)+1 + " " + (p2<0?p1 : p2)+1 + "\t" + c[0] + " " + c[1]);
			return c;
		}finally{
			iterator.close(); //chr21 82170491 82170991
		}
	}	
	
	public static void generateForCis(String cisOut, String mOut, String cisCsv, 
			String[] cBams, String[] tBams,	String[] cKeys, String[] tKeys, int maxSpan, int fcThreshold, Hashtable<String, Double> medians, int nt){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(cisCsv);
			
		try {		
			ArrayList<Hashtable<String, StringBuffer>> mStringMaps = new ArrayList<Hashtable<String, StringBuffer>>();
			for(int n=0;n<nt;n++){
				mStringMaps.add(new Hashtable<String, StringBuffer>());
			}
			ArrayList<Thread> threads = new ArrayList<Thread>();
			ArrayList<ScoredPosition> sps = parser.getPositions();
			for(int n=0;n<nt;n++){
				SamReader[] creaders1 = new SamReader[cBams.length];
				SamReader[] treaders1 = new SamReader[tBams.length];
				
				for(int i=0;i<creaders1.length;i++){
					creaders1[i] = SamReaderFactory.makeDefault().open(new File(cBams[i]));
				}	
				for(int i=0;i<treaders1.length;i++){
					treaders1[i] = SamReaderFactory.makeDefault().open(new File(tBams[i]));
				}	
				RunnerForCis runner = new RunnerForCis(cisOut+"."+n+".tmp", sps, cKeys, tKeys, creaders1, 
						treaders1, mStringMaps.get(n), new Hashtable<String, Double>(medians), fcThreshold, maxSpan, nt, n);
				Thread thread = new Thread(runner, "CISRunner " + n); //Thread created       
				thread.start();
				threads.add(thread);
			}
			
			for(Thread thread : threads){
				try {
					thread.join();
				} catch (InterruptedException e) {						
					e.printStackTrace();
				}
			}
						
			PrintStream out = new PrintStream(cisOut);
			
			out.print(parser.getHeader());
			for(int i=0; i<cKeys.length; i++){
				String ckey = cKeys[i];
				out.print("\t" + ckey + "_5p\t" + ckey + "_3p");
				for(int j=0; i + j*cKeys.length<tKeys.length; j++){
					String tkey = tKeys[i + j*cKeys.length];
					out.print("\t" + tkey + "_5p\t" + tkey + "_3p\t" + tkey + "_FC");
				}
			}
			out.println();
			for(int n=0;n<nt;n++){
				BufferedLineReader in = new BufferedLineReader(cisOut+"."+n+".tmp");
				String s;
				while((s=in.readLine())!=null){
					out.println(s);
				}
				in.close();
				new File(cisOut+"."+n+".tmp").delete();
			}
			
			Hashtable<String, StringBuffer> mStringMap = new Hashtable<String, StringBuffer>();
			
			for(Hashtable<String, StringBuffer> ms : mStringMaps){
				for(String k : ms.keySet()){
					if(!mStringMap.containsKey(k)) mStringMap.put(k, new StringBuffer());
					mStringMap.get(k).append(ms.get(k));
				}
			}
						
			PrintStream outm = new PrintStream(mOut);
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("]';");
				//outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=20 & " + mkey + "(:,2)>= 20),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			outm.close();
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	
	public static void generateForTrans(String transOut, String mOut, String transCsv, 
		  String cisCsv, String[] cKeys, String[] tKeys,  int maxSpan, boolean filtered, boolean control, int fcThreshold, Hashtable<String, Double> medians){
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
			String key5 = sp.getContig() + ";" + sp.getThreePPosition() + ";" + sp.isPlusStrand();
			String key3 = sp.getContig() + ";" + sp.getFivePPosition() + ";" + sp.isPlusStrand();
			
			cntrTable.put(key3, new Hashtable<String, Integer>());
			cntrTable.put(key5, new Hashtable<String, Integer>());
			
		//	System.out.println(sp.getMiscInfo());
			
			String[] misc = sp.getMiscInfo().split("\t");
			
			for(String ck : cKeys){
				int i3 = headerIndex.get(ck+"_3p");
				int r3 = Integer.parseInt(misc[i3]);
				int i5 = headerIndex.get(ck+"_5p");
				int r5 = Integer.parseInt(misc[i5]);
				cntrTable.get(key3).put(ck, r3);		
				cntrTable.get(key5).put(ck, r5);
			}
			
			for(String tk : tKeys){
				int i3 = headerIndex.get(tk+"_3p");
				int r3 = Integer.parseInt(misc[i3]);
				int i5 = headerIndex.get(tk+"_5p");
				int r5 = Integer.parseInt(misc[i5]);
				cntrTable.get(key3).put(tk, r3);		
				cntrTable.get(key5).put(tk, r5);
			}
		}
		try {	
			PrintStream out = new PrintStream(transOut);
			out.print(transParser.getHeader());
			for(int i=0; i<cKeys.length; i++){
				String ckey = cKeys[i];
				out.print("\t" + ckey + "_5p\t" + ckey + "_3p");
				for(int j=0; i + j*cKeys.length<tKeys.length; j++){
					String tkey = tKeys[i + j*cKeys.length];
					out.print("\t" + tkey + "_5p\t" + tkey + "_3p\t" + tkey + "_FC");
				}
			}
			out.println();		
			
			PrintStream outm = new PrintStream(mOut);
			Hashtable<String, StringBuffer> mStringMap = new Hashtable<String, StringBuffer>();
			//HashSet<ScoredPosition> positions = new HashSet<ScoredPositionOutputParser.ScoredPosition>();


				for(ScoredPair pair : transParser.getPairs()){
			//	if(pair.getOverHang() !=2) continue;
				
				out.print(pair);
				String ksuf = (control? "Control" : "") + pair.getClassification();
				
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
				
				
				for(int i=0; i<cKeys.length;i++){
					
					String cKey = (filtered? "TransFiltered" : "Trans") + ksuf + "_" + cKeys[i];				
					
					String key5 = pair.getThreePContig() + ";" + pair.getThreePPosition() + ";" + pair.getThreePStrand();
					String key3 = pair.getFivePContig() + ";" + pair.getFivePPosition() + ";" + pair.getFivePStrand();
						
					int c5 = cntrTable.get(key5).get(cKeys[i]);
					int c3 = cntrTable.get(key3).get(cKeys[i]);					
					
					out.print("\t" + c5 +  "\t" + c3);
					
					if(!mStringMap.containsKey(cKey)) mStringMap.put(cKey, new StringBuffer());
					StringBuffer sbc = mStringMap.get(cKey);
					sbc.append((c5));sbc.append(",");
					sbc.append((c3));sbc.append(",");
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
				
					for(int j=0; i + j*cKeys.length<tKeys.length; j++){
						String tKey = (filtered? "TransFiltered" : "Trans") + ksuf + "_" + tKeys[i + j*cKeys.length];						
						double median = medians.get(cKeys[i]+";"+tKeys[i + j*cKeys.length]);		
							
						int t5 = cntrTable.get(key5).get(tKeys[i + j*cKeys.length]);
						int t3 = cntrTable.get(key3).get(tKeys[i + j*cKeys.length]);					
						
						out.print("\t" + t5 + "\t" + t3 +  "\t" + ((t5+t3 + c5+c3 > fcThreshold)?  Math.log(((double)t5+t3)/((double)c5+c3))/Math.log(2) - median: " ") );
						if(!mStringMap.containsKey(tKey)) mStringMap.put(tKey, new StringBuffer());
						//mStringMap.get(tKey).append((t5 + t3) + ";");
						StringBuffer sbt = mStringMap.get(tKey);
						sbt.append((t5));sbt.append(",");
						sbt.append((t3));sbt.append(",");
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
				
				out.println();				
			}
		
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("]';"); 
			}
			outm.close();
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}	
	}	
	
	
	static private ArrayList<String> _generateRandomizedPositionsFromCsv(String csv, AnnotationFileParser annotationParser, ZeroBasedFastaParser fastaParser){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		ArrayList<String> ret = new ArrayList<String>();
		HashSet<String> accessions = new HashSet<String>();
		for(ScoredPosition sp : parser.getPositions()){
			if(sp.getContainingGeneAccessions() != null)
				accessions.addAll(sp.getContainingGeneAccessions());
		}			
		
		
		for(ScoredPosition sp : parser.getPositions()){
			for(int n=0;n<50;n++){
				int rd1 = 0;
				int rd2 = 0;
				int r1 = 0;
				int r2 = 0;
				while(true){
					ArrayList<AnnotatedGene> genes = annotationParser.getAnnotatedGenes(sp.getContig());
					int rdi = new Random().nextInt(genes.size()); // select random gene
					AnnotatedGene gene = genes.get(rdi);
					if(accessions.contains(gene.getAccession())) continue;
					
					rd1 = gene.getTxStart() + new Random().nextInt(gene.getTxEnd() - gene.getTxStart());	
					rd2 = rd1 + FCLIP_Scorer.getMinReadDiff() + new Random().nextInt(FCLIP_Scorer.getMaxReadDiff() - FCLIP_Scorer.getMinReadDiff());
									
					String region1 = annotationParser.getGenomicRegionNameAndFrameShift(sp.getContig(), sp.isPlusStrand(), rd1, gene, true).get(0);
					r1 = region1 == null? 0 : (region1.endsWith("Intron")? 1 : 2);
					if(r1 == 0) continue;
					
					String region2 = annotationParser.getGenomicRegionNameAndFrameShift(sp.getContig(), sp.isPlusStrand(), rd2, gene, true).get(0);
					r2 = region2 == null? 0 : (region2.endsWith("Intron")? 1 : 2);
					if(r2 == 0) continue;
					
					break;
						//if(rdi%2==0){
						//	if(r == r3) break;
						//}else{
						//	if(r == r5) break;
						//}
				//	}
				}
				
				String fstring = sp.getContig() + "\t" + rd1 + "\t" + rd2 + "\t" + r1 + "\t" + r2;
				ret.add(fstring);				
			}
		}
		
		return ret;		
	}
	
	static private ArrayList<String> generateRandomizedPositionsFromCsv(String csv, AnnotationFileParser annotationParser, ZeroBasedFastaParser fastaParser){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		ArrayList<String> ret = new ArrayList<String>();
		
		for(ScoredPosition sp : parser.getPositions()){
			for(int n=0;n<10;n++){
				int rd1 = 0;
				int rd2 = 0;
				int r1 = 0;
				int r2 = 0;
				while(true){				
					boolean sign = new Random().nextBoolean();
					rd1 = (sp.getFivePPosition() < 0? sp.getThreePPosition() : sp.getFivePPosition()) + (sign? 1 : -1) * (50000 + new Random().nextInt(20000)); 
					rd2 = rd1 + FCLIP_Scorer.getMinReadDiff() + new Random().nextInt(FCLIP_Scorer.getMaxReadDiff() - FCLIP_Scorer.getMinReadDiff());
							
					for(int i=0;i<2;i++){
						ArrayList<AnnotatedGene> genes = annotationParser.getContainingGenes(sp.getContig(), i==0, rd1);
						
						if(genes != null && !genes.isEmpty()){
							AnnotatedGene gene = genes.get(0);
							ArrayList<String> region1 = annotationParser.getGenomicRegionNameAndFrameShift(sp.getContig(), i==0, rd1, gene, true);							
							int tr1 = 0;
							if(region1 != null && !region1.isEmpty()){
								tr1 = 1;
								for(String r : region1){
									if(r.endsWith("ORF") || r.endsWith("UTR")){
										tr1 = 2;
										break;
									}
								}
							}
							r1 = Math.max(r1, tr1);							
							
							ArrayList<String> region2 = annotationParser.getGenomicRegionNameAndFrameShift(sp.getContig(), i==0, rd2, gene, true);
							int tr2 = 0;
							if(region2 != null && !region2.isEmpty()){
								tr2 = 1;
								for(String r : region2){
									if(r.endsWith("ORF") || r.endsWith("UTR")){
										tr2 = 2;
										break;
									}
								}
							}
							r2 = Math.max(r2, tr2);
						}						
					}
					break;
				}
				
				String fstring = sp.getContig() + "\t" + rd1 + "\t" + rd2 + "\t" + r1 + "\t" + r2;
				ret.add(fstring);				
			}
		}
		
		return ret;		
	}
	
	
	public static double getMedianFC(Vector<Double> c, Vector<Double> t, int th){
		Vector<Double> fcs = new Vector<Double>();
		for(int i=0;i<c.size();i++){
			double n = c.get(i);
			double m = t.get(i);
			if(n + m <= th) continue;
			
			fcs.add(Math.log((m+1)/(n+1))/Math.log(2));
		}
		Collections.sort(fcs);
		return fcs.get(fcs.size()/2);
	}
	
	
	
	
	public static void generateForBackground(String csv,  AnnotationFileParser annotationParser, ZeroBasedFastaParser fastaParser, String mOut, String[] bams, String[] keys, int cutThreshold, int maxSpan, int nt){
		try {	
			ArrayList<String> positions = generateRandomizedPositionsFromCsv(csv, annotationParser, fastaParser);
		
			ArrayList<Hashtable<String, StringBuffer>> mStringMaps = new ArrayList<Hashtable<String, StringBuffer>>();
			ArrayList<Hashtable<String, Vector<Double>>> ctms = new ArrayList<Hashtable<String, Vector<Double>>>();
					
			for(int n=0;n<nt;n++){
				mStringMaps.add(new Hashtable<String, StringBuffer>());
				ctms.add(new Hashtable<String, Vector<Double>>());
			}
			
			ArrayList<Thread> threads = new ArrayList<Thread>();
			for(int n=0;n<nt;n++){
				
				SamReader[] readers1 = new SamReader[bams.length];
				for(int i=0;i<readers1.length;i++){
					readers1[i] = SamReaderFactory.makeDefault().open(new File(bams[i]));
				}	
				RunnerForBackground runner = new RunnerForBackground(positions, keys, readers1, mStringMaps.get(n), ctms.get(n), maxSpan, nt, n, fastaParser);
				Thread thread = new Thread(runner, "BGRunner " + n); //Thread created       
				thread.start();
				threads.add(thread);
			}
			
			for(Thread thread : threads){
				try {
					thread.join();
				} catch (InterruptedException e) {						
					e.printStackTrace();
				}
			}			
			
			Hashtable<String, StringBuffer> mStringMap = new Hashtable<String, StringBuffer>();
			Hashtable<String, Vector<Double>> ctm = new Hashtable<String, Vector<Double>>();
			
			for(Hashtable<String, StringBuffer> ms : mStringMaps){
				for(String k : ms.keySet()){
					if(!mStringMap.containsKey(k)) mStringMap.put(k, new StringBuffer());
					mStringMap.get(k).append(ms.get(k));
				}
			}
			
			for(Hashtable<String, Vector<Double>> ct : ctms){
				for(String k : ct.keySet()){
					if(!ctm.containsKey(k)) ctm.put(k, new Vector<Double>());
					ctm.get(k).addAll(ct.get(k));
				}
			}
			PrintStream outm = new PrintStream(mOut);
			
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey.replace('-', '_') + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("]';");
				//outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=20 & " + mkey + "(:,2)>= 20),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			
			for(String key1 : ctm.keySet()){
				for(String key2 : ctm.keySet()){
					if(key1.equals(key2)) continue;					
					outm.print("%"+key1 + ";"+key2+";");
					outm.println(getMedianFC(ctm.get(key1), ctm.get(key2), cutThreshold));
				}
			}
			
			outm.close();
			
			
		}catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	static Hashtable<String, Double> getMedians(String file){
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
		//for(String k: medians.keySet()){
		//	System.out.println(k + " " + medians.get(k));
		//}
		
		return medians;
	}
	
	public static void main(String[] args) {
		if(args[args.length - 1].equals("CIS")){ // cis
			String cisOut = args[0];
			String cisMOut = args[1] + new File(cisOut).getName().replace('.', '_');
			cisMOut = cisMOut.substring(0, cisMOut.length() - 4) + ".m";			
			String cisCsv = args[2];			
			String[] cBams = args[3].split(",");
			String[] tBams = args[4].split(",");
			String[] cKeys = args[5].split(",");
			String[] tKeys = args[6].split(",");
			int maxSpan = Integer.parseInt(args[7]);	
			int fcThreshold = Integer.parseInt(args[8]);
			String bc = args[9];
		    int numThreads = Integer.parseInt(args[10]);
			generateForCis(cisOut, cisMOut, cisCsv, cBams, tBams, cKeys, tKeys, maxSpan, fcThreshold, getMedians(bc), numThreads);
		}
		if(args[args.length - 1].equals("TRANS")){
			String transOut = args[0];
			String transMOut =args[1] + new File(transOut).getName().replace('.', '_');
			transMOut = transMOut.substring(0, transMOut.length() - 4)  + ".m";
			String transCsv = args[2];
			String cisCsv = args[3];
			String[] cKeys = args[4].split(",");
			String[] tKeys = args[5].split(",");
			int maxSpan = Integer.parseInt(args[6]);		
			boolean filtered = args[7].equals("True");
			boolean control = args[8].equals("True");
			int fcThreshold = Integer.parseInt(args[9]);
			String bc = args[10];
			generateForTrans(transOut, transMOut, transCsv, cisCsv, cKeys, tKeys, maxSpan, filtered, control, fcThreshold, getMedians(bc));
		}
		
		if(args[args.length - 1].equals("BG")){
			String cisCsv = args[0];
			AnnotationFileParser annotationParser = new AnnotationFileParser(args[1]);
			ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser(args[2]);
			String mOut = args[3];
			String[] bams = args[4].split(",");
			String[] keys = args[5].split(",");	
			int maxSpan = Integer.parseInt(args[6]);		
			int cutThreshold = Integer.parseInt(args[7]);
		    int numThreads = Integer.parseInt(args[8]);
			generateForBackground(cisCsv, annotationParser, fastaParser, mOut, bams, keys, cutThreshold, maxSpan, numThreads);
		}		
	}

}























