package fCLIP.paper;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import fCLIP.FCLIP_Scorer;
import launcher.RNAcofoldLauncher;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.ZeroBasedFastaParser;

public class Figure1 {

	
	private static double[] normalize(double[] v, double threshold){
		double[] ret = new double[v.length];
		double sum = getSum(v);
		
		double ratioThreshold = 1;
		if(sum >= 5) ratioThreshold = 4.0/sum;
			
		if(sum > threshold){
			for(int i=0;i<v.length;i++){				
				ret[i] = v[i] / sum;	
				if(ret[i] < ratioThreshold) ret[i] = 0;
			}
		}
		return ret;
	}
	
	private static double getSum(double[] v){
		double sum = 0;
		for(double t : v){
			sum += t;
		}
		return sum;
	}
	
	private static double[] findNmaxValues(double[] v, int n){
		double[] ret = new double[n];
		double[] nv = Arrays.copyOf(v, v.length);
		Arrays.sort(nv);
		for(int i=0;i<n;i++){
			ret[i] = nv[nv.length-1-i];
		}
		return ret;
	}
	
	
	private static ArrayList<Integer> getIndicesExceedingThreshold(double[] v, int readSumThreshold, double threshold){ // v is not normalized
		ArrayList<Integer> indices = new ArrayList<Integer>();
			//if(sum < readSumThreshold) return indices;
		double[] nv = normalize(v, readSumThreshold);
		for(int i=0;i<nv.length;i++){
			if(nv[i] > threshold) indices.add(i);
		}
		return indices;
	}
	
	
	
	private static int getMaxIndex(double[] v){
		int index = 0;
		double max = 0;
		for(int i=0;i<v.length;i++){
			if(max < v[i]){
				max = v[i];
				index = i;
			}
		}
		return index;
	}
	
	private static ArrayList<Integer> getOffsets(ArrayList<Integer> p1s, ArrayList<Integer> p2s, int n){
		HashSet<Integer> offsetSet = new HashSet<Integer>();
		
		for(int p : p1s)
			offsetSet.add(p - (n - 1));
		
		for(int p : p2s)
			offsetSet.add(p - n);
			
		ArrayList<Integer> offsets = new ArrayList<Integer>(offsetSet);
		Collections.sort(offsets);
		return offsets;
	}
	
	private static ArrayList<ArrayList<Integer>> getGroupNumbers(int[] groupNumbers, double[][] v, int readSumThreshold, double readRatioThreshold, int flankingLength, int preLength){
		ArrayList<ArrayList<Integer>> offsets5p3p = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> indices = new ArrayList<ArrayList<Integer>>();
		for(int i=0;i<4;i++){
			indices.add(getIndicesExceedingThreshold(v[i], readSumThreshold, readRatioThreshold));
		}
		ArrayList<Integer> fp1s = indices.get(0);
		ArrayList<Integer> fp2s = indices.get(1);
		
		if(fp1s.isEmpty() && fp2s.isEmpty()){
			groupNumbers[0] = 0;
			offsets5p3p.add(new ArrayList<Integer>());
		}
		else{
			ArrayList<Integer> offsets = getOffsets(fp1s, fp2s, flankingLength);
			if(offsets.size() > 1) groupNumbers[0] = 10;
			else{
				if(!fp1s.isEmpty() && fp2s.isEmpty()){
					groupNumbers[0] = 1;
				}else if(!fp1s.isEmpty() && !fp2s.isEmpty()){
					groupNumbers[0] = 2;
				}else{
					groupNumbers[0] = 3;
				}
				if(offsets.get(0) >= 1) groupNumbers[0] += 6;
				if(offsets.get(0) <= -1) groupNumbers[0] += 3;
			}
			offsets5p3p.add(offsets);
		}
		
		ArrayList<Integer> tp2s = indices.get(2);
		ArrayList<Integer> tp1s = indices.get(3);
		
		if(tp1s.isEmpty() && tp2s.isEmpty()){
			groupNumbers[1] = 0;
			offsets5p3p.add(new ArrayList<Integer>());
		}
		else{
			ArrayList<Integer> offsets = getOffsets(tp2s, tp1s, preLength);
			if(offsets.size() > 1) groupNumbers[1] = 10;
			else{
				if(!tp1s.isEmpty() && tp2s.isEmpty()){
					groupNumbers[1] = 1;
				}else if(!tp1s.isEmpty() && !tp2s.isEmpty()){
					groupNumbers[1] = 2;
				}else{
					groupNumbers[1] = 3;
				}
				if(offsets.get(0) >= 1) groupNumbers[1] += 3;
				if(offsets.get(0) <= -1) groupNumbers[1] += 6;
			}
			offsets5p3p.add(offsets);
		}
		
		return offsets5p3p;
	}
	
	
	public static void main(String[] args) {
		CountDroshaCleavagesAroundmiRNAs.main(null);
		
		String input = "/media/kyowon/Data1/fCLIP/miRNA_ReadEnd_unpairedIncluded.txt";
		String output= "/media/kyowon/Data1/fCLIP/miRNA_ReadEnd_unpairedIncluded.m";
		String outcsv= "/media/kyowon/Data1/fCLIP/miRNA_ReadEnd_unpairedIncluded.csv"; 
		
		
		String parameterFileName = "/media/kyowon/Data1/fCLIP/samples/sample4/results/training/merged.param";
				
		int preLength = 10;
		int flankingLength = 10;
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/fCLIP/genomes/hg38.fa");
		
		int flankingNTNumber = 13;
		
		int readSumThreshold = 1;
		int totalReadSumThreshold = 10;
		double readRatioThreshold = .333333;
		
		
		String prevContig = null;
		Bed12Parser bedParser = null;
		
		try {
			BufferedLineReader in = new BufferedLineReader(input);
			PrintStream out = new PrintStream(output);
			PrintStream outc = new PrintStream(outcsv);
			outc.println("miRNA\tGroup5p\tGroup3p\tOverhang\tContig\tStrand\tmirBasePosition5p\tfCLIPPosition5p\tfCLIP5pScore\tDiff5p\tmirBasePosition3p\tfCLIPPosition3p\tfCLIP3pScore\tDiff3p\tTotalReadSum");
			
			StringBuffer f1 = new StringBuffer();  
			StringBuffer f2 = new StringBuffer();
			StringBuffer t1 = new StringBuffer();
			StringBuffer t2 = new StringBuffer();
			StringBuffer gp = new StringBuffer();
			StringBuffer oh = new StringBuffer();
			StringBuffer ll = new StringBuffer();
		//	StringBuffer mirP = new StringBuffer();
		//	StringBuffer fCLIPP = new StringBuffer();
			
			StringBuffer paires = new StringBuffer();
			f1.append("five1=[");
			f2.append("five2=[");
			t1.append("three1=[");
			t2.append("three2=[");
			gp.append("group=[");
			oh.append("overhang=[");
			ll.append("loopLength=[");
			paires.append("paired=[");
			//mirP.append("mirP=[");
			//fCLIPP.append("fCLIPP=[");
			String s;
			out.println("miRNAs={");
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				double[][] reads = new double[4][];
				double sum = 0;
				for(int i=6; i<token.length;i++){
					String[] stoken = token[i].split(",");
					reads[i-6] = new double[stoken.length];
					for(int j=0;j<stoken.length;j++){
						reads[i-6][j] = Double.parseDouble(stoken[j]);
						sum += reads[i-6][j];
					}					
				}
				if(sum < totalReadSumThreshold) continue;
				
				int[] groupNumbers = new int[2]; 
				ArrayList<ArrayList<Integer>> offsets = getGroupNumbers(groupNumbers, reads, readSumThreshold, readRatioThreshold, flankingLength , preLength);
				if(groupNumbers[0] == 0 && groupNumbers[1] == 0) continue;
				out.println("'" + token[0] + "'");
				
				Integer overhang = -1000;
				int loopLength = -1;
				
				ArrayList<Integer> fp5 = new ArrayList<Integer>();
				ArrayList<Integer> fp3 = new ArrayList<Integer>();
								
				int fu = 0;
				int fd = 0;
				int tu = 0;
				int td = 0;
				String chr = token[1];
				boolean isPlusStrand = token[5].equals("+");
				Integer p5 = null;
				Integer p3 = null;
				
				p5 = Integer.parseInt(token[3]);	
				p3 = Integer.parseInt(token[4]); // inclusive
				
				if(groupNumbers[0] >0){					
					ArrayList<Integer> offsets5p = offsets.get(0);
					for(int o : offsets5p){
						fp5.add(p5 + (isPlusStrand? o : -o));
					}					
				}
				
				
				if(groupNumbers[1] >0){
					ArrayList<Integer> offsets3p = offsets.get(1);
					for(int o : offsets3p){
						fp3.add(p3 + (isPlusStrand? o : -o));
					}				
				}
				
					
				if(fp5.size() == 1 && fp3.size() == 1){
					String seq = isPlusStrand? fastaParser.getSequence(chr, fp5.get(0) - flankingNTNumber, fp3.get(0) + flankingNTNumber + 1)
						    : ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(chr, 
						    		fp3.get(0)-flankingNTNumber, fp5.get(0) + flankingNTNumber + 1), true);
					RNAcofoldLauncher test = new RNAcofoldLauncher(seq , flankingNTNumber);
					
					overhang = test.getOverHang();	
					System.out.println(token[0] + " " + fp5.get(0) + " " + fp3.get(0) + " " + overhang + " "+ seq);
		//		groupMap.put(token[0], groupNumbers);
				}				
				
				for(double v : normalize(reads[0], readSumThreshold)){
					f1.append(v);
					f1.append(',');
				}
				f1.append('\n');
				
				for(double v : normalize(reads[1], readSumThreshold)){
					f2.append(v);
					f2.append(',');
				}
				f2.append('\n');
				
				for(double v : normalize(reads[2], readSumThreshold)){
					t2.append(v);
					t2.append(',');
				}
				t2.append('\n');
				
				for(double v : normalize(reads[3], readSumThreshold)){
					t1.append(v);
					t1.append(',');
				}
				t1.append('\n');				
				
				if(token[2].startsWith("false")) groupNumbers[0] = -groupNumbers[0];
				if(token[2].endsWith("false")) groupNumbers[1] = -groupNumbers[1];
					
				gp.append(groupNumbers[0]);gp.append(' ');gp.append(groupNumbers[1]); gp.append('\n');
				oh.append(overhang); oh.append('\n');
				ll.append(loopLength); ll.append('\n');
				paires.append(fd + "," + fu + "," + td + "," + tu + "\n");
				//mirP.append(p5 + "," + p3);
				//fCLIPP.append(fp5 + "," + fp3);
				if(prevContig == null || !prevContig.equals(chr)){
					bedParser = new Bed12Parser("/media/kyowon/Data1/fCLIP/samples/sample4/bed/merged.se.bed", 
						new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt"),	chr, true);
					prevContig = chr;
				}
				FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, null, null, parameterFileName);
				
			
				String fstr = "";
				
				for(int i=0 ;i<fp5.size();i++){
					int p = fp5.get(i);
					fstr += (p+1) + (i<fp5.size()-1 ? ";" : "");
				}
				fstr += "\t";
				
				for(int i=0 ;i<fp5.size();i++){
					int p = fp5.get(i);
					fstr += scorer.getScore(p, isPlusStrand, false) + (i<fp5.size()-1 ? ";" : "");
				}
				
				fstr += "\t";
				
				for(int i=0 ;i<fp5.size();i++){
					int p = fp5.get(i);
					fstr += (p5-p) + (i<fp5.size()-1 ? ";" : "");
				}
				
				
				String tstr = "";
				for(int i=0 ;i<fp3.size();i++){
					int p = fp3.get(i);
					tstr += (p+1) + (i<fp3.size()-1 ? ";" : "");
				}
				tstr += "\t";
				
				for(int i=0 ;i<fp3.size();i++){
					int p = fp3.get(i);
					tstr += scorer.getScore(p, isPlusStrand, true) + (i<fp3.size()-1 ? ";" : "");
				}
				
				tstr += "\t";
				
				for(int i=0 ;i<fp3.size();i++){
					int p = fp3.get(i);
					tstr += (p3-p) + (i<fp3.size()-1 ? ";" : "");
				}
				
					
				outc.println(token[0] + "\t" + groupNumbers[0] + "\t" + groupNumbers[1] + "\t" + overhang + "\t" + chr + "\t" + (isPlusStrand? '+' : '-') + "\t" +  (p5 != null? (p5+1) : "N/A") + "\t" + fstr + "\t" 
				        + (p3 != null? (p3+1) : "N/A") + "\t" + tstr + "\t" + sum);
				
			}	
			out.println("};");
			
			f1.append("];");
			f2.append("];");
			t1.append("];");
			t2.append("];");
			gp.append("];");
			oh.append("];");
			ll.append("];");
			paires.append("];");
		//	mirP.append("];");
		//	fCLIPP.append("];");
			out.println(f1);
			out.println(f2);
			out.println(t1);
			out.println(t2);
			out.println(gp);
			out.println(oh);
			out.println(ll);
			out.println(paires);
		//	out.println(mirP);
		//	out.println(fCLIPP);
			in.close();
			out.close();
			outc.close();
			
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
