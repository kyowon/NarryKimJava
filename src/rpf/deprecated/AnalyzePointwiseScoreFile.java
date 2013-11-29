package rpf.deprecated;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import net.sf.samtools.util.BufferedLineReader;

public class AnalyzePointwiseScoreFile {
	
	private ArrayList<String> lines;
	
	public AnalyzePointwiseScoreFile(){
		lines = new ArrayList<String>();
	}
	
	
	private void read(String fileName){
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(new FileInputStream(fileName));
			String s;
			while((s=in.readLine())!=null){			
				lines.add(s);
			}
			in.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void getQuantification(String refFlatMapped1, String refFlatMapped2, String outFile, double scoreThreshold){
		HashMap<String, String> t1 = new HashMap<String, String>();
		HashMap<String, String> t2 = new HashMap<String, String>();
		String s;
		try {
			BufferedLineReader refIn1 = new BufferedLineReader(new FileInputStream(refFlatMapped1));
			while((s=refIn1.readLine())!=null){
				String[] token = s.split("\t");
				if(!token[1].startsWith("NM")) continue; //TODO 
				if(token[3].equals("+")){
					String seq = token[token.length-1];
					if(!seq.equals("ATG") && !seq.equals("CTG")) continue;
				}
				if(token[3].equals("-")){
					String seq = token[token.length-1];
					if(!seq.equals("CAT") && !seq.equals("CAG")) continue;
				}
				StringBuffer key = new StringBuffer();
				key.append(token[1]); // gene name
				key.append(token[2]); // chr 
				key.append(token[9]); // pos
				key.append(token[10]); // strand
				String keyString = key.toString();
				//if(!t1.containsKey(keyString)) t1.put(keyString, new ArrayList<String>());
				t1.put(keyString, s);
			}
			refIn1.close();
			
			BufferedLineReader refIn2 = new BufferedLineReader(new FileInputStream(refFlatMapped2)); // yes, I am stupid.
			while((s=refIn2.readLine())!=null){
				String[] token = s.split("\t");
				if(!token[1].startsWith("NM")) continue;
				if(token[3].equals("+")){
					String seq = token[token.length-1];
					if(!seq.equals("ATG") && !seq.equals("CTG")) continue;
				}
				if(token[3].equals("-")){
					String seq = token[token.length-1];
					if(!seq.equals("CAT") && !seq.equals("CAG")) continue;
				}
				StringBuffer key = new StringBuffer();
				key.append(token[1]); // gene name
				key.append(token[2]); // chr
				key.append(token[9]); // pos
				key.append(token[10]); // strand
				String keyString = key.toString();
				//if(!t.containsKey(keyString)) t.put(keyString, new ArrayList<String>());
				t2.put(keyString, s);
			}			
			refIn2.close();
			
			HashMap<Double, ArrayList<String>> outLines = new HashMap<Double, ArrayList<String>>();
			ArrayList<Double> quanRatios = new ArrayList<Double>();
			ArrayList<Double> scoreRatios = new ArrayList<Double>();
			
			HashSet<String> keys = new HashSet<String>(t1.keySet());
			keys.addAll(t2.keySet());
			
			for(String key : keys){
				//ArrayList<String> v = t.get(key);
				//if(v.size() < 2) continue;
				String v1 = t1.get(key);
				String v2 = t2.get(key);
				
				double quan1, quan2, score1, score2;
				String[] token = null;
				
				if(v1 != null){
					String[] token1 = v1.split("\t");
					token = token1;
					score1 = Double.parseDouble(token1[11]);
					quan1 = Double.parseDouble(token1[12]);					
				}else{
					score1 = quan1 = 1e-6;
				}
				
				if(v2 != null){
					String[] token2 = v2.split("\t");
					token = token2;
					score2 = Double.parseDouble(token2[11]);
					quan2 = Double.parseDouble(token2[12]);					
				}else{
					score2 = quan2 = 1e-6;
				}
				//12
				double quanRatio = Math.log10(quan1/quan2);
				
				if(score1 < scoreThreshold && score2 < scoreThreshold) continue;
				
				double scoreRatio = Math.log10(score1/score2);
				if(!quanRatios.contains(quanRatio)) quanRatios.add(quanRatio);
				if(!scoreRatios.contains(scoreRatio)) scoreRatios.add(scoreRatio);
				
				StringBuffer outLine = new StringBuffer();
				for(int i=0;i<11;i++){
					outLine.append(token[i]);
					outLine.append("\t");
				}
				outLine.append(score1);outLine.append("\t");
				outLine.append(quan1);outLine.append("\t");
				outLine.append(score2);outLine.append("\t");
				outLine.append(quan2);outLine.append("\t");
				outLine.append(scoreRatio);outLine.append("\t");				
				outLine.append(quanRatio);outLine.append("\t");
				outLine.append(token[13]);   /// seq
				if(!outLines.containsKey(scoreRatio)) outLines.put(scoreRatio, new ArrayList<String>());
				outLines.get(scoreRatio).add(outLine.toString());
			}
			
			Collections.sort(scoreRatios);
			PrintStream out = new PrintStream(outFile);
			out.println("#GeneName\tName\tChrom\tStrand\ttxStart\ttxEnd\tcdsStart\tcsdEnd"
					+ "\t*\tDetectedPosition\tDetectedStrand\tScore1\tQuan1\tScore2\tQuan2\t"
					+ "\tScore1Score2Ratio\tQuan1Quan2Ratio\tSequence");
			for(double ratio : scoreRatios){
				for(String outLine : outLines.get(ratio)){
					out.println(outLine);
				}
			}		
			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private void generateRefFlatMappedOutput(String refFlatFile, String outFile, double scoreThreshold){
		HashMap<String, HashMap<Long, String>> map = new HashMap<String, HashMap<Long, String>>();
		HashMap<String, ArrayList<Long>> positions = new  HashMap<String, ArrayList<Long>>();
			
		for(String s : lines){
			String[] token = s.split("\t");
			if(Double.parseDouble(token[3]) < scoreThreshold) continue;
			String contig = token[0];
			long position = Long.parseLong(token[1]);			
			if(!map.containsKey(contig)){
				map.put(contig, new HashMap<Long, String>());
				positions.put(contig, new ArrayList<Long>());
			}
			StringBuffer toPut = new StringBuffer();
			for(int i=1;i<token.length;i++){
				toPut.append("\t");
				toPut.append(token[i]);				
			}
			
			map.get(contig).put(position, toPut.toString());
			positions.get(contig).add(position);
		}
		
		for(String contig : positions.keySet()){
			Collections.sort(positions.get(contig));			
		}		
		
		try {
			BufferedLineReader refIn = new BufferedLineReader(new FileInputStream(refFlatFile));
			PrintStream out = new PrintStream(outFile);
			String t;
			
			
			while((t=refIn.readLine())!=null){
				String[] token = t.split("\t");
				String contig = token[2];
				if(!positions.containsKey(contig)) continue;
				long position1 = Long.parseLong(token[4]);
				long position2 = Long.parseLong(token[5]);
				ArrayList<Long> pos = positions.get(contig);
				int index1 = Collections.binarySearch(pos, position1);
				int index2 = Collections.binarySearch(pos, position2);
				
				index1 = index1 < 0? -index1-1:index1;
				index2 = index2 < 0? -index2-1:index2;
				String prefix = token[0] + "\t" + token[1] + "\t" + token[2] + "\t" + token[3] + "\t" + token[4] + "\t" + token[5] + "\t" + token[6] + "\t" + token[7];
				for(int i=index1;i<=index2;i++){
					if(i >= pos.size()) continue;
					long cp = pos.get(i);
					if(cp >= position1 && cp <= position2){
						String tt = map.get(contig).get(cp);
						out.println(prefix + "\t*" + tt);
					}						
				}				
			}		
			out.close();
			refIn.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}

	private void printNAratios(double scoreThreshold, boolean isPlusStrand){
		
		HashMap<String, Integer> counter = new HashMap<String, Integer>();
		for(String s : lines)
		{
			String[] token = s.split("\t");		
			
			if(Double.parseDouble(token[3]) < scoreThreshold) continue;
			if(isPlusStrand && !token[2].equals("+")) continue;
			if(!isPlusStrand && !token[2].equals("-")) continue;
			
			if(!counter.containsKey(token[5])) counter.put(token[5], 0);
			counter.put(token[5], counter.get(token[5])+1);
		}
		
		float sum = 0;
		for(String tt : counter.keySet()){
			sum += counter.get(tt);
		}
		for(String tt : counter.keySet()){
			System.out.println(tt + "\t" + counter.get(tt) + "\t" + String.format("%.2f", 100 * counter.get(tt)/sum));
		}
		
		System.out.println();
	}
	
	public static void main(String[] args){
		String keyword = "NS_Harr10m";
		//if(args!=null) keyword = args[0];
		String covFileprefix = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/" + keyword + ".sorted";
		String refFlatFile = "/media/kyowon/Data1/RPF_Project/data/refFlatHuman.txt";
		String covFilePlus = covFileprefix + ".plus.out";
		String covFileMinus = covFileprefix + ".minus.out";
		AnalyzePointwiseScoreFile test = new AnalyzePointwiseScoreFile();
		test.read(covFilePlus);
		test.read(covFileMinus);
		//test.printNAratios(2, true);		
		//test.printNAratios(2, false);
		
		test.generateRefFlatMappedOutput(refFlatFile, covFileprefix + ".refFlatMapped.txt", -1);
		//System.exit(0);
		test.getQuantification("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NS_Harr10m.sorted.refFlatMapped.txt", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.refFlatMapped.txt", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NSvsThy.txt", 2);
		test.getQuantification("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NS_Harr10m.sorted.refFlatMapped.txt", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10m.sorted.refFlatMapped.txt", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NSvsNoco.txt", 2);
		test.getQuantification("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10m.sorted.refFlatMapped.txt", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.refFlatMapped.txt", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NocovsThy.txt", 2);
		
	}
	
	
}
