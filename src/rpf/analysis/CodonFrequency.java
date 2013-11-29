package rpf.analysis;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;

import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;

public class CodonFrequency {
	
	
	
	
	private static void writeCodonFrequencies(ScoringOutputParser parser, double scoreThreshold, PrintStream out){
		HashMap<String, Integer> freqMap = new HashMap<String, Integer>();
		int sum = 0;
		for(ScoredPosition position : parser.getPositions()){
			if(position.getScore() < scoreThreshold) continue;
			String codon = position.isPlusStrand()? position.getCodon() : getComplementaryCodon(position.getCodon());
			if(!freqMap.containsKey(codon)) freqMap.put(codon, 0);
			freqMap.put(codon, freqMap.get(codon) + 1);
			sum ++;
		}
		out.print(String.format("%.2f", scoreThreshold) + "\t");
		Integer num = freqMap.get("ATG");
		num = num == null? 0 : num;
		sum -= num;
		out.print(num + "\t");
		num = freqMap.get("CTG");
		num = num == null? 0 : num;
		sum -= num;
		out.print(num + "\t");
		num = freqMap.get("AAG");
		num = num == null? 0 : num;
		sum -= num;
		out.print(num + "\t");
		num = freqMap.get("CGC");
		num = num == null? 0 : num;
		sum -= num;
		out.print(num + "\t");
		num = freqMap.get("GTG");
		num = num == null? 0 : num;
		sum -= num;
		out.println(num + "\t" + sum);
		
	}	
	
	private static String getComplementaryCodon(String codon){
		StringBuffer cc = new StringBuffer();
		char[] nas = codon.toCharArray();
		for(int i = nas.length-1;i>=0;i--){
			char na = nas[i];
			cc.append(getComplementaryNA(na));
		}
		return cc.toString();
	}
	
	private static char getComplementaryNA(char na){
		if(na == 'A') return 'T';
		if(na == 'T') return 'A';
		if(na == 'C') return 'G';
		return 'C';
	}
	
	private static void writeHeader(PrintStream out){
		out.println("#ScoreThreshold\tATG\tCTG\tAAG\tCGC\tGTG\tOther");
	}
	
	public static void main(String[] args){
		String inFile = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10mNew.sorted.plus.cov.score.tsv";
		String outFile = inFile + ".freq";
		PrintStream out;
		ScoringOutputParser parser = new ScoringOutputParser(inFile);
		try {
			out = new PrintStream(outFile);		
			writeHeader(out);
			for(double scoreThreshold = 1.3; scoreThreshold < 2.5;scoreThreshold+=.05)
				CodonFrequency.writeCodonFrequencies(parser, scoreThreshold, out);
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
}
