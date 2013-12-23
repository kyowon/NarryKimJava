package rpf.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import parser.AnnotationFileParser;
import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;

public class GenomicRegionComposition {

	public static void main(String[] args) {
		double scoreThreshold = 2.0;
		ScoringOutputParser scoringOutputParser = new ScoringOutputParser("/media/kyowon/Data1/RPF_Project/samples/sample1/coverages5/Thy_Harr_10mHsum-uncollapsed.plus.cov.score.tsv");
		AnnotationFileParser annotationFileParser = new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/refFlatHuman.txt");
		HashMap<String, ArrayList<ScoredPosition>> positions = new HashMap<String, ArrayList<ScoredPosition>>();
		
		int sum = 0;
		for(ScoredPosition position : scoringOutputParser.getPositions()){
			if(position.getScore() < scoreThreshold) continue;
			String name = annotationFileParser.getGenomicRegionName(position.getContig(), position.isPlusStrand(), position.getPosition());
			if(!positions.containsKey(name)) positions.put(name, new ArrayList<ScoredPosition>());
			positions.get(name).add(position);
			sum ++;
		}
		
		ArrayList<String> keys = new ArrayList<String>(positions.keySet());		
		Collections.sort(keys);
		for(String key : keys){
			System.out.println(key + "\t" + String.format("%.1f", (double)positions.get(key).size()/sum*100));
			printFreq(getCodonFrequency(positions.get(key)));
		}		
	}
	
	private static void printFreq(HashMap<String, Double> freq){
		ArrayList<Double> n = new ArrayList<Double>(freq.values());
		Collections.sort(n, Collections.reverseOrder());
		HashSet<String> printedCodons = new HashSet<String>();
		for(int i=0; i<Math.min(4,  n.size());i++){
			for(String codon : freq.keySet()){
				if(freq.get(codon) == n.get(i)){
					if(!printedCodons.contains(codon)) System.out.print("\t" + codon + ": " + String.format("%.2f", n.get(i)));
					printedCodons.add(codon);
				}
			}
		}
		System.out.println();
	}
	
	private static HashMap<String, Double> getCodonFrequency(ArrayList<ScoredPosition> positions){
		HashMap<String, Double> freq = new HashMap<String, Double>();
		double sum = 0;
		for(ScoredPosition position : positions){
			String codon = position.getCodon();
			if(!freq.containsKey(codon)) freq.put(codon, 0.0);
			freq.put(codon, freq.get(codon)+1);
			sum++;
		}
		for(String codon : freq.keySet()){
			if(sum > 0)
				freq.put(codon, 100*freq.get(codon)/sum);
		}
		
		return freq;
	}
	
}
