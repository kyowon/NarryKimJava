package rpf;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BedCovFileParser;
import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;
import parser.ZeroBasedFastaParser;
import rpf.analysis.CodonFrequency;
import net.sf.samtools.util.BufferedLineReader;

public class Scorer {
	private int leftWindowSize = 30;
	private int rightWindowSize = 60;
	private int numberOfNonZeroElements = 5;
	//private int coverageThreshold = 1;
	
	private BedCovFileParser bedCovPlusFileParser;
	private BedCovFileParser bedCovMinusFileParser;
	private double[] filter;
	private double[] signal; // used for quantification
	private double filterNorm; // for speed up
	private PolynomialFunction likelihoodFunction;
	private AnnotationFileParser annotationFileParser = null; // if specified, only annotated start positions are considered
	
	public Scorer(String bedCovPlusFile, String bedCovMinusFile, String paramFile){
		bedCovPlusFileParser = new BedCovFileParser(bedCovPlusFile);
		bedCovMinusFileParser = new BedCovFileParser(bedCovMinusFile);		
		read(paramFile);		
	}
	
	public Scorer setAnnotationFileFile(String refFlat){
		annotationFileParser = new AnnotationFileParser(refFlat);
		return this;
	}
	
	public void scoreNWrite(double scoreThreshold, String fastaFile, String outFile){
		if(new File(outFile).exists()){
			System.out.println(outFile + " exists. Scoring is skipped.");
			return;
		}
		try {
			PrintStream out = new PrintStream(outFile);
			scoreNWrite(scoreThreshold, true, fastaFile, out);
			scoreNWrite(scoreThreshold, false, fastaFile, out);
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private void scoreNWrite(double scoreThreshold, boolean isPlusStrand, String fastaFile, PrintStream out){
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser(fastaFile);
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;	
		for(String contig : bedCovFileParser.getContigs()){
			Iterator<Integer> iterator;
			iterator = bedCovFileParser.getNonZeroCoveragePositionIterator(contig);
			int lastConsidered = 0;
			while(iterator.hasNext()){
				int position = iterator.next();
				AnnotatedGene gene = null;
					
				int	start = isPlusStrand? position - leftWindowSize - 1: position - rightWindowSize - 1;
				int	end = isPlusStrand? position + rightWindowSize + 1: position + leftWindowSize + 1;				
			
				for(int currentPosition=Math.max(lastConsidered + 1, start);currentPosition<end;currentPosition++){
					double[] cov = bedCovFileParser.getSqrtCoverages(contig, currentPosition, leftWindowSize, rightWindowSize, isPlusStrand);
					lastConsidered = currentPosition;
					if(numberOfNonZeroElements(cov) < numberOfNonZeroElements) continue;
					double score = getLRScore(getRawScore(cov));
					
					if(score > scoreThreshold){
						String codon;
						if(!fastaParser.containsContig(contig)){
							codon = "N/A";
						}else if(isPlusStrand){
							codon = fastaParser.getSequence(contig, currentPosition, currentPosition+3);
						}else{							
							codon = MakeProteinFasta.getComplementaryCodon(fastaParser.getSequence(contig, currentPosition-2, currentPosition+1));
						}						
						
						boolean isAnnotated = false;
						if(annotationFileParser != null){
							gene = annotationFileParser.getAnnotatedGene(contig, isPlusStrand, currentPosition);
							if(gene != null) isAnnotated = true;
							else gene = annotationFileParser.getContainingGene(contig, isPlusStrand, currentPosition);								
						}						
						ScoredPosition scoredPosition = ScoringOutputParser.getScoredPosition(contig, currentPosition, isPlusStrand, score, codon, gene, isAnnotated);
						out.println(scoredPosition);						
					}						
				}	
			}	
		}
	}
	
	
	
	void writeWindowFilteredOutput(String outFile, String outWindowFile, int window){
		if(new File(outWindowFile).exists()){
			return;
		}
		ScoringOutputParser outParser = new ScoringOutputParser(outFile);
		PrintStream out;
		try {
			out = new PrintStream(outWindowFile);
		
			for(String contig : outParser.getContigs()){
				ArrayList<ScoredPosition> positions = outParser.getPositions(contig);
				Collections.sort(positions);
				for(int i = 0; i < positions.size(); i++) {
			      int rank = 1;
			      
			      ScoredPosition thisPosition = positions.get(i);
			      
			      // move left
			      int prevIndex = i-1;
			      while(prevIndex >= 0) {
			    	ScoredPosition prevPosition = positions.get(prevIndex);
			        if(thisPosition.getPosition() - prevPosition.getPosition() > window)    break;
			        if(prevPosition.getScore() > thisPosition.getScore()) rank++;
			        prevIndex--;
			      }
		
			      // move right
			      int nextIndex = i+1;
			      while(nextIndex < positions.size()) {
			    	ScoredPosition nextPosition = positions.get(nextIndex);
			        if(nextPosition.getPosition() - thisPosition.getPosition() > window)    break;
			        if(nextPosition.getScore() > thisPosition.getScore()) rank++;
			        nextIndex++;
			      }
			    
			      if(rank <= 1) out.println(thisPosition);
				}
			}
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	private void read(String paramFile){
		BufferedLineReader in;
		double[] likelihoodFunctionCoefficients = null;
		try {
			in = new BufferedLineReader(new FileInputStream(paramFile));
			String s;
			int i = 0;
			int mod = 0;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				if(s.startsWith("#RIGHT")){
					rightWindowSize = Integer.parseInt(token[1]);
				}else if(s.startsWith("#LEFT")){
					leftWindowSize = Integer.parseInt(token[1]);
				}else if(s.startsWith("#NONZERO")){
					numberOfNonZeroElements = Integer.parseInt(token[1]);
				//}else if(s.startsWith("COVTHRESHOLD")){
				//	coverageThreshold = Integer.parseInt(token[1]);
				}else if(s.startsWith("#FILTER")){
					filter = new double[Integer.parseInt(token[1])];
					mod = 1;
					i = 0;
				}else if(s.startsWith("#SIGNAL")){
					signal = new double[Integer.parseInt(token[1])];
					mod = 2;
					i = 0;
				}else if(s.startsWith("#LIKELIHOOD")){
					likelihoodFunctionCoefficients = new double[Integer.parseInt(token[1])];
					mod = 3;
					i = 0;
				}else if(mod == 1){
					filter[i++] = Double.parseDouble(token[0]);
				}else if(mod == 2){
					signal[i++] = Double.parseDouble(token[0]);
				}else if(mod == 3){
					likelihoodFunctionCoefficients[i++] = Double.parseDouble(token[0]);
				}
			}
			if(likelihoodFunctionCoefficients!=null){
				likelihoodFunction = new PolynomialFunction(likelihoodFunctionCoefficients);
			}
			in.close();
			filterNorm = getNorm(filter);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	public double getLRScore(double score){ 
		if(likelihoodFunction == null) return score;
		return likelihoodFunction.value(score);
	}
	
	public static double getRawScore(double[] filter, double[] cov, int numberOfNonZeroElements){
		double filterNorm = getNorm(filter);
		return getRawScore(filter, cov, numberOfNonZeroElements, filterNorm);
	}
	
	
	public double getRawScore(double[] cov){
		return getRawScore(filter, cov, numberOfNonZeroElements, filterNorm);
	}
	
	private static double getRawScore(double[] filter, double[] cov, int numberOfNonZeroElements, double filterNorm){
		if(numberOfNonZeroElements(cov) < numberOfNonZeroElements) return 0;
	//	double[] sqrtCov = getSqrtVector(cov);
		double norm = getNorm(cov);
		if(norm <= 0) return 0;
		return getInnerProduct(filter, cov) / filterNorm / norm;
	}
	
	public int getLeftWindowSize() {
		return leftWindowSize;
	}

	public int getRightWindowSize() {
		return rightWindowSize;
	}
	
	//public int getCoverageThreshold(){
	//	return coverageThreshold;
	//}
	
	public double getQuantity(double[] cov){
		return getQuantity(cov, false);
	}
	
	public double getQuantity(double[] cov, boolean excludeStartPosition){
		double[] signalFilter;
		if(excludeStartPosition){
			signalFilter = new double[signal.length];
			for(int i=0;i<signal.length;i++){
				if(i>=leftWindowSize-12 && i<=leftWindowSize+12) continue; //TODO test
				signalFilter[i] = signal[i];
			}
		}else signalFilter = signal;
		//double[] sqrtCov = getSqrtVector(cov);
		return getInnerProduct(signalFilter, cov);
	}
	
	
	public static int numberOfNonZeroElements(double[] v){
		int c = 0;
		for(double e : v) if(e > 0) c++;
		return c;
	}
	
	public static double getNorm(double[] v){
		double sum = 0;
		for(double t : v)
			sum += t*t;
		return Math.sqrt(sum);
	}
	
	/*public static double[] getSqrtVector(double[] v){
		double[] ret = new double[v.length];
		for(int i=0; i<v.length;i++){
			ret[i] = Math.sqrt(v[i]);
		}
		return ret;
	}*/
	
	public static void normalize(double[] v){
		//double sum = 0;
		//for(int i=0; i<v.length;i++){
			//v[i] = v[i] * v[i];
		//	sum += v[i];
		//}
		double norm = getNorm(v);
		if(norm > 0){
			for(int i=0; i<v.length;i++){
				v[i] /= norm;
			}
		}
	}
	
	private static double getInnerProduct(double[] h, double[] s){
		double sum = 0;
		for(int i=0; i<Math.min(h.length, s.length);i++){
			sum += h[i] * s[i];
		}
		return sum;
	}
	
	static void writeCodonFrequencies(String outFile, String codonOutFile){
		HashMap<String, Double> freqMapPlus = new HashMap<String, Double>();
		HashMap<String, Double> freqMapMinus = new HashMap<String, Double>();
		ArrayList<Double> freqPlus = new ArrayList<Double>();
		ArrayList<Double> freqMinus = new ArrayList<Double>();
		
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(outFile));
			String s;
			HashMap<String, Double> freqMap;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				boolean isPlusStrand = token[2].equals("+");
				String codon = token[4];
				if(isPlusStrand) freqMap = freqMapPlus;
				else freqMap = freqMapMinus;
				
				if(!freqMap.containsKey(codon))freqMap.put(codon, 0.0);
				freqMap.put(codon, freqMap.get(codon) + 1.0);
			}
			in.close();
			PrintStream out = new PrintStream(codonOutFile);
			out.println("#Plus strand");
			double sum = 0;
			for(String codon : freqMapPlus.keySet()){
				sum += freqMapPlus.get(codon);
			}
			for(String codon : freqMapPlus.keySet()){
				double v = freqMapPlus.get(codon)/sum;
				freqMapPlus.put(codon, v);
				if(!freqPlus.contains(v)) freqPlus.add(v);
			}
			Collections.sort(freqPlus, Collections.reverseOrder());
			for(double v : freqPlus){
				for(String codon : freqMapPlus.keySet()){
					if(freqMapPlus.get(codon) == v){
						out.println(codon + "\t" + v);
					}
				}					
			}
			out.println("#Minus strand");
			sum = 0;
			for(String codon : freqMapMinus.keySet()){
				sum += freqMapMinus.get(codon);
			}
			for(String codon : freqMapMinus.keySet()){
				double v = freqMapMinus.get(codon)/sum;
				freqMapMinus.put(codon, v);
				if(!freqMinus.contains(v)) freqMinus.add(v);
			}
			Collections.sort(freqMinus, Collections.reverseOrder());
			for(double v : freqMinus){
				for(String codon : freqMapMinus.keySet()){
					if(freqMapMinus.get(codon) == v){
						out.println(codon + "\t" + v);
					}
				}					
			}
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args){
		String keyword = "Thy_Harr10m";
		String covFileprefix = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/" + keyword + ".sorted";
		String covFilePlus = covFileprefix + ".plus.cov";
		String covFileMinus = covFileprefix + ".minus.cov";
		String paramFile = covFileprefix + ".param";
		String outFile = covFileprefix + ".out";
		String fasta = "/media/kyowon/Data1/RPF_Project/data/hg19.fa";
		String refFlat = "/media/kyowon/Data1/RPF_Project/data/refFlatHuman.txt";
		
		//Scorer test = new Scorer(covFilePlus, covFileMinus, paramFile).setAnnotationFileFile(refFlat);
		//test.scoreNWrite(2, fasta, outFile);
		//test.writeWindowFilteredOutput(outFile, covFileprefix + ".windowed.out", 50);
		Scorer.writeCodonFrequencies("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.plus.cov.score.tsv2.windowed.tsv", 
				"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.plus.cov.score.tsv.windowed.tsv" + ".freq");
	}
	
}
