package rpf;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.Bed12Parser;
import parser.ZeroBasedFastaParser;
import rpf.parser.ScoringOutputParser;
import rpf.parser.ScoringOutputParser.ScoredPosition;
import net.sf.samtools.util.BufferedLineReader;

public class Scorer {
	private int leftWindowSize = 30;
	private int rightWindowSize = 60;
	private static int numberOfNonZeroElements = 5;
	//private int coverageThreshold = 1;
	
	private Bed12Parser bedParser;	
	private double[] startFilter;
	private double[] startSignal; // used for quantification
	private double startFilterNorm; // for speed up
	
	private double[] stopFilter;
	private double[] stopSignal; // used for quantification
	private double stopFilterNorm; // for speed up
//	private PolynomialFunction likelihoodFunction;
	private AnnotationFileParser annotationFileParser = null; // if specified, only annotated start positions are considered
	private ZeroBasedFastaParser fastaParser = null;
	
	public Scorer(Bed12Parser bedParser, String paramFile, AnnotationFileParser annotationFileParser, ZeroBasedFastaParser fastaParser){
		this.bedParser = bedParser;
		read(paramFile);		
		this.annotationFileParser = annotationFileParser;
		this.fastaParser = fastaParser;
	}
	
	public AnnotationFileParser getAnnotationFileParser(){
		return annotationFileParser;
	}
	
	public ArrayList<ScoredPosition> getScoredPositions(double scoreThreshold, HashSet<String> allowedCodons){
		ArrayList<ScoredPosition> positions = new ArrayList<ScoringOutputParser.ScoredPosition>();
		positions.addAll(getScoredPositions(scoreThreshold, false, allowedCodons));
		positions.addAll(getScoredPositions(scoreThreshold, true, allowedCodons));
		return positions;
	}
	
	private ArrayList<ScoredPosition> getScoredPositions(double scoreThreshold, boolean isPlusStrand, HashSet<String> allowedCodons){
		//contig = "chr1";
		System.out.println("Scoring for " + bedParser.getContig() + " " + (isPlusStrand? '+' : '-') + " strand");
		Iterator<Integer> iterator;
		iterator = bedParser.getNonZero5pPositionIterator(isPlusStrand);
		int lastConsidered = 0;
		ArrayList<ScoredPosition> positions = new ArrayList<ScoringOutputParser.ScoredPosition>();
		while(iterator.hasNext()){
			int position = iterator.next();
			int	start = isPlusStrand? position - leftWindowSize - 1: position - rightWindowSize - 1;
			int	end = isPlusStrand? position + rightWindowSize + 1: position + leftWindowSize + 1;				
			for(int currentPosition=Math.max(lastConsidered + 1, start);currentPosition<end;currentPosition++){
				//double[] cov = bedCovFileParser.getSqrtCoverages(contig, currentPosition, leftWindowSize, rightWindowSize, isPlusStrand);
				
				ArrayList<String> genomicRegionAndFrameShift = null;
				if(annotationFileParser != null){
					ArrayList<ArrayList<Integer>> coordinates = bedParser.getCoordinates(currentPosition, leftWindowSize, rightWindowSize, isPlusStrand);
					ArrayList<AnnotatedGene> containingGenes = annotationFileParser.getContainingGenes(bedParser.getContig(), isPlusStrand, currentPosition);
					
					for(ArrayList<Integer> coordinate : coordinates){
						
						double startScore = getStartScore(coordinate, isPlusStrand);
						if(startScore <= scoreThreshold) continue;
						//System.out.println((bedParser.getCoordinates("chr1", 4822442, 50, 50, true)));
						
					//	if(currentPosition == 4822442)
					//		System.out.println(currentPosition + " " + Bed12Parser.getSplices(coordinate));
						ArrayList<AnnotatedGene> matchingGenes = annotationFileParser.getMatchingGenes(containingGenes, isPlusStrand, coordinate);
						
						String seq = null;
						if(!fastaParser.containsContig(bedParser.getContig())){
							seq = "N/A";
						}else{
							seq = fastaParser.getSequence(bedParser.getContig(), coordinate.subList(leftWindowSize, coordinate.size()));
							if(!isPlusStrand)						
								seq = ZeroBasedFastaParser.getComplementarySequence(seq, false);
							//System.out.println(currentPosition + " " + coordinate.subList(leftWindowSize, leftWindowSize+3));
						}	
						String codon = seq.substring(0, Math.min(seq.length(), 3));
						
						if(matchingGenes == null){
							if(containingGenes == null){ // intergenic
								genomicRegionAndFrameShift = annotationFileParser.getGenomicRegionNameAndFrameShift(bedParser.getContig(), isPlusStrand, currentPosition, null, false);
								if(!genomicRegionAndFrameShift.get(0).equals("NM_3_UTR")){										
									if(!allowedCodons.contains(codon)) continue;
								}								
							
								ScoredPosition scoredPosition = ScoringOutputParser.getScoredPosition(bedParser.getContig(), currentPosition, coordinate, isPlusStrand, startScore, seq, null, false, genomicRegionAndFrameShift.get(0), genomicRegionAndFrameShift.get(1));
								positions.add(scoredPosition);						
							}else{
								for(AnnotatedGene gene : containingGenes){
									genomicRegionAndFrameShift = annotationFileParser.getGenomicRegionNameAndFrameShift(bedParser.getContig(), isPlusStrand, currentPosition, gene, false);
									if(!genomicRegionAndFrameShift.get(0).equals("NM_3_UTR")){										
										if(!allowedCodons.contains(codon)) continue;
									}								
								
									ScoredPosition scoredPosition = ScoringOutputParser.getScoredPosition(bedParser.getContig(), currentPosition, coordinate, isPlusStrand, startScore, seq, gene, false, genomicRegionAndFrameShift.get(0), genomicRegionAndFrameShift.get(1));
									positions.add(scoredPosition);	
								}
							}							
						}else{
							for(AnnotatedGene gene : matchingGenes){
								genomicRegionAndFrameShift = annotationFileParser.getGenomicRegionNameAndFrameShift(bedParser.getContig(), isPlusStrand, currentPosition, gene, true);
								if(!genomicRegionAndFrameShift.get(0).equals("NM_3_UTR")){										
									if(!allowedCodons.contains(codon)) continue;
								}								
							
								ScoredPosition scoredPosition = ScoringOutputParser.getScoredPosition(bedParser.getContig(), currentPosition, coordinate, isPlusStrand, startScore, seq, gene, gene.isAnnotated(currentPosition), genomicRegionAndFrameShift.get(0), genomicRegionAndFrameShift.get(1));
								positions.add(scoredPosition);						
							}							
						}
					}
				}	
				lastConsidered = currentPosition;				
			}
		}
		return positions;
	}

	
	
	
	
	
	private void read(String paramFile){
		BufferedLineReader in;
		//double[] likelihoodFunctionCoefficients = null;
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
				//}else if(s.startsWith("#NONZERO")){
				//	numberOfNonZeroElements = Integer.parseInt(token[1]);
				//}else if(s.startsWith("COVTHRESHOLD")){
				//	coverageThreshold = Integer.parseInt(token[1]);
				}else if(s.startsWith("#STARTFILTER")){
					startFilter = new double[Integer.parseInt(token[1])];
					mod = 1;
					i = 0;
				}else if(s.startsWith("#STARTSIGNAL")){
					startSignal = new double[Integer.parseInt(token[1])];
					mod = 2;
					i = 0;
					
				}else if(s.startsWith("#STOPFILTER")){
					stopFilter = new double[Integer.parseInt(token[1])];
					mod = 3;
					i = 0;
				}else if(s.startsWith("#STOPSIGNAL")){
					stopSignal = new double[Integer.parseInt(token[1])];
					mod = 4;
					i = 0;
				
			//}else if(s.startsWith("#LIKELIHOOD")){
				//	likelihoodFunctionCoefficients = new double[Integer.parseInt(token[1])];
				//	mod = 3;
			//		i = 0;
				}else if(mod == 1){
					startFilter[i++] = Double.parseDouble(token[0]);
				}else if(mod == 2){
					startSignal[i++] = Double.parseDouble(token[0]);
				}else if(mod == 3){
					stopFilter[i++] = Double.parseDouble(token[0]);
				}else if(mod == 4){
					stopSignal[i++] = Double.parseDouble(token[0]);
				
					//}else if(mod == 3){
				//	likelihoodFunctionCoefficients[i++] = Double.parseDouble(token[0]);
				}
			}
		//	if(likelihoodFunctionCoefficients!=null){
		//		likelihoodFunction = new PolynomialFunction(likelihoodFunctionCoefficients);
		//	}
			in.close();
			startFilterNorm = getNorm(startFilter);
			stopFilterNorm = getNorm(stopFilter);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	//public double getLRScore(double score){ 
	//	if(likelihoodFunction == null) return score;
	//	return likelihoodFunction.value(score);
	//}
	
	//public static double getRawScore(double[] filter, double[] cov, int numberOfNonZeroElements){
	//	double filterNorm = getNorm(filter);
	//	return getRawScore(filter, cov, numberOfNonZeroElements, filterNorm);
	//}
	
	
	private double getRawStartScore(double[] cov){
		return getRawScore(startFilter, cov, startFilterNorm);
	}
	
	private double getRawStopScore(double[] cov){
		return getRawScore(stopFilter, cov, stopFilterNorm);
	}
	
	public double getStartScore(ArrayList<Integer> coordinate, boolean isPlusStrand){
		double[] cov = bedParser.get5pCoverages(isPlusStrand, coordinate);
		return getRawStartScore(cov);
	}
	
	//stop position inclusive..
	public double getStopScore(AnnotatedGene gene, int stopPosition){
		//int stopPostion = -1;
		//ArrayList<ArrayList<Integer>> lps = gene.getLiftOverPositions(stopPosition, leftWindowSize, rightWindowSize, false);
		ArrayList<Integer> coordinate = gene.getLiftOverPositions(stopPosition, leftWindowSize, rightWindowSize, false);
		//for(ArrayList<Integer> lp : lps){
		//	coordinate.addAll(lp);
		//}
			//System.out.println(stopPostion);
		double[] cov = bedParser.get5pCoverages(gene.isPlusStrand(), coordinate);
	//	double score = 0;
	//	if(!considerNumberOfNonZeroElements || numberOfNonZeroElements(cov) >= numberOfNonZeroElements) 
		//	score = ();
		return getRawStopScore(cov);
		
	}
	
	
	//getSqrtCoveragesTillnextStopCodon(String contig, int position, int leftWindowSize, boolean isPlusStrand, int maxLength, ZeroBasedFastaParser fastaParser)
	
	/*public double getStartScore(String contig, int position, boolean isPlusStrand, int maxLength, boolean considerNumberOfNonZeorElements){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;	
		double[] cov = bedCovFileParser.getCoveragesTillnextStopCodon(contig, position, leftWindowSize, isPlusStrand, maxLength, fastaParser);
		double score = 0;
		//System.out.println("hh " + numberOfNonZeroElements(cov));
		if(!considerNumberOfNonZeorElements || numberOfNonZeroElements(cov) >= numberOfNonZeroElements){
			score = (getRawStartScore(cov));
			//
		}
		return score;
	}*/
	
	private static double getRawScore(double[] filter, double[] cov, double filterNorm){
		if(numberOfNonZeroElements(cov) < numberOfNonZeroElements) return 0;
	//	double[] sqrtCov = getSqrtVector(cov);
		//System.out.println("hh ");
		double norm = getNorm(cov);
		//System.out.println(norm);
		if(norm <= 0) return 0;
		
		double[] efilter = getLengthMatchedFilter(filter, cov);
		if(efilter != filter){
			filterNorm = getNorm(efilter);
		}
		
		double ip = getInnerProduct(efilter, cov);
		
		//System.out.println("hh " + ip + " " + filterNorm + " " + norm);
		
		return ip / filterNorm / norm;
	}
	
	private static double[] getLengthMatchedFilter(double[] h, double[] s){
		if(h.length == s.length) return h;
		
		double[] eh = new double[s.length];
		
		if(h.length > s.length){
			for(int i=0; i<eh.length;i++){
				eh[i] = h[i];
			}
			return eh;
		}
		
		for(int i=0; i<h.length;i++){
			eh[i] = h[i];
		}	
		
		double[] hp = new double[9];
		for(int i=0;i<hp.length;i++){
			hp[i] = h[h.length-hp.length+i];
		}
		
	//	hp[0] = h[h.length-3];
	//	hp[1] = h[h.length-2];
	//	hp[2] = h[h.length-1];
		
		double seed = 1.0001;
		double f = 1;
		for(int i=h.length; i<eh.length;i++){
			eh[i] = hp[(i-h.length)%hp.length];
			eh[i] /= f;
			f =seed * f;
		}
		return eh;
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
	/*
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
	*/


	
	public static int numberOfNonZeroElements(double[] v){
		int c = 0;
		for(double e : v) if(e != 0) c++;
		return c;
	}
	
	public static double sum(double[] v){
		double sum = 0;
		for(double s : v) sum += s;
		return sum;
	}
	
	public static double max(double[] v){
		double max = Double.NEGATIVE_INFINITY;
		for(double s : v) max = max < s? s : max;
		return max;
	}
	
	public static double getNorm(double[] v){
		double sum = 0;
		for(double t : v){
			sum += t*t;
			//System.out.print(t + " ");
		}
		//System.out.println();
		return Math.sqrt(sum);
	}
	
	/*public static double[] getSqrtVector(double[] v){
		double[] ret = new double[v.length];
		for(int i=0; i<v.length;i++){
			ret[i] = Math.sqrt(v[i]);
		}
		return ret;
	}*/
	
	private static void addRandomSmallValues(double[] s){
		double min = 1000000;
		for(int i=0;i<s.length;i++){
			if(s[i] <= 0)continue;
			min = min > s[i] ? s[i] : min;
		}
		for(int i=0;i<s.length;i++){
			s[i] += new Random().nextDouble() * min / 10;
		}
	}

	public static void addSmallValuesNNormalize(double[] v){
		addRandomSmallValues(v);
		normalize(v, 0);
	}
	
	public static void normalize(double[] v, double offset){
		//double sum = 0;
		//for(int i=0; i<v.length;i++){
			//v[i] = v[i] * v[i];
		//	sum += v[i];
		//}
		if(offset != 0){
			for(int i=0; i<v.length;i++){
				v[i] += offset;
			}
		}
		
		double norm = getNorm(v);
		if(norm > 0){
			for(int i=0; i<v.length;i++){
				v[i] /= norm;
			}
		}
	}
	
	public static void normalize(double[][] m){
		normalize(m, 0);
	}
	public static void normalize(double[][] m, double offset){
		//double sum = 0;
		//for(int i=0; i<v.length;i++){
			//v[i] = v[i] * v[i];
		//	sum += v[i];
		//}
		double[] v = new double[m.length];
		
		for(int i=0; i<v.length;i++){
			v[i] = m[i][0] + offset;
		}
		
		
		double norm = getNorm(v);
		if(norm > 0){
			for(int i=0; i<v.length;i++){
				v[i] /= norm;
			}
		}
	}
	
	public static double getInnerProduct(double[] h, double[] s){
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
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args){
		AnnotationFileParser annotationFileParser = new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.refFlat.txt");
		Bed12Parser bedParser = new Bed12Parser("/media/kyowon/Data1/RPF_Project/samples/sample3/bed/Harr_C-uncollapsed.bed", annotationFileParser, "chr1");
		
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.fa");
		String paramFile = "/media/kyowon/Data1/RPF_Project/samples/sample3/bed/Harr_C-uncollapsed.bed.param";
		Scorer test = new Scorer(bedParser, paramFile, annotationFileParser, fastaParser);
		HashSet allowedCodons = new HashSet<String>();
		allowedCodons.add("ATG");
		allowedCodons.add("CTG");
		
		
		//System.out.println((bedParser.getCoordinates("chr1", 4822442, 50, 50, true)));
		
		
		//test.scoreNWrite(0.3,allowedCodons, "/media/kyowon/Data1/RPF_Project/samples/sample3/bed/Harr_C-uncollapsed.bed.param.out", false);
		
	}
	
}
