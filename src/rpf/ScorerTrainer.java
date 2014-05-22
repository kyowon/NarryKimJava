package rpf;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.Bed12Parser;
import util.MatchedFilter;

public class ScorerTrainer {

	private int leftWindowSize = 30;
	private int rightWindowSize = 200;
	private int numberOfNonZeroElements = 7;
	//private int coverageThreshold = 1;
	
	private ArrayList<double[]> startSignal = null;
	private ArrayList<double[]> stopSignal = null;
	
	private ArrayList<double[]> startNoise = null;
	private ArrayList<double[]> stopNoise = null;
	
	private AnnotationFileParser annotationFileParser = null;
	//private BedCovFileParser bedCovPlusFileParser;
	//private BedCovFileParser bedCovMinusFileParser;
	private String bedFileName;
	private double[][] avgedStartSignal;
	private double[][] avgedStopSignal;
	
	private double[] startFilter;
	private double[] stopFilter;
	//private double[] likelihoodFunctionCoefficients;
	private String outParamFile;
	
	public ScorerTrainer(String bedFileName, AnnotationFileParser annotationFileParser, String outParamFile){
		this.bedFileName = bedFileName;
		this.annotationFileParser = annotationFileParser;	
		this.outParamFile = outParamFile;
	}
		
	public void train(int leftWindowSize, int rightWindowSize, int numberOfNonZeroElements){
		this.leftWindowSize = leftWindowSize;
		this.rightWindowSize = rightWindowSize;
		this.numberOfNonZeroElements = numberOfNonZeroElements;
		//this.coverageThreshold = coverageThreshold;
		startSignal = new ArrayList<double[]>();
		stopSignal = new ArrayList<double[]>();
		startNoise = new ArrayList<double[]>();
		stopNoise = new ArrayList<double[]>();
		
		for(String contig : annotationFileParser.getContigs()){
			getSignalNNoise(true, contig);
			getSignalNNoise(false, contig);
		}
		avgedStartSignal = MatchedFilter.getAvg(startSignal);
		startFilter = MatchedFilter.getFilter(avgedStartSignal, startNoise);
		
		avgedStopSignal = MatchedFilter.getAvg(stopSignal);
		stopFilter = MatchedFilter.getFilter(avgedStopSignal, stopNoise);
		//likelihoodFunctionCoefficients = calculateLikelihoodFunctionCoefficients(filter, signal, noise, numberOfNonZeroElements);
		write(outParamFile);
	}
	
	private void getSignalNNoise(boolean isPlusStrand, String contig){
		Iterator<AnnotatedGene> iterator = annotationFileParser.getAnnotatedGeneIterator(contig);
		Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationFileParser, contig);
		System.out.print("Training for " + contig + " " + (isPlusStrand? '+' : '-') + " strand");
		while(iterator.hasNext()){
			AnnotatedGene gene = iterator.next();
			if(isPlusStrand != gene.isPlusStrand()) continue;
			int startPosition = isPlusStrand? gene.getCdsStart() : gene.getCdsEnd() - 1;
			int stopPosition = isPlusStrand? gene.getCdsEnd() - 1: gene.getCdsStart();
					
			ArrayList<Integer> startCoordinate = gene.getLiftOverPositions(startPosition, leftWindowSize, rightWindowSize, false);
			double[] startCov = bedParser.get5pCoverages(isPlusStrand, startCoordinate);
			ArrayList<Integer> stopCoordinate = gene.getLiftOverPositions(stopPosition, leftWindowSize, rightWindowSize, false);		
			double[] stopCov = bedParser.get5pCoverages(isPlusStrand, stopCoordinate);
			
			if(Scorer.numberOfNonZeroElements(startCov) >= numberOfNonZeroElements){
				//double[] sqrtCov = Scorer.getSqrtVector(cov);				
				Scorer.addSmallValuesNNormalize(startCov);
				startSignal.add(startCov);				
			}
			
			if(Scorer.numberOfNonZeroElements(stopCov) >= numberOfNonZeroElements){
				//double[] sqrtCov = Scorer.getSqrtVector(cov);				
				Scorer.addSmallValuesNNormalize(stopCov);
				stopSignal.add(stopCov);				
			}
			
			HashSet<Integer> offsets = new HashSet<Integer>();
			for(int i=0;i<50;i++){
				int offset =  (new Random().nextInt(500)) + 150;
				if(offsets.contains(offset)) continue;
				offsets.add(offset);
				offset = isPlusStrand? offset : - offset;
				
				ArrayList<Integer> startNoiseCoordinate = gene.getLiftOverPositions(startPosition + offset, leftWindowSize, rightWindowSize, false);
				double[] startNoiseCov = bedParser.get5pCoverages(isPlusStrand, startNoiseCoordinate);
				
				if(Scorer.numberOfNonZeroElements(startNoiseCov) >= numberOfNonZeroElements){					
					Scorer.addSmallValuesNNormalize(startNoiseCov);
					startNoise.add(startNoiseCov);	
				}
				
				ArrayList<Integer> stopNoiseCoordinate = gene.getLiftOverPositions(startPosition + offset, leftWindowSize, rightWindowSize, false);
				double[] stopNoiseCov = bedParser.get5pCoverages(isPlusStrand, stopNoiseCoordinate);
				if(Scorer.numberOfNonZeroElements(stopNoiseCov) >= numberOfNonZeroElements){					
					Scorer.addSmallValuesNNormalize(stopNoiseCov);
					stopNoise.add(stopNoiseCov);	
				}
			}	
		}
		System.out.println("... done");
	}
	
	
	private void write(String outParamFile){
		try {
			PrintStream out = new PrintStream(outParamFile);
			
			out.println("#LEFT\t"+leftWindowSize);
			out.println("#RIGHT\t"+rightWindowSize);
			out.println("#NONZERO\t"+numberOfNonZeroElements);
			//out.println("#COVTHRESHOLD\t"+coverageThreshold);
			
			out.println("#STARTFILTER\t"+startFilter.length);
				
			for(int i=0;i<startFilter.length;i++){
				out.println(startFilter[i]);
			}
			
			out.println("#STARTSIGNAL\t"+avgedStartSignal.length);
			
			for(int i=0;i<avgedStartSignal.length;i++){
				out.println(avgedStartSignal[i][0]);
			}
			
			out.println("#STOPFILTER\t"+stopFilter.length);
			
			for(int i=0;i<stopFilter.length;i++){
				out.println(stopFilter[i]);
			}
			
			out.println("#STOPSIGNAL\t"+avgedStopSignal.length);
			
			for(int i=0;i<avgedStopSignal.length;i++){
				out.println(avgedStopSignal[i][0]);
			}
			
		//	out.println("#LIKELIHOOD\t" + likelihoodFunctionCoefficients.length);
		//	for(double c : likelihoodFunctionCoefficients){
		//		out.println(c);
		//	}
			
			out.close();
			
			} catch (IOException e) {
				e.printStackTrace();
			}
	}
	
	/*private static double[] calculateLikelihoodFunctionCoefficients(double[] filter, ArrayList<double[]> signal, ArrayList<double[]> noise, int numberOfNonZeroElements){
		double[] likelihoodFunctionCoefficients = new double[4];
		ArrayList<Double> signalScores = new ArrayList<Double>();
		ArrayList<Double> noiseScores = new ArrayList<Double>();
		CurveFitter<Parametric> fitter = new CurveFitter<Parametric>(new LevenbergMarquardtOptimizer());
		
		for(int i=0;i< signal.size();i++){
			signalScores.add(Scorer.getRawScore(filter, signal.get(i), numberOfNonZeroElements));
		}
		for(int i=0;i< noise.size();i++){
			noiseScores.add(Scorer.getRawScore(filter, noise.get(i), numberOfNonZeroElements));
		}		
		Collections.sort(signalScores);
		Collections.sort(noiseScores); 
		
		double min = signalScores.get(0);
		double max = signalScores.get(signalScores.size() - 1);
		double prevMax = -10000;
		for(int i = 0; i< likelihoodFunctionCoefficients.length; i++){
			double v = min + (max-min)*(i)/likelihoodFunctionCoefficients.length;
			int l = Collections.binarySearch(signalScores, v);
			l = l<0? -l+1 : l;
			l = signalScores.size() - l;
			l = l<=0? 1 : l;
			int k = Collections.binarySearch(noiseScores, v);
			k = k<0? -k+1 : k;
			k = noiseScores.size() - k;
			k = k<=0? 1 : k;
			double lr =  Math.log10(((double)l/signalScores.size())/((double)k/noiseScores.size()));
			if(prevMax > lr){
				lr = prevMax;
			}else{
				prevMax = lr;
			}				
			
			fitter.addObservedPoint(v, lr);
		
		}	
		likelihoodFunctionCoefficients = fitter.fit(new PolynomialFunction.Parametric(), likelihoodFunctionCoefficients);
			
		return likelihoodFunctionCoefficients;
	}
	*/
	public static void main(String[] args){
	//public ScorerTrainer(Bed12Parser bedParser, AnnotationFileParser annotationFileParser, String outParamFile){
	
		String annotationFile = "/media/kyowon/Data1/RPF_Project/genomes/hg19.refFlat.txt";
		ScorerTrainer test = new ScorerTrainer(annotationFile, new AnnotationFileParser(annotationFile), "/media/kyowon/Data1/RPF_Project/samples/sample1/bed/Noco_Harr_10mHsum-uncollapsed.bed.param");
		test.train(90, 300, 30);
	}
	
	
	
}
