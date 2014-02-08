package rpf;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction.Parametric;
import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BedCovFileParser;
import util.MC;

public class MatchedFilterTrainier {

	private int leftWindowSize = 30;
	private int rightWindowSize = 200;
	private int numberOfNonZeroElements = 7;
	//private int coverageThreshold = 1;
	
	private ArrayList<double[]> signal = null;
	private ArrayList<double[]> noise = null;
	private BedCovFileParser bedCovPlusFileParser;
	private BedCovFileParser bedCovMinusFileParser;
	private AnnotationFileParser annotationFileParser;
	private double[][] avgedSignal;
	private double[] filter;
	private double[] likelihoodFunctionCoefficients;
	private String outParamFile;
	
	public MatchedFilterTrainier(String bedCovPlusFile, String bedCovMinusFile, AnnotationFileParser annotationFileParser, String outParamFile){
		bedCovPlusFileParser = new BedCovFileParser(bedCovPlusFile, annotationFileParser);
		bedCovMinusFileParser = new BedCovFileParser(bedCovMinusFile, annotationFileParser);
		this.annotationFileParser = annotationFileParser;	
		this.outParamFile = outParamFile;
	}
		
	public void train(int leftWindowSize, int rightWindowSize, int numberOfNonZeroElements){
		this.leftWindowSize = leftWindowSize;
		this.rightWindowSize = rightWindowSize;
		this.numberOfNonZeroElements = numberOfNonZeroElements;
		//this.coverageThreshold = coverageThreshold;
		signal = new ArrayList<double[]>();
		noise = new ArrayList<double[]>();
		getSignalNNoise(true);
		getSignalNNoise(false);
		avgedSignal = getAvg(signal);
		filter = getFilter(avgedSignal, noise);
		likelihoodFunctionCoefficients = calculateLikelihoodFunctionCoefficients(filter, signal, noise, numberOfNonZeroElements);
		write(outParamFile);
	}
	
	private void getSignalNNoise(boolean isPlusStrand){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		Iterator<AnnotatedGene> iterator = annotationFileParser.getAnnotatedGeneIterator();
		while(iterator.hasNext()){
			AnnotatedGene gene = iterator.next();
			if(isPlusStrand != gene.isPlusStrand()) continue;
			int position = isPlusStrand? gene.getCdsStart() : gene.getCdsEnd() - 1;
			
			double[] cov = bedCovFileParser.getSqrtCoverages(gene.getContig(), position, leftWindowSize, rightWindowSize, isPlusStrand);
			
			if(Scorer.numberOfNonZeroElements(cov) >= numberOfNonZeroElements){
				//double[] sqrtCov = Scorer.getSqrtVector(cov);				
				Scorer.normalize(cov);
				signal.add(cov);				
			}
			HashSet<Integer> offsets = new HashSet<Integer>();
			for(int i=0;i<200;i++){
				int offset =  (new Random().nextInt(500)) + 150;
				if(offsets.contains(offset)) continue;
				offsets.add(offset);
				offset = isPlusStrand? offset : - offset;
				double[] noisyCov = bedCovFileParser.getSqrtCoverages(gene.getContig(), position + offset, leftWindowSize, rightWindowSize, isPlusStrand);
				if(Scorer.numberOfNonZeroElements(noisyCov) >= numberOfNonZeroElements){					
					Scorer.normalize(noisyCov);
					noise.add(noisyCov);	
				}
			}	
		}
	}
	
	private void write(String outParamFile){
		try {
			PrintStream out = new PrintStream(outParamFile);
			
			out.println("#LEFT\t"+leftWindowSize);
			out.println("#RIGHT\t"+rightWindowSize);
			out.println("#NONZERO\t"+numberOfNonZeroElements);
			//out.println("#COVTHRESHOLD\t"+coverageThreshold);
			
			out.println("#FILTER\t"+filter.length);
				
			for(int i=0;i<filter.length;i++){
				out.println(filter[i]);
			}
			
			out.println("#SIGNAL\t"+avgedSignal.length);
			
			for(int i=0;i<avgedSignal.length;i++){
				out.println(avgedSignal[i][0]);
			}
			
			out.println("#LIKELIHOOD\t" + likelihoodFunctionCoefficients.length);
			for(double c : likelihoodFunctionCoefficients){
				out.println(c);
			}
			
			out.close();
			
			} catch (IOException e) {
				e.printStackTrace();
			}
	}
	
	private static double[][] getAvg(ArrayList<double[]> observations){
		if(observations.isEmpty()) return null;
		double[] avg = new double[observations.get(0).length];
		for(double[] o : observations){
			for(int i=0;i<avg.length;i++){
				avg[i] += o[i];
			}
		}
		float sum = 0;
		for(double v : avg){
			sum += v;
		}
		for(int i=0;i<avg.length;i++){
			avg[i] /= sum;
		}
		
		double[][] ret = new double[avg.length][1];
		for(int i=0; i<avg.length;i++){
			ret[i][0] = avg[i];
		}
		
		return ret;
	}
	
	private static double[] getFilter(double[][] avgedSignal, ArrayList<double[]> noise){		
		double[][] R = getCovMatrix(noise);
		double[][] Rinv = MC.invert(R);
		double[][] Rinvs = MC.multiply(Rinv, avgedSignal);
		double div = Math.sqrt(MC.multiply(MC.transpose(avgedSignal),Rinvs)[0][0]);
		double[][] h = MC.multiply(Rinvs, 1/div);
		double[] filter = new double[h.length];
		for(int i=0;i<filter.length;i++)
			filter[i] = h[i][0];
		return filter;
	}
	
	
	private static double[][] getCovMatrix(ArrayList<double[]> noise){
		double num = 0;
		double[][] cov = new double[noise.get(0).length][noise.get(0).length];
		for(double[] v : noise){
			double[][] arg = new double[v.length][1];
			for(int i=0;i<v.length;i++){
				arg[i][0] = v[i];
			}
			num ++;
			cov = MC.sum(cov, MC.multiply(arg, MC.transpose(arg)),1);
		}
		cov = MC.multiply(cov, 1/num);
		return cov;
	}
	
	private static double[] calculateLikelihoodFunctionCoefficients(double[] filter, ArrayList<double[]> signal, ArrayList<double[]> noise, int numberOfNonZeroElements){
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
	
	/*public static void main(String[] args){
		String keyword =  "RPF6_NS_RPF_1";
		String annotationkey = "uORF";
		String covFileprefix = "/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/" + keyword + "-uncollapsed";
		String covFilePlus = covFileprefix + ".plus.cov";
		String covFileMinus = covFileprefix + ".minus.cov";
		String paramFile = covFileprefix + ".test.param";
		String annotationFile = "/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/"+keyword.substring(0, keyword.indexOf('_'))+"_"+ annotationkey + "1.5.txt";
		annotationFile = "/media/kyowon/Data1/RPF_Project/genomes/refFlatHuman.txt";
		
		MatchedFilterTrainier test = new MatchedFilterTrainier(covFilePlus, covFileMinus, annotationFile, annotationFile+ keyword + ".param");
		test.train(30, 50, 7);
	}*/
	
}
