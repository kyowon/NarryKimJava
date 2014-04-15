package fCLIP;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import fCLIP.MirGff3FileParser.MiRNA;
import parser.AnnotationFileParser;
import parser.BedCovFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import rpf.Scorer;
import rpf.ScorerTrainer;

public class fCLIPTrainer {
	private int leftWindowSize = 30;
	private int rightWindowSize = 200;
	private int numberOfNonZeroElements = 1;
	//private int coverageThreshold = 1;
	
	private ArrayList<double[]> signal = null;
	private ArrayList<double[]> noise = null;
	private double[][] avgedNoise = null;
	
	private BedCovFileParser bedCov5PFileParser;
	private BedCovFileParser bedCov3PFileParser;
	private MirGff3FileParser mirParser;
	private double[][] avgedSignal;
	private double[] filter;
	//private double[] likelihoodFunctionCoefficients;
	private String outParamFile;
	
	public fCLIPTrainer(String bedCov5PFileParser, String bedCov3PFileParser, MirGff3FileParser mirParser, AnnotationFileParser annotationParser, String outParamFile){
		this.bedCov5PFileParser = new BedCovFileParser(bedCov5PFileParser, annotationParser);
		this.bedCov3PFileParser = new BedCovFileParser(bedCov3PFileParser, annotationParser);
		this.mirParser = mirParser;	
		this.outParamFile = outParamFile;
	}
	
	public void train(int leftWindowSize, int rightWindowSize, int numberOfNonZeroElements){
		this.leftWindowSize = leftWindowSize;
		this.rightWindowSize = rightWindowSize;
		this.numberOfNonZeroElements = numberOfNonZeroElements;
		//this.coverageThreshold = coverageThreshold;
		signal = new ArrayList<double[]>();
		noise = new ArrayList<double[]>();
		
		getSignalNNoise();
		
		avgedSignal = ScorerTrainer.getAvg(signal);
		avgedNoise = ScorerTrainer.getAvg(noise);
		filter = ScorerTrainer.getFilter(avgedSignal, noise);
		
		write(outParamFile);
	}
	
	private int getMaxCenterOffset(double[] cov){
		int maxOffset = -1;
		double max = 0;
		for(int i=leftWindowSize;i<cov.length;i++){
			if(max < cov[i]){
				//max = cov[i];
				maxOffset = i;
			}
		}
		return maxOffset - leftWindowSize;
		
	}
	
	private void getSignalNNoise(){
		Iterator<MiRNA> iterator = mirParser.getMiRNAIterator();
		double valoffset = .001;
		while(iterator.hasNext()){
			MiRNA mi = iterator.next();
			boolean isPlusStrand = mi.isPlusStrand();
			//if(isPlusStrand != mi.isPlusStrand()) continue;
			int pos5p = mi.isPlusStrand() ? mi.getStart() : mi.getEnd() - 1;
			int pos3p = mi.isPlusStrand() ? mi.getEnd(): mi.getStart() - 1;
			
			double[] cov5p = bedCov3PFileParser.getCoverages(mi.getContig(), pos5p, leftWindowSize, rightWindowSize, isPlusStrand);
			int maxOffset5p = getMaxCenterOffset(cov5p);
			cov5p = bedCov3PFileParser.getCoverages(mi.getContig(), pos5p + (isPlusStrand? maxOffset5p : -maxOffset5p), leftWindowSize, rightWindowSize, isPlusStrand);
			
			double[] cov3p = bedCov5PFileParser.getCoverages(mi.getContig(), pos3p, leftWindowSize, rightWindowSize, !isPlusStrand);
			int maxOffset3p = getMaxCenterOffset(cov3p);
			cov3p = bedCov5PFileParser.getCoverages(mi.getContig(), pos3p + (isPlusStrand? -maxOffset3p : maxOffset3p), leftWindowSize, rightWindowSize, !isPlusStrand);
			
			if(maxOffset5p >= 0 && Scorer.numberOfNonZeroElements(cov5p) >= numberOfNonZeroElements){				
				Scorer.normalize(cov5p, valoffset);
				signal.add(cov5p);				
			}
			if(maxOffset3p >= 0 && Scorer.numberOfNonZeroElements(cov3p) >= numberOfNonZeroElements){				
				Scorer.normalize(cov3p, valoffset);
				signal.add(cov3p);				
			}
			
			HashSet<Integer> offsets = new HashSet<Integer>();
			for(int i=0;i<200;i++){
				int offset =  (new Random().nextInt(50000)) + 1500;
				if(offsets.contains(offset)) continue;
				offsets.add(offset);
				offset = isPlusStrand? offset : - offset;
				double[] noiseCov5p = bedCov3PFileParser.getCoverages(mi.getContig(), pos5p + offset, leftWindowSize, rightWindowSize, isPlusStrand);
				//noiseCov5p = bedCov5PFileParser.getCoverages(mi.getContig(), pos5p + offset + getMaxCenterOffset(noiseCov5p), leftWindowSize, rightWindowSize, isPlusStrand);
				if(Scorer.numberOfNonZeroElements(noiseCov5p) > numberOfNonZeroElements){					
					Scorer.normalize(noiseCov5p, valoffset);
					noise.add(noiseCov5p);
				}
				
				double[] noiseCov3p = bedCov5PFileParser.getCoverages(mi.getContig(), pos3p + offset, leftWindowSize, rightWindowSize, !isPlusStrand);
				//noiseCov3p = bedCov5PFileParser.getCoverages(mi.getContig(), pos3p + offset + getMaxCenterOffset(noiseCov3p), leftWindowSize, rightWindowSize, !isPlusStrand);
				if(Scorer.numberOfNonZeroElements(noiseCov3p) > numberOfNonZeroElements){					
					Scorer.normalize(noiseCov3p, valoffset);
					noise.add(noiseCov3p);
				}
			}	
		}
		System.out.println(signal.size() + " " + noise.size());
	}
	
	
	
	private void write(String outParamFile){
		try {
			PrintStream out = new PrintStream(outParamFile);
			
			out.println("#LEFT\t"+leftWindowSize);
			out.println("#RIGHT\t"+rightWindowSize);
			out.println("#NONZERO\t"+numberOfNonZeroElements);
			//out.println("#COVTHRESHOLD\t"+coverageThreshold);
			
			out.println("#SIGNAL\t"+avgedSignal.length);
			
			for(int i=0;i<avgedSignal.length;i++){
				out.println(avgedSignal[i][0]);
			}
			
			out.println("#NOISE\t"+avgedNoise.length);
			
			for(int i=0;i<avgedNoise.length;i++){
				out.println(avgedNoise[i][0]);
			}
			
			out.println("#FILTER\t"+ filter.length);
			
			for(int i=0;i<filter.length;i++){
				out.println(filter[i]);
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
	
	public static void main(String[] args){
		fCLIPTrainer test = new fCLIPTrainer("/media/kyowon/Data1/fCLIP/Data/Drosha2-5pends.bedgraph", "/media/kyowon/Data1/fCLIP/Data/Drosha2-3pends.bedgraph", 
				new MirGff3FileParser("/media/kyowon/Data1/fCLIP/Genome/hsa.gff3"), 
				new AnnotationFileParser("/media/kyowon/Data1/fCLIP/Genome/hg19.refFlat.txt"), "/media/kyowon/Data1/fCLIP/Genome/test.param");
		test.train(40, 40, 1);
	}
	
}
