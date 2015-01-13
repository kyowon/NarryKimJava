package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import fCLIP.MirGff3FileParser.MiRNA;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import rpf.Scorer;
import util.MatchedFilter;

public class FCLIP_ScorerTrainer {
	private int leftWindowSize = 30;
	private int rightWindowSize = 60;
	private int maxDepthThreshold = 1;
	//private int coverageThreshold = 1;
	
	private ArrayList<double[]> signal3p = null;
	private ArrayList<double[]> noise3p = null;
	private double[][] avgedNoise3p = null;
	private double[][] avgedSignal3p;
	private double[] filter3p;
	
	private ArrayList<double[]> signal5p = null;
	private ArrayList<double[]> noise5p = null;
	private double[][] avgedNoise5p = null;
	private double[][] avgedSignal5p;
	private double[] filter5p;
	
	private String bedfileName;
	private MirGff3FileParser mirParser;
	private  AnnotationFileParser annotationParser;
	private String outParamFile;
	
	public FCLIP_ScorerTrainer(String bedfileName, MirGff3FileParser mirParser, AnnotationFileParser annotationParser, String outParamFile){
		this.bedfileName = bedfileName;
		this.mirParser = mirParser;	
		this.annotationParser = annotationParser;
		this.outParamFile = outParamFile;
	}
	
	public void train(int leftWindowSize, int rightWindowSize, int numberOfNonZeroElements){
		this.leftWindowSize = leftWindowSize;
		this.rightWindowSize = rightWindowSize;
		this.maxDepthThreshold = numberOfNonZeroElements;
		//this.coverageThreshold = coverageThreshold;
		signal5p = new ArrayList<double[]>();
		noise5p = new ArrayList<double[]>();
		
		signal3p = new ArrayList<double[]>();
		noise3p = new ArrayList<double[]>();
		
		for(String contig : annotationParser.getContigs())
			getSignalNNoise(contig);
		
		avgedSignal5p = MatchedFilter.getAvg(signal5p);
		avgedNoise5p = MatchedFilter.getAvg(noise5p);
		filter5p = MatchedFilter.getFilter(avgedSignal5p, noise5p);
		
		avgedSignal3p = MatchedFilter.getAvg(signal3p);
		avgedNoise3p = MatchedFilter.getAvg(noise3p);
		filter3p = MatchedFilter.getFilter(avgedSignal3p, noise3p);
		
		write(outParamFile);
	}
	
	private int getMaxCenterOffset(double[] cov, int window, boolean leftSide){
		int maxOffset = window;
		double max = 0;
		for(int i=0;i<cov.length;i++){
			if(leftSide){
				if(i>window) break;
			}else{
				if(i<=window) continue;
			}
			if(max < cov[i]){
				max = cov[i];
				maxOffset = i;
			}
		}
		return maxOffset - window;
		
	}
	
	private double[] invert(double[] c){
		double[] ic = new double[c.length];
		for(int i=0;i<ic.length;i++){
			ic[i] = c[c.length - i - 1];
		}		
		return ic;
	}
	
	private void getSignalNNoise(String contig){
	//	System.out.println(contig);
		Iterator<MiRNA> iterator = mirParser.getMiRNAIterator(contig);
		if(iterator == null) return;
		double valoffset = .0;
		Bed12Parser bedParser = new Bed12Parser(bedfileName, annotationParser, contig, true);
		while(iterator.hasNext()){
			MiRNA mi = iterator.next();
			//System.out.println(mi);
			boolean isPlusStrand = mi.isPlusStrand();		
			double[] cov5p = bedParser.get5pCoverages(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize));
			
			int maxOffset5p = getMaxCenterOffset(cov5p, leftWindowSize, true);
			cov5p = bedParser.get5pCoverages(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p));
			double[] _cov3p = bedParser.get3pCoverages(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p));
			
			double[] dep5p = bedParser.get5pSignalForfCLIP(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p));
					
			double[] cov3p = bedParser.get3pCoverages(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize));
			int maxOffset3p = getMaxCenterOffset(cov3p, rightWindowSize, false) + 1;
			cov3p = bedParser.get3pCoverages(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize, maxOffset3p));
			double[] _cov5p = bedParser.get5pCoverages(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize, maxOffset3p));
			
			double[] dep3p = bedParser.get3pSignalForfCLIP(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize, maxOffset3p));
			
			//System.out.println(mi.get3pCoordinate(leftWindowSize, rightWindowSize));
			//System.out.println("*" + mi.get5pCoordinate(rightWindowSize, leftWindowSize));

		/*	if(Scorer.sum(cov5p) >= readCountThreshold && Scorer.sum(cov3p) < readCountThreshold){
				System.out.println("* " + mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p));
				for(double c : cov5p) System.out.print(c + " ");
				System.out.println();
				System.out.println("# " + mi.get5pCoordinate(leftWindowSize, rightWindowSize, maxOffset3p));
				for(double c : cov3p) System.out.print(c + " ");
				System.out.println();
			}
		*/	
			
			
			boolean recorded = false;
			if(Scorer.max(cov5p) + Scorer.max(_cov3p) >= maxDepthThreshold){		
				Scorer.normalize(dep5p, (new Random().nextDouble() - .5) * valoffset);
				signal5p.add(dep5p);
				recorded = true;
			}
			if(Scorer.max(cov3p) + Scorer.max(_cov5p) >= maxDepthThreshold){
				Scorer.normalize(dep3p, (new Random().nextDouble() - .5) * valoffset);
				signal3p.add(dep3p);
				recorded = true;
			}
			
			
			if(!recorded) continue;
			
			//System.out.println("signal found");
			HashSet<Integer> offsets = new HashSet<Integer>();
			for(int i=0;i<3;){
				int offset =  (new Random().nextInt(leftWindowSize)) - leftWindowSize/2;
				if(offsets.size() >= 50) break;
				if(offsets.contains(offset)) continue;
				offsets.add(offset);				
				if(offset ==0) continue;
				//System.out.println(offset);
				//offset = isPlusStrand? offset : - offset;
				double[] noiseCov5p =  bedParser.get5pCoverages(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p + offset));
				double[] _noiseCov3p =  bedParser.get3pCoverages(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p + offset));
				
				
				//noiseCov5p = bedCov5PFileParser.getCoverages(mi.getContig(), pos5p + offset + getMaxCenterOffset(noiseCov5p), leftWindowSize, rightWindowSize, isPlusStrand);
				if(Scorer.max(noiseCov5p) + Scorer.max(_noiseCov3p) >= maxDepthThreshold){	
					double[] noiseDep5p = bedParser.get5pSignalForfCLIP(isPlusStrand, mi.get3pCoordinate(leftWindowSize, rightWindowSize, maxOffset5p + offset));
					Scorer.normalize(noiseDep5p, (new Random().nextDouble() - .5) * valoffset);
					//noiseDep5p = MC.multiply(noiseDep5p, 0.5);
					noise5p.add(noiseDep5p);
					i++;
				}
				
				double[] noiseCov3p = bedParser.get3pCoverages(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize, maxOffset3p + offset));
				double[] _noiseCov5p = bedParser.get5pCoverages(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize, maxOffset3p + offset));
				//noiseCov3p = invert(noiseCov3p);
				//noiseCov3p = bedCov5PFileParser.getCoverages(mi.getContig(), pos3p + offset + getMaxCenterOffset(noiseCov3p), leftWindowSize, rightWindowSize, !isPlusStrand);
				if(Scorer.max(noiseCov3p) + Scorer.max(_noiseCov5p) >= maxDepthThreshold){	
					double[] noiseDep3p = bedParser.get3pSignalForfCLIP(isPlusStrand, mi.get5pCoordinate(rightWindowSize, leftWindowSize,maxOffset3p + offset));
					Scorer.normalize(noiseDep3p, (new Random().nextDouble() - .5) * valoffset);
					//noiseDep3p = MC.multiply(noiseDep3p, 0.5);
					noise3p.add(noiseDep3p);
					i++;
				}
				
			//	System.out.println(Scorer.max(noiseCov5p));
			//	System.out.println(Scorer.max(_noiseCov3p));
			//	System.out.println(Scorer.max(noiseCov3p));
			//	System.out.println(Scorer.max(_noiseCov5p));
				
				//if(i>0) System.out.println(i);
			}
			//System.out.println("noise found");
		}
		System.out.println(contig + " " + signal3p.size() + " " + signal5p.size());
//		System.out.println(contig + " " + signal.size() + " " + noise.size());
	}
	
	private void write(String outParamFile){
		try {
			PrintStream out = new PrintStream(outParamFile);
			
			out.println("#LEFT\t"+leftWindowSize);
			out.println("#RIGHT\t"+rightWindowSize);
			out.println("#MAXDEPTHTHRESHOLD\t"+maxDepthThreshold);
			//out.println("#COVTHRESHOLD\t"+coverageThreshold);
			
			out.println("#SIGNAL5P\t"+avgedSignal5p.length);
			
			for(int i=0;i<avgedSignal5p.length;i++){
				out.println(avgedSignal5p[i][0]);
			}
			
			out.println("#SIGNAL3P\t"+avgedSignal3p.length);
			
			for(int i=0;i<avgedSignal3p.length;i++){
				out.println(avgedSignal3p[i][0]);
			}
			
			out.println("#NOISE5P\t"+avgedNoise5p.length);
			
			for(int i=0;i<avgedNoise5p.length;i++){
				out.println(avgedNoise5p[i][0]);
			}

			out.println("#NOISE3P\t"+avgedNoise3p.length);
			
			for(int i=0;i<avgedNoise3p.length;i++){
				out.println(avgedNoise3p[i][0]);
			}

			out.println("#FILTER5P\t"+ filter5p.length);
			
			for(int i=0;i<filter5p.length;i++){
				out.println(filter5p[i]);
			}
			
			out.println("#FILTER3P\t"+ filter3p.length);
			
			for(int i=0;i<filter3p.length;i++){
				out.println(filter3p[i]);
			}
			
			out.close();
			
			} catch (IOException e) {
				e.printStackTrace();
			}
	}
	
	public static void train(String outParam, String bed, MirGff3FileParser mirParser, AnnotationFileParser annotationParser){
		if(new File(outParam).exists()) return;
		FCLIP_ScorerTrainer test = new FCLIP_ScorerTrainer(bed , mirParser, annotationParser, outParam);
		test.train(30, 20, 5);
	}
	
	public static void main(String[] args){
	//	System.out.println(args.length);
	//	for(String arg : args)
	//		System.out.println(arg);
		
		String outParam;
		String bed;
		MirGff3FileParser mirParser;
		AnnotationFileParser annotationParser;
		int leftWindow; // 30 
		int rightWindow; // 20
		int numNonzero; // 5
		
		if(args.length > 0){
		//	System.out.println(args.length);
		//		for(String arg : args)
		//			System.out.println(arg);
			bed = args[1];
			mirParser = new MirGff3FileParser(args[2]);
			annotationParser = new AnnotationFileParser(args[3]);
			outParam = args[0];
			leftWindow = Integer.parseInt(args[4]);;
			rightWindow = Integer.parseInt(args[5]);;
			numNonzero = Integer.parseInt(args[6]);;
		}else{
			bed = "/media/kyowon/Data1/fCLIP/samples/sample4/bed/merged.se.bed";
			mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3");
			annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt");
			outParam = "/media/kyowon/Data1/fCLIP/samples/sample4/erase.txt";
			leftWindow = 30;
			rightWindow = 20;
			numNonzero = 5;
		}
		FCLIP_ScorerTrainer test = new FCLIP_ScorerTrainer(bed , mirParser, annotationParser, outParam);
		test.train(leftWindow, rightWindow, numNonzero);
	}
	
}
