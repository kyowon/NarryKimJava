package rpf.deprecated;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import util.MC;
import net.sf.samtools.util.BufferedLineReader;

public class TrainingDataSet {
	
	private static int leftWindowSize = 60;
	private static int nonZeroPointThreshold = 5;
	private static int rightWindowSize = 120; // to be written in file
	
	private ArrayList<double[]> signal = null;
	private ArrayList<double[]> shiftedSignal = null;
	private ArrayList<double[]> noise = null;
	private String inCovFilePlus;
	private String inCovFileMinus;
	private String paramFile;
	private String annotationFile;
	private HashMap<String, HashMap<Long, Integer>> coverage;
	private double[][] filter;
	private double[][] lr;
	
	public static int getLeftWindowSize() {return leftWindowSize;}
	public static int getRightWindowSize() {return rightWindowSize;}
	public static void setLeftWindowSize(int w){leftWindowSize = w;}
	public static void setRightWindowSize(int w){rightWindowSize = w;}
	public static void setNonZeroPointThreshold(int w){nonZeroPointThreshold = w;}
	
	public static void main(String[] args){
		//TrainingDataSet test = new TrainingDataSet("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NS_Harr10m.sorted.5.cov", 
		//		"/media/kyowon/Data1/RPF_Project/data/refFlatHuman.txt", "/media/kyowon/Data1/RPF_Project/results/param.txt");
		//TrainingDataSet test = new TrainingDataSet("/media/kyowon/Data1/RPF_Project/data/RNA-10_mm9.sorted.cov", "", true);

	}
	
	public TrainingDataSet(String inCovFilePlus, String inCovFileMinus, String paramFile, String annotationFile){
		this.inCovFilePlus = inCovFilePlus;
		this.inCovFileMinus = inCovFileMinus;
		this.paramFile = paramFile;
		this.annotationFile = annotationFile;		
		signal = new ArrayList<double[]>();
		noise = new ArrayList<double[]>();
	}
	
	
	public void train(){ 
		train(true);
		//writeSignal("/media/kyowon/Data1/RPF_Project/results/signal5.m");
		train(false);
		//writeSignal("/media/kyowon/Data1/RPF_Project/results/signal3.m");
		filter = getFilter(signal, noise);
		lr = calculateLikelihood(filter, signal, noise);
		
		write(); 
	}
	
	private void train(boolean isPlusStrand){
		coverage = new HashMap<String, HashMap<Long, Integer>>();	
		String inCovFile = isPlusStrand? inCovFilePlus : inCovFileMinus;
		getAllCoverages(coverage, inCovFile);		
		getSignalNNoise(isPlusStrand);
		
	}
	
	private static double[][] calculateLikelihood(double[][] filter, ArrayList<double[]> signal, ArrayList<double[]> noise){
		double[][] likelihoodScores = new double[100][2];
		ArrayList<Double> signalScores = new ArrayList<Double>();
		ArrayList<Double> noiseScores = new ArrayList<Double>();
		
		for(int i=0;i< signal.size();i++){
			signalScores.add(getFilteredValue(filter, signal.get(i)));
		}
		for(int i=0;i< noise.size();i++){
			noiseScores.add(getFilteredValue(filter, noise.get(i)));
		}		
		Collections.sort(signalScores);
		Collections.sort(noiseScores); 
		
		double min = signalScores.get(0);
		double max = signalScores.get(signalScores.size() - 1);
		double prevlr = -100;
		
		for(int i = 0; i< likelihoodScores.length; i++){
			double v = min + (max-min)*(i)/likelihoodScores.length;
			int l = Collections.binarySearch(signalScores, v);
			l = l<0? -l+1 : l;
			l = signalScores.size() - l;
			l = l<=0? 1 : l;
			int k = Collections.binarySearch(noiseScores, v);
			k = k<0? -k+1 : k;
			k = noiseScores.size() - k;
			k = k<=0? 1 : k;
			double lr =  Math.log10(((double)l/signalScores.size())/((double)k/noiseScores.size()));
			lr = Math.max(lr, prevlr);
			
			prevlr = lr;
			likelihoodScores[i][0] = v;
			likelihoodScores[i][1] = lr;
		}
		return likelihoodScores;
	}
	
	private void writeSignal(String outFile){
		PrintStream out;
		try {
			out = new PrintStream(outFile);
		
			out.println("filter=[");
			for(double[] v : filter){
				out.println(v[0]);
			}
			out.println("];");
			
			out.println("avg=[");
			for(double s[] : getAvg(signal)){
				out.println(s[0]);
			}
			
			out.println("];");
			
			out.println("avgn=[");
			for(double s[] : getAvg(noise)){
				out.println(s[0]);
			}
			
			out.println("];");
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	private void write(){
		try {
			PrintStream out = new PrintStream(paramFile);
			
			out.println("#LEFT\t"+leftWindowSize);
			out.println("#RIGHT\t"+rightWindowSize);
			out.println("#NONZERO\t"+nonZeroPointThreshold);
			
			out.println("#FILTER\t"+filter.length);
				
			for(int i=0;i<filter.length;i++){
				out.println(filter[i][0]);
			}
			
			out.println("#LIKELIHOOD\t" + lr.length);
			for(double[] l : lr){
				out.println(l[0]+"\t"+ l[1]);
			}
			
			out.close();
			
			} catch (IOException e) {
				e.printStackTrace();
			}
	}
	
	public static double getFilteredValue(double[][] h, double[] s){
		return MC.multiply(MC.transpose(h), s)[0];
	}
	
	private static double[][] getFilter(ArrayList<double[]> signal, ArrayList<double[]> noise){
		double[][] s = getAvg(signal);
		double[][] R = getCovMatrix(noise);
		double[][] Rinv = MC.invert(R);
		double[][] Rinvs = MC.multiply(Rinv, s);
		double div = Math.sqrt(MC.multiply(MC.transpose(s),Rinvs)[0][0]);
		double[][] h = MC.multiply(Rinvs, 1/div);
				
		return h;
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
	
	public static void getAllCoverages(HashMap<String, HashMap<Long, Integer>> coverage, String inCovFile){		
		try {
			if(inCovFile == null) return;
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(inCovFile));
			String s;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				if(!coverage.containsKey(token[0]))
					coverage.put(token[0], new HashMap<Long, Integer>());
				HashMap<Long, Integer> cov = coverage.get(token[0]);
				cov.put(Long.parseLong(token[1]), (int)Double.parseDouble(token[2]));
			}
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public static boolean isValidStartPosition(long position, HashMap<Long, Integer> covmap, boolean isPlusStrand){
		int nonZeroCntr = 0;		
		for(int i=0;i<rightWindowSize + leftWindowSize;i++){
			if(covmap.containsKey(position + (isPlusStrand? i : -i)))
				nonZeroCntr++;			
		}						
		return nonZeroCntr > nonZeroPointThreshold;
	}
	
	
	public static ArrayList<Long> getAnnotatedStartPositions(HashMap<Long, Integer> covmap, String refFlat, String contig, boolean isSignal, boolean isPlusStrand) throws IOException{
		ArrayList<Long> validStartPositions = new ArrayList<Long>();
		HashSet<Long> positions = new HashSet<Long> ();
	
		BufferedLineReader in = new BufferedLineReader(new FileInputStream(refFlat));
		String key = isPlusStrand? "+" : "-";
		int loc = isPlusStrand? 6 : 7;
		String s;
		while((s=in.readLine())!=null){
			String[] token = s.split("\t");
			if(!token[3].equals(key) || !token[2].equals(contig)) continue; // fix later
			positions.add(Long.parseLong(token[loc]));			
		}
		
		in.close();
	
		for(long position : positions){
			int repeat = isSignal? 1 : 20;
				for(int i=0;i<repeat;i++){
					int offset = (isSignal? 0 : 150 + new Random().nextInt(60)) - leftWindowSize;			
					if(!isPlusStrand) offset = - (offset + 1);
					if(isValidStartPosition(position + offset, covmap, isPlusStrand)){ // 
						validStartPositions.add(position + offset);				
					}		
			}		
		}
		Collections.sort(validStartPositions);
		return validStartPositions;
	}
	
	private static void normalize(double[] v){
		float sum = 0;
		for(int i=0; i<v.length;i++){
			//v[i] = v[i] * v[i];
			sum += v[i];
		}
		if(sum > 0){
			for(int i=0; i<v.length;i++){
				v[i] /= sum;
			}
		}
	}
	
	private static void getCoverages(ArrayList<Long> positions, HashMap<Long, Integer> covmap, ArrayList<double[]> observations, boolean is5isPlusStrandprime){
		for(long position : positions){
			observations.add(getCoveragesAtOnePoint(position, covmap, is5isPlusStrandprime, true));
		}
	}
	
	public static double[] getCoveragesAtOnePoint(long position, HashMap<Long, Integer> covmap, boolean isPlusStrand, boolean normalize){
		double[] v = new double[rightWindowSize + leftWindowSize + 1];
		for(int i=0; i<v.length;i++){
			double d = 0;
			long pos = position + (isPlusStrand? i : - i);
			if(covmap.containsKey(pos)){
				d = covmap.get(pos);
			}				
			v[i] = d;					
		}
		if(normalize)normalize(v);	
		return v;
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
	
	private void getSignalNNoise(boolean isPlusStrand){
		try {
		//	validIntervalMutiplyFactor = 5;
			//fasta = new IndexedFastaSequenceFile(new File(inFastaFile));
			
		//	fasta.
			//ReferenceSequence t = null;
			//while((t=fasta.nextSequence())!=null)
			//	System.out.println(t);
			
			for(String contig : coverage.keySet()){
				ArrayList<Long> startPositions = new ArrayList<Long>();
				//ArrayList<Long> shiftedStartPositions = new ArrayList<Long>();
				ArrayList<Long> randomPositions = new ArrayList<Long>();
				
				HashMap<Long, Integer> covmap = coverage.get(contig);
				
				for(long start : getAnnotatedStartPositions(covmap, annotationFile, contig, true, isPlusStrand)){ 
					startPositions.add(start); //
					//shiftedStartPositions.add(start + (isPlusStrand? 30 : -30));
				}
				
				for(long start : getAnnotatedStartPositions(covmap, annotationFile, contig, false, isPlusStrand)){		
					randomPositions.add(start); // 		
				}
				//System.out.println(startPositions.size() + "\t" + randomPositions.size());
				getCoverages(startPositions, covmap, signal, isPlusStrand);
				//getCoverages(shiftedStartPositions, covmap, shiftedSignal, isPlusStrand);
				
				getCoverages(randomPositions, covmap, noise, isPlusStrand);				
			}
			//fasta.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
	}
	//
	
	
}
