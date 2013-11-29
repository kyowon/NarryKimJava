package rpf.deprecated;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.BufferedLineReader;

public class GeneratePositionwiseScore {
	
	private double[][] filter;
	private double[][] lr;
	private String inCovFilePlus;
	private String inCovFileMinus;
	private String paramFile;
	private String fastaFile;
	
	public GeneratePositionwiseScore(String inCovFilePlus, String inCovFileMinus, String paramFile, String fastaFile){
		this.inCovFilePlus = inCovFilePlus;
		this.inCovFileMinus = inCovFileMinus;
		this.paramFile = paramFile;
		this.fastaFile = fastaFile;
		read();
	}
	
	public void read(){
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(new FileInputStream(paramFile));
			String s;
			int i = 0;
			int mod = 0;
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				//System.out.println(s);
				if(s.startsWith("#RIGHT")){
					TrainingDataSet.setRightWindowSize(Integer.parseInt(token[1]));
				}else if(s.startsWith("#LEFT")){
					TrainingDataSet.setLeftWindowSize(Integer.parseInt(token[1]));
				}else if(s.startsWith("#NONZERO")){
					TrainingDataSet.setNonZeroPointThreshold(Integer.parseInt(token[1]));
				}else if(s.startsWith("#FILTER")){
					filter = new double[Integer.parseInt(token[1])][1];
					mod = 1;
					i = 0;
				}else if(s.startsWith("#LIKELIHOOD")){
					lr = new double[Integer.parseInt(token[1])][2];
					mod = 4;
					i = 0;
				}else if(mod == 1){
					filter[i++][0] = Double.parseDouble(token[0]);
				}else if(mod == 4){
					lr[i][0] = Double.parseDouble(token[0]);
					lr[i++][1] = Double.parseDouble(token[1]);
				}
			}
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	public void run(double scoreThreshold, String outFile, boolean isPlusStrand){
		HashMap<String, HashMap<Long, Integer>> coverage = new HashMap<String, HashMap<Long, Integer>>();		
		try {
			IndexedFastaSequenceFile f = new IndexedFastaSequenceFile(new File(fastaFile));
			PrintStream out = new PrintStream(outFile);
			TrainingDataSet.getAllCoverages(coverage, (isPlusStrand? inCovFilePlus : inCovFileMinus));			
			ReferenceSequence r;
			while((r = f.nextSequence())!=null){	
				String contig = r.getName();
				if(!coverage.containsKey(contig)) continue;
				HashMap<Long, Double> scores = new HashMap<Long, Double>();
				HashMap<Long, Double> filteredValue = new HashMap<Long, Double>();
				
				HashMap<Long, Integer> covmap = coverage.get(contig);
				
				HashSet<Long> positionsToConsider = new HashSet<Long>();
				
				for(long position : covmap.keySet()){
					for(long offset = - TrainingDataSet.getLeftWindowSize();offset<=0;offset++){
						long pos = position + (isPlusStrand ? offset : - offset) ;
						positionsToConsider.add(pos);
					}
				}
				
				for(long position : positionsToConsider){ // TODO fix later too much memory..
					double score = getLRScore(getRawScore(position, covmap, filter, isPlusStrand));
					if(score > scoreThreshold){
						scores.put(position, score);
						filteredValue.put(position, getUnnormalizedFilteredValue(position, covmap, filter, isPlusStrand));
					}
				}				
				write(contig, f, scores, filteredValue, out, isPlusStrand);
			}
			f.close();
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static HashMap<Long, Double> windowFilter(HashMap<Long, Double> subScores, int lwindow, int rwindow){
		HashMap<Long, Double> ret = new HashMap<Long, Double>();
		ArrayList<Long> positions = new ArrayList<Long>(subScores.keySet());
		Collections.sort(positions);
		
	    for(int i = 0; i < positions.size(); i++) {
	      int rank = 1;
	      
	      long thisPosition = positions.get(i);
	      double thisScore = subScores.get(thisPosition);
	      
	      // move left
	      int prevIndex = i-1;
	      while(prevIndex >= 0) {
	        long prevPosition = positions.get(prevIndex);
	        if(thisPosition - prevPosition > lwindow)    break;
	        if(subScores.get(prevPosition) > thisScore)            rank++;
	        prevIndex--;
	      }

	      // move right
	      int nextIndex = i+1;
	      while(nextIndex < positions.size()) {
	        long nextPosition = positions.get(nextIndex);
	        if(nextPosition - thisPosition > rwindow)    break;
	        if(subScores.get(nextPosition) > thisScore)            rank++;
	        nextIndex++;
	      }
	    
	      if(rank <= 1) ret.put(thisPosition, thisScore);
	    }
	    return ret;
	}
	
	
	
	private static double getRawScore(long position, HashMap<Long, Integer> covmap, double[][] filter, boolean isPlusStrand){ // positionwise score
		if(TrainingDataSet.isValidStartPosition(position, covmap, isPlusStrand)){
			double[] s = TrainingDataSet.getCoveragesAtOnePoint(position, covmap, isPlusStrand, true);		
			return TrainingDataSet.getFilteredValue(filter, s);
		}else return -10;
	}
	
	public static double getUnnormalizedFilteredValue(long position, HashMap<Long, Integer> covmap, double[][] filter, boolean isPlusStrand){ // positionwise score
		double[] s = TrainingDataSet.getCoveragesAtOnePoint(position, covmap, isPlusStrand, false);		
		return TrainingDataSet.getFilteredValue(filter, s);
	}
	
	
	private double getLRScore(double score){ // TODO stupid and slow..
		for(int i=0; i<lr.length;i++){
			if(lr[i][0] > score){
				return lr[i][1];
			}
		}
		return lr[lr.length-1][1];
	}
	
	private void write(String contig, IndexedFastaSequenceFile f, HashMap<Long, Double> scores, HashMap<Long, Double> filteredValue, PrintStream out, boolean isPlusStrand){
		//f = new IndexedFastaSequenceFile(new File(fastaFile));
		long offset = isPlusStrand? TrainingDataSet.getLeftWindowSize() : - TrainingDataSet.getLeftWindowSize() + 1;
		long length = f.getSequence(contig).getBases().length;
		for(long position : scores.keySet()){
			long start = position + offset;				
			if(length <= start + 3 || start < 3) continue;				
			
			StringBuffer nac = new StringBuffer();
			for(byte b : f.getSubsequenceAt(contig, start + (isPlusStrand?1:-2), start + (isPlusStrand?3:0)).getBases())
				nac.append((char)b);
			
			String nacString = nac.toString();
				
			out.println(contig + "\t" + start + "\t" + (isPlusStrand? "+" : "-") + "\t" + scores.get(position) + "\t" +  filteredValue.get(position) + "\t" + nacString);
		}	
		out.close();
	}
	
	
	public static void main(String[] args) {
		String keyword = "Thy_Harr10m";
		if(args!=null) keyword = args[0];
		String covFileprefix = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/" + keyword + ".sorted";
		String covFilePlus = covFileprefix + ".plus.cov";
		String covFileMinus = covFileprefix + ".minus.cov";
		String paramFile = covFileprefix + ".param";
		
		TrainingDataSet training = new TrainingDataSet(covFilePlus, covFileMinus, paramFile, "/media/kyowon/Data1/RPF_Project/data/refFlatHuman.txt");
		training.train();
		
		System.out.println("Training done..");
		
		GeneratePositionwiseScore test = new GeneratePositionwiseScore(covFilePlus, covFileMinus, paramFile, "/media/kyowon/Data1/RPF_Project/data/hg19.fa");
		test.run(-3, covFileprefix + ".plus.out", true);
		test.run(-3, covFileprefix + ".minus.out", false);
			
	}
}
