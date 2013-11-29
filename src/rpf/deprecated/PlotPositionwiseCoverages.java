package rpf.deprecated;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;

import net.sf.samtools.util.BufferedLineReader;

public class PlotPositionwiseCoverages {

	private HashMap<String, HashMap<Long, Integer>> coveragesNocoPlus;
	private HashMap<String, HashMap<Long, Integer>> coveragesNSPlus;
	private HashMap<String, HashMap<Long, Integer>> coveragesThyPlus;
	private HashMap<String, HashMap<Long, Integer>> coveragesNocoMinus;
	private HashMap<String, HashMap<Long, Integer>> coveragesNSMinus;
	private HashMap<String, HashMap<Long, Integer>> coveragesThyMinus;
	private HashMap<String, Long> positions; // key format - chr:gene_name:+/-
	
	
	public PlotPositionwiseCoverages(String covNocoPlus, String covNSPlus, String covThyPlus, String covNocoMinus, String covNSMinus, String covThyMinus){
		coveragesNocoPlus = new HashMap<String, HashMap<Long, Integer>>();
		coveragesNSPlus = new HashMap<String, HashMap<Long, Integer>>();
		coveragesThyPlus = new HashMap<String, HashMap<Long, Integer>>();
		coveragesNocoMinus = new HashMap<String, HashMap<Long, Integer>>();
		coveragesNSMinus = new HashMap<String, HashMap<Long, Integer>>();
		coveragesThyMinus = new HashMap<String, HashMap<Long, Integer>>();
		positions = new HashMap<String, Long>();
		
		readCovs(covNocoPlus, covNSPlus, covThyPlus, covNocoMinus, covNSMinus, covThyMinus);
		
	}
	
	private void readCovs(String covNocoPlus, String covNSPlus, String covThyPlus, String covNocoMinus, String covNSMinus, String covThyMinus){
		TrainingDataSet.getAllCoverages(coveragesNocoPlus, covNocoPlus);
		TrainingDataSet.getAllCoverages(coveragesNSPlus, covNSPlus);
		TrainingDataSet.getAllCoverages(coveragesThyPlus, covThyPlus);
		TrainingDataSet.getAllCoverages(coveragesNocoMinus, covNocoMinus);
		TrainingDataSet.getAllCoverages(coveragesNSMinus, covNSMinus);
		TrainingDataSet.getAllCoverages(coveragesThyMinus, covThyMinus);
	}
	
	private void readPositions(String pos){
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(new FileInputStream(pos));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				StringBuffer key = new StringBuffer();
				key.append(token[2]); key.append(':');
				key.append(token[1]); key.append(':');
				key.append(token[10]); 
				//System.out.println(key);
				long position = Long.parseLong(token[9]);
				positions.put(key.toString(), position);
			}
			
			in.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	private double[] getCoveragesAtOnePosition(HashMap<String, HashMap<Long, Integer>> coverages, String contig, boolean isPlusStrand, long position){
		HashMap<Long, Integer> covmap = coverages.get(contig);
		if(covmap == null) return null;
		return TrainingDataSet.getCoveragesAtOnePoint(position + (isPlusStrand? -30 : + 30), covmap, isPlusStrand, false);
	}
	
	private void write(String outFile, String key){
		try {
			PrintStream out = new PrintStream(outFile);
			
			for(String s : positions.keySet()){
				long position = positions.get(s);
				String[] token = s.split(":");
				String contig = token[0];
				String gene = token[1];
				boolean isPlusStrand = token[2].equals("+");
				double[] noco = getCoveragesAtOnePosition((isPlusStrand? coveragesNocoPlus : coveragesNocoMinus), contig, isPlusStrand, position);
				//double[] ns = getCoveragesAtOnePosition((isPlusStrand? coveragesNSPlus : coveragesNSMinus), contig, isPlusStrand, position);
				double[] thy = getCoveragesAtOnePosition((isPlusStrand? coveragesThyPlus : coveragesThyMinus), contig, isPlusStrand, position);
				
				if(noco == null) continue;
				//if(gene.equals("NM_000019")) System.out.println(contig + " " + position);
				out.println(gene + "_" +  key + "=[");
				for(int i=0;i<Math.min(120, noco.length);i++){
					out.println(noco[i]  + "\t" + thy[i]); 
				}
				out.println("];");
			}			
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		String key = "2_RPF";
		String covNocoPlus = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_"+key+".sorted.plus.cov";
		//String covNSPlus = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NS_"+key+".sorted.plus.cov";
		String covThyPlus = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_"+key+".sorted.plus.cov";
		String covNocoMinus = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_"+key+".sorted.minus.cov";
		//String covNSMinus = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NS_"+key+".sorted.minus.cov";
		String covThyMinus = "/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_"+key+".sorted.minus.cov";
		
		PlotPositionwiseCoverages test = new PlotPositionwiseCoverages(covNocoPlus, null, covThyPlus, covNocoMinus, null, covThyMinus);
		test.readPositions("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NocovsThy.txt");
		//test.readPositions("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NSvsNoco.txt");
		//test.readPositions("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/NSvsThy.txt");
		
		test.write("/media/kyowon/Data1/RPF_Project/results/point"+key+".m", key);
		
	}

}
