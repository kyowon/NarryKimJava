package rpf.analysis;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import parser.AnnotationFileParser;
import rpf.parser.MergedFileParser;
import rpf.parser.MergedFileParser.MergedResult;
import rpf.parser.ScoringOutputParser.ScoredPosition;
import util.Codon;

public class SmORFselector {

	public static void main(String[] args) throws FileNotFoundException {
		String sample = "sample5";
		MergedFileParser test = new MergedFileParser("/media/kyowon/Data1/RPF_Project/samples/" + sample + "/results/out_new_0.3.csv");
		PrintStream out = new PrintStream("/media/kyowon/Data1/RPF_Project/samples/" + sample + "/results/out_0.3_smORF.csv");
		PrintStream outfasta = new PrintStream("/media/kyowon/Data1/RPF_Project/samples/" + sample + "/results/out_0.3_smORF.fasta");
		AnnotationFileParser anntationParser = new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.refFlat.txt");
		//ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.fa");
	//	ArrayList<MergedResult> list = new ArrayList<MergedFileParser.MergedResult>();
		HashSet<ScoredPosition> pp = new HashSet<ScoredPosition>();
		boolean start = true;
		int j=2;
		int num = 0;
		HashMap<String, Double> codonUsage = new HashMap<String, Double>();
		
		for(MergedResult result: test.getList()){
			if(start){
				out.println(result.getSimpleHeader());
				start = false;
			}
					
			//if(result.getGene()==null || !result.getGene().getGeneName().contains("LINC")) continue;
			//String region = result.getGenomicRegion();
			//if(region.equals("InterGenic")) System.out.println(result.getStopPosition() + " " + result.getLength());
			
			if(result.getAnnotatedClass().equals("TS") || result.getAnnotatedClass().equals("T") || result.getAnnotatedClass().equals("_")) continue; 
			//if(region.contains("ORF") && !region.contains("Intron")) continue;
			String[] predicted = result.getPredictedClasses();
			double[] predictionScores = result.getPredictedClassesScores();
			
			if(result.getStopPosition() < 20 || result.getStopPosition() > 450 && result.getLength() > 450) continue;
			boolean keep = true;
			//if(pp!=null && result.getGene().equals(pp.getGene())) 
			//	System.out.println(pp.getPosition() + "\n* " + result.getPosition() +"\n " + pp.equals(result.getScoredPosition()));
			if(pp.contains(result.getScoredPosition())){
				continue;
			}
			//boolean containingTS = false;
			int n = 0;
			int m = 0;
			for(int i=0; i< predicted.length;i++){
				if(result.getRpfStartScores()[i]<.25) continue;
				//if(n++ > 3) break;
				String p = predicted[i];
				//if(p.equals("TS")) containingTS = true;
				if(p.toUpperCase().equals("TS") || p.toUpperCase().equals("T")){
					n++;
				//	continue;
				}
				//if(m++>3) break;
			//	break;
				//keep = false;
				//break;
			}
			
			if(n<1) keep = false;
			//if(((float)n) / predicted.length < .3) keep = false;
			//if(n != m) keep = false;
			if(result.getScoredPosition().isPlusStrand()){
			//	if(result.getGene() != null && result.getStopPosition() > result.getGene().getCdsStart() ) 
					//keep = false;
			}else{
			//	if(result.getGene() != null && result.getStopPosition() <= result.getGene().getCdsEnd()) 
					//keep = false;
			}
				
			if(!keep) continue;
			
			String frame = anntationParser.getGenomicRegionNameAndFrameShift(result.getContig(), result.isPlusStrand(), result.getPosition(), result.getGene(), true).get(1);
			result.setFrameShift(frame);
			pp.add(result.getScoredPosition());
			
		//	list.add(result);
			//if(r.getGene() != null){
				
				//	if(Math.abs(r.getPosition() - r.getGene().getCdsStart())<30) continue;
				//	if(Math.abs(r.getPosition() - r.getGene().getCdsEnd())<30) continue;
				//}
			out.print(result.toSimpleString());
					
			String seq = result.getSequence();
			
			double gc = 0;
			for(int k=0;k<seq.toUpperCase().length();k++){
				if(seq.charAt(k) == 'G' || seq.charAt(k) == 'C') gc ++; 
			}
			if(Codon.getAminoAcids(seq).startsWith("L")) continue; // TODO 
			//System.out.println(Codon.getAminoAcids(seq));
			
			for(int k=0;k<seq.length()-2;k+=3){
				String codon = seq.toUpperCase().substring(k, k+3);
				if(!codonUsage.containsKey(codon)) codonUsage.put(codon, .0);
				codonUsage.put(codon, codonUsage.get(codon)+1.0);
			}
			
			double sum = 0;
			for(double v : codonUsage.values()) sum += v;
			for(String codon : codonUsage.keySet()){
				codonUsage.put(codon, codonUsage.get(codon) / sum);
			}
			num++;
			out.println("\t" + (gc/seq.length()) + "\t"+ Codon.getAminoAcids(seq));
			outfasta.println(">e_" + j++);
			outfasta.println(Codon.getAminoAcids(seq));
			
			if(Codon.getAminoAcids(seq).contains("SSPVFQVPKDDTTELDSLGLASPPK")) System.out.println(seq);
		//	if(Codon.getAminoAcids(seq).contains("WDYPE")) System.out.println("seq");
		}
		//for(String codon : codonUsage.keySet()){
		//	System.out.println(codon + "\t" + codonUsage.get(codon));
		//}
		System.out.println(num);
		outfasta.close();
		out.close();
	}

}
