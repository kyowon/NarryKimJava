package rpf.analysis;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import parser.AnnotationFileParser;
import rpf.parser.MergedFileParser;
import rpf.parser.MergedFileParser.MergedResult;
import rpf.parser.ScoringOutputParser.ScoredPosition;
import util.Codon;

public class SmORFselector {

	public static void main(String[] args) throws FileNotFoundException {
		String sample = "sample3"; // sample4 e, sample3, 5 h h2
		String fastaPrefix = "h";
		MergedFileParser test = new MergedFileParser("/media/kyowon/Data1/RPF_Project/samples/" + sample + "/results/out_ssh_0.3.csv");
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
		HashSet<String> seqs = new HashSet<String>();
		HashMap<String, Integer> regionCountMap = new HashMap<String, Integer>();
		
		for(MergedResult result: test.getList()){
			String region = "";
			
			if(start){
				out.println(result.getSimpleHeader());
				start = false;
			}
			if(result.getGene()!=null && !result.getGenomicRegion().startsWith("NR")){
				if(result.getAnnotatedClass().equals("TS") || result.getAnnotatedClass().equals("T")){
					if(result.getFrameShift().equals("0") && result.getGene().isPlusStrand() == result.isPlusStrand()){
						ArrayList<String> gr = anntationParser.getGenomicRegionNameAndFrameShift(result.getContig(), result.isPlusStrand(), result.getStopPosition()+ (result.isPlusStrand()?-2:2), result.getGene(), true);
						
						//if(result.getGene().getGeneName().equals("Mterf1a")){
						//	System.out.println(gr);
						//}
						
						if(gr.get(0).equals("NM_ORF") && gr.get(1).equals("0")){
							region = "ORF";
							continue;	
						}
						else region = "ExomSpanning/NovelIsoform";
					}else if(!result.getFrameShift().equals("0")){
						region = "FrameShifted";
					}else{
						region = "StrandReversed";
					}
				}else{
					if(result.getGenomicRegion().contains("5_UTR")) region = "5UTR";
					else if(result.getGenomicRegion().contains("3_UTR")) region = "3UTR";
					else if(result.getGenomicRegion().contains("Intron")){
						System.out.println(result.getGenomicRegion());
						region = "Intron";
					}
					else{
						//System.out.println(result);
						region = "ExomSpanning/NovelIsoform";
					}
				}
			}else if(result.getGene() == null){
				region = "Intergenic";
			}else if(result.getGenomicRegion().startsWith("NR")){
				region = "NR";
			}else{
				
			}
			
			if(!regionCountMap.containsKey(region)) regionCountMap.put(region, 0);
			regionCountMap.put(region, regionCountMap.get(region)+1);
			
			//
			String[] predicted = result.getPredictedClasses();
			
			//System.out.println(result.getStopPosition() + " " + result.getLength() );
			//if(result.getStopPosition() < 20 || (result.getStopPosition() > 450 && result.getLength() > 450)) continue;
			if(result.getLength() > 600) continue;
			boolean keep = true;
			//if(pp!=null && result.getGene().equals(pp.getGene())) 
			//	System.out.println(pp.getPosition() + "\n* " + result.getPosition() +"\n " + pp.equals(result.getScoredPosition()));
			
			if(pp.contains(result.getScoredPosition())){
				continue;
			}
			
			//boolean containingTS = false;
			int n = 0;
			for(int i=0; i< predicted.length;i++){
				String p = predicted[i];
				if(!p.toUpperCase().equals("NC")){
					n++;
				}
			}
			
			if(n<1) keep = false;			
			
			//if(!keep) continue;
			
			String frame = anntationParser.getGenomicRegionNameAndFrameShift(result.getContig(), result.isPlusStrand(), result.getPosition(), result.getGene(), true).get(1);
			result.setFrameShift(frame);
			pp.add(result.getScoredPosition());
					
			String seq = result.getSequence();
			
			double gc = 0;
			for(int k=0;k<seq.toUpperCase().length();k++){
				if(seq.charAt(k) == 'G' || seq.charAt(k) == 'C') gc ++; 
			}
			
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
			String translated = Codon.getAminoAcids(seq);
			if(seqs.contains(translated)) continue;
			if(translated.length() < 7) continue;
			seqs.add(translated);
			num++;
			out.print(result.toSimpleString());
			
			out.println("\t" + (gc/seq.length()) + "\t"+ translated);
			outfasta.println(">"+fastaPrefix+"_" + j++);
			outfasta.println(translated);
			
		}
		for(String region : regionCountMap.keySet()){
			System.out.println(region+ " " + regionCountMap.get(region));
		}
		
		System.out.println(num);
		outfasta.close();
		out.close();
	}

}
