package fCLIP;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import fCLIP.MirGff3FileParser.MiRNA;
import parser.AnnotationFileParser;
import parser.BedCovFileParser;
import parser.ScoringOutputParser;
import parser.ZeroBasedFastaParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;

public class fCLIPScorer {
	
	static private int leftWindowSize = 40;
	static private int rightWindowSize = 40;
	
	static private HashSet<String> scoreNWrite(boolean isPlusStrand, boolean is5prime, BedCovFileParser bedCovFileParser, ZeroBasedFastaParser fastaParser, MirGff3FileParser mirParser, PrintWriter out){
		HashSet<String> miRNAs = new HashSet<String>();
		for(String contig : bedCovFileParser.getContigs()){
			if(!fastaParser.containsContig(contig)){
				continue;
			}
			System.out.println("Scoring for " + contig + " " + (isPlusStrand? "+" : "-") + " strand");
			Iterator<Integer> iterator;
			iterator = bedCovFileParser.getNonZeroCoveragePositionIterator(contig);
			//int lastConsidered = 0;
			int positionOffset =  (is5prime && isPlusStrand || !is5prime && !isPlusStrand? - 30: 30);
			while(iterator.hasNext()){
				int position = iterator.next();
				double maxCov = 0;
				double[] covs = bedCovFileParser.getCoverages(contig, position, leftWindowSize, rightWindowSize, isPlusStrand);
				if(covs[leftWindowSize] < 5) continue;
				
				for(double cov : covs){
					maxCov = maxCov > cov? maxCov : cov;
				}
				if(maxCov > covs[leftWindowSize]) continue;
				
				position += positionOffset;
				
			//	int	start = isPlusStrand? position - leftWindowSize - 1: position - rightWindowSize - 1;
			//	int	end = isPlusStrand? position + rightWindowSize + 1: position + leftWindowSize + 1;				
			
				//for(int currentPosition=Math.max(lastConsidered + 1, start);currentPosition<end;currentPosition++){
					
				String seq = null;
				if(isPlusStrand){
					seq = fastaParser.getSequence(contig, position - leftWindowSize, position + rightWindowSize);
				}else{							
					seq = ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(contig, position- rightWindowSize + 1, position+ leftWindowSize + 1), true);
				}	
				
				ArrayList<Double> ret = RunRNAfold.run(seq);
				double energy = ret.get(0);
				double depth = ret.get(1);
				MiRNA miRNA = mirParser.getContainingMiRNA(contig, isPlusStrand, position);
				
				out.println(energy + "\t" + depth + "\t" + (miRNA == null ? 0 : 1));		
				if(miRNA !=null){
					//System.out.println(miRNA  + "\t" + position);
					miRNAs.add(miRNA.getName());
				}
				
				//}
			}
		}
		return miRNAs;
	}
	
	static private void scoreNWrite(MirGff3FileParser mirParser, ZeroBasedFastaParser fastaParser, PrintWriter out){
		Iterator<MiRNA> iterator = mirParser.getMiRNAIterator();
		while(iterator.hasNext()){
			MiRNA mi = iterator.next();
			String contig = mi.getContig();
			
			if(!fastaParser.containsContig(contig)){
				continue;
			}
			//System.out.println(contig);
			//System.out.println("Scoring for " + contig + " " + (isPlusStrand? "+" : "-") + " strand");
			boolean isPlusStrand = mi.isPlusStrand();
			int	start = mi.getStart();
			int	end = mi.getEnd();
			int currentPosition = (start + end)/2;
			
			String seq = null;
			if(isPlusStrand){
				seq = fastaParser.getSequence(contig, currentPosition - leftWindowSize, currentPosition + rightWindowSize);
			}else{							
				seq = ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(contig, currentPosition- rightWindowSize + 1, currentPosition+ leftWindowSize + 1), true);
			}	
			
			ArrayList<Double> ret = RunRNAfold.run(seq);
			double energy = ret.get(0);
			double depth = ret.get(1);
			out.println(energy + "\t" + depth);
			//if(depth < 25) System.out.println(mi.getName() + " " + seq);
		//	lastConsidered = currentPosition;			
		
		
		}
	}

	static public void main(String[] args) throws IOException{
		// -20 20
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		PrintWriter out = new PrintWriter("/media/kyowon/Data1/fCLIP/scores.m");
		MirGff3FileParser mirParser = new MirGff3FileParser("/media/kyowon/Data1/fCLIP/Genome/hsa.gff3");
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/Genome/hg19.refFlat.txt");
		BedCovFileParser bedCovFilePlus5pParser = new BedCovFileParser("/media/kyowon/Data1/fCLIP/Data/Drosha2_paired_Plus_5.bedgraph", annotationParser);
		BedCovFileParser bedCovFilePlus3pParser = new BedCovFileParser("/media/kyowon/Data1/fCLIP/Data/Drosha2_paired_Plus_3.bedgraph", annotationParser);
		BedCovFileParser bedCovFileMinus5pParser = new BedCovFileParser("/media/kyowon/Data1/fCLIP/Data/Drosha2_paired_Minus_5.bedgraph", annotationParser);
		BedCovFileParser bedCovFileMinus3pParser = new BedCovFileParser("/media/kyowon/Data1/fCLIP/Data/Drosha2_paired_Minus_3.bedgraph", annotationParser);
		
		out.println("s1=[");
		fCLIPScorer.scoreNWrite(mirParser, fastaParser, out);
		out.println("];");
		out.close();
		
		HashSet<String> miRNAs = new HashSet<String>();
		
		PrintWriter out2 = new PrintWriter("/media/kyowon/Data1/fCLIP/scores2.m");
		out2.println("s2=[");
		miRNAs.addAll(fCLIPScorer.scoreNWrite(true, true,  bedCovFilePlus5pParser, fastaParser, mirParser, out2));
		miRNAs.addAll(fCLIPScorer.scoreNWrite(true, false,  bedCovFilePlus3pParser, fastaParser, mirParser, out2));
		miRNAs.addAll(fCLIPScorer.scoreNWrite(false, true,  bedCovFileMinus5pParser, fastaParser, mirParser, out2));
		miRNAs.addAll(fCLIPScorer.scoreNWrite(false, false,  bedCovFileMinus3pParser, fastaParser, mirParser, out2));
		
		out2.println("];");
		out2.close();
		System.out.println("#Micro RNAs : " + miRNAs.size());
	}
	
}
