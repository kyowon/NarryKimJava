package fCLIP.paper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.MirGff3FileParser;
import parser.MirGff3FileParser.MiRNA;
import fCLIP.FCLIP_Scorer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountDroshaCleavagesAroundmiRNAs {

	public static void main(String[] args) {
		
		String fCLIPBam = "/media/kyowon/Data1/fCLIP/samples/sample4/alignments/merged.sorted.bam";
		String miRNAgff3 = "/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3";
		MirGff3FileParser mirParser = new MirGff3FileParser(miRNAgff3);
		String outFile = "/media/kyowon/Data1/fCLIP/miRNA_ReadEnd_unpairedIncluded.txt";
		String miRNAsorted = "/media/kyowon/Data1/fCLIP/pri-miRNAs_sorted.txt"; // retreive ... paired or unpaired
		int flankingLength = 10;
		int preLength = 10;
		String fCLIPbed = "/media/kyowon/Data1/fCLIP/samples/sample4/bed/merged.se.bed";
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt");
		String parameterFileName = "/media/kyowon/Data1/fCLIP/samples/sample4/results/training/merged.param";
		String prevContig = null;
		Bed12Parser bedParser = null;
		
		//flankingLength--;
		SamReader reader = SamReaderFactory.makeDefault().open(new File(fCLIPBam));
		try {
			BufferedLineReader oin = new BufferedLineReader(miRNAsorted);
			HashSet<String> miRNAs = new HashSet<String>();
			String s;
			while((s=oin.readLine())!=null){
				String[] token = s.split("\t");
				if(!token[1].equals(".")) continue;
				miRNAs.add(token[0].trim());
			}
			oin.close();
			PrintStream out = new PrintStream(outFile);
			//BufferedLineReader in = new BufferedLineReader(miRNABed);
			
			Iterator<MiRNA> miIterator = mirParser.getMiRNAIterator();
			
			while(miIterator.hasNext()){
				MiRNA mi = miIterator.next();
				if(!miRNAs.contains(mi.getName())) continue;
				boolean isPlusStrand = mi.isPlusStrand();
				String contig = mi.getContig();
				
				Integer f = mi.get5p();
				Integer t = mi.get3p();
				
				out.print(mi.getName() + "\t" + contig + "\t" + (f != null) + ";" + (t != null) + "\t");
				
				
				int[][] fsignal  = new int[2][preLength+flankingLength];
				int[][] tsignal  = new int[2][preLength+flankingLength];
				
				if(f!=null){
					int tmpf = f;
					if(isPlusStrand) tmpf -= 1;
					fsignal = getRead5p3pSignal(reader, contig, tmpf - flankingLength + 1, tmpf + preLength, isPlusStrand);
				}else{ // get five location..
					if(prevContig == null || !prevContig.equals(contig)){
						bedParser = new Bed12Parser(fCLIPbed, annotationParser,	contig, true);
						prevContig = contig;
					}
					FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, null, null, parameterFileName);
					double maxScore = -100;
					int maxf = 0;
					for(int i=40;i<150;i++){
						int tf = t - (isPlusStrand? i:-i);
						double score = scorer.getScore(tf, isPlusStrand, false);
						if(maxScore < score){
							maxScore = score;
							maxf = tf;
						}
					}
					f = maxf;
					if(maxf > 0){
						if(isPlusStrand) maxf -= 1;
						fsignal = getRead5p3pSignal(reader, contig, maxf - flankingLength + 1, maxf + preLength, isPlusStrand);
					}
				}
				
				if(t!=null){
					int tmpt = t;
					if(!isPlusStrand) tmpt -= 1;
					tsignal = getRead5p3pSignal(reader, contig, tmpt - preLength + 1, tmpt + flankingLength, isPlusStrand);
				}else{
					if(prevContig == null || !prevContig.equals(contig)){
						bedParser = new Bed12Parser(fCLIPbed, annotationParser,	contig, true);
						prevContig = contig;
					}
					//bedParser = new Bed12Parser(fCLIPbed, annotationParser,	contig, true);
					FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, null, null, parameterFileName);
					double maxScore = -100;
					int maxt = 0;
					for(int i=40;i<150;i++){
						int tt = f + (isPlusStrand? i:-i);
						double score = scorer.getScore(tt, isPlusStrand, true);
						if(maxScore < score){
							maxScore = score;
							maxt = tt;
						}
					}
					t = maxt;
					if(maxt > 0){
						if(!isPlusStrand) maxt -= 1;
						tsignal = getRead5p3pSignal(reader, contig, maxt - preLength + 1, maxt + flankingLength, isPlusStrand);
					}
				}
				
				out.print(f + "\t" + t + "\t" + (isPlusStrand? "+" : "-") + "\t");
				
				for(int p : fsignal[1]){
					out.print(p+",");
				}
				out.print("\t");
				for(int p : fsignal[0]){
					out.print(p+",");
				}
				out.print("\t");
				for(int p : tsignal[1]){
					out.print(p+",");
				}
				out.print("\t");
				for(int p : tsignal[0]){
					out.print(p+",");
				}
				out.println();
				
			}
			
			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	
	
	private static int[][] getRead5p3pSignal(SamReader reader, String contig, int start, int end, boolean isPlusStrand){ // inclusive 0-based begin end
		int[][] signal = new int[2][end-start+1];
		SAMRecordIterator iterator = reader.query(contig, start+1, end+1, false);
		
		while(iterator.hasNext()){
			SAMRecord read = iterator.next();
			if(isPlusStrand == read.getReadNegativeStrandFlag()) continue;
			if(isPlusStrand){
				int index = read.getAlignmentStart() - start - 1;
				if(index >=0 && index < signal[0].length) signal[0][index] ++;
				
				index = read.getAlignmentEnd() - start - 1;
				if(index >=0 && index < signal[1].length) signal[1][index] ++;
			}else{
				int index = end - read.getAlignmentEnd() + 1;
				if(index >=0 && index < signal[0].length) signal[0][index] ++;
				
				index = end - read.getAlignmentStart() + 1;
				if(index >=0 && index < signal[1].length) signal[1][index] ++;
			}
		}
		iterator.close();
		return signal;
	}
	
}
