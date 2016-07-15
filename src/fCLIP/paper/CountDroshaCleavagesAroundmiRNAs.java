package fCLIP.paper;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.MirGff3FileParser;
import parser.MirGff3FileParser.MiRNA;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class CountDroshaCleavagesAroundmiRNAs {

	public static void run(String fCLIPbed, String miRNAgff3, String miRNAsorted, String outFile, String annotationFile, int flankingLength, int preLength, String mirGeneDB){
		/*String fCLIPbed = args[0];//"/media/kyowon/Data1/fCLIP/samples/sample13_HelaH/bed/Hela.se.bed";//"/media/kyowon/Data1/fCLIP/samples/sample8/bed/bamMerged.se.bed";
		String miRNAgff3 = args[1];//"/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3";
		String miRNAsorted = args[2];//"/media/kyowon/Data1/fCLIP/pri-miRNAs_sorted.txt"; // retreive ... paired or unpaired
		String outFile = args[3];// "/media/kyowon/Data1/fCLIP/NoFCLIPCalculatedRNA_ReadEnd_unpairedIncludedHela.txt";
		//"/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt");
		
		int flankingLength = Integer.parseInt(args[5]);//10;
		int preLength = Integer.parseInt(args[6]);//10;
		*/
		AnnotationFileParser annotationParser = new AnnotationFileParser(annotationFile);
		MirGff3FileParser mirParser = new MirGff3FileParser(miRNAgff3);
		if(mirGeneDB!=null) mirParser.correctAnnotationUsignMirGeneDB(mirGeneDB);
		//"/media/kyowon/Data1/fCLIP/genomes/mirgenedb.txt"
		
		
		
		
		String prevContig = null;
		Bed12Parser bedParser = null;
		boolean calculateForUnannotated = true;
		
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
				if(f == null && t == null) continue;
				
				out.print(mi.getName() + "\t" + contig + "\t" + (f != null) + ";" + (t != null) + "\t");
				
				
				int[][] fsignal  = new int[2][preLength+flankingLength];
				int[][] tsignal  = new int[2][preLength+flankingLength];
				if(prevContig == null || !prevContig.equals(contig)){
					bedParser = new Bed12Parser(fCLIPbed, annotationParser,	contig, true);
					prevContig = contig;
				}
				if(f!=null){
					int tmpf = f;
					if(isPlusStrand) tmpf -= 1;
					fsignal = getRead5p3pSignal(bedParser, tmpf - flankingLength + 1, tmpf + preLength, isPlusStrand);//getRead5p3pSignal(reader, contig, tmpf - flankingLength + 1, tmpf + preLength, isPlusStrand);					
				}else if(calculateForUnannotated){ // get five location..					
					int maxf = -1;
					int maxDepth = -1;
					for(int i=40;i<=160;i++){
						int tmpf = t - (isPlusStrand? i : -i);						
						int depth = (bedParser.get5pDepth(isPlusStrand, tmpf) +1 ) * (bedParser.get3pDepth(isPlusStrand, tmpf - (isPlusStrand? 1 : -1))+1);
						if(maxDepth < depth){
							maxDepth = depth;
							maxf = tmpf;
						}
					}
					f = maxf;
					if(isPlusStrand) maxf -= 1;
					fsignal = getRead5p3pSignal(bedParser, maxf - flankingLength + 1, maxf + preLength, isPlusStrand);//getRead5p3pSignal(reader, contig, tmpf - flankingLength + 1, tmpf + preLength, isPlusStrand);
					
				}
				
				if(t!=null){
					int tmpt = t;
					if(!isPlusStrand) tmpt -= 1;
					tsignal = getRead5p3pSignal(bedParser, tmpt - preLength + 1, tmpt + flankingLength, isPlusStrand);
				}else if(calculateForUnannotated){	
					
					int maxt = -1;
					int maxDepth = -1;
					for(int i=40;i<=160;i++){
						int tmpt = f + (isPlusStrand? i : -i);						
						int depth = (bedParser.get3pDepth(isPlusStrand, tmpt)+1 )* (bedParser.get5pDepth(isPlusStrand, tmpt + (isPlusStrand? 1 : -1))+1);
						if(maxDepth < depth){
							maxDepth = depth;
							maxt = tmpt;
						}
					}
					t = maxt;
					if(!isPlusStrand) maxt -= 1;
					tsignal = getRead5p3pSignal(bedParser, maxt - preLength + 1, maxt + flankingLength, isPlusStrand);				
					
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

	
	
	public static int[][] getRead5p3pSignal(SamReader reader, String contig, int start, int end, boolean isPlusStrand){ // inclusive 0-based begin end
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
	
	
	public static int[][] getRead5p3pSignal(Bed12Parser parser, int start, int end, boolean isPlusStrand){ // inclusive 0-based begin end
		int[][] signal = new int[2][end-start+1];
		
		for(int l=start;l<=end;l++){
			signal[0][isPlusStrand? l-start : end - l] = parser.get5pDepth(isPlusStrand, l);
			signal[1][isPlusStrand? l-start : end - l] = parser.get3pDepth(isPlusStrand, l);
		}	
		
		return signal;
	}
	
	
}
