package fCLIP.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import parser.AnnotationFileParser;
import parser.Bed12Parser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class GenerateDepthsForEncodeDataSets {

	public static void generate(String infile, String outfile, AnnotationFileParser annotationParser ) {
		//String infile = args[0];//"/media/kyowon/Data1/Dropbox/h19x2.sorted.out.csv";
		//String outfile = args[1];//"/media/kyowon/Data1/Dropbox/h19x2.sorted.out.encode.csv";
		String bedfolder = "/media/kyowon/Data1/fCLIP/genomes/ENCODE";
	//	AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg19.refFlat.txt");
		
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(infile);
		HashMap<ScoredPosition, int[]> depthMap3p = new HashMap<ScoredPosition, int[]>();
		HashMap<ScoredPosition, int[]> depthMap5p = new HashMap<ScoredPosition, int[]>();
		HashMap<ScoredPosition, int[]> depthMapPre = new HashMap<ScoredPosition, int[]>();
		
		HashSet<String> contigs = new HashSet<String>();
		ArrayList<ScoredPosition> positions = new ArrayList<ScoredPosition>();
		int numfile = 0;
		
		for(File bedfile : new File(bedfolder).listFiles()){
			if(!bedfile.getName().endsWith(".bed")) continue;
			numfile++;
		}
		
		for(ScoredPosition position : parser.getPositions()){
			if(!position.getClassification().toUpperCase().equals("M")) continue;				
			depthMap3p.put(position, new int[numfile]);
			depthMap5p.put(position, new int[numfile]);
			depthMapPre.put(position, new int[numfile]);
			contigs.add(position.getContig());
			positions.add(position);
		}
		
		System.out.println(positions.size());
		
		int i = 0;
		String bedfilenames = "";
		for(File bedfile : new File(bedfolder).listFiles()){
			if(!bedfile.getName().endsWith(".bed")) continue;
			int j = 0;
			System.out.println(bedfile.getName());
			for(String contig : contigs){
				System.out.println(contig);
				Bed12Parser bedparser = new Bed12Parser(bedfile.getAbsolutePath(), annotationParser, contig, true);
				for(ScoredPosition position : positions){
					if(!position.getContig().equals(contig)) continue;
					int sign = position.isPlusStrand()? 1 : -1;
					int d3 = bedparser.get5pDepth(position.isPlusStrand(), position.getThreePposition()); // TODO test
					d3 += bedparser.get5pDepth(position.isPlusStrand(), position.getThreePposition() + sign); // TODO test
					d3 += bedparser.get5pDepth(position.isPlusStrand(), position.getThreePposition() + sign * 2); // TODO test
					
					int d5 = bedparser.get3pDepth(position.isPlusStrand(), position.getFivePposition()); // TODO test
					d5 += bedparser.get3pDepth(position.isPlusStrand(), position.getFivePposition() - sign); // TODO test
					d5 += bedparser.get3pDepth(position.isPlusStrand(), position.getFivePposition() - sign * 2); // TODO test
					
					int dp = bedparser.getReadDepth(position.isPlusStrand(), position.getThreePposition(), position.getFivePposition());
					dp += bedparser.getReadDepth(position.isPlusStrand(), position.getThreePposition() + sign, position.getFivePposition() - sign);
					dp += bedparser.getReadDepth(position.isPlusStrand(), position.getThreePposition() + sign * 2, position.getFivePposition() - sign);
					dp += bedparser.getReadDepth(position.isPlusStrand(), position.getThreePposition() + sign, position.getFivePposition() - sign * 2);
					dp += bedparser.getReadDepth(position.isPlusStrand(), position.getThreePposition() + sign * 2, position.getFivePposition() - sign * 2);
					
					depthMap3p.get(position)[i] += d3;
					depthMap5p.get(position)[i] += d5;
					depthMapPre.get(position)[i] += dp;
					j++;
					System.out.println(bedfile.getName() + " " + j + " " + positions.size());
				}
			}
			bedfilenames += "\t" + bedfile.getName() + "5p\t" + bedfile.getName() + "3p\t" + bedfile.getName() + "pre";
			i++;
		}
			
		bedfilenames += "\t5pSum\t3pSum\tpreSum";
		
		try {
			PrintStream outStream = new PrintStream(outfile);
			outStream.print(parser.getHeader());
			outStream.println(bedfilenames);
			
			for(ScoredPosition position : positions){
				outStream.print(position);
				int sum5 = 0;
				int sum3 = 0;
				int sumPre = 0;
				int[] d3s = depthMap3p.get(position);
				int[] d5s = depthMap5p.get(position);
				int[] dps = depthMapPre.get(position);
				
				for(int n=0;n<numfile;n++){
					outStream.print("\t" + d3s[n] + "\t" + d5s[n] + "\t" + dps[n]);
					sum5 += d3s[n];
					sum3 += d5s[n];
					sumPre += dps[n];
				}
				outStream.println("\t" + sum5 + "\t" + sum3 + "\t" + sumPre);				
			}
			
			
			outStream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
