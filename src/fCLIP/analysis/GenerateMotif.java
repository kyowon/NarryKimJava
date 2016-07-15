package fCLIP.analysis;

import parser.ZeroBasedFastaParser;
import fCLIP.parser.ScoredPositionOutputParser;

public class GenerateMotif {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String csv = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.csv";
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser("/media/kyowon/Data1/fCLIP/genomes/hg38.fa");
		int flankingNTNumber = 15;
		String fasta5p = "/media/kyowon/Data1/fCLIP/motifFigure/ghg5p.fasta";
		String fasta3p = "/media/kyowon/Data1/fCLIP/motifFigure/ghg3p.fasta";
		String mfile = "/media/kyowon/Data1/fCLIP/motifFigure/ghg.m";
		String classification = "miRNA";
		int option = 3;
		int aluOption = 1;
		
		ScoredPositionOutputParser.generateFastaForMotif(csv,  fastaParser, flankingNTNumber, fasta5p, fasta3p, 
				mfile, classification, option, aluOption);
	//	if(args.length == 4) ScoredPositionOutputParser.generateFastaForMotif(csv, out5p, out3p, out5p + ".m" , args[3]);
	//	else if(args.length == 3) ScoredPairOutputParser.generateFastaForMotif(csv, out5p, out3p);
		
	}

}
