package fCLIP.analysis;

import parser.AnnotationFileParser;
import parser.ZeroBasedFastaParser;
import fCLIP.parser.ScoredPositionOutputParser;

public class GenerateCoverageBedfile {

	public static void main(String[] args) {
		if(args[args.length-1].equals("KD"))
			ScoredPositionOutputParser.generateBedFromCsv(args[0], args[1], args[2], false);
		if(args[args.length-1].equals("BG")){
			ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser(args[3]);
			AnnotationFileParser annotationParser = new AnnotationFileParser(args[2]);
			ScoredPositionOutputParser.generateRandomizedBedFromCsv(args[0], args[1], annotationParser, fastaParser);
		}
	}

}
