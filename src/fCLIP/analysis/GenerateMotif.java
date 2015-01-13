package fCLIP.analysis;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;

public class GenerateMotif {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String out5p = args[0];
		String out3p = args[1];
		String csv = args[2];
		
		if(args.length == 4) ScoredPositionOutputParser.generateFastaForMotif(csv, out5p, out3p, args[3]);
		else if(args.length == 3) ScoredPairOutputParser.generateFastaForMotif(csv, out5p, out3p);
		
	}

}
