package fCLIP.analysis;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;

public class OverlapWithBam {
	public static void main(String[] args){		
		String out = args[0];
		String in = args[1];
		String bam = args[2];
		if(args[args.length-1].equals("CIS")) ScoredPositionOutputParser.intersectWithBam(out, in, bam, true);
		if(args[args.length-1].equals("TRANS")) ScoredPairOutputParser.intersectWithBam(out, in, bam, true);
	}
}
