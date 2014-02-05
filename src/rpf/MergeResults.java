package rpf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import parser.BedCovFileParser;
import parser.RPFMergedFileParser;
import parser.RPFMergedFileParser.MergedResult;
import parser.ScoringOutputParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;
import parser.ZeroBasedFastaParser;

public class MergeResults {
	private ScoringOutputParser[] scoringOutputParsers;
	private Scorer[] scorers = null; 
	private Quantifier[] rpfQuantifiers = null;
	private Quantifier[] rnaQuantifiers = null;
	
	public MergeResults(String[] scoreOutFiles, 
			String[] harrCovPlusFiles, String[] harrCovMinusFiles, 
			String[] rpfCovPlusFiles, String[] rpfCovMinusFiles,
			String[] rnaCovPlusFiles, String[] rnaCovMinusFiles,
			String[] paramFiles,
			String annotationFile,
			String fastaFile){
		scoringOutputParsers = new ScoringOutputParser[scoreOutFiles.length];
		scorers = new Scorer[scoreOutFiles.length];
		if(rpfCovPlusFiles!=null)
			rpfQuantifiers = new Quantifier[rpfCovPlusFiles.length];
		if(rnaCovPlusFiles!=null)
			rnaQuantifiers = new Quantifier[rnaCovPlusFiles.length];
		
		for(int i=0;i<scoringOutputParsers.length;i++){
			scoringOutputParsers[i] = new ScoringOutputParser(scoreOutFiles[i]);
			scorers[i] = new Scorer(harrCovPlusFiles[i], harrCovMinusFiles[i], paramFiles[i], annotationFile);
		}
		if(rpfQuantifiers!=null){
			for(int i=0;i<rpfQuantifiers.length;i++){
				rpfQuantifiers[i] = new Quantifier(rpfCovPlusFiles[i], rpfCovMinusFiles[i], annotationFile, fastaFile);
			}
		}
		if(rnaQuantifiers!=null){
			for(int i=0;i<rnaQuantifiers.length;i++){
				rnaQuantifiers[i] = new Quantifier(rnaCovPlusFiles[i], rnaCovMinusFiles[i], annotationFile, fastaFile);
			}
		}
	}
	
	void merge(String outFile, double scoreThreshold){
		try {
			PrintStream out = new PrintStream(outFile);
			PrintStream mout = new PrintStream(outFile.replace('.', '_')+".m");
			ArrayList<MergedResult> mergedResults = new ArrayList<MergedResult>();
			HashSet<String> contigs = new HashSet<String>();
			for(int i=0;i<scoringOutputParsers.length;i++){
				contigs.addAll(scoringOutputParsers[i].getContigs());
			}
			for(String contig : contigs){
				for(ScoredPosition position : ScoringOutputParser.getUnionPositions(scoringOutputParsers, contig, scoreThreshold)){
					mergedResults.add(new RPFMergedFileParser().new MergedResult(position, scorers, rpfQuantifiers, rnaQuantifiers));	
					writeMfile(mout, position);	
				}
			}		
			out.println(mergedResults.get(0).getHeader());
			for(MergedResult m : mergedResults){
				out.println(m);
			}
			out.close();
			mout.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}	
	
	private void writeMfile(PrintStream mout, ScoredPosition position){
		String gname;
		if(!position.isAnnotated()) gname = position.getContig() +"_"+position.getPosition()+"_"+(position.isPlusStrand()? "P" : "M");
		else gname = position.getGene().getGeneName();
		mout.print(gname);
		//a{3}=[ 4 5 ; 7 8 ]'
		for(int i=0;i<scorers.length;i++){
			mout.println("{"+(i+1)+"}=[");
			BedCovFileParser sbp =  position.isPlusStrand()? scorers[i].getBedCovPlusFileParser() : scorers[i].getBedCovMinusFileParser();
			subWriteMfile(mout, position, sbp);
			if(rpfQuantifiers != null && rpfQuantifiers[i] != null){
				BedCovFileParser rpfbp =  position.isPlusStrand()? rpfQuantifiers[i].getBedCovPlusFileParser() : rpfQuantifiers[i].getBedCovMinusFileParser();
				subWriteMfile(mout, position, rpfbp);
			}
			if(rnaQuantifiers != null && rnaQuantifiers[i] != null){
				BedCovFileParser rnabp =  position.isPlusStrand()? rnaQuantifiers[i].getBedCovPlusFileParser() : rnaQuantifiers[i].getBedCovMinusFileParser();
				subWriteMfile(mout, position, rnabp);
			}
			mout.println("]';");	
		}				
					
	}
	
	private void subWriteMfile(PrintStream mout, ScoredPosition position, BedCovFileParser bp){
		double[] cov = bp.getCoverages(position.getContig(), position.getPosition(), 30, 100, position.isPlusStrand());
		for(Double c : cov){
			mout.append(c.toString());
			mout.append('\t');
		}
		mout.append('\n');
	}
		
	public static void main(String[] args){ // quantities can be from different files!!!!  param and bedCovPlusFiles should be added for quan.
		String[] scoreOutFiles = {	"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10m.sorted.out",
									"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.out"};
		String[][] bedCovPlusFiles = {
				{"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10m.sorted.plus.cov",
		"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.plus.cov"},
		{"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_2_RPF.sorted.plus.cov",
			"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_2_RPF.sorted.plus.cov"}
		};
		String[][] bedCovMinusFiles = {
				{"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10m.sorted.minus.cov",
		"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.minus.cov"},
		{"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_2_RPF.sorted.minus.cov",
		"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_2_RPF.sorted.minus.cov"},
		};
		String[][] paramFiles = {
				{"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10m.sorted.param",
		"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_Harr10m.sorted.param"},
		{"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_2_RPF.sorted.param",
			"/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Thy_2_RPF.sorted.param"}
		};
				
		//MergeResults test = new MergeResults(scoreOutFiles, bedCovPlusFiles, bedCovMinusFiles, paramFiles);
		//test.merge("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/out.txt", 2);
	}	
}