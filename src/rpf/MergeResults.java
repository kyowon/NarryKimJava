package rpf;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.AnnotationFileParser;
import parser.BedCovFileParser;
import parser.MergedFileParser;
import parser.MergedFileParser.MergedResult;
import parser.MafParser;
import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;
import parser.ZeroBasedFastaParser;

public class MergeResults {
	private ScoringOutputParser[] scoringOutputParsers;
	private Scorer[] harrScorers = null; 
	private Scorer[] rpfScorers = null; 
	private Quantifier[] rpfQuantifiers = null;
	private Quantifier[] rnaQuantifiers = null;
	private MafParser mafParser = null;
	private AnnotationFileParser annotationFileParser = null;
	public MergeResults(String[] scoreOutFiles, 
			String[] harrCovPlusFiles, String[] harrCovMinusFiles, 
			String[] rpfCovPlusFiles, String[] rpfCovMinusFiles,
			String[] rnaCovPlusFiles, String[] rnaCovMinusFiles,
			String[] harrParamFiles,
			String[] rpfParamFiles,
			AnnotationFileParser annotationFileParser,
			ZeroBasedFastaParser fastaFileParser,
			MafParser mafParser){
		scoringOutputParsers = new ScoringOutputParser[scoreOutFiles.length];
		harrScorers = new Scorer[harrParamFiles.length];
		rpfScorers = new Scorer[rpfParamFiles.length];
		if(rpfCovPlusFiles!=null)
			rpfQuantifiers = new Quantifier[rpfCovPlusFiles.length];
		if(rnaCovPlusFiles!=null)
			rnaQuantifiers = new Quantifier[rnaCovPlusFiles.length];
		
		for(int i=0;i<scoringOutputParsers.length;i++){
			scoringOutputParsers[i] = new ScoringOutputParser(scoreOutFiles[i]);
			harrScorers[i] = new Scorer(harrCovPlusFiles[i], harrCovMinusFiles[i], harrParamFiles[i], annotationFileParser, fastaFileParser);
		}
		
		for(int i=0;i<rpfScorers.length;i++){			
			rpfScorers[i] = new Scorer(rpfCovPlusFiles[i], rpfCovMinusFiles[i], rpfParamFiles[i], annotationFileParser, fastaFileParser);
		}
		
		if(rpfQuantifiers!=null){
			for(int i=0;i<rpfQuantifiers.length;i++){
				rpfQuantifiers[i] = new Quantifier(rpfCovPlusFiles[i], rpfCovMinusFiles[i], annotationFileParser, fastaFileParser);
			}
		}
		if(rnaQuantifiers!=null){
			for(int i=0;i<rnaQuantifiers.length;i++){
				rnaQuantifiers[i] = new Quantifier(rnaCovPlusFiles[i], rnaCovMinusFiles[i], annotationFileParser, fastaFileParser);
			}
		}
		this.mafParser = mafParser;
		this.annotationFileParser = annotationFileParser;
	}
	
	public void merge(String outFile, double scoreThreshold, double rpfRPKMThreshold, double rnaRPKMThreshold, int positionQuantityChangeLength, int positionQuantityOffset, int maxLengthUntilStopcodon){
		try {
			PrintStream out = new PrintStream(outFile);
			PrintStream mout = new PrintStream(outFile.replace('.', '_')+".m");
			PrintStream tout = new PrintStream(outFile + ".train.csv");
			//ArrayList<MergedResult> mergedResults = new ArrayList<MergedResult>();
			HashSet<String> contigs = new HashSet<String>();
			for(int i=0;i<scoringOutputParsers.length;i++){
				contigs.addAll(scoringOutputParsers[i].getContigs());
			}
			boolean start = true;
			
			for(String contig : contigs){
				System.out.println("Merging " + contig);
				for(ScoredPosition position : ScoringOutputParser.getUnionPositions(scoringOutputParsers, rpfQuantifiers, rnaQuantifiers, contig, scoreThreshold, rpfRPKMThreshold, rnaRPKMThreshold)){
					MergedResult mr = new MergedFileParser().new MergedResult(position, harrScorers, rpfScorers, rpfQuantifiers, rnaQuantifiers, mafParser, positionQuantityChangeLength, positionQuantityOffset, maxLengthUntilStopcodon);
					if(start){
						out.println(mr.getHeader());
						tout.println(mr.GetTrainingSetHeader());
						start = false;
					}
					out.println(mr);
					
					for(int index : position.getConditionIndices()){
						String ts = mr.ToTrainingSetString(annotationFileParser, index);
						if(ts != null) tout.println(ts);
					}
					
					writeMfile(mout, position);	
				}
			}
			out.close();
			mout.close();	
			tout.close();
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
		for(int i=0;i<harrScorers.length;i++){
			mout.println("{"+(i+1)+"}=[");
			BedCovFileParser sbp =  position.isPlusStrand()? harrScorers[i].getBedCovPlusFileParser() : harrScorers[i].getBedCovMinusFileParser();
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