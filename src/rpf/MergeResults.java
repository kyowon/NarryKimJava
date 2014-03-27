package rpf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
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
	private ScoringOutputParser[] rpfScoringOutputParsers;
	private Scorer[] harrScorers = null; 
	private Scorer[] rpfScorers = null; 
	private Quantifier[] rpfQuantifiers = null;
	private Quantifier[] rnaQuantifiers = null;
	private MafParser mafParser = null;
	private AnnotationFileParser annotationFileParser = null;
	public MergeResults(String[] rpfScoreOutFiles, 
			String[] harrCovPlusFiles, String[] harrCovMinusFiles, 
			String[] rpfCovPlusFiles, String[] rpfCovMinusFiles,
			String[] rnaCovPlusFiles, String[] rnaCovMinusFiles,
			String[] harrParamFiles,
			String[] rpfParamFiles,
			AnnotationFileParser annotationFileParser,
			ZeroBasedFastaParser fastaFileParser,
			MafParser mafParser){
		rpfScoringOutputParsers = new ScoringOutputParser[rpfScoreOutFiles.length];
		harrScorers = new Scorer[harrParamFiles.length];
		rpfScorers = new Scorer[rpfParamFiles.length];
		if(rpfCovPlusFiles!=null)
			rpfQuantifiers = new Quantifier[rpfCovPlusFiles.length];
		if(rnaCovPlusFiles!=null)
			rnaQuantifiers = new Quantifier[rnaCovPlusFiles.length];
		
		for(int i=0;i<rpfScoringOutputParsers.length;i++){
			rpfScoringOutputParsers[i] = new ScoringOutputParser(rpfScoreOutFiles[i]);			
		}
		
		for(int i=0;i<harrScorers.length;i++){		
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
	
	public void merge(String outFile, double scoreThreshold, int positionQuantityChangeLength, int positionQuantityOffset, int maxLengthUntilStopcodon, HashSet<String> allowedCodons){
		try {
			String dir = new File(outFile).getParent();
			new File(dir +System.getProperty("file.separator") +  new File(outFile).getName() + "_classification").mkdir();
		//	new File(dir +System.getProperty("file.separator") + new File(outFile).getName() + "_test").mkdir();
		//	new File(dir +System.getProperty("file.separator") +  new File(outFile).getName() + "_m").mkdir();
			
			PrintStream[] trainOut = new PrintStream[rpfScorers.length];//(dir + System.getProperty("file.separator") + new File(outFile).getName() + "_train" + System.getProperty("file.separator") + "train.arff");
			PrintStream[] testOut = new PrintStream[rpfScorers.length]; 
			String[] trainString = new String[rpfScorers.length];
			String[] testString = new String[rpfScorers.length];
			
			PrintStream mout = new PrintStream(dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "a.m");
			PrintStream mgout = new PrintStream(dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "g.m");
			
			for(int i=0;i<testOut.length;i++){
				testString[i] = dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "test_" + i + ".arff";
				testOut[i] = new PrintStream(testString[i]);
				trainString[i] = dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "train_" + i + ".arff";
				trainOut[i] = new PrintStream(trainString[i]);
			}
			
			
			//ArrayList<MergedResult> mergedResults = new ArrayList<MergedResult>();
			HashSet<String> contigs = new HashSet<String>();
			for(int i=0;i<rpfScoringOutputParsers.length;i++){
				contigs.addAll(rpfScoringOutputParsers[i].getContigs());
			}
			boolean start = true;
			boolean harrRPFsynced = harrScorers.length == rpfQuantifiers.length;
			
			ArrayList<MergedResult> mrs = new ArrayList<MergedResult>();
			
			for(String contig : contigs){
				System.out.println("Merging " + contig);
				for(ScoredPosition position : ScoringOutputParser.getUnionPositions(rpfScoringOutputParsers, annotationFileParser, contig, scoreThreshold, positionQuantityChangeLength, positionQuantityOffset)){
					//if(position.getGenomicRegion().equals("NM_3_UTR")){
						
					//}else if(!position.getCodon().equals("ATG") && !position.getCodon().equals("CTG")) continue;
					
					//if(position.getCodon() != )
					MergedResult mr = new MergedFileParser().new MergedResult(position, harrScorers, rpfScorers, rpfQuantifiers, rnaQuantifiers, mafParser, positionQuantityChangeLength, positionQuantityOffset, maxLengthUntilStopcodon);
					if(start){						
						for(int i=0;i<testOut.length;i++) trainOut[i].println(mr.GetArffSetHeader());
						for(int i=0;i<testOut.length;i++) testOut[i].println(mr.GetArffSetHeader());		
						mout.println(mr.GetArffSetMfileHeader());
						mgout.println("group={");
						start = false;
					}
					
					if(Math.abs(mr.getStopPostion() - position.getPosition()) > 25){
					//	if((position.getCodon().equals("ATG") || position.getCodon().equals("CTG"))) {
						if(allowedCodons.contains(position.getCodon())){
							//out.println(mr);
							mrs.add(mr);
							for(int i=0;i<rpfScorers.length;i++){
								String test = mr.ToTestSetString(i, harrRPFsynced);
								testOut[i].println(test);
							}
						}
						//}					
						
						for(int rpfrnaIndex : position.getRPFRNAIndices()){
							String ts = mr.ToTrainingSetString(rpfrnaIndex, harrRPFsynced);
							if(ts != null){
								trainOut[rpfrnaIndex].println(ts);
								mout.println(ts.replace("?", "NaN").substring(0, ts.lastIndexOf(',')));
								mgout.println("'" + position.getGenomicRegion().replace('_', ' ') + (position.isAnnotated()? "(Translation Start)" : "")+ "'");
							}
							
						}
						
					//	writeMfile(mout, position);
					}	
				}
			}
			
			mout.println("];");
			mgout.println("};");
			
			System.out.println("Merged positions : " + mrs.size());
			for(int i=0;i<testOut.length;i++){
				Classifier.evaluate(trainString[i], true);
				Classifier.classify(trainString[i], testString[i], mrs, i, true);
			}
			PrintStream out = new PrintStream(outFile);
			
			//out.println(mr.getHeader()); 
			for(int i=0; i<mrs.size();i++){
				MergedResult mr = mrs.get(i);
				if(i == 0) out.println(mr.getHeader());
				out.println(mr);
			}
			
			
			out.close();
			mout.close();	
			mgout.close();
			for(int i=0;i<testOut.length;i++){
				testOut[i].close();
				trainOut[i].close();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}	
	
	
	private void writeMfile(PrintStream mout, ScoredPosition position){
		
		
		/*String gname;
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
			*/		
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