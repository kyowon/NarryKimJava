package rpf;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import jspp.SearchMatrix;
import jspp.SignalPeptidePredictor;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.MafParser;
import parser.ZeroBasedFastaParser;
import rpf.parser.ScoringOutputParser;
import rpf.parser.MergedFileParser.MergedResult;
import rpf.parser.ScoringOutputParser.ScoredPosition;

public class MergeResults {
	private ScoringOutputParser[] rpfScoringOutputParsers;
	private Scorer[] harrScorers = null; 
	private Scorer[] rpfScorers = null; 
	private Quantifier[] rpfQuantifiers = null;
	private Quantifier[] rnaQuantifiers = null;
	private MafParser mafParser = null;
	private AnnotationFileParser annotationFileParser = null;
	private SignalPeptidePredictor pd = null;
	private ZeroBasedFastaParser fastaFileParser = null;
	private int[][] groups;
	private String[] harrBedFiles;
	private String[] rpfBedFiles;
	private String[] rnaBedFiles;
	private String[] harrParamFiles;
	private String[] rpfParamFiles;
	
	public MergeResults(String[] rpfScoreOutFiles, 
			String[] harrBedFiles, 
			String[] rpfBedFiles,
			String[] rnaBedFiles,
			String[] harrParamFiles,
			String[] rpfParamFiles,
			int[][] groups,
			AnnotationFileParser annotationFileParser,
			ZeroBasedFastaParser fastaFileParser,
			MafParser mafParser){
		this.groups = groups;
		this.mafParser = mafParser;
		this.annotationFileParser = annotationFileParser;
		
		this.harrBedFiles = harrBedFiles;
		this.rpfBedFiles = rpfBedFiles;
		this.rnaBedFiles = rnaBedFiles;
		this.harrParamFiles = harrParamFiles;
		this.rpfParamFiles = rpfParamFiles;
		
		rpfScoringOutputParsers = new ScoringOutputParser[rpfScoreOutFiles.length];
		
		for(int i=0;i<rpfScoringOutputParsers.length;i++){
			rpfScoringOutputParsers[i] = new ScoringOutputParser(rpfScoreOutFiles[i]);			
		}
		///media/kyowon/Data1/RPF_Project/tools/matrices(SearchMatrix-objects)/eukarya.smx
		File pmatrix=new File("./../tools/matrices(SearchMatrix-objects)/eukarya.smx");
	      ObjectInputStream in;
		try {
			in = new ObjectInputStream(new FileInputStream(pmatrix));
			SearchMatrix smp;
			try {
				smp = (SearchMatrix) in.readObject();
				pd = new SignalPeptidePredictor(smp);
			} catch (ClassNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		    in.close();    
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.fastaFileParser = fastaFileParser;
	}
	
	private void initiate(String contig){
		if(harrBedFiles != null){
			harrScorers = new Scorer[harrParamFiles.length];			
			for(int i=0;i<harrScorers.length;i++){		
				Bed12Parser harrParser = new Bed12Parser(harrBedFiles[i], annotationFileParser, contig);
				harrScorers[i] = new Scorer(harrParser, harrParamFiles[i], annotationFileParser, fastaFileParser);
			}
		}
		
		rpfScorers = new Scorer[rpfBedFiles.length];
		rpfQuantifiers = new Quantifier[rpfBedFiles.length];
		if(rnaBedFiles!=null)
			rnaQuantifiers = new Quantifier[rnaBedFiles.length];
		
		for(int i=0;i<rpfBedFiles.length;i++){
			Bed12Parser rpfParser = new Bed12Parser(rpfBedFiles[i], annotationFileParser, contig);
			rpfScorers[i] = new Scorer(rpfParser, rpfParamFiles[i], annotationFileParser, fastaFileParser);
		//	for(int i=0;i<rpfQuantifiers.length;i++){
			rpfQuantifiers[i] = new Quantifier(rpfParser, annotationFileParser, fastaFileParser);
		//	}
		
		}
		
		
		if(rnaQuantifiers!=null){
			for(int i=0;i<rnaBedFiles.length;i++){
				Bed12Parser rnaParser = new Bed12Parser(rnaBedFiles[i], annotationFileParser, contig);	
				rnaQuantifiers[i] = new Quantifier(rnaParser, annotationFileParser, fastaFileParser);
			}
		}
	
	}
	
	public void merge(String outFile, double scoreThreshold, int positionQuantityChangeLength, int positionQuantityOffset, int maxLengthUntilStopcodon, HashSet<String> allowedCodons){
		try {
			String dir = new File(outFile).getParent();
			new File(dir +System.getProperty("file.separator") +  new File(outFile).getName() + "_classification").mkdir();
		//	new File(dir +System.getProperty("file.separator") + new File(outFile).getName() + "_test").mkdir();
		//	new File(dir +System.getProperty("file.separator") +  new File(outFile).getName() + "_m").mkdir();
			
			PrintStream[] trainOut = new PrintStream[rpfBedFiles.length];//(dir + System.getProperty("file.separator") + new File(outFile).getName() + "_train" + System.getProperty("file.separator") + "train.arff");
			PrintStream[] testOut = new PrintStream[rpfBedFiles.length]; 
			String[] trainString = new String[rpfBedFiles.length];
			String[] testString = new String[rpfBedFiles.length];
			
			PrintStream mout = new PrintStream(dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "a.m");
			PrintStream mgout = new PrintStream(dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "g.m");
			
			for(int i=0;i<testOut.length;i++){
				testString[i] = dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "test_" + i + ".arff";
				testOut[i] = new PrintStream(testString[i]);
				trainString[i] = dir + System.getProperty("file.separator") + new File(outFile).getName()  + "_classification" + System.getProperty("file.separator") + "train_" + i + ".arff";
				trainOut[i] = new PrintStream(trainString[i]);
			}
						
			//ArrayList<MergedResult> mergedResults = new ArrayList<MergedResult>();
			HashSet<String> contigs = new HashSet<String>(annotationFileParser.getContigs());
			
			boolean start = true;
			
			ArrayList<MergedResult> mrs = new ArrayList<MergedResult>();
			
			for(String contig : contigs){
				//contig = "chr5";
				System.out.println("Merging " + contig);
				initiate(contig);
				boolean harrRPFsynced = harrScorers != null && harrScorers.length == rpfQuantifiers.length;
				
				ScoredPosition prevPosition = null;
				for(ScoredPosition position : ScoringOutputParser.getUnionPositions(rpfScoringOutputParsers, annotationFileParser, contig, scoreThreshold, positionQuantityChangeLength, positionQuantityOffset)){
					//if(position.getGenomicRegion().equals("NM_3_UTR")){
					if(prevPosition !=null && !prevPosition.equals(position)){
						if(position.getGene() == null){
							if(prevPosition.getGene() == null) continue;
						}else{
							if(position.getGene().equals(prevPosition.getGene())) continue;
						}
						//(position.getGene() == null && prevPosition.getGene() == null || position.getGene().equals(prevPosition.getGene()))) continue;	
					}
					//}else if(!position.getCodon().equals("ATG") && !position.getCodon().equals("CTG")) continue;
					
					//if(position.getCodon() != )
					MergedResult mr = new MergedResult(position, harrScorers, rpfScorers, rpfQuantifiers, rnaQuantifiers, groups, fastaFileParser, mafParser, pd, positionQuantityChangeLength, positionQuantityOffset, maxLengthUntilStopcodon);
					
					if(mr.getLength() < 25) continue;
					
					if(start){						
						for(int i=0;i<testOut.length;i++) trainOut[i].println(mr.GetArffSetHeader());
						for(int i=0;i<testOut.length;i++) testOut[i].println(mr.GetArffSetHeader());		
						mout.println(mr.GetArffSetMfileHeader());
						mgout.println("group={");
						start = false;
					}
					
					//if(Math.abs(mr.getStopPostion() - position.getPosition()) > 25){
					//	if((position.getCodon().equals("ATG") || position.getCodon().equals("CTG"))) {
						if(allowedCodons.contains(position.getSequence().substring(0, Math.min(position.getSequence().length(), 3)))){
							//out.println(mr);
							//if(!mrs.isEmpty() && mrs.get(mrs.size()-1).get)
							mrs.add(mr);
							for(int i=0;i<rpfScorers.length;i++){
								String test = mr.ToTestSetString(i, harrRPFsynced);
								testOut[i].println(test);
							}
						}
						//}					
						
						if(prevPosition !=null && !prevPosition.equals(position)){
							for(int rpfrnaIndex : position.getRPFRNAIndices()){
								String ts = mr.ToTrainingSetString(rpfrnaIndex, harrRPFsynced);
								if(ts != null){
									trainOut[rpfrnaIndex].println(ts);
									mout.println(ts.replace("?", "NaN").substring(0, ts.lastIndexOf(',')));
									//mgout.println("'" + position.getGenomicRegion().replace('_', ' ') + (position.isAnnotated()? "(Translation Start)" : "")+ "'");
									mgout.println("'" + ts.substring(ts.lastIndexOf(',')+1)  + "'");
								}	
							}
						}
						
					//	writeMfile(mout, position);
					//}	
					
					prevPosition = position;
				}
			}
			
			mout.println("];");
			mgout.println("};");
			//TODO make scorethreshold..
			double[] predictionThresholds = {0.25, 0.5, 0.75, 0.9};
			
			PrintStream out = new PrintStream(outFile);
			PrintStream out_statistics = new PrintStream(outFile.substring(0, outFile.lastIndexOf(".")) + ".statistics.txt");
			
			PrintStream[] out_diff_high = new PrintStream[predictionThresholds.length]; 
			
			for(int i=0;i<testOut.length;i++){
				Classifier.classify(trainString[i], testString[i], mrs, i, true, .5);
			}
				
			
			for(int j=0;j<predictionThresholds.length;j++){
				System.out.println("Merged positions : " + mrs.size());
				out_statistics.println("Prediction score threshold: " + predictionThresholds[j]);
				System.out.println("Prediction score threshold: " + predictionThresholds[j]);
				for(int i=0;i<testOut.length;i++){
					String stat = Classifier.evaluate(trainString[i], true, predictionThresholds[j]);
					out_statistics.println(stat);
					System.out.println(stat);
				//	Classifier.classify(trainString[i], testString[i], mrs, i, true, .5);
				}
			
				
				out_diff_high[j] = new PrintStream(outFile.substring(0, outFile.lastIndexOf("."))+".diff_higher_than_" + predictionThresholds[j] + ".csv");				
				
				
				for(int i=0; i<mrs.size();i++){
					MergedResult mr = mrs.get(i);
					if(j == 0){
						if(i == 0){
							out.println(mr.getHeader());					
						}
						out.println(mr);
					}
					if(i == 0){
						out_diff_high[j].println(mr.getHeader());		
					}
					if(!mr.isPredictionCorrect()){
						boolean high = false;
						double[] scores = mr.getPredictedClassesScores();
						for(double score : scores){
							if(score >= predictionThresholds[j]){
								high = true;
								break;
							}
						}
						if(high) out_diff_high[j].println(mr);
						//else lowScoredResults.add(mr);
					}
				}				
				out_diff_high[j].close();
			}
			out_statistics.close();
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
	

	

}