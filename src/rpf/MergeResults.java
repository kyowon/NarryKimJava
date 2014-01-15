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

public class MergeResults {
	private BedCovFileParser[][] bedCovPlusFileParsers; // [][0] : harr [][1] : RPF TODO make "sample number" and replace it with '2'
	private BedCovFileParser[][] bedCovMinusFileParsers;
	private ScoringOutputParser[] scoringOutputParsers;
	private Scorer[][] scorers;
	
	public MergeResults(String[] scoreOutFiles, String[][] bedCovPlusFiles, String[][] bedCovMinusFiles, String[][] paramFiles, String annotationFile){
		int n = scoreOutFiles.length;
		scoringOutputParsers = new ScoringOutputParser[n];
		bedCovPlusFileParsers = new BedCovFileParser[n][2];
		bedCovMinusFileParsers = new BedCovFileParser[n][2];
		scorers = new Scorer[n][2];
		
		for(int i=0; i<n;i++){
			scoringOutputParsers[i] = new ScoringOutputParser(scoreOutFiles[i]);
			for(int j=0;j<2;j++){
				bedCovPlusFileParsers[i][j] = new BedCovFileParser(bedCovPlusFiles[i][j], annotationFile);
				bedCovMinusFileParsers[i][j] = new BedCovFileParser(bedCovMinusFiles[i][j], annotationFile);
				scorers[i][j] = new Scorer(bedCovPlusFiles[i][j], bedCovMinusFiles[i][j], paramFiles[i][j], annotationFile);
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
					String gname;;
					if(!position.isAnnotated()) gname = contig+"_"+position.getPosition()+"_"+(position.isPlusStrand()? "P" : "M");
					else gname = position.getGene().getGeneName();
					mout.println(gname + "=[");
					mergedResults.add(getSinglePositionInformation(position, mout));
					mout.println("]';");					
				}
			}
			
			ArrayList<Double> values = new ArrayList<Double>();
			for(MergedResult m : mergedResults){
				values.add(m.getScoreRatio());
			}
			for(MergedResult m : mergedResults){
				m.setScorePvalue(getPvalue(values, m.getScoreRatio()));
			}
			values.clear();
			for(MergedResult m : mergedResults){
				values.add(m.getQuantityRatio());
			}
			for(MergedResult m : mergedResults){
				m.setQuantityPvalue(getPvalue(values, m.getQuantityRatio()));
			}
			out.println(getHeader());
			for(MergedResult m : mergedResults){
				out.println(m);
			}
			out.close();
			appendPvalues(outFile,6,7);
			appendPvalues(outFile,8,9);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
		
	private void appendPvalues(String outFile, int valueColumn, int pvalueColumn){
		ArrayList<Double> values = new ArrayList<Double>();
		try {
            BufferedReader br = new BufferedReader(new FileReader(outFile));
            String s;
            while ((s = br.readLine()) != null) {
            	if(s.startsWith("#")) continue;
            	values.add(Double.parseDouble(s.split("\t")[valueColumn]));
            }
            br.close();            
            Collections.sort(values);
            br = new BufferedReader(new FileReader(outFile));
            ArrayList<String> lines = new ArrayList<String>();
            while ((s = br.readLine()) != null) {
            	if(s.startsWith("#")){
            		lines.add(s + "\n");
            		continue;
            	}
            	String[] token = s.split("\t");
            	double value = Double.parseDouble(token[valueColumn]);
            	Double pval = getPvalue(values, value);
            	token[pvalueColumn] = pval.toString();
            	StringBuffer line = new StringBuffer();
            	for(int i=0;i<token.length-1;i++){
            		line.append(token[i]);line.append("\t");
            	}
            	line.append(token[token.length-1]);line.append("\n");
            	lines.add(line.toString());
            }
            br.close();
            BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
            for(String line : lines){
            	bw.write(line);
            }
            bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private double getPvalue(ArrayList<Double> values, double value){
		int index = Collections.binarySearch(values, value);
		if(index < 0) index = - index - 1;
		double pval = (double) index / values.size();
		return pval;//Math.min(pval, 1-pval);
	}
	
	private String getHeader(){
		String ret = "#Contig\tPosition\tStrand\tCodon\t";
		for(int i=0;i<scorers.length;i++){
			ret+= "Scores" + (i+1); 
			if(i<scorers.length-1) ret+=";";
			else ret+="\t";
		}
		/*for(int i=0;i<scorers.length;i++){
			ret+= "HarrQuantities" + (i+1);
			if(i<scorers.length-1) ret+=",";
			else ret+="\t";
		}*/
		for(int i=0;i<scorers.length;i++){
			ret+= "RPFQuantities" + (i+1);
			if(i<scorers.length-1) ret+=";";
			else ret+="\t";
		}
		ret += "mostExtremeHarrScoreRatio\tmostExtremeHarrScoreRatioPvalue\tmostExtremeRPFQuantityRatio\tmostExtremeRPFQuantityRatioPvalue";
		ret += "\tgenomicRegion\tframeShift\tIsAnnotated\t" + AnnotatedGene.getHeader();
		return ret;
	}
		
	private MergedResult getSinglePositionInformation(ScoredPosition position, PrintStream outMFileStream){
		BedCovFileParser[][] bedCovFileParsers = position.isPlusStrand()? bedCovPlusFileParsers : bedCovMinusFileParsers; 
		double[] scores = new double[scorers.length];
		double[] quantities = new double[scorers.length];
		
		for(int i=0;i<scorers.length;i++){
			double[][] cov = new double[2][];
			for(int j=0;j<2;j++){
				cov[j] = bedCovFileParsers[i][j].getSqrtCoverages(position.getContig(), position.getPosition(), scorers[i][j].getLeftWindowSize(), scorers[i][j].getRightWindowSize(), position.isPlusStrand());
				double[] covForMfile = bedCovFileParsers[i][j].getCoverages(position.getContig(), position.getPosition(), 30, 100, position.isPlusStrand());
				for(Double c : covForMfile){
					outMFileStream.append(c.toString());
					outMFileStream.append('\t');
				}
				outMFileStream.append('\n');		
			}		
			scores[i] = scorers[i][0].getLRScore(scorers[i][0].getRawScore(cov[0]));		
			quantities[i] = scorers[i][1].getQuantity(cov[1], true);
		}
		
		return new RPFMergedFileParser().new MergedResult(position, scores, quantities);
				
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