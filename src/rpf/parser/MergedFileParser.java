package rpf.parser;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import jspp.SignalPeptidePredictor;
import launcher.PhyloPLauncher;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.MafParser;
import parser.ZeroBasedFastaParser;
import parser.AnnotationFileParser.AnnotatedGene;
import rpf.Quantifier;
import rpf.Scorer;
import rpf.parser.ScoringOutputParser.ScoredPosition;
import util.Codon;
import util.DnDsCalculator;

public class MergedFileParser {
	private ArrayList<MergedResult> list;
	public ArrayList<MergedResult> getList(){ return list; }
	public static class MergedResult{
		//	IsAnnotated	ContainingGeneName	ContainingGBGeneName	txStart	txEnd	cdsStart	cdsEnd	genomicRegion	frameShift
		//private int positionQuantityChangeLength;
		//private int positionQuantityOffset;
		//private int maxLengthUntilStopcodon;
		
		private String contig;
		private ScoredPosition position;
		private boolean isPlusStrand;
		private String seq;
		private String annotatedClass;
		private String[] predictedClasses = null;
		private double[] predictedClassesScores;
		private double[] harrStartScores;
		private double[] harrStopScores;
		private double[] rpfStartScores;		
		private double[] rpfStopScores;
		private double[] rpfCDSQuantities;//
		private double[] cdsTE;//
		private double[] rpfPositionQuantities;//
		private double[] positionTE;//	
		private double[] rpfPositionQuantityChanges;//
		private double[] rpfCDSRPKMchanges;
		private double[] releaseScore;//
		private int stopPosition = -1;
		private double[] rnaCDSQuantities;//
		private double[] rnaPositionQuantities;//
		private double[] rnaPositionQuantityChanges;//
		private String genomicRegion;
		private String frameShift;
		private double dndsRatioAfter = -1;
		private double dndsRatioBefore = -1;
		private double phyloP = 1;
		private double signalPeptideScore = -1;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		//private int numSample = 0;
		private int[][] groups;
		
		public String getContig() {
			return contig;
		}

		public int getPosition() {
			return position.getPosition();
		}

		public String getAnnotatedClass(){
			return annotatedClass;
		}
		public ScoredPosition getScoredPosition(){
			return position;
		}
		
		public boolean isPlusStrand() {
			return isPlusStrand;
		}

		public String getSequence() {
			return seq;
		}

		public String[] getPredictedClasses() {
			return predictedClasses;
		}

		public double[] getPredictedClassesScores() {
			return predictedClassesScores;
		}

		public double[] getHarrStartScores() {
			return harrStartScores;
		}

		public double[] getHarrStopScores() {
			return harrStopScores;
		}

		public double[] getRpfStartScores() {
			return rpfStartScores;
		}

		public double[] getRpfStopScores() {
			return rpfStopScores;
		}

		public double[] getRpfCDSQuantities() {
			return rpfCDSQuantities;
		}

		public double[] getCdsTE() {
			return cdsTE;
		}

		public double[] getRpfPositionQuantities() {
			return rpfPositionQuantities;
		}

		public double[] getPositionTE() {
			return positionTE;
		}

		public double[] getRpfPositionQuantityChanges() {
			return rpfPositionQuantityChanges;
		}

		public double[] getReleaseScore() {
			return releaseScore;
		}

		public int getStopPosition() {
			return stopPosition;
		}

		public double[] getRnaCDSQuantities() {
			return rnaCDSQuantities;
		}

		public double[] getRnaPositionQuantities() {
			return rnaPositionQuantities;
		}

		public double[] getRnaPositionQuantityChanges() {
			return rnaPositionQuantityChanges;
		}

		public String getGenomicRegion() {
			return genomicRegion;
		}

		public String getFrameShift() {
			return frameShift;
		}

		public void setFrameShift(String f){
			frameShift = f;
		}
		public double getDnDsRatioAfter() {
			return dndsRatioAfter;
		}

		public double getPhyloP(){
			return phyloP;
		}
		
		public double getDnDsRatioBefore() {
			return dndsRatioBefore;
		}

		public boolean isAnnotated() {
			return isAnnotated;
		}

		public AnnotatedGene getGene() {
			return gene;
		}
		
		public int getLength(){
			return seq.length();
			//if(stopPosition < 0) return -1;
			//return Math.abs(position.getPosition() - stopPosition) + 1;
		}
		public boolean isPredictionCorrect(){
			boolean correct = true;
			for(String prediction : predictedClasses){
				if(prediction.toUpperCase().equals(annotatedClass)) continue;
				correct = false;
				break;
			}
			return correct;
		}
		
		public MergedResult(ScoredPosition scoredPosition, Scorer[] harrScorers, Scorer[] rpfScorers, 
				Quantifier[] rpfQuantifiers, Quantifier[] rnaQuantifiers, int[][] groups, ZeroBasedFastaParser fastaFileParser, MafParser mafParser,
				SignalPeptidePredictor pd, int positionQuantityChangeLength, int positionQuantityOffset, int maxLengthUntilStopcodon){
			this.contig = scoredPosition.getContig();
			this.position = scoredPosition;
			this.isPlusStrand = scoredPosition.isPlusStrand();
			this.seq = scoredPosition.getSequence();
			this.genomicRegion = scoredPosition.getGenomicRegion();
			this.frameShift = scoredPosition.getFrameShift();
			this.isAnnotated = scoredPosition.isAnnotated();
			this.gene = scoredPosition.getGene();
			this.groups = groups;
			
			stopPosition = -1;
			if(gene!=null && !genomicRegion.endsWith("ISOFORM")){
				ArrayList<Integer> lp = gene.getLiftOverPositionsTillNextStopCodon(position.getPosition(), 2, maxLengthUntilStopcodon, fastaFileParser);
				if(lp.size() > 0) stopPosition = lp.get(lp.size()-1);
				else{
					stopPosition = this.position.getPosition();// + (isPlusStrand? 3 : -3);
					return;
				}
			}else{
				//stopPosition = position.getPosition();				
				
				for(int i=0;i<seq.length()-2;i+=3){
					String codon = seq.substring(i, i+3);
					if(codon.equals("TAG") || codon.equals("TAA") || codon.equals("TGA")){		
						//stopPosition += (isPlusStrand? +3 : -3);
						break;
					}
					stopPosition = position.getCoordinate().get(i+2);//(isPlusStrand? +3 : -3);
				}
				//position.getCoordinate()
			}
			
			
			
			if(harrScorers != null){
				this.harrStartScores = new double[harrScorers.length];
					
				for(int i=0;i<this.harrStartScores.length;i++){
					this.harrStartScores[i] = harrScorers[i].getStartScore(this.position.getCoordinate(), isPlusStrand);//.getStartScore(this.contig, this.position, this.isPlusStrand, false);
				}
				
				if(gene!=null && !this.genomicRegion.endsWith("Intron") && !this.genomicRegion.endsWith("ISOFORM")){	
					this.harrStopScores = new double[harrScorers.length];
					for(int i=0;i<this.harrStopScores.length;i++){
						this.harrStopScores[i] = harrScorers[i].getStopScore(gene, stopPosition);
						//getStopScore(this.contig, this.position, this.isPlusStrand, false, maxLengthUntilStopcodon);
					}
				}
			}
			//if(rpfScorers !=null){
			//this.numSample = rpfScorers.length;
			this.rpfStartScores = new double[rpfScorers.length];
			this.predictedClasses = new String[rpfScorers.length];
			this.predictedClassesScores = new double[rpfScorers.length];
			for(int i=0;i<this.rpfStartScores.length;i++){
				this.rpfStartScores[i] = rpfScorers[i].getStartScore(this.position.getCoordinate(), this.isPlusStrand);
			}
			if(gene!=null && !this.genomicRegion.endsWith("Intron") && !this.genomicRegion.endsWith("ISOFORM")){
				this.rpfStopScores = new double[rpfScorers.length];				
				for(int i=0;i<this.rpfStopScores.length;i++){
					this.rpfStopScores[i] = rpfScorers[i].getStopScore(gene, stopPosition);
				}
			}
			//}
			
			
			String sequence = null;
			if(!fastaFileParser.containsContig(contig)){
				sequence = "N/A";
			}else if(isPlusStrand){
				sequence = fastaFileParser.getSequence(contig, position.getPosition(), position.getPosition() + 70*3);
			}else{							
				sequence = ZeroBasedFastaParser.getComplementarySequence(fastaFileParser.getSequence(contig, position.getPosition() - 70*3 + 1, position.getPosition() + 1), true);
			}	
			
			sequence =  Codon.getAminoAcids(sequence);
			int cpos = pd.predictEnhancedPosition(sequence);
			//pd.isSignalPeptide();
			signalPeptideScore = pd.getScore();			
			//
			if(this.isAnnotated){
				//System.out.println(sequence + " " + signalPeptideScore + " " + pd.isSignalPeptide());
			}			
	        
			if(rpfQuantifiers!=null){
				if(gene != null && !genomicRegion.startsWith("NR") && !this.genomicRegion.endsWith("ISOFORM")){				
					rpfCDSQuantities = new double[rpfQuantifiers.length];
					for(int i=0; i<rpfQuantifiers.length;i++){
						if(rpfQuantifiers[i]==null) continue;
						rpfCDSQuantities[i] = rpfQuantifiers[i].getCDSRPKM(gene);
					}
				}
				rpfPositionQuantities = new double[rpfQuantifiers.length];
				rpfPositionQuantityChanges = new double[rpfQuantifiers.length];
				if(gene!=null && !this.genomicRegion.endsWith("Intron") && !this.genomicRegion.endsWith("ISOFORM")) releaseScore = new double[rpfQuantifiers.length];
				if(gene!=null && this.genomicRegion.equals("NM_ORF") && !this.isAnnotated) rpfCDSRPKMchanges = new double[rpfQuantifiers.length];
					
				for(int i=0; i<rpfPositionQuantities.length;i++){
					if(rpfQuantifiers[i]==null) continue;
					rpfPositionQuantities[i] = rpfQuantifiers[i].getPositionRPKM(position.getContig(), isPlusStrand, position.getCoordinate());
					//.getPositionRPKM(contig, this.position, isPlusStrand, maxLengthUntilStopcodon, positionQuantityOffset);
					rpfPositionQuantityChanges[i] = rpfQuantifiers[i].getPositionQuantityChangeRatio(position, positionQuantityChangeLength, positionQuantityOffset);
									//.getPositionQuantityChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength, positionQuantityOffset);
					//rpfPositionQuantityChanges[i] = Math.log(rpfPositionQuantityChanges[i]+.0001)/Math.log(2);
					if(gene!=null && !this.genomicRegion.endsWith("Intron") && !this.genomicRegion.endsWith("ISOFORM")){
						releaseScore[i] = rpfQuantifiers[i].getReleaseScore(gene, stopPosition, positionQuantityChangeLength);
						//.getNextStopCodonQuantityChangeRatioNStopPosition(contig, this.position, isPlusStrand, positionQuantityChangeLength, maxLengthUntilStopcodon);
						
					}//else rpfAroundStopCodonQuantityChanges = null;
					if(rpfCDSRPKMchanges!=null){
						rpfCDSRPKMchanges[i] = rpfQuantifiers[i].getCDSRPKMChangeRatio(gene, position.getPosition());//
						//getgetCDSRPKMChangeRatio(gene, this.position, isPlusStrand, positionQuantityOffset);
					}
					
				}
				
				if(mafParser!=null){ 
					//this.dsdnRatioAfter = DsDnCalculator.calculate(mafParser.getSeqs(contig, position + (isPlusStrand? -30:30), isPlusStrand, Math.min(60, 30 + Math.abs(position - stopPosition))));
				//	int offset = 0;//cpos*3;
					String[] after = mafParser.getSeqs(contig, position.getCoordinate(), position.getPosition(), isPlusStrand, positionQuantityChangeLength);
					String[] before = mafParser.getSeqs(contig, position.getCoordinate(), position.getPosition() + (isPlusStrand? - positionQuantityChangeLength : positionQuantityChangeLength), isPlusStrand, positionQuantityChangeLength);
				//	String[] after = mafParser.getSeqs(contig, position + (isPlusStrand? offset : -offset), isPlusStrand , positionQuantityChangeLength);
				//	String[] before = mafParser.getSeqs(contig, position + (isPlusStrand? -positionQuantityChangeLength:positionQuantityChangeLength), isPlusStrand , positionQuantityChangeLength);
					this.dndsRatioAfter = DnDsCalculator.calculate(after); 
					this.dndsRatioBefore = DnDsCalculator.calculate(before);
					
					String mafString = mafParser.getSeqsInMafFormat(contig, position.getCoordinate(), position.getPosition(), isPlusStrand, positionQuantityChangeLength);
					this.phyloP = new PhyloPLauncher(mafString).getPvalConservation();
					//this.dsdnRatio = Math.log(dsdnRatio + .0001)/Math.log(2);
				}
				//rpfAroundStopCodonQuantityChanges
			}
			
			
			if(rnaQuantifiers!=null){
				if(gene != null && !genomicRegion.startsWith("NR")  && !this.genomicRegion.endsWith("ISOFORM")){				
					rnaCDSQuantities = new double[rnaQuantifiers.length];
					for(int i=0; i<rnaQuantifiers.length;i++){
						if(rnaQuantifiers[i]==null) continue;
						rnaCDSQuantities[i] = rnaQuantifiers[i].getCDSRPKM(gene);
					}
				}
				rnaPositionQuantities = new double[rnaQuantifiers.length];
				rnaPositionQuantityChanges = new double[rnaQuantifiers.length];
				for(int i=0; i<rnaPositionQuantities.length;i++){
					if(rnaQuantifiers[i]==null) continue;
					rnaPositionQuantities[i] = rnaQuantifiers[i].getPositionRPKM(position.getContig(), isPlusStrand, position.getCoordinate());//getPositionRPKM(contig, this.position, isPlusStrand, maxLengthUntilStopcodon, 0);
					rnaPositionQuantityChanges[i] = rnaQuantifiers[i].getPositionQuantityChangeRatio(position, positionQuantityChangeLength, positionQuantityOffset);
					//.getPositionQuantityChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength, 0);
					//rnaPositionQuantityChanges[i] = Math.log(rnaPositionQuantityChanges[i]+.0001)/Math.log(2);
				}
			}
			if(rpfQuantifiers!=null && rnaQuantifiers!=null){
				if(gene != null && !genomicRegion.startsWith("NR") && !this.genomicRegion.endsWith("Intron") && !this.genomicRegion.endsWith("ISOFORM")){
					cdsTE = new double[rpfQuantifiers.length];
					for(int i=0; i<cdsTE.length;i++){
						cdsTE[i] = rnaCDSQuantities[i] == 0?  0 : rpfCDSQuantities[i] / rnaCDSQuantities[i];
					//	rpfCDSQuantitiesNormalized[i] = Math.log(rpfCDSQuantitiesNormalized[i]+.0001)/Math.log(2);
					}
				}
				positionTE = new double[rpfQuantifiers.length];
				for(int i=0; i<positionTE.length;i++){
					positionTE[i] = rnaPositionQuantities[i] == 0? 0 : rpfPositionQuantities[i] / rnaPositionQuantities[i];
					
				//	rpfPositionQuantitiesNormalized[i] = Math.log(rpfPositionQuantitiesNormalized[i]+.0001)/Math.log(2);
				}					
			}
						
			if(isAnnotated)	annotatedClass = "TS";
			else if(genomicRegion.equals("NM_ORF")) annotatedClass = "T";
			else if(genomicRegion.equals("NM_5_UTR")) annotatedClass = "5U";
			else if(genomicRegion.equals("NM_3_UTR") || genomicRegion.endsWith("Intron")) annotatedClass = "NC";
			else annotatedClass = "_";
		}
		
		
		public static int[][] getGroups(String header){
			String s = header.substring(header.indexOf(':') + 1);
			String[] token = s.split("\t");
			int[][] groups = new int[token.length][];
			for(int i=0;i<token.length;i++){
				String[] t = token[i].split(";");
				groups[i] = new int[t.length];
				for(int j=0;j<t.length;j++){
					groups[i][j] = Integer.parseInt(t[j].trim());
				}
			}
			
			return groups;
		}
		
		public MergedResult(String s, int[][] groups){ 
			String[] token = s.split("\t");
			this.groups = groups;
			int i = 0;
			this.contig = token[i++];
			int position = Integer.parseInt(token[i++]);
			
			String[] mrStarts = token[i++].split(",");
			String[] mrEnds = token[i++].split(",");
			ArrayList<ArrayList<Integer>> mappedRegionStartsEnds = new ArrayList<ArrayList<Integer>>();
			for(int j=0;j<mrStarts.length;j++){
				ArrayList<Integer> splice = new ArrayList<Integer>();
				splice.add(Integer.parseInt(mrStarts[j]));
				splice.add(Integer.parseInt(mrEnds[j]));
				mappedRegionStartsEnds.add(splice);
			}
			
			this.isPlusStrand = token[i++].equals("+");
			
			
			
			this.seq = token[i++];
			StringBuilder geneString = new StringBuilder();
			geneString.append(token[i++]);
			geneString.append('\t');
			geneString.append(token[i++]);
			geneString.append('\t');
			geneString.append(this.contig);
			geneString.append('\t');
			
			this.genomicRegion = token[i++];
			
			this.annotatedClass = token[i];
			this.isAnnotated = token[i++].equals("TS");			
			
			this.predictedClasses = token[i++].split(";");
			i++;// div
			this.predictedClassesScores = subParse(token[i++]);
			this.frameShift = token[i++];
			this.dndsRatioBefore = Double.parseDouble(token[i++]);
			this.dndsRatioAfter = token[i].equals("_") ? -1 : Double.parseDouble(token[i]);
			i++;
			
			this.phyloP = Double.parseDouble(token[i++]);
			
			this.signalPeptideScore = token[i].equals("_") ? -1 : Double.parseDouble(token[i]);
			i++;
			this.stopPosition = token[i].equals("_") ? -1 : Integer.parseInt(token[i]);
			i++;
			i++;//div
			this.harrStartScores = subParse(token[i++]);
			i++;//div
			this.harrStopScores = subParse(token[i++]);
			i++;//div
			this.rpfStartScores = subParse(token[i++]);
			i++;//div
			this.rpfStopScores = subParse(token[i++]);
			i++;//div
			this.rpfPositionQuantities = subParse(token[i++]);
			i++;//div
			this.rnaPositionQuantities = subParse(token[i++]);
			i++;//div
			this.positionTE = subParse(token[i++]);
			i++;//div
			this.rpfPositionQuantityChanges = subParse(token[i++]);
			i++;//div
			this.rpfCDSRPKMchanges = subParse(token[i++]);
			i++;//div
			this.releaseScore = subParse(token[i++]);
			i++;//div
			this.rnaPositionQuantityChanges = subParse(token[i++]);
			i++;//div
			this.rpfCDSQuantities = subParse(token[i++]);
			i++;//div
			this.rnaCDSQuantities = subParse(token[i++]);
			i++;//div
			this.cdsTE = subParse(token[i++]);
			i++;//div
			i++;//te diff
			
			boolean isGenic = true;
			for(;i<token.length;i++){
				if(token[i].equals("_")){
					isGenic = false;
					break;
				}
				geneString.append(token[i]);
				if(i<token.length-1) geneString.append('\t');
			}
			//System.out.println(geneString);
			if(isGenic) this.gene = new AnnotatedGene(geneString.toString());
			
			this.position = 
					new ScoringOutputParser().new ScoredPosition(contig, position, Bed12Parser.getCoordinate(mappedRegionStartsEnds, isPlusStrand), isPlusStrand, 1, seq, gene, isAnnotated, genomicRegion, frameShift);	
		
		}
		
		public void setPredictedClasses(String predictedClass, double predictedClassesScore, int index){
			this.predictedClasses[index] = predictedClass; 
			this.predictedClassesScores[index] = predictedClassesScore;
		}
		
		
		private double[] subParse(String s){
			String[] token = s.split(";");
			if(token.length == 0 || token[0].equals("_")) return null;
			double[] ret = new double[token.length];
			for(int i=0;i<ret.length;i++){
				if(token[i].equals("_")) continue;
				ret[i] = Double.parseDouble(token[i]);
			}
			return ret;
		}
		
		public String getHeader(){
			String groupHeader = "@Group: ";
			
			for(int i=0; i<groups.length;i++){
				int[] group = groups[i];
				for(int j=0;j<group.length;j++){
					groupHeader += group[j];
					if(j<group.length-1) groupHeader += ";";
				}
				if(i<groups.length-1) groupHeader += "\t";
			}
			
			String header =  groupHeader + "\nContig\tPosition\tReadStarts\tReadEnds\tStrand\tSequence\t" + AnnotatedGene.getSimpleHeader() + "\tGenomicRegion\tAnnotated\t";
			header += subGetHeader("Predicted", groups) + "\t";
			header += subGetHeader("PredictionScore", groups, false) + "\t";
					//"Predicted\tPredictionScore\t";
			header += "Frame\tDnDsBefore\tDnDsAfter\tPhyloP\tSignalPeptideScore\tStopCodonPosition\tLength\t";
			//header += "";
			header += subGetHeader("HarrMatchScore", harrStartScores == null? null : new int[harrStartScores.length][1]) + "\t";
			header += subGetHeader("HarrStopCodonMatchScore", harrStartScores == null? null : new int[harrStartScores.length][1]) + "\t";
			header += subGetHeader("RPFMatchScore", groups) + "\t";
			header += subGetHeader("RPFStopCodonMatchScore", groups) + "\t";
			
			header += subGetHeader("RPFRPKM", groups) + "\t";
			header += subGetHeader("RNARPKM", groups) + "\t";
			header += subGetHeader("TE", groups) + "\t"; 
			
			header += subGetHeader("RPFFoldChange", groups) + "\t";
			header += subGetHeader("RPFFoldChange(a/b position)", groups) + "\t";
			header += subGetHeader("Release", groups) + "\t";			
		
			header += subGetHeader("RNAFoldChange", groups) + "\t";
			
			header += subGetHeader("RPFCDSRPKMs", groups) + "\t";
			header += subGetHeader("RNACDSRPKMs", groups) + "\t";
			header += subGetHeader("CDSTE", groups) + "\t";
			
			header += "TEDiff(Pos vs. CDS)\t";
			
			header += AnnotatedGene.getDetailedInfoHeader();
			
			return header;
		}
		
		public String getSimpleHeader(){
			String header = "Contig\tPosition\tStopPosition\tStrand\tCodon\t" + AnnotatedGene.getSimpleHeader() + "\tGenomicRegion\tFrame\tAnnotated\tDnDs\tPhyloP\t";
			header += subGetHeader("Predicted", groups) + "\t";
			
			header += "Length\t";
			header += subGetHeader("HarrMatchScore", harrStartScores == null? null : new int[harrStartScores.length][1]) + "\t";
			header += subGetHeader("HarrStopCodonMatchScore", harrStartScores == null? null : new int[harrStartScores.length][1]) + "\t";
			header += subGetHeader("RPFMatchScore", groups) + "\t";
			header += subGetHeader("RPFStopCodonMatchScore", groups) + "\t";
			
			header += subGetHeader("RPFRPKM", groups) + "\t";
			header += subGetHeader("TE", groups) + "\t"; 
			
			header += subGetHeader("RPFFoldChange", groups) + "\t";
			header += subGetHeader("Release", groups) + "\t";			
			
			return header;
		}
		
		private String subGetHeader(String prefix, int[][] groups){
			return subGetHeader(prefix, groups, true);
		}
		
		private String subGetHeader(String prefix, int[][] groups, boolean appendUnuniformity){
			String header = "";
			
			if(groups == null){
				header += prefix;
			}else{
				HashMap<Integer, String> suffixMap = new HashMap<Integer, String>();
				int sm = 0;
				for(int i=0;i<groups.length;i++){
					for(int j=0;j<groups[i].length;j++){
						suffixMap.put(groups[i][j], (i+1) + "_" + (j+1));
						sm ++;
					}
				}
				
				for(int i=0;i<sm;i++){
					header += prefix + suffixMap.get(i);
					if(i == sm - 1) break;
					header += ";";
				
				}
			}
			if(appendUnuniformity) header += "\t" + prefix + "_Diversity";
			
			return header;
		}
		
		
		public String toSimpleString(){
			/*String header = "Contig\tPosition\tStrand\tCodon\t" + AnnotatedGene.getSimpleHeader() + "\tGenomicRegion\tAnnotated\t";
			header += subGetHeader("Predicted", groups) + "\t";
			
			header += "Length\t";
			header += subGetHeader("HarrMatchScore", harrStartScores == null? null : new int[harrStartScores.length][1]) + "\t";
			header += subGetHeader("HarrStopCodonMatchScore", harrStartScores == null? null : new int[harrStartScores.length][1]) + "\t";
			header += subGetHeader("RPFMatchScore", groups) + "\t";
			header += subGetHeader("RPFStopCodonMatchScore", groups) + "\t";
			
			header += subGetHeader("RPFRPKM", groups) + "\t";
			header += subGetHeader("TE", groups) + "\t"; 
			
			header += subGetHeader("RPFFoldChange", groups) + "\t";
			header += subGetHeader("Release", groups) + "\t";			
			*/
			StringBuffer sb = new StringBuffer();
			sb.append(contig);sb.append('\t');
			sb.append(position.getPosition());sb.append('\t');
			sb.append(stopPosition);sb.append('\t');
			sb.append(isPlusStrand? '+' : '-');sb.append('\t');
			sb.append(seq.substring(0, 3));sb.append('\t');
			
			sb.append(gene == null ? AnnotatedGene.getEmptySimpleString() : gene.toSimpleString());
			sb.append('\t');
			sb.append(genomicRegion);sb.append('\t');
			sb.append(frameShift);sb.append('\t');
			//TODO predicted class make this as a function..
			sb.append(annotatedClass);sb.append('\t');
			sb.append(dndsRatioAfter);sb.append('\t');
			sb.append(phyloP);sb.append('\t');
			
			subToString(predictedClasses, sb, groups);
			sb.append('\t');
			
			sb.append(getLength()>=0? getLength() : "_");
			sb.append('\t');
			
			subToString(harrStartScores, sb, harrStartScores == null? 0 : harrStartScores.length);
			sb.append('\t');
			
			subToString(harrStopScores, sb, harrStartScores == null? 0 : harrStartScores.length);
			sb.append('\t');
			
			subToString(rpfStartScores, sb, groups);
			sb.append('\t');
			
			subToString(rpfStopScores, sb, groups);
			sb.append('\t');
			
			subToString(rpfPositionQuantities, sb, groups);
			sb.append('\t');
			
			
			subToString(positionTE, sb, groups);
			sb.append('\t');
			
			subToString(rpfPositionQuantityChanges, sb, groups);
			sb.append('\t');			
			
			subToString(releaseScore, sb, groups);
			
			return sb.toString();
		}
		
		@Override
		public String toString(){			
			StringBuffer sb = new StringBuffer();
			sb.append(contig);sb.append('\t');
			sb.append(position.getPosition());sb.append('\t');
			position.toMappedRegionString(sb);
			sb.append('\t');
			sb.append(isPlusStrand? '+' : '-');sb.append('\t');
			sb.append(seq);sb.append('\t');
			
			sb.append(gene == null ? AnnotatedGene.getEmptySimpleString() : gene.toSimpleString());
			sb.append('\t');
			sb.append(genomicRegion);sb.append('\t');
			//TODO predicted class make this as a function..
			sb.append(annotatedClass);sb.append('\t');
			
			subToString(predictedClasses, sb, groups);
			sb.append('\t');
			
			
			subToString(predictedClassesScores, sb, groups, false);
			//subToString(, sb, false);
			sb.append('\t');
			
			
			sb.append(frameShift);sb.append('\t');
			sb.append(dndsRatioBefore>=0?dndsRatioBefore : '_');sb.append('\t');
			sb.append(dndsRatioAfter>=0?dndsRatioAfter : '_');sb.append('\t');
			sb.append(phyloP);sb.append('\t');
			
			sb.append(signalPeptideScore>=0?signalPeptideScore : '_');sb.append('\t');
			
			
			
			//sb.append(gene == null ? '_' : (isAnnotated? 'T' : 'F'));sb.append('\t');
			
			sb.append(stopPosition>=0? stopPosition : "_");
			sb.append('\t');
			sb.append(seq.length());
			//sb.append(stopPosition>=0? Math.abs(position.getPosition() - stopPosition) : "_");
			sb.append('\t');
			
			
			subToString(harrStartScores, sb, harrStartScores == null? 0 : harrStartScores.length);
			sb.append('\t');
			
			subToString(harrStopScores, sb, harrStartScores == null? 0 : harrStartScores.length);
			sb.append('\t');
			
			
			subToString(rpfStartScores, sb, groups);
			sb.append('\t');
			
			subToString(rpfStopScores, sb, groups);
			sb.append('\t');
			
			subToString(rpfPositionQuantities, sb, groups);
			sb.append('\t');
			
			subToString(rnaPositionQuantities, sb, groups);
			sb.append('\t');
			
			subToString(positionTE, sb, groups);
			sb.append('\t');
			
			subToString(rpfPositionQuantityChanges, sb, groups);
			sb.append('\t');
			
			subToString(rpfCDSRPKMchanges, sb, groups);
			sb.append('\t');
						
			
			subToString(releaseScore, sb, groups);
			sb.append('\t');	
			
			subToString(rnaPositionQuantityChanges, sb, groups);
			sb.append('\t');
			
			subToString(rpfCDSQuantities, sb, groups);
			sb.append('\t');
			
			subToString(rnaCDSQuantities, sb, groups);
			sb.append('\t');
			
			subToString(cdsTE, sb, groups);
			sb.append('\t');
			
			if(cdsTE != null && positionTE != null)
				sb.append(String.format("%.3e", getKL(positionTE, cdsTE, groups)));
			else sb.append("_");
			sb.append('\t');
			
			sb.append(gene == null ? AnnotatedGene.getEmptyInfoString() : gene.toDetailedInfoString());
			return sb.toString();
		}
		
		
		public String GetArffSetHeader(){
			return "@relation w\n\n@attribute DnDsA numeric\n@attribute DnDsB numeric\n@attribute DnDsD numeric\n@attribute SigPep numeric\n@attribute HStart numeric\n@attribute HStop numeric\n@attribute RStart numeric"
			+ "\n@attribute RStop numeric\n@attribute RPFChange numeric\n@attribute RNAChange numeric"
			+ "\n@attribute RPFRPKM numeric\n@attribute RNARPKM numeric\n@attribute TE numeric\n@attribute RELEASE numeric\n@attribute LEN numeric"
			+ "\n@attribute Class {5U,T,TS,NC}\n\n@data";
			//return "DsDn,HStart,HStop,RStart, RStop, RPFChange,RNAChange,RPFRPKM,RNARPKM,TE,RELEASE,Class";
		}
		
		
		public String GetArffSetMfileHeader(){
			return "att={'DnDs After';'DnDs Before';'DnDs Diff'; 'SigPep Score'; 'Harr Match Score'; 'Harr Stop Codon Match Score'; 'RPF Match Score'; 'RPF Stop Codon Match Score'; 'RPF Fold Change'; 'RNA Fold Change'; 'RPF RPKM';"
					+ "'RNA RPKM'; 'Tranlation Efficiency'; 'Release'; 'Length'; 'Class'};\n value=[";
			//return "DsDn,HStart,HStop,RStart, RStop, RPFChange,RNAChange,RPFRPKM,RNARPKM,TE,RELEASE,Class";
		}
		
		public String ToTestSetString(int rpfrnaIndex, boolean harrRPFsynced){
			StringBuffer sb = new StringBuffer();
			
			sb.append(Math.log(dndsRatioAfter + .0001)/Math.log(2));sb.append(',');
			sb.append(Math.log(dndsRatioBefore + .0001)/Math.log(2));sb.append(',');
			sb.append(Math.log((dndsRatioBefore + .0001)/(dndsRatioAfter+ .0001))/Math.log(2));sb.append(',');
			sb.append(signalPeptideScore<0? '?' :  signalPeptideScore);sb.append(',');
			if(harrStartScores!=null){
				if(harrRPFsynced){					
					sb.append(harrStartScores[rpfrnaIndex] == 0.0? "?" : harrStartScores[rpfrnaIndex]);sb.append(',');			
				}else{
					double max = -10000;
					int maxIndex = 0;
					for(int i=0;i<harrStartScores.length;i++){
						if(max < harrStartScores[i]){
							maxIndex = i;
							max = harrStartScores[i];
						}
					}						
					sb.append(harrStartScores[maxIndex] == 0.0? "?" : harrStartScores[maxIndex]);sb.append(',');
				}
			}else sb.append("?,");
				
			if(harrStopScores!=null){
				if(harrRPFsynced){					
					sb.append(harrStopScores[rpfrnaIndex] == 0.0? "?" : harrStopScores[rpfrnaIndex]);sb.append(',');		
				}else{
					double max = -10000;
					int maxIndex = 0;
					for(int i=0;i<harrStopScores.length;i++){
						if(max < harrStopScores[i]){
							maxIndex = i;
							max = harrStopScores[i];
						}
					}						
					sb.append(harrStopScores[maxIndex] == 0.0? "?" : harrStopScores[maxIndex]);sb.append(',');
				}
			}else sb.append("?,");
			
	
			if(rpfStartScores!=null && rpfStartScores.length > rpfrnaIndex){
				sb.append(rpfStartScores[rpfrnaIndex] == 0.0 ? "?" : rpfStartScores[rpfrnaIndex]);sb.append(',');
			}else sb.append("?,");
			
			if(rpfStopScores!=null && rpfStopScores.length > rpfrnaIndex){
				sb.append(rpfStopScores[rpfrnaIndex] == 0.0 ? "?" : rpfStopScores[rpfrnaIndex]);sb.append(',');
			}else sb.append("?,");
			
			if(rpfPositionQuantityChanges!=null && rpfPositionQuantityChanges.length > rpfrnaIndex){
				sb.append(rpfPositionQuantityChanges[rpfrnaIndex]<0? "?" : Math.log(rpfPositionQuantityChanges[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
			}else sb.append("?,");
			
			if(rnaPositionQuantityChanges!=null && rnaPositionQuantityChanges.length > rpfrnaIndex){
				sb.append(rnaPositionQuantityChanges[rpfrnaIndex]<0? "?" : Math.log(rnaPositionQuantityChanges[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
			}else sb.append("?,");
			
			if(rpfPositionQuantities!=null && rpfPositionQuantities.length > rpfrnaIndex){
				sb.append(Math.log(rpfPositionQuantities[rpfrnaIndex]+0.0001)/Math.log(2));sb.append(',');
			}else sb.append("?,");
			
			if(rnaPositionQuantities!=null && rnaPositionQuantities.length > rpfrnaIndex){
				sb.append(Math.log(rnaPositionQuantities[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
			}else sb.append("?,");
			
			
			if(positionTE!=null && positionTE.length > rpfrnaIndex && positionTE[rpfrnaIndex]!= 0.0){
				sb.append(Math.log(positionTE[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
			}else sb.append("?,");
			
			//System.out.println(rpfrnaIndex + " " + rpfPositionQuantities[rpfrnaIndex] + " " + rnaPositionQuantities[rpfrnaIndex] + " " + positionTE[rpfrnaIndex]);
			//System.out.println("** " + Math.log(rpfPositionQuantities[rpfrnaIndex]+0.0001)/Math.log(2) + " " + Math.log(rnaPositionQuantities[rpfrnaIndex]+.0001)/Math.log(2) + " " + Math.log(positionTE[rpfrnaIndex]+0.0001)/Math.log(2));
				
			if(releaseScore!=null && releaseScore.length > rpfrnaIndex){
				sb.append(releaseScore[rpfrnaIndex]<0? "?" : Math.log(releaseScore[rpfrnaIndex]+0.0001)/Math.log(2));sb.append(',');
			}else sb.append("?,");
			
			if(getLength() >=0){
				sb.append(Math.log(getLength()+0.0001)/Math.log(2));sb.append(',');
			} else sb.append("?,");
			
			sb.append("?");
			
			return sb.toString();
		}		
		
		public String ToTrainingSetString(int rpfrnaIndex, boolean harrRPFsynced){
			
		//	if(rnaPositionQuantities!=null && rnaPositionQuantities.length <= index) return null;
		//	if(rpfPositionQuantities!=null && rpfPositionQuantities.length <= index) return null;
		//	if(rnaPositionQuantities!=null && rnaPositionQuantities.length > index && rnaPositionQuantities[index]<=rnaRPKMThreshold) return null;
		//	if(rpfPositionQuantities!=null && rpfPositionQuantities.length > index && rpfPositionQuantities[index]<=rpfRPKMThreshold) return null;
			//if(Math.abs(position - stopPosition) < 120) return null;
			int c = 0 ;
			
			if(isAnnotated && genomicRegion.equals("NM_ORF")){
				c = 1;
			}else if(genomicRegion.equals("NM_ORF") && frameShift.equals("0") && seq.startsWith("ATG")){
				c = 2;
			}else if(genomicRegion.equals("NM_3_UTR")){// && parser.getContainingGene(contig, isPlusStrand, stopPosition) == null){
				c = 3;
			}else if(genomicRegion.equals("NM_5_UTR") && seq.startsWith("ATG")){
				c = 4;
			}//else if(genomicRegion.endsWith("Intron") && (codon.equals("ATG") || codon.equals("CTG"))){
			//	c = 5;
			//}
			
			if(c>0){
				StringBuffer sb = new StringBuffer();
				sb.append(Math.log(dndsRatioAfter + .0001)/Math.log(2));sb.append(',');
				sb.append(Math.log(dndsRatioBefore + .0001)/Math.log(2));sb.append(',');
				sb.append(Math.log((dndsRatioBefore + .0001)/(dndsRatioAfter+ .0001))/Math.log(2));sb.append(',');
				sb.append(signalPeptideScore<0? '?' :  signalPeptideScore);sb.append(',');
				
				if(harrStartScores!=null){
					if(harrRPFsynced){
						//if(harrStartScores.length > rpfrnaIndex){
						sb.append(harrStartScores[rpfrnaIndex] == 0.0? "?" : harrStartScores[rpfrnaIndex]);sb.append(',');
						//}else return null;
					}else{
						double max = -10000;
						int maxIndex = 0;
						for(int i=0;i<harrStartScores.length;i++){
							if(max < harrStartScores[i]){
								maxIndex = i;
								max = harrStartScores[i];
							}
						}						
						sb.append(harrStartScores[maxIndex] == 0.0? "?" : harrStartScores[maxIndex]);sb.append(',');
					}
				}else sb.append("?,");
					
				if(harrStopScores!=null){
					if(harrRPFsynced){
						//if(harrStopScores.length > rpfrnaIndex){
						sb.append(harrStopScores[rpfrnaIndex] == 0.0? "?" : harrStopScores[rpfrnaIndex]);sb.append(',');
						//}else return null;
					}else{
						double max = -10000;
						int maxIndex = 0;
						for(int i=0;i<harrStopScores.length;i++){
							if(max < harrStopScores[i]){
								maxIndex = i;
								max = harrStopScores[i];
							}
						}						
						sb.append(harrStopScores[maxIndex] == 0.0? "?" : harrStopScores[maxIndex]);sb.append(',');
					}
				}else sb.append("?,");
				
		
				if(rpfStartScores!=null && rpfStartScores.length > rpfrnaIndex){
					sb.append(rpfStartScores[rpfrnaIndex] == 0.0 ? "?" : rpfStartScores[rpfrnaIndex]);sb.append(',');
				}else sb.append("?,");
				
				if(rpfStopScores!=null && rpfStopScores.length > rpfrnaIndex){
					sb.append(rpfStopScores[rpfrnaIndex] == 0.0 ? "?" : rpfStopScores[rpfrnaIndex]);sb.append(',');
				}else sb.append("?,");
				
				if(rpfPositionQuantityChanges!=null && rpfPositionQuantityChanges.length > rpfrnaIndex){
					sb.append(rpfPositionQuantityChanges[rpfrnaIndex]<0? "?" : Math.log(rpfPositionQuantityChanges[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
				}else sb.append("?,");
				
				if(rnaPositionQuantityChanges!=null && rnaPositionQuantityChanges.length > rpfrnaIndex){
					sb.append(rnaPositionQuantityChanges[rpfrnaIndex]<0? "?" : Math.log(rnaPositionQuantityChanges[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
				}else sb.append("?,");
				
				if(rpfPositionQuantities!=null && rpfPositionQuantities.length > rpfrnaIndex){
					sb.append(Math.log(rpfPositionQuantities[rpfrnaIndex]+0.0001)/Math.log(2));sb.append(',');
				}else sb.append("?,");
				
				if(rnaPositionQuantities!=null && rnaPositionQuantities.length > rpfrnaIndex){
					sb.append(Math.log(rnaPositionQuantities[rpfrnaIndex]+0.0001)/Math.log(2));sb.append(',');
				}else sb.append("?,");
				
				if(positionTE!=null && positionTE.length > rpfrnaIndex && positionTE[rpfrnaIndex]!= 0.0){
					sb.append(Math.log(positionTE[rpfrnaIndex]+.0001)/Math.log(2));sb.append(',');
				}else sb.append("?,");
					
				if(releaseScore!=null && releaseScore.length > rpfrnaIndex){
					sb.append(releaseScore[rpfrnaIndex]<0? "?" : Math.log(releaseScore[rpfrnaIndex]+0.0001)/Math.log(2));sb.append(',');
				}else sb.append("?,");
				
				if(getLength() >=0){
					sb.append(Math.log(getLength()+0.0001)/Math.log(2));sb.append(',');
				} else sb.append("?,");
				
				if(c==1){
					sb.append("TS");
				}else if(c==2){
					sb.append("T");
				}else if(c==3){
					sb.append("NC");
				//	if(this.getRpfStartScores()[0] > 0.87) System.out.println(this.getPosition());
				}else if(c==4){
					sb.append("5U");
				}//else{
				//	sb.append("IN");
				//}
				return sb.toString();
			}
			
			return null;
			
			// dsdn, stop, te, positionchange, harrscore, rpfscore // per sample. or for sample 0 (or for samples with those properties)
			
			// translation start (T) :  for sample 0, take annotated points (already rpkm is large)
			// non translating (N) : for sample, take intergenic region (ATG or CTG)
			// inner translation start (I) : for sample 0, take unannotated but NM_ORF points (non ATG or CTG) 
		}

		private void subToString(double[] quantities, StringBuffer sb, int[][] groups){
			subToString(quantities, sb, groups, true);
		}

		private void subToString(String[] predictions, StringBuffer sb, int[][] groups){
			double[] quantities = new double[predictions.length];
			for(int i=0;i<quantities.length;i++){
				if(predictions[i].toUpperCase().equals("TS")){
					quantities[i] = 4;
				}else if(predictions[i].toUpperCase().equals("T")){
					quantities[i] = 3;
				}else if(predictions[i].toUpperCase().equals("5U")){
					quantities[i] = 2;
				}else quantities[i] = 1;
			}
			for(int i=0; i<predictions.length-1; i++){
				sb.append(predictions[i]); sb.append(';');
			}
			sb.append(predictions[predictions.length-1]); 
			
			sb.append('\t');
			if(quantities.length > 1){
				sb.append(String.format("%.3e", getNonuniformity(quantities, groups)));
			}else sb.append("_");
		}
		
		private void subToString(double[] quantities, StringBuffer sb, int numSample){
			subToString(quantities, sb, new int[numSample][1], true);
		}
		
		private void subToString(double[] quantities, StringBuffer sb, int[][] groups, boolean appendUnuniformity){ 
			if(quantities!=null && quantities.length>0){
				for(int i=0; i<quantities.length-1; i++){
					sb.append(String.format("%.3e", quantities[i])); sb.append(';');
				}
				sb.append(String.format("%.3e", quantities[quantities.length-1])); 
				if(appendUnuniformity){
				sb.append('\t');
					if(quantities.length > 1){
						sb.append(String.format("%.3e", getNonuniformity(quantities, groups)));
					}else sb.append("_");
				}
			}else{
				//sb.append("_");
				
				for(int i=0; i<groups.length; i++){
					for(int j=0;j<groups[i].length;j++){
						if(i == groups.length-1 && j == groups[i].length-1) break;
						sb.append("_;");
					}
				}
				sb.append("_");
				if(appendUnuniformity) sb.append("\t_");
			}
		}
		
		
		// check if non negative!!
		private double getNonuniformity(double[] dist, int[][] groups){
			
			double[] gdist = getGroupDistribution(dist, groups, 1.0);
			
			double kl = 0;
			for(double v : gdist){
				if(v > 0) kl += v * Math.log(v*gdist.length);
			}
			return kl;		
		}
		
		private double[] getGroupDistribution(double[] dist, int[][] groups, double offset){
			double[] gdist = new double[groups.length];
			for(int i=0; i<gdist.length;i++){
				for(int j : groups[i]){
					gdist[i] += dist[j];
				}
				gdist[i] += offset;
			}
			
			double sum = 0;
			
			for(double v : gdist){
				sum += v;
			}
			
			if(sum != 0){
				for(int i=0;i<gdist.length;i++){
					gdist[i] /= sum;
				}
			}
			return gdist;
		}
		
		private double getKL(double[] dist1, double[] dist2, int[][] groups){// s
			double sum = 0;
			
			double[] dist1cp = new double[dist1.length];
			double[] dist2cp = new double[dist2.length];
			
			for(double v : dist1){
				sum += v;
			}
			for(int i=0;i<dist1.length;i++){
				dist1cp[i] = dist1[i] / sum;
			}
			sum = 0;
			double min = 1;
			for(int i=0;i<dist2.length;i++){
				if(dist2[i] == 0) continue;
				min = min > dist2[i]? dist2[i] : min; 
			}
			for(int i=0;i<dist2.length;i++){
				if(dist2[i] == 0) sum += min *.1;
				else sum += dist2[i];
			}
		
			for(int i=0;i<dist2.length;i++){
				dist2cp[i] = dist2[i] / sum;
			}
			double kl = 0;
			for(int i=0;i<dist1cp.length;i++){
				double v= dist1cp[i];
				if(v > 0) kl += v * Math.log(v/(dist2cp[i] == 0? min * .1 : dist2cp[i]));
			}
			return kl;					
		}
		/*
		private double getNonuniformity(double[] dist){
			double lambda = 0;
			for(double v : dist){
				lambda += v;
			}
			lambda /= dist.length; 
			
			double ret = 0;
			for(double v : dist){
				double p = Math.exp(-lambda)
						
				if(v < lambda){
					
				}else{
					
				}
				
			}
			
		}
		*/
		
		
	}
	
	public MergedFileParser(String file){
		list = new ArrayList<MergedFileParser.MergedResult>();
		int[][] groups = null;
		try {
			BufferedLineReader in = new BufferedLineReader(file);
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("@")){
					groups = MergedResult.getGroups(s);
					continue;
				}else if(s.startsWith("Contig")){
					continue;
				}
				list.add(new MergedResult(s, groups));
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void write(ArrayList<MergedResult> list, String file){
		try {
			PrintStream out = new PrintStream(file);
			for(int i=0;i<list.size();i++){
				if(i == 0){
					out.println(list.get(i).getHeader());
				}
				out.println(list.get(i));
			}
			
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private static double getAvg(double[] v){
		double sum = 0;
		for(double t : v){
			sum += t;
		}
		return sum / v.length;
	}
	
	public static void main(String[] args) throws IOException{
		MergedFileParser test = new MergedFileParser("/media/kyowon/Data1/RPF_Project/samples/sample4/results/out_ssh_0.3.csv");
	//	ZeroBasedFastaParser fasta = new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		ArrayList<MergedResult> list = test.getList();
		HashMap<String, double[]> terpkmMap = new HashMap<String, double[]>();
		for(MergedResult mr : list){
			if( mr.getGene() == null) continue;
			String gn = mr.getGene().getGeneName();
			double[] te = mr.getCdsTE();
			double[] rpkm = mr.getRpfCDSQuantities();
			if(te == null || rpkm == null) continue;
			double[] v = new double[2];
			v[0] = getAvg(te);
			v[1] = getAvg(rpkm);
			terpkmMap.put(gn, v);
		}
		
		BufferedLineReader in = new  BufferedLineReader("/home/kyowon/protein.csv");
		String s;
		PrintStream out = new PrintStream("/home/kyowon/protein.out.csv");
		HashSet<String> gns = new HashSet<String>(); 
		while((s=in.readLine())!=null){
			if(!s.startsWith(">")){
				out.println("Gene\tTE\tRPKM\tSpecInt");
				continue;
			}
			String[] token = s.split("\t")[0].split(";");
			
			for(String t : token){
				if(t.indexOf("GN=")<0 || t.indexOf(" PE=")<0){					
					continue;
				}
				String gn = t.substring(t.indexOf("GN=") + 3, t.indexOf(" PE="));
				if(gns.contains(gn)) continue;
				gns.add(gn);
				double[] v = terpkmMap.get(gn);
				if(v == null) continue;
				out.println(gn+"\t"+v[0]+"\t"+v[1] + "\t" + s.split("\t")[1].replace(".E", "e"));				
			}
			
		}
		out.close();
		in.close();
	}
	
}
