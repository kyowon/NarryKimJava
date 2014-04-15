package parser;

import java.util.ArrayList;
import java.util.HashMap;

import jspp.SignalPeptidePredictor;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;
import rpf.Quantifier;
import rpf.Scorer;
import util.Codon;
import util.DnDsCalculator;

public class MergedFileParser {
	public class MergedResult{
		//	IsAnnotated	ContainingGeneName	ContainingGBGeneName	txStart	txEnd	cdsStart	cdsEnd	genomicRegion	frameShift
		//private int positionQuantityChangeLength;
		//private int positionQuantityOffset;
		//private int maxLengthUntilStopcodon;
		
		private String contig;
		private int position;
		private boolean isPlusStrand;
		private String codon;
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
		private double signalPeptideScore = -1;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		//private int numSample = 0;
		private int[][] groups;
		
		public String getContig() {
			return contig;
		}

		public int getPosition() {
			return position;
		}

		public boolean isPlusStrand() {
			return isPlusStrand;
		}

		public String getCodon() {
			return codon;
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

		public double getDnDsRatioAfter() {
			return dndsRatioAfter;
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
		
		public boolean isPredictionCorrect(){
			boolean correct = true;
			for(String prediction : predictedClasses){
				if(prediction.toUpperCase().equals(annotatedClass)) continue;
				correct = false;
				break;
			}
			return correct;
		}
		
		public MergedResult(ScoredPosition scoredPosition, Scorer[] harrScorers, Scorer[] rpfScorers, Quantifier[] rpfQuantifiers, Quantifier[] rnaQuantifiers, int[][] groups, ZeroBasedFastaParser fastaFileParser, MafParser mafParser,
				SignalPeptidePredictor pd, int positionQuantityChangeLength, int positionQuantityOffset, int maxLengthUntilStopcodon){
			this.contig = scoredPosition.getContig();
			this.position = scoredPosition.getPosition();
			this.isPlusStrand = scoredPosition.isPlusStrand();
			this.codon = scoredPosition.getCodon();
			this.genomicRegion = scoredPosition.getGenomicRegion();
			this.frameShift = scoredPosition.getFrameShift();
			this.isAnnotated = scoredPosition.isAnnotated();
			this.gene = scoredPosition.getGene();
			this.groups = groups;
			
			if(harrScorers != null){
				this.harrStartScores = new double[harrScorers.length];
					
				for(int i=0;i<this.harrStartScores.length;i++){
					this.harrStartScores[i] = harrScorers[i].getStartScore(this.contig, this.position, this.isPlusStrand, false);
				}
				
				if(gene!=null && !this.genomicRegion.endsWith("Intron")){	
					this.harrStopScores = new double[harrScorers.length];
					for(int i=0;i<this.harrStopScores.length;i++){
						this.harrStopScores[i] = harrScorers[i].getStopScore(this.contig, this.position, this.isPlusStrand, false, maxLengthUntilStopcodon);
					
					}
				}
			}
			//if(rpfScorers !=null){
			//this.numSample = rpfScorers.length;
			this.rpfStartScores = new double[rpfScorers.length];
			this.predictedClasses = new String[rpfScorers.length];
			this.predictedClassesScores = new double[rpfScorers.length];
			for(int i=0;i<this.rpfStartScores.length;i++){
				this.rpfStartScores[i] = rpfScorers[i].getStartScore(this.contig, this.position, this.isPlusStrand, false);
			}
			if(gene!=null && !this.genomicRegion.endsWith("Intron")){
				this.rpfStopScores = new double[rpfScorers.length];				
				for(int i=0;i<this.rpfStopScores.length;i++){
					this.rpfStopScores[i] = rpfScorers[i].getStopScore(this.contig, this.position, this.isPlusStrand, false, maxLengthUntilStopcodon);
				}
			}
			//}
			
			
			String sequence = null;
			if(!fastaFileParser.containsContig(contig)){
				sequence = "N/A";
			}else if(isPlusStrand){
				sequence = fastaFileParser.getSequence(contig, position, position + 70*3);
			}else{							
				sequence = ZeroBasedFastaParser.getComplementarySequence(fastaFileParser.getSequence(contig, position - 70*3 + 1, position + 1), true);
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
				if(gene != null){				
					rpfCDSQuantities = new double[rpfQuantifiers.length];
					for(int i=0; i<rpfQuantifiers.length;i++){
						rpfCDSQuantities[i] = rpfQuantifiers[i].getCDSRPKM(gene);
					}
				}
				rpfPositionQuantities = new double[rpfQuantifiers.length];
				rpfPositionQuantityChanges = new double[rpfQuantifiers.length];
				if(gene!=null && !this.genomicRegion.endsWith("Intron")) releaseScore = new double[rpfQuantifiers.length];
				if(gene!=null && this.genomicRegion.equals("NM_ORF") && !this.isAnnotated) rpfCDSRPKMchanges = new double[rpfQuantifiers.length];
					
				for(int i=0; i<rpfPositionQuantities.length;i++){
					rpfPositionQuantities[i] = rpfQuantifiers[i].getPositionRPKM(contig, this.position, isPlusStrand, maxLengthUntilStopcodon, positionQuantityOffset);
					rpfPositionQuantityChanges[i] = rpfQuantifiers[i].getPositionQuantityChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength, positionQuantityOffset);
					//rpfPositionQuantityChanges[i] = Math.log(rpfPositionQuantityChanges[i]+.0001)/Math.log(2);
					if(gene!=null && !this.genomicRegion.endsWith("Intron")){
						ArrayList<Double> qs = rpfQuantifiers[i].getNextStopCodonQuantityChangeRatioNStopPosition(contig, this.position, isPlusStrand, positionQuantityChangeLength, maxLengthUntilStopcodon);
						if(qs!=null && !this.genomicRegion.endsWith("Intron")){
							releaseScore[i] = qs.get(0);
							//rpfAroundStopCodonQuantityChanges[i] = Math.log(rpfAroundStopCodonQuantityChanges[i]+.0001)/Math.log(2);
							stopPosition = (int)((double)qs.get(1));
						}
					}//else rpfAroundStopCodonQuantityChanges = null;
					if(rpfCDSRPKMchanges!=null){
						rpfCDSRPKMchanges[i] = rpfQuantifiers[i].getCDSRPKMChangeRatio(gene, this.position, isPlusStrand, positionQuantityOffset);
					}
					
				}
				
				if(mafParser!=null){ 
					//this.dsdnRatioAfter = DsDnCalculator.calculate(mafParser.getSeqs(contig, position + (isPlusStrand? -30:30), isPlusStrand, Math.min(60, 30 + Math.abs(position - stopPosition))));
					int offset = 0;//cpos*3;
					String[] after = mafParser.getSeqs(contig, position + (isPlusStrand? offset : -offset), isPlusStrand , positionQuantityChangeLength);
					String[] before = mafParser.getSeqs(contig, position + (isPlusStrand? -positionQuantityChangeLength:positionQuantityChangeLength), isPlusStrand , positionQuantityChangeLength);
					this.dndsRatioAfter = DnDsCalculator.calculate(after); 
					this.dndsRatioBefore = DnDsCalculator.calculate(before);
					//this.dsdnRatio = Math.log(dsdnRatio + .0001)/Math.log(2);
				}
				//rpfAroundStopCodonQuantityChanges
			}
			
			
			if(rnaQuantifiers!=null){
				if(gene != null && !genomicRegion.startsWith("NR")){				
					rnaCDSQuantities = new double[rnaQuantifiers.length];
					for(int i=0; i<rnaQuantifiers.length;i++){
						rnaCDSQuantities[i] = rnaQuantifiers[i].getCDSRPKM(gene);
					}
				}
				rnaPositionQuantities = new double[rnaQuantifiers.length];
				rnaPositionQuantityChanges = new double[rnaQuantifiers.length];
				for(int i=0; i<rnaPositionQuantities.length;i++){
					rnaPositionQuantities[i] = rnaQuantifiers[i].getPositionRPKM(contig, this.position, isPlusStrand, maxLengthUntilStopcodon, 0);
					rnaPositionQuantityChanges[i] = rnaQuantifiers[i].getPositionQuantityChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength, 0);
					//rnaPositionQuantityChanges[i] = Math.log(rnaPositionQuantityChanges[i]+.0001)/Math.log(2);
				}
			}
			if(rpfQuantifiers!=null && rnaQuantifiers!=null){
				if(gene != null && !genomicRegion.startsWith("NR") && !this.genomicRegion.endsWith("Intron")){
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
		
		/*public MergedResult(String s){ //TODO fix next time- consider "_"
			int i = 0;
			String[] token = s.split("\t");
			contig = token[i++];
			position = Integer.parseInt(token[i++]);
			isPlusStrand = token[i++].equals("+");
			codon = token[i++];
			genomicRegion = token[i++];
			
			//TODO predicted, score
			frameShift = token[i++];
			dsdnRatioAfter = token[i].equals('_') ? -1 : Double.parseDouble(token[i]);
			i++;
			if(!token[i].equals("_"))
				isAnnotated = token[i].equals("T");
			i++;
			
			stopPosition = Integer.parseInt(token[i++]); //TODO
			i++;
			
			harrStartScores = subParse(token[i++]);
			harrStopScores = subParse(token[i++]);
			rpfStartScores = subParse(token[i++]);
			rpfStopScores = subParse(token[i++]);
			rpfPositionQuantities = subParse(token[i++]);
			rnaPositionQuantities = subParse(token[i++]);
			positionTE = subParse(token[i++]);
			rpfPositionQuantityChanges = subParse(token[i++]);
			releaseScore = subParse(token[i++]);
			
			rnaPositionQuantityChanges = subParse(token[i++]);
			rpfCDSQuantities = subParse(token[i++]);
			rnaCDSQuantities = subParse(token[i++]);
			cdsTE = subParse(token[i++]);
			//i++;
			
			if(!token[i].equals("_")){
				StringBuffer gs = new StringBuffer();
				for(;i<token.length;i++){
					gs.append(token[i]); gs.append('\t');
				}
				gene = new AnnotatedGene(gs.toString());
			}
		}*/
		
		public void setPredictedClasses(String predictedClass, double predictedClassesScore, int index){
			this.predictedClasses[index] = predictedClass; 
			this.predictedClassesScores[index] = predictedClassesScore;
		}
		
		public int getStopPostion() { return stopPosition; }
		
		private double[] subParse(String s){
			String[] token = s.split(";");
			double[] ret = new double[token.length];
			for(int i=0;i<ret.length;i++){
				ret[i] = Double.parseDouble(token[i]);
			}
			return ret;
		}
		
		public String getHeader(){
			String header =  "Contig\tPosition\tStrand\tCodon\t" + AnnotatedGene.getSimpleHeader() + "\tGenomicRegion\tAnnotated\t";
			header += subGetHeader("Predicted", groups) + "\t";
			header += subGetHeader("PredictionScore", groups, false) + "\t";
					//"Predicted\tPredictionScore\t";
			header += "Frame\tDnDsBefore\tDnDsAfter\tSignalPeptideScore\tStopCodonPosition\tLength\t";
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
		
		
		@Override
		public String toString(){			
			StringBuffer sb = new StringBuffer();
			sb.append(contig);sb.append('\t');
			sb.append(position);sb.append('\t');
			sb.append(isPlusStrand? '+' : '-');sb.append('\t');
			sb.append(codon);sb.append('\t');
			
			sb.append(gene == null ? AnnotatedGene.getEmptySimpleString() : gene.toSimpleString());
			sb.append('\t');
			sb.append(genomicRegion);sb.append('\t');
			//TODO predicted class make this as a function..
			sb.append(annotatedClass);sb.append('\t');
			
			/*boolean isConsistent = true;
			for(int i=0; i<predictedClasses.length-1;i++){
				sb.append(predictedClasses[i]);
				sb.append(";");
				if(!predictedClasses[i].toUpperCase().equals(predictedClasses[i+1].toUpperCase())) isConsistent = false;
			}
			sb.append(predictedClasses[predictedClasses.length-1] + "\t");
			sb.append(isConsistent);
			
			sb.append('\t');
			*/
			
			subToString(predictedClasses, sb, groups);
			sb.append('\t');
			
			
			subToString(predictedClassesScores, sb, groups, false);
			//subToString(, sb, false);
			sb.append('\t');
			
			
			sb.append(frameShift);sb.append('\t');
			sb.append(dndsRatioBefore>=0?dndsRatioBefore : '_');sb.append('\t');
			sb.append(dndsRatioAfter>=0?dndsRatioAfter : '_');sb.append('\t');
			sb.append(signalPeptideScore>=0?signalPeptideScore : '_');sb.append('\t');
			
			
			
			//sb.append(gene == null ? '_' : (isAnnotated? 'T' : 'F'));sb.append('\t');
			
			sb.append(stopPosition>=0? stopPosition : "_");
			sb.append('\t');
			
			sb.append(stopPosition>=0? Math.abs(position - stopPosition) : "_");
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
			
			if(stopPosition >=0){
				sb.append(Math.log(Math.abs(position-stopPosition)+0.0001)/Math.log(2));sb.append(',');
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
			}else if(genomicRegion.equals("NM_ORF") && frameShift.equals("0") && codon.equals("ATG")){
				c = 2;
			}else if(genomicRegion.equals("NM_3_UTR")){// && parser.getContainingGene(contig, isPlusStrand, stopPosition) == null){
				c = 3;
			}else if(genomicRegion.equals("NM_5_UTR") && codon.equals("ATG")){
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
				
				if(stopPosition >=0){
					sb.append(Math.log(Math.abs(position-stopPosition)+0.0001)/Math.log(2));sb.append(',');
				} else sb.append("?,");
				
				if(c==1){
					sb.append("TS");
				}else if(c==2){
					sb.append("T");
				}else if(c==3){
					sb.append("NC");
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
				if(predictions[i].equals("TS")){
					quantities[i] = 3;
				}else if(predictions[i].equals("T")){
					quantities[i] = 2;
				}else if(predictions[i].equals("5U")){
					quantities[i] = 1;
				}				
			}
			
			subToString(quantities, sb, groups, true);
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
}
