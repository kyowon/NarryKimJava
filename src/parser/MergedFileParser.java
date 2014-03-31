package parser;

import java.util.ArrayList;

import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;
import rpf.Quantifier;
import rpf.Scorer;
import util.DsDnCalculator;

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
		private double[] releaseScore;//
		private int stopPosition = -1;
		private double[] rnaCDSQuantities;//
		private double[] rnaPositionQuantities;//
		private double[] rnaPositionQuantityChanges;//
		private String genomicRegion;
		private String frameShift;
		private double dsdnRatioAfter = -1;
		private double dsdnRatioBefore = -1;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		
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

		public double getDsdnRatioAfter() {
			return dsdnRatioAfter;
		}

		public double getDsdnRatioBefore() {
			return dsdnRatioBefore;
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
		
		public MergedResult(ScoredPosition scoredPosition, Scorer[] harrScorers, Scorer[] rpfScorers, Quantifier[] rpfQuantifiers, Quantifier[] rnaQuantifiers, MafParser mafParser,
				int positionQuantityChangeLength, int positionQuantityOffset, int maxLengthUntilStopcodon){
			this.contig = scoredPosition.getContig();
			this.position = scoredPosition.getPosition();
			this.isPlusStrand = scoredPosition.isPlusStrand();
			this.codon = scoredPosition.getCodon();
			this.genomicRegion = scoredPosition.getGenomicRegion();
			this.frameShift = scoredPosition.getFrameShift();
			this.isAnnotated = scoredPosition.isAnnotated();
			this.gene = scoredPosition.getGene();
			
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
				}
				
				if(mafParser!=null){ 
					//this.dsdnRatioAfter = DsDnCalculator.calculate(mafParser.getSeqs(contig, position + (isPlusStrand? -30:30), isPlusStrand, Math.min(60, 30 + Math.abs(position - stopPosition))));
					String[] after = mafParser.getSeqs(contig, position, isPlusStrand , positionQuantityChangeLength);
					String[] before = mafParser.getSeqs(contig, position + (isPlusStrand? -positionQuantityChangeLength:positionQuantityChangeLength), isPlusStrand , positionQuantityChangeLength);
					this.dsdnRatioAfter = DsDnCalculator.calculate(after, positionQuantityChangeLength);
					this.dsdnRatioBefore = DsDnCalculator.calculate(before, positionQuantityChangeLength);
					//this.dsdnRatio = Math.log(dsdnRatio + .0001)/Math.log(2);
				}
				//rpfAroundStopCodonQuantityChanges
			}
			
			if(rnaQuantifiers!=null){
				if(gene != null){				
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
				if(gene != null && !this.genomicRegion.endsWith("Intron")){
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
			String header =  "Contig\tPosition\tStrand\tCodon\tGenomicRegion\tAnnotated\t";
			header += subGetHeader("Predicted", predictedClasses) + "\t";
			header += subGetHeader("PredictionScore", predictedClassesScores, false) + "\t";
					//"Predicted\tPredictionScore\t";
			header += "Frame\tDnDsratio\tStopCodonPosition\tLength\t";
			//header += "";
			header += subGetHeader("HarrMatchScore", harrStartScores) + "\t";
			header += subGetHeader("HarrStopCodonMatchScore", harrStopScores) + "\t";
			header += subGetHeader("RPFMatchScore", rpfStartScores) + "\t";
			header += subGetHeader("RPFStopCodonMatchScore", rpfStopScores) + "\t";
			
			header += subGetHeader("RPFRPKM", rpfPositionQuantities) + "\t";
			header += subGetHeader("RNARPKM", rnaPositionQuantities) + "\t";
			header += subGetHeader("TE", positionTE) + "\t"; 
			
			header += subGetHeader("RPFFoldChange", rpfPositionQuantityChanges) + "\t";
			header += subGetHeader("Release", releaseScore) + "\t";			
		
			header += subGetHeader("RNAFoldChange", rnaPositionQuantityChanges) + "\t";
			
			header += subGetHeader("RPFCDSRPKMs", rpfCDSQuantities) + "\t";
			header += subGetHeader("RNACDSRPKMs", rnaCDSQuantities) + "\t";
			header += subGetHeader("CDSTE", cdsTE) + "\t";
			
			header += "TEDiff(Pos vs. CDS)\t";
			
			header += AnnotatedGene.getHeader();
			
			return header;
		}
		
		private String subGetHeader(String prefix, double[] quantities){
			return subGetHeader(prefix, quantities, true);
		}
		
		private String subGetHeader(String prefix, double[] quantities, boolean appendUnuniformity){
			String header = "";
			
			if(quantities != null){
				for(int i=0;i<quantities.length-1;i++){
					header += prefix + (i+1) + ";";
				}
				header += prefix + quantities.length;
			}else{
				header += prefix;
			}
			if(appendUnuniformity) header += "\t" + prefix + "_Ununiformity";
			
			return header;
		}
		
		private String subGetHeader(String prefix, String[] quantities){
			String header = "";
			
			if(quantities != null){
				for(int i=0;i<quantities.length-1;i++){
					header += prefix + (i+1) + ";";
				}
				header += prefix + quantities.length;
			}else{
				header += prefix;
			}
		//	header += prefix + "Ununiformity";
			
			return header;
		}
		
		@Override
		public String toString(){			
			StringBuffer sb = new StringBuffer();
			sb.append(contig);sb.append('\t');
			sb.append(position);sb.append('\t');
			sb.append(isPlusStrand? '+' : '-');sb.append('\t');
			sb.append(codon);sb.append('\t');
			
			sb.append(genomicRegion);sb.append('\t');
			//TODO predicted class make this as a function..
			sb.append(annotatedClass);sb.append('\t');
			
			if(predictedClasses == null){
				sb.append("_");
			}else{
				for(int i=0; i<predictedClasses.length-1;i++){
					sb.append(predictedClasses[i]);
					sb.append(";");
				}
				sb.append(predictedClasses[predictedClasses.length-1]);
			}
			sb.append('\t');
			
			subToString(predictedClassesScores, sb, false);
			sb.append('\t');
			
			
			sb.append(frameShift);sb.append('\t');
			sb.append(dsdnRatioAfter>=0?dsdnRatioAfter : '_');sb.append('\t');
			//sb.append(gene == null ? '_' : (isAnnotated? 'T' : 'F'));sb.append('\t');
			
			sb.append(stopPosition>=0? stopPosition : "_");
			sb.append('\t');
			
			sb.append(stopPosition>=0? Math.abs(position - stopPosition) : "_");
			sb.append('\t');
			
			
			subToString(harrStartScores, sb);
			sb.append('\t');
			
			subToString(harrStopScores, sb, harrStartScores == null? 0 : harrStartScores.length);
			sb.append('\t');
			
			
			subToString(rpfStartScores, sb);
			sb.append('\t');
			
			subToString(rpfStopScores, sb, rpfStartScores == null? 0 : rpfStartScores.length);
			sb.append('\t');
			
			subToString(rpfPositionQuantities, sb);
			sb.append('\t');
			
			subToString(rnaPositionQuantities, sb);
			sb.append('\t');
			
			subToString(positionTE, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantityChanges, sb);
			sb.append('\t');
			
			subToString(releaseScore, sb, rpfPositionQuantities == null? 0 : rpfPositionQuantities.length);
			sb.append('\t');	
			
			subToString(rnaPositionQuantityChanges, sb);
			sb.append('\t');
			
			subToString(rpfCDSQuantities, sb, rpfPositionQuantities == null? 0 : rpfPositionQuantities.length);
			sb.append('\t');
			
			subToString(rnaCDSQuantities, sb, rnaPositionQuantities == null? 0 : rnaPositionQuantities.length);
			sb.append('\t');
			
			subToString(cdsTE, sb, cdsTE == null? 0 : cdsTE.length);
			sb.append('\t');
			
			if(cdsTE != null && positionTE != null)
				sb.append(String.format("%.3e", getKL(positionTE, cdsTE)));
			else sb.append("_");
			sb.append('\t');
			
			sb.append(gene == null ? AnnotatedGene.getEmptyGeneString() : gene.toString());
			return sb.toString();
		}
		
		
		public String GetArffSetHeader(){
			return "@relation w\n\n@attribute DsDnAfter numeric\n@attribute DsDnBefore numeric\n@attribute DsDnDiff numeric\n@attribute HStart numeric\n@attribute HStop numeric\n@attribute RStart numeric"
			+ "\n@attribute RStop numeric\n@attribute RPFChange numeric\n@attribute RNAChange numeric"
			+ "\n@attribute RPFRPKM numeric\n@attribute RNARPKM numeric\n@attribute TE numeric\n@attribute RELEASE numeric\n@attribute LEN numeric"
			+ "\n@attribute Class {5U,T,TS,NC}\n\n@data";
			//return "DsDn,HStart,HStop,RStart, RStop, RPFChange,RNAChange,RPFRPKM,RNARPKM,TE,RELEASE,Class";
		}
		
		
		public String GetArffSetMfileHeader(){
			return "att={'DnDs After';'DnDs Before';'DnDs Change'; 'Harr Match Score'; 'Harr Stop Codon Match Score'; 'RPF Match Score'; 'RPF Stop Codon Match Score'; 'RPF Fold Change'; 'RNA Fold Change'; 'RPF RPKM';"
					+ "'RNA RPKM'; 'Tranlation Efficiency'; 'Release Score'; 'Length'; 'Class'};\n value=[";
			//return "DsDn,HStart,HStop,RStart, RStop, RPFChange,RNAChange,RPFRPKM,RNARPKM,TE,RELEASE,Class";
		}
		
		public String ToTestSetString(int rpfrnaIndex, boolean harrRPFsynced){
			StringBuffer sb = new StringBuffer();
			
			sb.append(Math.log(dsdnRatioAfter + .0001)/Math.log(2));sb.append(',');
			sb.append(Math.log(dsdnRatioBefore + .0001)/Math.log(2));sb.append(',');
			sb.append(Math.log((dsdnRatioBefore + .0001)/(dsdnRatioAfter+ .0001))/Math.log(2));sb.append(',');

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
				sb.append(Math.log(dsdnRatioAfter + .0001)/Math.log(2));sb.append(',');
				sb.append(Math.log(dsdnRatioBefore + .0001)/Math.log(2));sb.append(',');
				sb.append(Math.log((dsdnRatioBefore + .0001)/(dsdnRatioAfter+ .0001))/Math.log(2));sb.append(',');

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
		
		private void subToString(double[] quantities, StringBuffer sb){
			subToString(quantities, sb, true);
		}
		
		private void subToString(double[] quantities, StringBuffer sb, boolean appendUnuniformity){
			subToString(quantities, sb, quantities == null? 0 : quantities.length, appendUnuniformity);
		}
		
		private void subToString(double[] quantities, StringBuffer sb, int numSample){
			subToString(quantities, sb, numSample, true);
		}
		
		private void subToString(double[] quantities, StringBuffer sb, int numSample, boolean appendUnuniformity){ 
			if(quantities!=null && quantities.length>0){
				for(int i=0; i<quantities.length-1; i++){
					sb.append(String.format("%.3e", quantities[i])); sb.append(';');
				}
				sb.append(String.format("%.3e", quantities[quantities.length-1])); 
				if(appendUnuniformity){
				sb.append('\t');
					if(quantities.length > 1){
						sb.append(String.format("%.3e", getNonuniformity(quantities)));
					}else sb.append("_");
				}
			}else{
				//sb.append("_");
				for(int i=0; i<numSample-1; i++){
					sb.append("_;");
				}
				sb.append("_");
				if(appendUnuniformity) sb.append("\t_");
			}
		}
		
		private double getNonuniformity(double[] dist){
			double sum = 0;
			for(double v : dist){
				sum += v + 1;
			}
			if(sum == 0) return 0;
			double kl = 0;
			for(double v : dist){
				double vn = (v + 1)/sum;
				if(vn > 0) kl += vn * Math.log(vn*dist.length);
			}
			return kl;		
		}
		
		private double getKL(double[] dist1, double[] dist2){// s
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
