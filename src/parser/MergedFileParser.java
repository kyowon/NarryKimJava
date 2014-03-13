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
		private double[] harrStartScores;
		private double[] harrStopScores;
		private double[] rpfStartScores;		
		private double[] rpfStopScores;
		
		private double[] rpfCDSQuantities;//
		private double[] rpfCDSQuantitiesNormalized;//
		
		private double[] rpfPositionQuantities;//
		private double[] rpfPositionQuantitiesNormalized;//	
		private double[] rpfPositionQuantityChanges;//
		private double[] rpfAroundStopCodonQuantityChanges;//
		private int stopPosition = -1;
		
		private double[] rnaCDSQuantities;//
		
		private double[] rnaPositionQuantities;//
		private double[] rnaPositionQuantityChanges;//
		
		private String genomicRegion;
		private String frameShift;
		private double dsdnRatio = -1;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		
		
		
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
			
			if(mafParser!=null){ 
				this.dsdnRatio = DsDnCalculator.calculate(mafParser.getSeqs(contig, position, isPlusStrand, 150));
				this.dsdnRatio = Math.log(dsdnRatio + .0001)/Math.log(2);
			}
			
			
			this.harrStartScores = new double[harrScorers.length];
			for(int i=0;i<this.harrStartScores.length;i++){
				this.harrStartScores[i] = harrScorers[i].getStartScore(this.contig, this.position, this.isPlusStrand, false);
			}
			
			this.harrStopScores = new double[harrScorers.length];
			for(int i=0;i<this.harrStopScores.length;i++){
				this.harrStopScores[i] = harrScorers[i].getStopScore(this.contig, this.position, this.isPlusStrand, false, maxLengthUntilStopcodon);
			}
			
			if(rpfScorers !=null){
				this.rpfStartScores = new double[rpfScorers.length];
				for(int i=0;i<this.rpfStartScores.length;i++){
					this.rpfStartScores[i] = rpfScorers[i].getStartScore(this.contig, this.position, this.isPlusStrand, maxLengthUntilStopcodon, false);
				}
				this.rpfStopScores = new double[rpfScorers.length];
				for(int i=0;i<this.rpfStopScores.length;i++){
					this.rpfStopScores[i] = rpfScorers[i].getStopScore(this.contig, this.position, this.isPlusStrand, false, maxLengthUntilStopcodon);
				}
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
				rpfAroundStopCodonQuantityChanges = new double[rpfQuantifiers.length];
				for(int i=0; i<rpfPositionQuantities.length;i++){
					rpfPositionQuantities[i] = rpfQuantifiers[i].getPositionRPKM(contig, this.position, isPlusStrand, maxLengthUntilStopcodon);
					rpfPositionQuantityChanges[i] = rpfQuantifiers[i].getPositionQuantityChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength, positionQuantityOffset);
					rpfPositionQuantityChanges[i] = Math.log(rpfPositionQuantityChanges[i]+.0001)/Math.log(2);
					ArrayList<Double> qs = rpfQuantifiers[i].getNextStopCodonQuantityChangeRatioNStopPosition(contig, this.position, isPlusStrand, positionQuantityChangeLength, maxLengthUntilStopcodon);
					if(qs!=null){
						rpfAroundStopCodonQuantityChanges[i] = qs.get(0);
						rpfAroundStopCodonQuantityChanges[i] = Math.log(rpfAroundStopCodonQuantityChanges[i]+.0001)/Math.log(2);
						stopPosition = (int)((double)qs.get(1));
					}//else rpfAroundStopCodonQuantityChanges = null;
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
					rnaPositionQuantities[i] = rnaQuantifiers[i].getPositionRPKM(contig, this.position, isPlusStrand, maxLengthUntilStopcodon);
					rnaPositionQuantityChanges[i] = rnaQuantifiers[i].getPositionQuantityChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength, positionQuantityOffset);
					rnaPositionQuantityChanges[i] = Math.log(rnaPositionQuantityChanges[i]+.0001)/Math.log(2);
				}
			}
			if(rpfQuantifiers!=null && rnaQuantifiers!=null){
				if(gene != null){
					rpfCDSQuantitiesNormalized = new double[rpfQuantifiers.length];
					for(int i=0; i<rpfCDSQuantitiesNormalized.length;i++){
						rpfCDSQuantitiesNormalized[i] = rpfCDSQuantities[i] == 0? rpfCDSQuantities[i] : rpfCDSQuantities[i] / rnaCDSQuantities[i];
						rpfCDSQuantitiesNormalized[i] = Math.log(rpfCDSQuantitiesNormalized[i]+.0001)/Math.log(2);
					}
				}
				rpfPositionQuantitiesNormalized = new double[rpfQuantifiers.length];
				for(int i=0; i<rpfPositionQuantitiesNormalized.length;i++){
					rpfPositionQuantitiesNormalized[i] = rpfPositionQuantities[i] == 0? rpfPositionQuantities[i] : rpfPositionQuantities[i] / rnaPositionQuantities[i];
					rpfPositionQuantitiesNormalized[i] = Math.log(rpfPositionQuantitiesNormalized[i]+.0001)/Math.log(2);
				}					
			}
		}
		
		public MergedResult(String s){ //TODO fix next time- consider "_"
			int i = 0;
			String[] token = s.split("\t");
			contig = token[i++];
			position = Integer.parseInt(token[i++]);
			isPlusStrand = token[i++].equals("+");
			codon = token[i++];
			genomicRegion = token[i++];
			frameShift = token[i++];
			dsdnRatio = token[i].equals('_') ? -1 : Double.parseDouble(token[i]);
			i++;
			if(!token[i].equals("_"))
				isAnnotated = token[i].equals("T");
			i++;
			harrStartScores = subParse(token[i++]);
			harrStopScores = subParse(token[i++]);
			rpfStartScores = subParse(token[i++]);
			rpfStopScores = subParse(token[i++]);
			rpfPositionQuantities = subParse(token[i++]);
			rnaPositionQuantities = subParse(token[i++]);
			rpfPositionQuantitiesNormalized = subParse(token[i++]);
			rpfPositionQuantityChanges = subParse(token[i++]);
			rpfAroundStopCodonQuantityChanges = subParse(token[i++]);
			stopPosition = Integer.parseInt(token[i++]); //TODO
			rnaPositionQuantityChanges = subParse(token[i++]);
			rpfCDSQuantities = subParse(token[i++]);
			rnaCDSQuantities = subParse(token[i++]);
			rpfCDSQuantitiesNormalized = subParse(token[i++]);
			//i++;
			
			if(!token[i].equals("_")){
				StringBuffer gs = new StringBuffer();
				for(;i<token.length;i++){
					gs.append(token[i]); gs.append('\t');
				}
				gene = new AnnotatedGene(gs.toString());
			}
		}
		
		private double[] subParse(String s){
			String[] token = s.split(";");
			double[] ret = new double[token.length];
			for(int i=0;i<ret.length;i++){
				ret[i] = Double.parseDouble(token[i]);
			}
			return ret;
		}
		
		public String getHeader(){
			String header =  "Contig\tPosition\tStrand\tCodon\tGenomicRegion\tFrameShift\tDsDnratio(log2)\tisAnnotated\t";
			header += subGetHeader("Harr Start Score", harrStartScores) + "\t";
			header += subGetHeader("Harr Stop Score", harrStopScores) + "\t";
			header += subGetHeader("RPF Start Score", rpfStartScores) + "\t";
			header += subGetHeader("RPF Stop Score", rpfStopScores) + "\t";
			
			header += subGetHeader("RPF_Pos_RPKMs", rpfPositionQuantities) + "\t";
			header += subGetHeader("RNA_Pos_RPKMs", rnaPositionQuantities) + "\t";
			header += subGetHeader("RPF/RNA_Pos_RPKMs(log2)", rpfPositionQuantitiesNormalized) + "\t";
			
			header += subGetHeader("RPF_QuanChange(log2)", rpfPositionQuantityChanges) + "\t";
			header += subGetHeader("RPF_StopCodonQuanChange(log2)", rpfAroundStopCodonQuantityChanges) + "\t";			
			header += "StopCodonPosition\t";
			header += subGetHeader("RNA_QuanChange(log2)", rnaPositionQuantityChanges) + "\t";
			
			header += subGetHeader("RPF_CDS_RPKMs", rpfCDSQuantities) + "\t";
			header += subGetHeader("RNA_CDS_RPKMs", rnaCDSQuantities) + "\t";
			header += subGetHeader("RPF/RNA_CDS_RPKMs", rpfCDSQuantitiesNormalized) + "\t";
			
			header += "DistBetweenPos&CDS(RPF/RNA)\t";
			
			header += AnnotatedGene.getHeader();
			
			return header;
		}
		
		private String subGetHeader(String prefix, double[] quantities){
			String header = "";
			
			if(quantities != null){
				for(int i=0;i<quantities.length-1;i++){
					header += prefix + (i+1) + ";";
				}
				header += prefix + quantities.length + "\t";
			}else{
				header += prefix + "\t";
			}
			header += prefix + "Ununiformity";
			
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
			sb.append(frameShift);sb.append('\t');
			sb.append(dsdnRatio>=0?dsdnRatio : '_');sb.append('\t');
			sb.append(gene == null ? '_' : (isAnnotated? 'T' : 'F'));sb.append('\t');
			
			subToString(harrStartScores, sb);
			sb.append('\t');
			
			subToString(harrStopScores, sb);
			sb.append('\t');
			
			
			subToString(rpfStartScores, sb);
			sb.append('\t');
			
			subToString(rpfStopScores, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantities, sb);
			sb.append('\t');
			
			subToString(rnaPositionQuantities, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantitiesNormalized, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantityChanges, sb);
			sb.append('\t');
			
			subToString(rpfAroundStopCodonQuantityChanges, sb);
			sb.append('\t');
			
			sb.append(stopPosition>=0? stopPosition : "_");
			sb.append('\t');
			
			subToString(rnaPositionQuantityChanges, sb);
			sb.append('\t');
			
			subToString(rpfCDSQuantities, sb, rpfPositionQuantities == null? 0 : rpfPositionQuantities.length);
			sb.append('\t');
			
			subToString(rnaCDSQuantities, sb, rnaPositionQuantities == null? 0 : rnaPositionQuantities.length);
			sb.append('\t');
			
			subToString(rpfCDSQuantitiesNormalized, sb, rpfPositionQuantitiesNormalized == null? 0 : rpfPositionQuantitiesNormalized.length);
			sb.append('\t');
			
			if(rpfCDSQuantitiesNormalized != null && rpfPositionQuantitiesNormalized != null)
				sb.append(String.format("%.3e", getKL(rpfPositionQuantitiesNormalized, rpfCDSQuantitiesNormalized)));
			else sb.append("_");
			sb.append('\t');
			
			sb.append(gene == null ? AnnotatedGene.getEmptyGeneString() : gene.toString());
			return sb.toString();
		}
		
		
		public String GetTrainingSetHeader(){
			return "DsDn,HStart,HStop,RStart, RStop, RPFChange,RNAChange,RPFRPKM,RNARPKM,TE,RELEASE,Class";
		}
		
		
		public String ToTrainingSetString(AnnotationFileParser parser, int index){
			
		//	if(rnaPositionQuantities!=null && rnaPositionQuantities.length <= index) return null;
		//	if(rpfPositionQuantities!=null && rpfPositionQuantities.length <= index) return null;
		//	if(rnaPositionQuantities!=null && rnaPositionQuantities.length > index && rnaPositionQuantities[index]<=rnaRPKMThreshold) return null;
		//	if(rpfPositionQuantities!=null && rpfPositionQuantities.length > index && rpfPositionQuantities[index]<=rpfRPKMThreshold) return null;
			
			int c = 0 ;
			
			if(isAnnotated && genomicRegion.equals("NM_ORF")){
				c = 1;
			}else if(genomicRegion.equals("NM_ORF") && frameShift.equals("0") && (codon.equals("ATG") || codon.equals("CTG"))){
				c = 2;
			}else if(genomicRegion.equals("NM_3_UTR")){// && parser.getContainingGene(contig, isPlusStrand, stopPosition) == null){
				c = 3;
			}else if(genomicRegion.equals("NM_5_UTR") && (codon.equals("ATG") || codon.equals("CTG"))){
				c = 4;
			}//else if(genomicRegion.endsWith("Intron") && (codon.equals("ATG") || codon.equals("CTG"))){
			//	c = 5;
			//}
			
			if(c>0){
				StringBuffer sb = new StringBuffer();
				sb.append(dsdnRatio);sb.append(',');

				if(harrStartScores!=null && harrStartScores.length > index){
					sb.append(harrStartScores[index]);sb.append(',');
				}else return null;

				if(harrStopScores!=null && harrStopScores.length > index){
					sb.append(harrStopScores[index]);sb.append(',');
				}else return null;
				if(rpfStartScores!=null && rpfStartScores.length > index){
					sb.append(rpfStartScores[index]);sb.append(',');
				}else return null;
				
				if(rpfStopScores!=null && rpfStopScores.length > index){
					sb.append(rpfStopScores[index]);sb.append(',');
				}else return null;
				
				if(rpfPositionQuantityChanges!=null && rpfPositionQuantityChanges.length > index){
					sb.append(rpfPositionQuantityChanges[index]);sb.append(',');
				}else return null;
				
				if(rnaPositionQuantityChanges!=null && rnaPositionQuantityChanges.length > index){
					sb.append(rnaPositionQuantityChanges[index]);sb.append(',');
				}else return null;
				
				if(rpfPositionQuantities!=null && rpfPositionQuantities.length > index){
					sb.append(Math.log(rpfPositionQuantities[index]+0.0001)/Math.log(2));sb.append(',');
				}else return null;
				
				if(rnaPositionQuantities!=null && rnaPositionQuantities.length > index){
					sb.append(Math.log(rnaPositionQuantities[index]+.0001)/Math.log(2));sb.append(',');
				}else return null;
				
				
				if(rpfPositionQuantitiesNormalized!=null && rpfPositionQuantitiesNormalized.length > index){
					sb.append(rpfPositionQuantitiesNormalized[index]);sb.append(',');
				}else return null;
	
				if(rpfAroundStopCodonQuantityChanges!=null && rpfAroundStopCodonQuantityChanges.length > index){
					sb.append(rpfAroundStopCodonQuantityChanges[index]);sb.append(',');
				}else return null;
				if(c==1){
					sb.append("TS");
				}else if(c==2){
					sb.append("T");
				}else if(c==3){
					sb.append("3U");
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
			subToString(quantities, sb, quantities == null? 0 : quantities.length);
		}
		
		private void subToString(double[] quantities, StringBuffer sb, int numSample){ 
			if(quantities!=null && quantities.length>0){
				for(int i=0; i<quantities.length-1; i++){
					sb.append(String.format("%.3e", quantities[i])); sb.append(';');
				}
				sb.append(String.format("%.3e", quantities[quantities.length-1])); 
				sb.append('\t');
				if(quantities.length > 1){
					sb.append(String.format("%.3e", getNonuniformity(quantities)));
				}
			}else{
				//sb.append("_");
				for(int i=0; i<numSample-1; i++){
					sb.append("_;");
				}
				sb.append("_\t_");
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
			for(double v : dist1){
				sum += v;
			}
			for(int i=0;i<dist1.length;i++){
				dist1[i] /= sum;
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
				dist2[i] /= sum;
			}
			double kl = 0;
			for(int i=0;i<dist1.length;i++){
				double v= dist1[i];
				if(v > 0) kl += v * Math.log(v/(dist2[i] == 0? min * .1 : dist2[i]));
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
