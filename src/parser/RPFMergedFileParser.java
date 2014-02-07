package parser;

import java.util.ArrayList;

import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;
import rpf.Quantifier;
import rpf.Scorer;

public class RPFMergedFileParser {
	public class MergedResult{
		//	IsAnnotated	ContainingGeneName	ContainingGBGeneName	txStart	txEnd	cdsStart	cdsEnd	genomicRegion	frameShift
		private int positionQuantityChangeLength = 300;
		private int positionQuantityMaxLength = 500000;
		
		private String contig;
		private int position;
		private boolean isPlusStrand;
		private String codon;
		private double[] scores;
		
		private double[] rpfCDSQuantities;//
		private double[] rpfCDSQuantitiesNormalized;//
		
		private double[] rpfPositionQuantities;//
		private double[] rpfPositionQuantitiesNormalized;//	
		private double[] rpfPositionQuantityChanges;//
		
		private double[] rnaCDSQuantities;//
		
		private double[] rnaPositionQuantities;//
		private double[] rnaPositionQuantityChanges;//
		
		private String genomicRegion;
		private String frameShift;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		
		public MergedResult(ScoredPosition scoredPosition, Scorer[] scorers, Quantifier[] rpfQuantifiers, Quantifier[] rnaQuantifiers){
			this.contig = scoredPosition.getContig();
			this.position = scoredPosition.getPosition();
			this.isPlusStrand = scoredPosition.isPlusStrand();
			this.codon = scoredPosition.getCodon();
			this.genomicRegion = scoredPosition.getGenomicRegion();
			this.frameShift = scoredPosition.getFrameShift();
			this.isAnnotated = scoredPosition.isAnnotated();
			this.gene = scoredPosition.getGene();
			
			this.scores = new double[scorers.length];
			for(int i=0;i<this.scores.length;i++){
				this.scores[i] = scorers[i].getScore(this.contig, this.position, this.isPlusStrand, false);
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
				for(int i=0; i<rpfPositionQuantities.length;i++){
					rpfPositionQuantities[i] = rpfQuantifiers[i].getPositionRPKM(contig, this.position, isPlusStrand, positionQuantityMaxLength);
					rpfPositionQuantityChanges[i] = rpfQuantifiers[i].getPositionQuantatyChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength);
				}
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
					rnaPositionQuantities[i] = rnaQuantifiers[i].getPositionRPKM(contig, this.position, isPlusStrand, positionQuantityMaxLength);
					rnaPositionQuantityChanges[i] = rnaQuantifiers[i].getPositionQuantatyChangeRatio(contig, this.position, isPlusStrand, positionQuantityChangeLength);
				}
			}
			if(rpfQuantifiers!=null && rnaQuantifiers!=null){
				if(gene != null){
					rpfCDSQuantitiesNormalized = new double[rpfQuantifiers.length];
					for(int i=0; i<rpfCDSQuantitiesNormalized.length;i++){
						rpfCDSQuantitiesNormalized[i] = rpfCDSQuantities[i] == 0? rpfCDSQuantities[i] : rpfCDSQuantities[i] / rnaCDSQuantities[i];
					}
				}
				rpfPositionQuantitiesNormalized = new double[rpfQuantifiers.length];
				for(int i=0; i<rpfPositionQuantitiesNormalized.length;i++){
					rpfPositionQuantitiesNormalized[i] = rpfPositionQuantities[i] == 0? rpfPositionQuantities[i] : rpfPositionQuantities[i] / rnaPositionQuantities[i];
				}					
			}
		}
		
		public MergedResult(String s){
			int i = 0;
			String[] token = s.split("\t");
			contig = token[i++];
			position = Integer.parseInt(token[i++]);
			isPlusStrand = token[i++].equals("+");
			codon = token[i++];
			scores = subParse(token[i++]);
			rpfPositionQuantities = subParse(token[i++]);
			rnaPositionQuantities = subParse(token[i++]);
			rpfPositionQuantitiesNormalized = subParse(token[i++]);
			rpfPositionQuantityChanges = subParse(token[i++]);
			rnaPositionQuantityChanges = subParse(token[i++]);
			rpfCDSQuantities = subParse(token[i++]);
			rnaCDSQuantities = subParse(token[i++]);
			rpfCDSQuantitiesNormalized = subParse(token[i++]);
			i++;
			genomicRegion = token[i++];
			frameShift = token[i++];
			if(!token[i].equals("_")){
				isAnnotated = token[i++].equals("T");
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
			String header =  "Contig\tPosition\tStrand\tCodon\t";
			header += subGetHeader("Score", scores) + "\t";
			header += subGetHeader("RPF_Pos_RPKMs", rpfPositionQuantities) + "\t";
			header += subGetHeader("RNA_Pos_RPKMs", rnaPositionQuantities) + "\t";
			header += subGetHeader("RPF/RNA_Pos_RPKMs", rpfPositionQuantitiesNormalized) + "\t";
			
			header += subGetHeader("RPF_QuanChange", rpfPositionQuantityChanges) + "\t";
			header += subGetHeader("RNA_QuanChange", rnaPositionQuantityChanges) + "\t";
			
			header += subGetHeader("RPF_CDS_RPKMs", rpfCDSQuantities) + "\t";
			header += subGetHeader("RNA_CDS_RPKMs", rnaCDSQuantities) + "\t";
			header += subGetHeader("RPF/RNA_CDS_RPKMs", rpfCDSQuantitiesNormalized) + "\t";
			
			header += "DistBetweenPos&CDS(RPF/RNA)\t";
			
			header += "GenomicRegion\tFrameShift\tisAnnotated\t" + AnnotatedGene.getHeader();
			
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
		
		public String toString(){			
			StringBuffer sb = new StringBuffer();
			sb.append(contig);sb.append('\t');
			sb.append(position);sb.append('\t');
			sb.append(isPlusStrand? '+' : '-');sb.append('\t');
			sb.append(codon);sb.append('\t');
			
			subToString(scores, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantities, sb);
			sb.append('\t');
			
			subToString(rnaPositionQuantities, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantitiesNormalized, sb);
			sb.append('\t');
			
			subToString(rpfPositionQuantityChanges, sb);
			sb.append('\t');
			
			subToString(rnaPositionQuantityChanges, sb);
			sb.append('\t');
			
			subToString(rpfCDSQuantities, sb);
			sb.append('\t');
			
			subToString(rnaCDSQuantities, sb);
			sb.append('\t');
			
			subToString(rpfCDSQuantitiesNormalized, sb);
			sb.append('\t');
			
			if(rpfCDSQuantitiesNormalized != null && rpfPositionQuantitiesNormalized != null)
				sb.append(getKL(rpfCDSQuantitiesNormalized, rpfPositionQuantitiesNormalized));
			else sb.append("N/A");
			sb.append('\t');
			
			sb.append(genomicRegion);sb.append('\t');
			sb.append(frameShift);sb.append('\t');
			sb.append(gene == null ? '_' : (isAnnotated? 'T' : 'F'));sb.append('\t');
			sb.append(gene == null ? AnnotatedGene.getEmptyGeneString() : gene.toString());
			return sb.toString();
		}
		
		private void subToString(double[] quantities, StringBuffer sb){
			if(quantities!=null && quantities.length>0){
				for(int i=0; i<quantities.length-1; i++){
					sb.append(quantities[i]); sb.append(';');
				}
				sb.append(quantities[quantities.length-1]); 
				sb.append('\t');
				if(quantities.length > 1){
					sb.append(getNonuniformity(quantities));
				}
			}else{
				sb.append("N/A\tN/A");
			}
		}
		
		
		private double getNonuniformity(double[] dist){
			double sum = 0;
			for(double v : dist){
				sum += v + 10;
			}
			if(sum == 0) return 0;
			double kl = 0;
			for(double v : dist){
				double vn = (v + 10)/sum;
				if(vn > 0) kl += vn * Math.log(vn*dist.length);
			}
			return kl;		
		}
		
		private double getKL(double[] dist1, double[] dist2){
			double sum = 0;
			for(double v : dist1){
				sum += v;
			}
			for(int i=0;i<dist1.length;i++){
				dist1[i] /= sum;
			}
			sum = 0;
			for(double v : dist2){
				sum += v;
			}
			for(int i=0;i<dist2.length;i++){
				dist2[i] /= sum;
			}
			double kl = 0;
			for(int i=0;i<dist1.length;i++){
				double v= dist1[i];
				if(v > 0 && dist2[i] > 0) kl += v * Math.log(v/dist2[i]);
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
