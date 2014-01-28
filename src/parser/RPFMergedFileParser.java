package parser;

import java.util.ArrayList;

import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;
import rpf.Quantifier;

public class RPFMergedFileParser {
	public class MergedResult{
		//	IsAnnotated	ContainingGeneName	ContainingGBGeneName	txStart	txEnd	cdsStart	cdsEnd	genomicRegion	frameShift
		private int quantityChangeLength = 300;
		
		private String contig;
		private int position;
		private boolean isPlusStrand;
		private String codon;
		private double[] scores;
		private double scoreNonuniformity;//
		
		private double[] rpfCDSQuantities;//
		private double[] rpfCDSQuantitiesNormalized;//
		private double rpfCDSNonuniformity;	//
		
		private double[] rpfPositionQuantities;//
		private double[] rpfPositionQuantitiesNormalized;//	
		private double rpfPositionNonuniformity;	//
		private double[] rpfPositionQuantityChanges;//
		
		private double[] rnaCDSQuantities;//
		private double rnaCDSNonuniformity;//
		
		private double[] rnaPositionQuantities;//
		private double rnaPositionNonuniformity;//
		private double[] rnaPositionQuantityChanges;//
		
		private String genomicRegion;
		private String frameShift;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		
		public MergedResult(ArrayList<ScoredPosition> positions, Quantifier[] rpfQuantifiers, Quantifier[] rnaQuantifiers){
			this.contig = positions.get(0).getContig();
			this.position = positions.get(0).getPosition();
			this.isPlusStrand = positions.get(0).isPlusStrand();
			this.codon = positions.get(0).getCodon();
			this.genomicRegion = positions.get(0).getGenomicRegion();
			this.frameShift = positions.get(0).getFrameShift();
			this.isAnnotated = positions.get(0).isAnnotated();
			this.gene = positions.get(0).getGene();
			
			this.scores = new double[positions.size()];
			for(int i=0;i<this.scores.length;i++){
				this.scores[i] = positions.get(i).getScore();
			}
			scoreNonuniformity = getNonuniformity(scores);
			
			if(rpfQuantifiers!=null){
				if(gene != null){				
					rpfCDSQuantities = new double[rpfQuantifiers.length];
					for(int i=0; i<rpfQuantifiers.length;i++){
						rpfCDSQuantities[i] = rpfQuantifiers[i].getCDSQuantity(gene);
					}
				}
				rpfPositionQuantities = new double[rpfQuantifiers.length];
				rpfPositionQuantityChanges = new double[rpfQuantifiers.length];
				for(int i=0; i<rpfPositionQuantities.length;i++){
					rpfPositionQuantities[i] = rpfQuantifiers[i].getPositionQuantity(contig, position, isPlusStrand);
					rpfPositionQuantityChanges[i] = rpfQuantifiers[i].getPositionQuantatyChangeRatio(contig, position, isPlusStrand, quantityChangeLength);
				}
			}
			
			if(rnaQuantifiers!=null){
				if(gene != null){				
					rnaCDSQuantities = new double[rnaQuantifiers.length];
					for(int i=0; i<rnaQuantifiers.length;i++){
						rnaCDSQuantities[i] = rnaQuantifiers[i].getCDSQuantity(gene);
					}
					rnaCDSNonuniformity = getNonuniformity(rnaCDSQuantities);
				}
				rnaPositionQuantities = new double[rnaQuantifiers.length];
				rnaPositionQuantityChanges = new double[rnaQuantifiers.length];
				for(int i=0; i<rnaPositionQuantities.length;i++){
					rnaPositionQuantities[i] = rnaQuantifiers[i].getPositionQuantity(contig, position, isPlusStrand);
					rnaPositionQuantityChanges[i] = rnaQuantifiers[i].getPositionQuantatyChangeRatio(contig, position, isPlusStrand, quantityChangeLength);
				}
				rnaPositionNonuniformity = getNonuniformity(rnaPositionQuantities);
			}
			if(rpfQuantifiers!=null && rnaQuantifiers!=null){
				if(gene != null){
					rpfCDSQuantitiesNormalized = new double[rpfQuantifiers.length];
					for(int i=0; i<rpfCDSQuantitiesNormalized.length;i++){
						rpfCDSQuantitiesNormalized[i] = rpfCDSQuantities[i] == 0? rpfCDSQuantities[i] : rpfCDSQuantities[i] / rnaCDSQuantities[i];
					}
					rpfCDSNonuniformity = getNonuniformity(rpfCDSQuantitiesNormalized);
				}
				rpfPositionQuantitiesNormalized = new double[rpfQuantifiers.length];
				for(int i=0; i<rpfPositionQuantitiesNormalized.length;i++){
					rpfPositionQuantitiesNormalized[i] = rpfPositionQuantitiesNormalized[i] == 0? rpfPositionQuantitiesNormalized[i] : rpfPositionQuantitiesNormalized[i] / rnaPositionQuantities[i];
				}	
				rpfPositionNonuniformity = getNonuniformity(rpfPositionQuantitiesNormalized);
			}else{
				if(gene != null){
					rpfCDSNonuniformity = getNonuniformity(rpfCDSQuantities);
				}
				rpfPositionNonuniformity = getNonuniformity(rpfPositionQuantities);
			}
		}
		
		public MergedResult(String s){
			String[] token = s.split("\t");
			contig = token[0];
			position = Integer.parseInt(token[1]);
			isPlusStrand = token[2].equals("+");
			codon = token[3];
			String[] st = token[4].split(";");
			String[] qt = token[5].split(";");
			scores = new double[st.length-1];
			rpfORFQuantities = new double[st.length-1];
			for(int i=0;i<st.length-1;i++){
				scores[i] = Double.parseDouble(st[i]);
				rpfORFQuantities[i] = Double.parseDouble(qt[i]);
			}
			scoreRatio = Double.parseDouble(token[6]);
			if(!token[7].equals("*")) scorePvalue = Double.parseDouble(token[7]);
			quantityRatio = Double.parseDouble(token[8]);
			if(!token[9].equals("*")) quantityPvalue = Double.parseDouble(token[9]);
			genomicRegion = token[10];
			frameShift = token[11];
			if(!token[12].equals("_")){
				isAnnotated = token[12].equals("T");
				StringBuffer gs = new StringBuffer();
				for(int i=13;i<token.length;i++){
					gs.append(token[i]); gs.append('\t');
				}
				gene = new AnnotatedGene(gs.toString());
			}
		}
		
		public String toString(){
			StringBuffer sb = new StringBuffer();
			sb.append(contig);sb.append('\t');
			sb.append(position);sb.append('\t');
			sb.append(isPlusStrand? '+' : '-');sb.append('\t');
			sb.append(codon);sb.append('\t');
			for(int i=0; i<scores.length - 1;i++){
				sb.append(scores[i]); sb.append(';');
			}
			sb.append(scores[scores.length-1]);
			sb.append('\t');
			for(int i=0; i<rpfORFQuantities.length-1; i++){
				sb.append(rpfORFQuantities[i]); sb.append(';');
			}
			sb.append(rpfORFQuantities[rpfORFQuantities.length-1]);
			sb.append('\t');
			sb.append(scoreRatio);sb.append('\t');
			sb.append(scorePvalue<0? '*' : scorePvalue);sb.append('\t');
			sb.append(quantityRatio);sb.append('\t');
			sb.append(quantityPvalue<0? '*' : quantityPvalue);sb.append('\t');
			sb.append(genomicRegion);sb.append('\t');
			sb.append(frameShift);sb.append('\t');
			sb.append(gene == null ? '_' : (isAnnotated? 'T' : 'F'));sb.append('\t');
			sb.append(gene == null ? AnnotatedGene.getEmptyGeneString() : gene.toString());
			return sb.toString();
		}
		
		private double getNonuniformity(double[] dist){
			double sum = 0;
			for(double v : dist){
				sum += v;
			}
			if(sum == 0) return 0;
			double kl = 0;
			for(double v : dist){
				double vn = v/sum;
				if(vn > 0) kl += vn * Math.log(vn*dist.length);
			}
			return kl;		
		}
		
	}
}
