package parser;

import java.util.ArrayList;

import parser.AnnotationFileParser.AnnotatedGene;
import parser.ScoringOutputParser.ScoredPosition;

public class RPFMergedFileParser {
	public class MergedResult{
		//	IsAnnotated	ContainingGeneName	ContainingGBGeneName	txStart	txEnd	cdsStart	cdsEnd	genomicRegion	frameShift

		private String contig;
		private int position;
		private boolean isPlusStrand;
		private String codon;
		private double[] scores;
		private double[] rpfORFQuantities;
		private double[] rpfORFQuantitiesNormalized;
		
		private double[] rpfQuantityChanges;
		private double[] rpfuORFQuantities;
		private double[] rpfuORFQuantitiesNormalized;
		
		private double[] rnaORFQuantities;
		//private double[] rnaQuantityChanges;
		//private double[] rnauORFQuantities;
		
		private double scoreRatio;
		private double scorePvalue = -1;		
		private double rpfORFQuantitiesNormalizedRatio;
		private double rpfORFQuantitiesNormalizedPvalue = -1;
		
		private double rpfuORFQuantitiesNormalizedRatio;
		private double rpfuORFQuantitiesNormalizedPvalue = -1;
		
		private String genomicRegion;
		private String frameShift;
		private boolean isAnnotated;
		private AnnotatedGene gene;
		
		public MergedResult(ArrayList<ScoredPosition> positions){
			this.contig = positions.get(0).getContig();
			this.position = positions.get(0).getPosition();
			this.isPlusStrand = positions.get(0).isPlusStrand();
			this.codon = positions.get(0).getCodon();
			this.genomicRegion = positions.get(0).getGenomicRegion();
			this.frameShift = positions.get(0).getFrameShift();
			this.isAnnotated = positions.get(0).isAnnotated();
			this.gene = positions.get(0).getGene();
			
			this.scores = new double[positions.size()];
			
			
			double[] mostExtremeRatios = new double[2];
			int[] signs = new int[2];		
			for(int i=0;i<scores.length;i++){
				for(int j=i+1;j<scores.length;j++){	
					for(int k=0;k<mostExtremeRatios.length;k++){
						double v = k==0? Math.log10(scores[i]/scores[j]) : Math.log10(rpfQuantities[i]/rpfQuantities[j]);
						if(mostExtremeRatios[k] < Math.abs(v)){
							mostExtremeRatios[k] = Math.abs(v);
							if(v<0) signs[k] = -1;
							else signs[k] = 1;
						}					
					}				
				}
			}
			
			for(int k=0;k<mostExtremeRatios.length;k++){
				mostExtremeRatios[k] *= signs[k];
			}
			scoreRatio = mostExtremeRatios[0];
			quantityRatio = mostExtremeRatios[1];
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
		
		public void setScorePvalue(double scorePvalue) {
			this.scorePvalue = scorePvalue;
		}

		public void setQuantityPvalue(double quantityPvalue) {
			this.quantityPvalue = quantityPvalue;
		}

		public double getScoreRatio() {
			return scoreRatio;
		}

		public double getQuantityRatio() {
			return quantityRatio;
		}

		
	}
}
