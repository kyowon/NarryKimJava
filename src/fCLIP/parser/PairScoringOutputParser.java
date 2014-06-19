package fCLIP.parser;

import parser.MafParser;
import util.DnDsCalculator;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import fCLIP.RNAfoldLauncher;
import fCLIP.RNAzLauncher;
import fCLIP.Scorer;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;

public class PairScoringOutputParser {
	
	public static class ScoredPair{
		private ScoredPosition[] sps;
		private int depth;
		private double energy;
		private int hairpinNum;
		private int leftPaired;
		private int rightPaired;
		private double sci = -1;
		private double zscore = 10000;
		private double dnds = -1;
		private int overHang;
		private String seq;
		private String classification;
		private double predictionScore;
		
		public void getClassification(String classification, double predictionsScore){
			this.classification = classification;
			this.predictionScore = predictionScore;
		}
		
		public int getDepth() {
			return depth;
		}

		public double getEnergy() {
			return energy;
		}

		public int getHairpinNum() {
			return hairpinNum;
		}

		public int getLeftPaired() {
			return leftPaired;
		}

		public int getRightPaired() {
			return rightPaired;
		}

		public int getOverHang() {
			return overHang;
		}

		public String getSeq() {
			return seq;
		}

		public String getClassification() {
			return classification;
		}

		public double getPredictionScore() {
			return predictionScore;
		}

		
		static public String getHeader(){
			return "Contig1\tStrand1\t3pPosition1\t5pPosition1\tAccession1\tGeneName1\tGenomicRegion1\tContig2\tStrand2\t3pPosition2\t5pPosition2\tAccession2\tGeneName2"
					+ "\tGenomicRegion2\tClass\tPredictionScore\tSCI\tz-score\tDnDs\tDepth\tEnergy"
					+ "\tHairpinNumber\tLeftPaired\tRightPaired\tOverhang\tSequence";
		}
		
		public ScoredPair(ScoredPosition p1, ScoredPosition p2, int seqLength){
			sps = new ScoredPosition[2];
			sps[0] = p1;
			sps[1] = p2;
			String seq1 = getTruncatedSequence(p1, seqLength);
			String seq2 = getTruncatedSequence(p2, seqLength);
			seq = p1.is3pScored()? seq1 + seq2 : seq2 + seq1;
			RNAfoldLauncher fold = new RNAfoldLauncher(seq, Scorer.flankingNTNumber);
			this.depth = fold.getDepth();
			this.energy = fold.getEnergy();
			this.leftPaired = fold.getLeftPaired();
			this.rightPaired = fold.getRightPaired();
			this.overHang = fold.getOverHang();
			this.hairpinNum = fold.getNumberOfHairPins();
		}
		
		private String getTruncatedSequence(ScoredPosition p, int seqLength){
			String seq = p.getSeq();
			if(p.is3pScored()) seq = seq.substring(0, Math.min(seq.length(), seqLength/2));
			else seq = seq.substring(Math.max(0, seq.length() - seqLength/2));
			return seq;
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<2;i++){
				sb.append(sps[i].getContig()); sb.append('\t');
				sb.append(sps[i].isPlusStrand()? '+': '-'); sb.append('\t');
				sb.append(sps[i].getThreePposition()); sb.append('\t');
				sb.append(sps[i].getFivePposition()); sb.append('\t');
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getContainingGeneAccessions());
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getContainingGeneNames());
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getGenomicRegions());
			}
			
			sb.append(classification); sb.append('\t');
			sb.append(predictionScore); sb.append('\t');
			sb.append(sci); sb.append('\t');
			sb.append(zscore); sb.append('\t');
			sb.append(dnds); sb.append('\t');
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNum); sb.append('\t');
			sb.append(leftPaired); sb.append('\t');
			sb.append(rightPaired); sb.append('\t');
			sb.append(overHang); sb.append('\t');
			sb.append(seq);
			return sb.toString();
		}	
			
		public Instance toInstance(Instances insts){
			Instance inst = new DenseInstance(6);
			inst.setDataset(insts);
			inst.setValue(0, getEnergy());
			inst.setValue(1, hairpinNum);
			if(sci<0) inst.setMissing(2);
			else inst.setValue(2, sci);
			if(zscore == 10000)inst.setMissing(3);
			else inst.setValue(3, zscore); 
			inst.setValue(4, overHang);
			inst.setClassMissing();
			return inst;
		}
		
		
		
		public String toArffString(){
			StringBuilder sb = new StringBuilder();
			//sb.append("1");sb.append(','); // always paired
		//	sb.append(depth);sb.append(',');
			sb.append(energy); sb.append(',');
			sb.append(hairpinNum); sb.append(',');
			sb.append(sci<0? "?" : sci); sb.append(',');
			sb.append(zscore==10000? "?" : zscore); sb.append(',');
			//	sb.append(leftPaired); sb.append(',');
		//	sb.append(rightPaired); sb.append(',');
			sb.append(overHang); sb.append(',');
			sb.append("?");
			return sb.toString();
		}

		public void setClassification(String c, double predictionScore) {
			classification = c;
			this.predictionScore = predictionScore;	
		}
	
		
		public void setRNAzScores(MafParser mafParser){
			RNAzLauncher rna = getRNAzLauncher(mafParser);
			this.sci = rna.getSCI();
			this.zscore = rna.getZScore();
		}
		
		
		public void setDnDs(MafParser mafParser){
			String[] contigs = new String[2];
			int[] positions = new int[2];
			boolean[] strands = new boolean[2];
			
			int length = seq.length() / 2;
			
			for(int i=0;i<2;i++){
				contigs[i] = sps[i].getContig();
				strands[i] = sps[i].isPlusStrand();
				positions[i] = 
						sps[i].is3pScored() ? sps[i].getThreePposition() + (strands[i]? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
								: sps[i].getFivePposition() + (strands[i]? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
						
			}
			
			String[] seqs = mafParser.getSeqs(contigs, positions, strands, length);
			//for(String seq : seqs) System.out.println(seq);
			
			this.dnds = DnDsCalculator.calculate(seqs, false);
		}
				
		private String getMafString(MafParser mafParser){
			String[] contigs = new String[2];
			int[] positions = new int[2];
			boolean[] strands = new boolean[2];
			
			int length = seq.length() / 2;
			
			for(int i=0;i<2;i++){
				contigs[i] = sps[i].getContig();
				strands[i] = sps[i].isPlusStrand();
				positions[i] = 
						sps[i].is3pScored() ? sps[i].getThreePposition() + (strands[i]? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
								: sps[i].getFivePposition() + (strands[i]? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
						
			}			
			
			return mafParser.getSeqsInMafFormat(contigs, positions, strands, length);
		}
		
		
		public RNAzLauncher getRNAzLauncher(MafParser mafParser){
			String maf = getMafString(mafParser);
			//System.out.println(maf);
			return new RNAzLauncher(maf);
		}
		
		static public void main(String[] args){
			MafParser mafParser = new MafParser("/media/kyowon/Data1/RPF_Project/genomes/hg19/maf/");
			mafParser.readIndexFile();
			ScoredPosition sp1 = new ScoredPosition("chr5	+	10437334	10437404	T	U	0.8928715112	_	-1	10000	-1	NM_001270660,NM_001270661,NM_005885	MARCH6	_	10437275	10437375	10437364	10437464	0.0840371969	0.0862782024	14	-0.1273684211	2	9	3	28	0.5894736842	GGGTTGAGTACTCTAAGCTCTCTACTTACTGTGAGAAGTTTTCTACATTGTAAAATTAAAAGATTATATTTAAAATACTTCTTTAGGATGTTATT");
			ScoredPosition sp2 = new ScoredPosition("chr6	+	7590461	7590527	F	U	0.5722689076	_	0.01	-0.86	1.8486003876	NM_152551	SNRNP48	NM_ORF			7590487	7590587	-1	0.0903463068	34	-0.5857142857	1	12	8	6	0.5934065934	CTGCAGTTCGGCCGCTTCCTCTTGGCGGGTGGGCTGCAGCTATGGAGGGCGAGCCTCCACCTGTGGAGGAGCGGCGGCGGCTGCAGGAGGA");
			
			ScoredPair test = new ScoredPair(sp1, sp2, 50);
			
			RNAzLauncher rna = test.getRNAzLauncher(mafParser);
			test.setDnDs(mafParser);
			System.out.println(rna.getSCI() + " " + rna.getZScore());
			
			
		}
		
	}
}
