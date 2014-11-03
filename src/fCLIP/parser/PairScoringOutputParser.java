package fCLIP.parser;

import launcher.RNAfoldLauncher;
import launcher.RNAzLauncher;
import parser.MafParser;
import util.DnDsCalculator;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
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
		private int blatHits;
		private double seqEntropy = 0;
		private double structureEntropy = 0;
		private double hetero = 0;
		
		public int getDepth() {
			return depth;
		}

		public double getEnergy() {
			return energy;
		}
		
		public double getSCI(){
			return sci;
		}
		
		public double getSeqEntropy(){
			return seqEntropy;
		}
		
		public double getStructureEntropy(){
			return structureEntropy;
		}
		
		public double getHetero(){
			return hetero;
		}
				
		public double getZScore(){
			return zscore;
		}

		public int getHairpinNumber() {
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
		
		public double getThreePScore(){
			return sps[0].getThreePScore()>=0?sps[0].getThreePScore() : sps[0].getFivePScore();
		}
		
		public double getFivePScore(){
			return sps[1].getThreePScore()>=0?sps[1].getThreePScore() : sps[1].getFivePScore();
		}
		
		static public String getHeader(){
			return "Contig1\tStrand1\t5pPosition1\t3pPosition1\t5p5preads1\t5p3preads1\t3p5preads1\t3p3preads1\tAccession1\tGeneName1\tGenomicRegion5p1\tGenomicRegion3p1\tContig2\tStrand2\t5pPosition2\t3pPosition2\t5p5preads2\t5p3preads2\t3p5preads2\t3p3preads2\tAccession2\tGeneName2"
					+ "\tGenomicRegion5p2\tGenomicRegion3p2\tClass\tPredictionScore\tBlatHits\t%SCI\tz-score\tDnDs\tSeedConservation\tDepth\tEnergy"
					+ "\tHairpinNumber\tLeftPaired\tRightPaired\tOverhang\t5PScore\t3pScore\tDepth1\tDepth2\tEnergy1\tEnergy2\tHairpinNumber1\tHairpinNumber2\tblastHit1\tblastHit2\tSeqEntropy\tStructureEntropy\tHetero"
					+ "\tSequence\tSequence1\tSequence2";
		}
		//5p5preads\t5p3preads\t3p5preads\t3p3preads\t
		public ScoredPair(ScoredPosition p1, ScoredPosition p2, int seqLength){
			sps = new ScoredPosition[2];
			sps[0] = p1;
			sps[1] = p2;
			String seq1 = getTruncatedSequence(p1, seqLength);
			String seq2 = getTruncatedSequence(p2, seqLength);
			seq = p1.is3pScored()? seq1 + seq2 : seq2 + seq1;
			RNAfoldLauncher fold = new RNAfoldLauncher(seq, Scorer.flankingNTNumber);
			this.depth = fold.getDepth();
			this.energy = fold.getEnergyPerNT();
			this.leftPaired = fold.getLeftPaired();
			this.rightPaired = fold.getRightPaired();
			this.overHang = fold.getOverHang();
			this.hairpinNum = fold.getNumberOfHairPins();
			this.seqEntropy = fold.getSeqEntropy();
			this.structureEntropy = fold.getStructureEntropy();
		//	int r1 = sps[0].getFivePReads()[0] + sps[0].getFivePReads()[1] + sps[0].getThreePReads()[0] + sps[0].getThreePReads()[1];
		//	int r2 = sps[1].getFivePReads()[0] + sps[1].getFivePReads()[1] + sps[1].getThreePReads()[0] + sps[1].getThreePReads()[1];
			this.hetero = -1;//(sps[0].getHetero() * r1 + sps[1].getHetero() * r2) / (r1 + r2);
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
				
				sb.append(sps[i].getThreePReads() == null? '_' : sps[i].getThreePReads()[0]); sb.append('\t');
				sb.append(sps[i].getThreePReads() == null? '_' : sps[i].getThreePReads()[1]); sb.append('\t');
				
				sb.append(sps[i].getFivePReads() == null? '_' : sps[i].getFivePReads()[0]); sb.append('\t');
				sb.append(sps[i].getFivePReads() == null? '_' : sps[i].getFivePReads()[1]); sb.append('\t');
				
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getContainingGeneAccessions());
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getContainingGeneNames());
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getGenomicRegions3p());
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getGenomicRegions5p());
			}
			
			sb.append(classification); sb.append('\t');
			sb.append(predictionScore); sb.append('\t');
			sb.append(blatHits); sb.append('\t');
			
			sb.append(sci); sb.append('\t');
			sb.append(zscore); sb.append('\t');
			sb.append(dnds); sb.append('\t');
			sb.append(sps[0].getThreePScore()>=0?sps[0].getSeedConservation() : sps[1].getSeedConservation()); sb.append('\t');
			
			
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNum); sb.append('\t');
			sb.append(leftPaired); sb.append('\t');
			sb.append(rightPaired); sb.append('\t');
			sb.append(overHang); sb.append('\t');
			sb.append(sps[0].getThreePScore()>=0?sps[0].getThreePScore() : sps[1].getThreePScore()); sb.append('\t');
			sb.append(sps[0].getFivePScore()>=0?sps[0].getFivePScore() : sps[1].getFivePScore()); sb.append('\t');
			
			
			sb.append(sps[0].getDepth()); sb.append('\t');
			sb.append(sps[1].getDepth()); sb.append('\t');			
			
			sb.append(sps[0].getEnergy()); sb.append('\t');
			sb.append(sps[1].getEnergy()); sb.append('\t');			
			sb.append(sps[0].getHairpinNumber()); sb.append('\t');
			sb.append(sps[1].getHairpinNumber()); sb.append('\t');

			sb.append(sps[0].getBlatHits()); sb.append('\t');
			sb.append(sps[1].getBlatHits()); sb.append('\t');
	//		sb.append(sps[0].getReadLengthVariance()); sb.append('\t');
	//		sb.append(sps[1].getReadLengthVariance()); sb.append('\t');		
			sb.append(seqEntropy); sb.append('\t');
			sb.append(structureEntropy); sb.append('\t');
			sb.append(hetero); sb.append('\t');
			//sb.append(fivePhetero); sb.append('\t');
			
			sb.append(ScoredPosition.getCleavedSequence(seq, Scorer.flankingNTNumber)); sb.append('\t');
			sb.append(ScoredPosition.getCleavedSequence(sps[0].getSeq(), Scorer.flankingNTNumber)); 
			sb.append('\t');
			sb.append(ScoredPosition.getCleavedSequence(sps[1].getSeq(), Scorer.flankingNTNumber)); 
			
			//\tDepth1\tDepth2\tEnergy1\tEnergy2\tHairpinNumber1\tHairpinNumber2\tSequence1\tSequence2
			return sb.toString();
		}	
			
		/*	
	    public Instance toInstance(Instances insts){
			Instance inst = new DenseInstance(6);
			inst.setDataset(insts);
			inst.setValue(0, getEnergy());
			inst.setValue(1, hairpinNum);
			if(sci<0) inst.setMissing(2);
			else inst.setValue(2, sci);
			if(zscore == 10000) inst.setMissing(3);
			else inst.setValue(3, zscore); 
			//inst.setMissing(4);
			inst.setValue(4, overHang);
			inst.setClassMissing();
			return inst;
		}*/
		
		
		
		public String toArffString(){
			
			StringBuilder sb = new StringBuilder();
			sb.append(energy); sb.append(',');
			sb.append(hairpinNum); sb.append(',');
			sb.append(sci<0? "?" : sci); sb.append(',');
			sb.append(zscore==10000? "?" : zscore); sb.append(',');
			//sb.append(sps[0].getThreePScore()>=0?sps[0].getThreePScore() : sps[1].getThreePScore()); sb.append(',');
			//sb.append(sps[0].getFivePScore()>=0?sps[0].getFivePScore() : sps[1].getFivePScore()); sb.append(',');
		//	sb.append(seqEntropy); sb.append(',');
			sb.append(structureEntropy); sb.append(',');
			sb.append("?"); sb.append(',');
			//sb.append(fivePhetero); sb.append(',');
			
			sb.append(classification == null? "?" : classification.toUpperCase());
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
		
		public void setBlatHits(int n){
			blatHits = n;
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
						sps[i].is3pScored() ? sps[i].getThreePposition() // + (strands[i]? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
								: sps[i].getFivePposition() + (strands[i]? - length + 1 : length - 1);			
						
			}
			
			String[] seqs = mafParser.getSeqs(contigs, positions, strands, length);
			
			this.dnds = DnDsCalculator.calculate(seqs);
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
						sps[i].is3pScored() ? sps[i].getThreePposition()// + (strands[i]? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
								: sps[i].getFivePposition() + (strands[i]?  - length + 1: length - 1);			
						
			}			
		//	System.out.println(ScoredPosition.getCleavedSequence(this.seq, Scorer.flankingNTNumber));
		//	System.out.println(mafParser.getSeqsInMafFormat(contigs, positions, strands, length));			
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
