package fCLIP.parser;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import launcher.RNAcofoldLauncher;
import parser.BufferedLineReader;
import parser.MafParser;
import fCLIP.FCLIP_Scorer;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class ScoredPairOutputParser {
	
	
	
	public static class ScoredPair implements Comparable<ScoredPair>{
		private ScoredPosition[] sps;
		private int depth;
		private double energy;
		private int hairpinNum;
		private int leftPaired;
		private int rightPaired;
		private boolean isFeasibleFold;
	//	private double sci = -1;
	//	private double zscore = 10000;
	//	private double[] phyloPConservationPvalues = new double[2];
	//	private double[] phyloPAccPvalues = new double[2];
		//private double seedConservation = 0;
		private double threePScore = 0;
		private double fivePScore = 0;
	//	private int[] blatHits;
		private String[] seqs;
		
//		private double dnds = -1;
		private int overHang;
		private String seq;
		private String classification = "";
		private double predictionScore;
		private int[] blatHits;
		private double seqEntropy = 0;
		private double structureEntropy = 0;
		private double hetero = 0;
		private int[] matchedNum = new int[2];
		private String misc = "";
		
		public int compareTo(ScoredPair o) {
			return Double.compare(getEnergy(), o.getEnergy());
		}
		
		public int getDepth() {
			return depth;
		}

		public double getEnergy() {
			return energy;
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
		
		public String[] getOriginalSeqs(){
			return seqs;
		}

		public String getClassification() {
			return classification;
		}

		public double getPredictionScore() {
			return predictionScore;
		}
		
		public double getThreePScore(){
			return threePScore;
		}
		
		public double getFivePScore(){
			return fivePScore;
		}
		
		public String getThreePContig(){
			return sps[0].getContig();
		}
		
		public String getFivePContig(){
			return sps[1].getContig();
		}
		
		public boolean getThreePStrand(){
			return sps[0].isPlusStrand();
		}
		
		public boolean getFivePStrand(){
			return sps[1].isPlusStrand();
		}
		
		public int getThreePPosition(){
			return sps[0].getThreePposition();
		}
		
		public int getFivePPosition(){
			return sps[1].getFivePposition();
		}
		
		public double getThreePEnergy(){
			return sps[0].getEnergy();
		}
		
		public double getFivePEnergy(){
			return sps[1].getEnergy();
		}
		
		public String getMisc(){
			return misc;
		}
		
		static public String getHeader(){
			return "Contig1\tStrand1\t5pPosition1\t3pPosition1\t5p5preads1\t5p3preads1\t3p5preads1\t3p3preads1\tAccession1\tGeneName1\tGenomicRegion5p1\tGenomicRegion3p1\tDistFromExon5p1\tDistFromExon3p1\tDepth1\tEnergy1\tHairpinNumber1\tOverhang1\t"
					+ "Contig2\tStrand2\t5pPosition2\t3pPosition2\t5p5preads2\t5p3preads2\t3p5preads2\t3p3preads2\tAccession2\tGeneName2\tGenomicRegion5p2\tGenomicRegion3p2\tDistFromExon5p2\tDistFromExon3p2\tDepth2\tEnergy2\tHairpinNumber2\tOverhang2\t"
					+ "Class\tPredictionScore\tNumMatched3p\tNumMatched5p\tDepth\tEnergy\tHairpinNumber\tLeftPaired\tRightPaired\tOverhang\t5pScore\t3pScore\tblastHit1\tblastHit2"
					+ "\tSeqEntropy\tStructureEntropy\tHetero\tSeq1\tSeq2\tSeq";
		}
		//5p5preads\t5p3preads\t3p5preads\t3p3preads\t
		public ScoredPair(ScoredPosition p1, ScoredPosition p2, int seqLength){ 
			sps = new ScoredPosition[2];
			sps[0] = p1;
			sps[1] = p2;
			String seq1 = getTruncatedSequence(p1, seqLength, true);
			String seq2 = getTruncatedSequence(p2, seqLength, false);
			seq = seq1 + seq2;
			RNAcofoldLauncher fold = new RNAcofoldLauncher(seq, FCLIP_Scorer.flankingNTNumber);
			isFeasibleFold = fold.isFeasibleFold();
			this.depth = fold.getDepth();
			this.energy = fold.getEnergyPerNT();
			this.leftPaired = fold.getLeftPaired();
			this.rightPaired = fold.getRightPaired();
			this.overHang = fold.getOverHang();
			this.hairpinNum = fold.getNumberOfHairPins();
			this.seqEntropy = fold.getSeqEntropy();
			this.structureEntropy = fold.getStructureEntropy();
		//	this.seedConservation = sps[0].getThreePScore()>=0?sps[0].getSeedConservation() : sps[1].getSeedConservation();
			this.threePScore = sps[0].getThreePScore()>=0?sps[0].getThreePScore() : sps[1].getThreePScore();
			this.fivePScore = sps[0].getFivePScore()>=0?sps[0].getFivePScore() : sps[1].getFivePScore();
			this.blatHits = new int[2];
			blatHits[0] = sps[0].getBlatHits();
			blatHits[1] = sps[1].getBlatHits();
			this.seqs = new String[2];
			seqs[0] = sps[0].getSeq();
			seqs[1] = sps[1].getSeq();
			
		//	int r1 = sps[0].getFivePReads()[0] + sps[0].getFivePReads()[1] + sps[0].getThreePReads()[0] + sps[0].getThreePReads()[1];
		//	int r2 = sps[1].getFivePReads()[0] + sps[1].getFivePReads()[1] + sps[1].getThreePReads()[0] + sps[1].getThreePReads()[1];
			this.hetero = -1;//(sps[0].getHetero() * r1 + sps[1].getHetero() * r2) / (r1 + r2);
		}
		
		public boolean isFeasibleFold() { return isFeasibleFold; }
		
		public ScoredPair(String s){//16
			String[] token = s.split("\t");
			String[] ps = new String[2];
			ps[0] = new String();
			ps[1] = new String();
			int i = 0;
			for(;i<18;i++){
				ps[0] += token[i] + "\t";
			}
			for(;i<36;i++){
				ps[1] += token[i] + "\t";
			}
			this.sps = new ScoredPosition[2];
			for(int j=0;j<2;j++){
				sps[j] = ScoredPosition.getScoredPositionForScoredPair(ps[j]);
			}
			this.classification = token[i++];
			this.predictionScore = Double.parseDouble(token[i++]);
			this.matchedNum = new int[2];
			this.matchedNum[0] =  Integer.parseInt(token[i++]);
			this.matchedNum[1] =  Integer.parseInt(token[i++]);
			//this.sci = Double.parseDouble(token[i++]);
			//this.zscore = Double.parseDouble(token[i++]);
			//this.dnds = Double.parseDouble(token[i++]);
			//this.phyloPConservationPvalues[0] = Double.parseDouble(token[i++]);
			//this.phyloPConservationPvalues[1] = Double.parseDouble(token[i++]);
			//this.phyloPAccPvalues[0] = Double.parseDouble(token[i++]);
			//this.phyloPAccPvalues[1] = Double.parseDouble(token[i++]);
			//this.seedConservation = Double.parseDouble(token[i++]);
			this.depth = Integer.parseInt(token[i++]);
			this.energy = Double.parseDouble(token[i++]);
			this.hairpinNum = Integer.parseInt(token[i++]);
			this.leftPaired = Integer.parseInt(token[i++]);
			this.rightPaired = Integer.parseInt(token[i++]);
			this.overHang = Integer.parseInt(token[i++]);
			this.threePScore = Double.parseDouble(token[i++]);
			this.fivePScore = Double.parseDouble(token[i++]);
			this.blatHits = new int[2];
			this.blatHits[0] = Integer.parseInt(token[i++]);
			this.blatHits[1] = Integer.parseInt(token[i++]);
			this.seqEntropy = Double.parseDouble(token[i++]);
			this.structureEntropy = Double.parseDouble(token[i++]);
			this.hetero = Double.parseDouble(token[i++]);
			this.seqs = new String[2];
			this.seqs[0] = token[i++].replaceAll(" ", "");
			sps[0].setSeq(this.seqs[0]);
			this.seqs[1] = token[i++].replaceAll(" ", "");
			sps[1].setSeq(this.seqs[1]);
			this.seq = token[i++].replaceAll(" ", "");
			
			for(;i<token.length;i++){
				this.misc += token[i] + (i<token.length-1 ? "\t" : "");
			}
		}
		
		public ScoredPosition[] getPairedScoredPositions(){
			return sps;
		}
		
		public void setMatched3pNum(int n){
			matchedNum[0] = n;
		}
		
		public void setMatched5pNum(int n){
			matchedNum[1] = n;
		}
		
		public int getMatched3pNum() { return matchedNum[0]; }
		public int getMatched5pNum() { return matchedNum[1]; }
		
		private String getTruncatedSequence(ScoredPosition p, int seqLength, boolean isThreePside){
			String seq = p.getSeq();
			if(isThreePside) seq = seq.substring(0, Math.min(seq.length(), seqLength/2));
			else seq = seq.substring(Math.max(0, seq.length() - seqLength/2));
			return seq;
		}
		
		public ArrayList<String> getGenomicRegions3p(){
			return sps[0].getGenomicRegions3p();
		}
		
		public ArrayList<String> getGenomicRegions5p(){
			return sps[1].getGenomicRegions5p();
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<2;i++){
				sb.append(sps[i].toSimpleString());
				/*sb.append(sps[i].getContig()); sb.append('\t');
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
				ScoringOutputParser.ScoredPosition.writeStringArray(sb, sps[i].getGenomicRegions5p());*/
				sb.append('\t');
			}
			sb.append(classification); sb.append('\t');
			sb.append(predictionScore); sb.append('\t');
			sb.append(matchedNum[0]); sb.append('\t');
			sb.append(matchedNum[1]); sb.append('\t');
	//		sb.append(sci); sb.append('\t');
	//		sb.append(zscore); sb.append('\t');
	//		sb.append(dnds); sb.append('\t');
	//		sb.append(phyloPConservationPvalues[0]); sb.append('\t');
	//		sb.append(phyloPConservationPvalues[1]); sb.append('\t');
	//		sb.append(phyloPAccPvalues[0]); sb.append('\t');
	//		sb.append(phyloPAccPvalues[1]); sb.append('\t');
			
			
	//		sb.append(seedConservation); sb.append('\t'); 
			
			
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNum); sb.append('\t');
			sb.append(leftPaired); sb.append('\t');
			sb.append(rightPaired); sb.append('\t');
			sb.append(overHang); sb.append('\t');
			sb.append(threePScore); sb.append('\t');
			sb.append(fivePScore); sb.append('\t'); 
			
			sb.append(blatHits[0]); sb.append('\t');
			sb.append(blatHits[1]); sb.append('\t');
	//		sb.append(sps[0].getReadLengthVariance()); sb.append('\t');
	//		sb.append(sps[1].getReadLengthVariance()); sb.append('\t');		
			sb.append(seqEntropy); sb.append('\t');
			sb.append(structureEntropy); sb.append('\t');
			sb.append(hetero); sb.append('\t');
			//sb.append(fivePhetero); sb.append('\t');
			sb.append(ScoredPosition.getCleavedSequence(seqs[0], FCLIP_Scorer.flankingNTNumber));  sb.append('\t');
			sb.append(ScoredPosition.getCleavedSequence(seqs[1], FCLIP_Scorer.flankingNTNumber));  sb.append('\t');
			
			sb.append(ScoredPosition.getCleavedSequence(seq, FCLIP_Scorer.flankingNTNumber)); 
			
			if(!misc.isEmpty()){
				sb.append('\t');
				sb.append(misc);
			}
			
			//\tDepth1\tDepth2\tEnergy1\tEnergy2\tHairpinNumber1\tHairpinNumber2\tSequence1\tSequence2
			return sb.toString();
		}	
		
		public void setClassification(String c, double predictionScore) {
			classification = c;
			this.predictionScore = predictionScore;	
		}
	
		
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Arrays.hashCode(sps);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			ScoredPair other = (ScoredPair) obj;
			if (!Arrays.equals(sps, other.sps))
				return false;
			return true;
		}

	/*	public void setRNAzScores(MafParser mafParser){
			RNAzLauncher rna = getRNAzLauncher(mafParser);
			this.sci = rna.getSCI();
			this.zscore = rna.getZScore();
		}
		
		public void setDnDs(MafParser mafParser){
			String[] contigs = new String[2];
			int[] positions = new int[2];
			boolean[] strands = new boolean[2];
			
			int length = 20;//seq.length() / 2 - Scorer.flankingNTNumber;
			
			for(int i=0;i<2;i++){
				contigs[i] = sps[i].getContig();
				strands[i] = sps[i].isPlusStrand();
				positions[i] = 
					//	i == 0? sps[i].getThreePposition() // + (strands[i]? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
					//			: sps[i].getFivePposition() + (strands[i]? - length + 1 : length - 1);	
						i == 0 ? sps[i].getThreePposition() + (strands[i]?  - length/2: length/2)
								: sps[i].getFivePposition() + (strands[i]?  - length/2 + 1: length/2 - 1);
						
			}
			
			String[] seqs = mafParser.getSeqs(contigs, positions, strands, length);
	//		for(String se : seqs){
	//			System.out.println(se);
	//		}
			this.dnds = DnDsCalculator.calculate(seqs);
		}
		
		public void setPhyloPPvalues(MafParser mafParser, String modFolder){
			int length = 20;//seq.length() / 2 - Scorer.flankingNTNumber;
			phyloPConservationPvalues = new double[2];
			phyloPAccPvalues = new double[2];
			//System.out.println(ScoredPosition.getCleavedSequence(seq, Scorer.flankingNTNumber));
			for(int i=0;i<2;i++){
				String contig = sps[i].getContig();
				boolean isPlusStrand = sps[i].isPlusStrand();
				int position =	i == 0 ? sps[i].getThreePposition() + (isPlusStrand?  - length/2: length/2)
						: sps[i].getFivePposition() + (isPlusStrand?  - length/2 + 1: length/2 - 1);			
				String s = mafParser.getSeqsInMafFormat(contig, position, isPlusStrand, length);
				PhyloPLauncher.setModFile( modFolder + contig + ".maf.mod");
				//if(i == 0 && !isPlusStrand) System.out.println(s);
				//for(String tt : mafParser.getSeqs(contig, position, isPlusStrand, length))
				//	System.out.println(tt);
				PhyloPLauncher pp = new PhyloPLauncher(s);
				phyloPConservationPvalues[i] = pp.getPvalConservation(); 
				phyloPAccPvalues[i] = pp.getPvalAcceleration();
			}	
		}*/
				
		private String getMafString(MafParser mafParser){
			String[] contigs = new String[2];
			int[] positions = new int[2];
			boolean[] strands = new boolean[2];
			
			int length = seq.length() / 2 - FCLIP_Scorer.flankingNTNumber;
			
			for(int i=0;i<2;i++){
				contigs[i] = sps[i].getContig();
				strands[i] = sps[i].isPlusStrand();
				positions[i] = 
						i == 0 ? sps[i].getThreePposition()// + (strands[i]? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
								: sps[i].getFivePposition() + (strands[i]?  - length + 1: length - 1);			
						
			}			
		//	System.out.println(ScoredPosition.getCleavedSequence(this.seq, Scorer.flankingNTNumber));
		//	System.out.println(mafParser.getSeqsInMafFormat(contigs, positions, strands, length));			
			return mafParser.getSeqsInMafFormat(contigs, positions, strands, length);
		}
	

			
		
	}
	
	
	public ScoredPairOutputParser(String inFile){
		read(inFile);
	}
	
	private ArrayList<ScoredPair> pairs;
	private String header;
	
	private void read(String outFile){
		pairs = new ArrayList<ScoredPair>();
		try {
			BufferedLineReader in  = new BufferedLineReader(outFile);
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("Contig")){
					header = s;
					continue;
				}
				ScoredPair position = new ScoredPair(s);				
				pairs.add(position);
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public String getHeader() { return header; }
	
	public ArrayList<ScoredPair> getPairs(){
		return pairs;
	}
	
	static public void generateFastaForMotif(String csv, String fasta5p, String fasta3p){
		PrintStream out5p, out3p;
		ScoredPairOutputParser parser = new ScoredPairOutputParser(csv);
		try {
			out5p = new PrintStream(fasta5p);
			out3p = new PrintStream(fasta3p);
			int i = 0;
			HashSet<String> fseqs = new HashSet<String>();
			HashSet<String> tseqs = new HashSet<String>();
			//int k = 0;
			for(ScoredPair pair : parser.getPairs()){
				String seq = ScoredPosition.getCleavedSequence(pair.getSeq(), FCLIP_Scorer.flankingNTNumber);
				//if(pair.getClassification().equals("M")){
					String fseq = seq.substring(0, seq.indexOf(' '));
					if(!fseqs.contains(pair.getOriginalSeqs()[0])){
						fseqs.add(pair.getOriginalSeqs()[0]);
						out5p.println(">"+ i);
						out5p.println(fseq);
						i++;
					}
					
					String tseq = seq.substring(seq.lastIndexOf(' ') + 1);					
					if(!tseqs.contains(pair.getOriginalSeqs()[1])){
						tseqs.add(pair.getOriginalSeqs()[1]);
						out3p.println(">"+ i);
						out3p.println(tseq);
						i++;
					}
				//}
			}
			out5p.close();
			out3p.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	static public void generateMfileForScatterHist(String csv, String mFile){
		PrintStream outM;
		ScoredPairOutputParser parser = new ScoredPairOutputParser(csv);
		try {
			outM = new PrintStream(mFile);
			outM.print("value=[");
			for(ScoredPair pair : parser.getPairs()){
				outM.println(pair.getOverHang() + " " + pair.getEnergy() + " " + pair.getThreePEnergy() + " " + pair.getFivePEnergy());
			}
			outM.println("];");
			
			outM.print("group={");
			for(ScoredPair pair : parser.getPairs()){
				outM.println("'" + pair.getClassification() + "'");
			}
			outM.println("};");
			
			outM.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public static void main(String[] args){
		String file = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/merged.trans.control.M.csv";
		ScoredPairOutputParser parser = new ScoredPairOutputParser(file);
		HashMap<ScoredPosition, Integer> map = new HashMap<ScoredPosition, Integer>();
		for(ScoredPair pair : parser.getPairs()){
			map.put(pair.getPairedScoredPositions()[0], pair.getMatched3pNum());
		}
		for(int v : map.values()){
			System.out.println(v);
		}
	}
}
