package fCLIP.parser;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.Deflater;

import launcher.RNAcofoldLauncher;
import launcher.RNAfoldLauncher;
import launcher.ShellLauncher;
import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.MirGff3FileParser.MiRNA;
import parser.Bed12Parser;
import parser.BufferedLineReader;
import parser.MafParser;
import parser.ZeroBasedFastaParser;
import fCLIP.Classifier;
import fCLIP.FCLIP_Scorer;

public class ScoredPositionOutputParser {
	public static class ScoredPosition implements Comparable<ScoredPosition>{
		private String contig;		
		private boolean isPlusStrand;
		private int fivePposition = -1;
		private int threePposition = -1;
		//private ArrayList<Integer> fivePcoordinate;
		//private ArrayList<Integer> threePcoordinate;
		
		//private double sci = 0;
		//private double zscore = 1;
		//private double seedConservation = 1;
		//private double dnds = -1;
		private boolean isFeasibleFold = true;
		private double fivePScore = -1;
		private double threePScore = -1;
		private int depth;
		private double energy;
		private int hairpinNumInStem;
		private int leftPaired;
		private int rightPaired;
		private int overHang;
		private int preLength;
		private String seq;
		private ArrayList<String> miRNAs = null;
		private ArrayList<String> containingGeneAccessions;
		private ArrayList<String> containingGeneNames;
		private ArrayList<String> genomicRegions5p;
		private ArrayList<String> genomicRegions3p;
		private ArrayList<Integer> distToNextExon5p;
		private ArrayList<Integer> distToNextExon3p;
		//private double complexity = 0;
		private String classification;
		private double predictionScore;
		private int[] threePreads;
		private int[] fivePreads;
		private int blatHits = 0;
		private double seqEntropy = 0;
		private double structureEntropy = 0;
		private double hetero = 0;
		private String misc = "";
		private boolean isALU3p = false;
		private boolean isALU5p = false;
		private boolean isRepeat3p = false;
		private boolean isRepeat5p = false;
		
		private int longestGeneIndex = -1;
		
		static int rmsk3pHeaderColumnNumber = -1;
		static int rmsk5pHeaderColumnNumber = -1;
		
		public void setLongestGeneIndex(AnnotationFileParser aparser){
			if(getContainingGeneAccessions() == null) return;
			int maxSize = 0;
			for(int i=0;i<getContainingGeneAccessions().size();i++){
				String acc = getContainingGeneAccessions().get(i);
				AnnotatedGene gene = aparser.getGeneByAccession(acc);
				int size = gene.getExonSize();
				if(maxSize < size){
					maxSize = size;
					longestGeneIndex = i;
				}
			}			
		}
			
		public int getLongestGeneIndex(){ return longestGeneIndex; }
		
		public static ScoredPosition get3p5pSwappedPosition(ScoredPosition p, int flankingNTNumber){
			ScoredPosition sp = new ScoredPosition();
			sp.contig = p.contig;
			sp.isPlusStrand = p.isPlusStrand;
			sp.fivePposition = p.threePposition;
			sp.threePposition = p.fivePposition;
		//	sp.sci = p.sci;
		//	sp.zscore = p.zscore;
		//	sp.seedConservation = p.seedConservation;
		//	sp.dnds = p.dnds;
			sp.fivePScore = p.threePScore;
			sp.threePScore = p.fivePScore;
			sp.fivePreads = p.threePreads;
			sp.threePreads = p.fivePreads;
			String seqw = "";
			for(char s : p.seq.toCharArray()){
				seqw = s + seqw;
			}
			sp.seq = seqw;
			RNAcofoldLauncher fold = new RNAcofoldLauncher(seqw, flankingNTNumber);
			sp.set(fold);		
			sp.miRNAs = p.miRNAs;
			sp.containingGeneAccessions = p.containingGeneAccessions;
			sp.containingGeneNames = p.containingGeneNames;
			sp.genomicRegions3p = p.genomicRegions5p;
			sp.genomicRegions5p = p.genomicRegions3p;
			sp.distToNextExon5p = p.distToNextExon5p;
			sp.distToNextExon3p = p.distToNextExon3p;
			
			sp.classification = p.classification;
			sp.predictionScore = 0;
			sp.blatHits = p.blatHits;
			sp.hetero = p.hetero;
			return sp;
		}
		//private int readCount = 0;
		//private double readLengthVariance = -1;
		
		public static String getHeader(){
			return "Contig\t5p\t3p\tStrand\tPreLength\tPaired\tClass\tPredictionScore\t5p5preads\t5p3preads\t3p5preads\t"
					+ "3p3preads\tMatching miRNAs\tBlatHits\tAccession\tGeneName\tGenomicRegion5p\tGenomicRegion3p\tDistToNextExon5p"
					+ "\tDistToNextExon3p\t5pSocre\t3pScore\tDepth\tEnergy\tHairpinNumInStem\tLeftPaired\tRightPaired\tOverhang\tSeqEntropy\tStructureEntropy\tHetero\tSeq";
		}
		
		public ArrayList<String> getContainingGeneAccessions() {
			return containingGeneAccessions == null? new ArrayList<String>() : containingGeneAccessions;
		}

		public ArrayList<String> getContainingGeneNames() {
			return containingGeneNames;
		}
		
		public boolean isALU3p() {
			return isALU3p;
		}

		public boolean isALU5p() {
			return isALU5p;
		}
		
		public boolean isRepeat3p(){
			return isRepeat3p;
		}
		
		public boolean isRepeat5p(){
			return isRepeat5p;
		}
		
		public double getSeqEntropy(){
			return seqEntropy;
		}
		
		public double getStructureEntropy(){
			return structureEntropy;
		}
		
		public int getFivePPosition() {
			return fivePposition;
		}

		public int getThreePPosition() {
			return threePposition;
		}

		public double getHetero(){
			return hetero;
		}
				
		public int getOverHang(){
			return overHang;
		}
		
		public int getPreLength(){
			return preLength;
		}
		
		public ArrayList<String> getGenomicRegions5p(){
			return genomicRegions5p;
		}
		
		public ArrayList<String> getGenomicRegions3p(){
			return genomicRegions3p;
		}
		
		public ArrayList<Integer> getDistToNextExon5p(){
			return distToNextExon5p;
		}
		
		public ArrayList<Integer> getDistToNextExon3p(){
			return distToNextExon3p;
		}
		//public double getReadLengthVariance(){
		//	return readLengthVariance;
		//}
		
		public int getDepth() {
			return depth;
		}
			
		public double getEnergy() {
			return energy;
		}
		
		public double getPredictionScore(){
			return predictionScore;
		}
		
		public int getHairpinNumberInStem(){
			return hairpinNumInStem;
		}
		
		public ArrayList<String> getMiRNAs(){
			return miRNAs;
		}
		
		public boolean hasMatchingMiRNA(){
			return miRNAs != null && !miRNAs.isEmpty();
		}
		
		public boolean isPlusStrand(){ return isPlusStrand;}
		
		public boolean is3pScored(){
			return threePScore >= 0;
		}
	
		public boolean is5pScored(){
			return fivePScore >= 0;
		}
		
		public boolean isPaired(){
			return is3pScored() && is5pScored();
		}
		
		public String getClassification(){
			return classification;
		}
		
		public int getBlatHits(){
			return blatHits;
		}
		
		public String getMiscInfo(){
			return misc;
		}
		//public int getReadCount(){
		//	return readCount;
		//}
		private ScoredPosition(){
			
		}
		
		
		
		public ScoredPosition(String s){
			String[] token = s.split("\t");
			int i = 0;
			this.contig = token[i++];
			this.threePposition = Integer.parseInt(token[i++]);
			this.fivePposition = Integer.parseInt(token[i++]);
			this.isPlusStrand = token[i++].equals("+");
			this.preLength = Integer.parseInt(token[i++]);			
			i++;
			this.classification = token[i++];
			this.predictionScore = Double.parseDouble(token[i++]);
			this.threePreads = new int[2];
			this.threePreads[0] = Integer.parseInt(token[i++]);
			this.threePreads[1] = Integer.parseInt(token[i++]);
			this.fivePreads = new int[2];
			this.fivePreads[0] = Integer.parseInt(token[i++]);
			this.fivePreads[1] = Integer.parseInt(token[i++]);
			
			//this.readCount = Integer.parseInt(token[i++]);
		//	this.readLengthVariance = Double.parseDouble(token[i++]);
			if(!token[i].startsWith("_")){
				String[] miRNAsString = token[i].split(";");
				this.miRNAs = new ArrayList<String>();
				for(String st : miRNAsString){
					miRNAs.add(st);
				}
			}
			i++;
			this.blatHits = Integer.parseInt(token[i++]);
		//	this.sci = Double.parseDouble(token[i++]);
		//	this.zscore = Double.parseDouble(token[i++]);
		//	this.dnds = Double.parseDouble(token[i++]);
		//	this.seedConservation  = Double.parseDouble(token[i++]);
			
			if(!token[i].startsWith("_")){
				String[] accessions = token[i].split(";");
				this.containingGeneAccessions = new ArrayList<String>();
				for(String st : accessions){
					containingGeneAccessions.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] names = token[i].split(";");
				this.containingGeneNames = new ArrayList<String>();
				for(String st : names){
					containingGeneNames.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] regions = token[i].split(";");
				this.genomicRegions3p = new ArrayList<String>();
				for(String st : regions){
					genomicRegions3p.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] regions = token[i].split(";");
				this.genomicRegions5p = new ArrayList<String>();
				for(String st : regions){
					genomicRegions5p.add(st);
				}
			}
			i++;
			if(!token[i].startsWith(" ")){
				String[] dists = token[i].split(";");
				this.distToNextExon3p = new ArrayList<Integer>();
				for(String d : dists){
					if(d.equals(" ")) distToNextExon3p.add(null);
					else distToNextExon3p.add(Integer.parseInt(d));
				}
			}
			i++;
			if(!token[i].startsWith(" ")){
				String[] dists = token[i].split(";");
				this.distToNextExon5p = new ArrayList<Integer>();
				for(String d : dists){
					if(d.equals(" ")) distToNextExon5p.add(null);
					else distToNextExon5p.add(Integer.parseInt(d));
				}
			}
			i++;
			
			
			/*
			if(!token[i].isEmpty()){
				String[] mrStarts = token[i++].split(",");
				String[] mrEnds = token[i++].split(",");
				ArrayList<ArrayList<Integer>> mappedRegionStartsEnds = new ArrayList<ArrayList<Integer>>();
				for(int j=0;j<mrStarts.length;j++){
					ArrayList<Integer> splice = new ArrayList<Integer>();
					splice.add(Integer.parseInt(mrStarts[j]));
					splice.add(Integer.parseInt(mrEnds[j]));
					mappedRegionStartsEnds.add(splice);
				}
				this.threePcoordinate = Bed12Parser.getCoordinate(mappedRegionStartsEnds, isPlusStrand);
			}else i+=2;
			
			if(!token[i].isEmpty()){
				String[] mrStarts = token[i++].split(",");
				String[] mrEnds = token[i++].split(",");
				ArrayList<ArrayList<Integer>> mappedRegionStartsEnds = new ArrayList<ArrayList<Integer>>();
				for(int j=0;j<mrStarts.length;j++){
					ArrayList<Integer> splice = new ArrayList<Integer>();
					splice.add(Integer.parseInt(mrStarts[j]));
					splice.add(Integer.parseInt(mrEnds[j]));
					mappedRegionStartsEnds.add(splice);
				}
				this.fivePcoordinate = Bed12Parser.getCoordinate(mappedRegionStartsEnds, isPlusStrand);
			}else i+=2;
			*/
			
			this.threePScore = Double.parseDouble(token[i++]);
			this.fivePScore = Double.parseDouble(token[i++]);
			this.depth = Integer.parseInt(token[i++]);
			this.energy = Double.parseDouble(token[i++]);
			this.hairpinNumInStem = Integer.parseInt(token[i++]);
			this.leftPaired = Integer.parseInt(token[i++]);
			this.rightPaired = Integer.parseInt(token[i++]);
			this.overHang = Integer.parseInt(token[i++]);
			this.seqEntropy = Double.parseDouble(token[i++]);
			this.structureEntropy = Double.parseDouble(token[i++]);
			this.hetero = Double.parseDouble(token[i++]);
			this.seq = token[i++].replaceAll(" ", "");
			
			for(;i<token.length;i++){
				this.misc += token[i] + (i<token.length-1 ? "\t" : "");
			}
			
			if(rmsk3pHeaderColumnNumber> 0){
				this.isALU3p = token[rmsk3pHeaderColumnNumber].toUpperCase().contains("ALU");
				this.isRepeat3p = !token[rmsk3pHeaderColumnNumber].equals("_");// toUpperCase().contains("ALU");
			}
			
			if(rmsk5pHeaderColumnNumber> 0){
				this.isALU5p = token[rmsk5pHeaderColumnNumber].toUpperCase().contains("ALU");
				this.isRepeat5p = !token[rmsk5pHeaderColumnNumber].equals("_");
			}
		}
		
		public static void writeStringArray(StringBuilder sb, ArrayList<?> array, String emptyCharacter){
			if(array != null && !array.isEmpty()){
				//System.out.println(array);
				for(int i=0; i<array.size();i++){
					sb.append(array.get(i) == null ? emptyCharacter : array.get(i).toString()); 
					sb.append(i == array.size()-1 ? "\t" : ";");					
				}			
			}else sb.append(emptyCharacter + "\t");
		}
		
		
//		public RNAzLauncher getRNAzLauncher(MafParser mafParser){
//			String mafSeq = getMafString(mafParser);
//			return new RNAzLauncher(mafSeq);
//		}
		
//		public void setRNAzScoresNSeedConservation(MafParser mafParser){
//		//	String mafSeq = getMafString(mafParser);
//			RNAzLauncher rna = getRNAzLauncher(mafParser);
//			this.sci = rna.getSCI();
//			this.zscore = rna.getZScore();
//			this.seedConservation = new PhyloPLauncher(getSeedMafString(mafParser)).getPvalConservation();
//		}
		
	//	public void setSeedConvervation(MafParser mafParser){
	//		String mafString = mafParser.getSeqsInMafFormat(contig, position.getCoordinate(), position.getPosition(), isPlusStrand, positionQuantityChangeLength);
	//		this.phyloP = new PhyloPLauncher(mafString).getPvalConservation();
	//	}
		
		
	/*	public void setDnDs(MafParser mafParser){
			int length = seq.length() - 2 * Scorer.flankingNTNumber;
			int position = this.threePposition; 
					//+ (isPlusStrand? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) 
				//	: 
				//this.fivePposition + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
			if(position < 0 || this.fivePposition < 0) return;
			String[] seqs = mafParser.getSeqs(contig, position, isPlusStrand, length);
			if(seqs.length < 2) return;
			//
			//for(String seq : seqs) System.out.println(seq);
			
			this.dnds = DnDsCalculator.calculate(seqs);
			
		}
		*/
		private String getMafString(MafParser mafParser){
			int length = seq.length() - 2 * FCLIP_Scorer.getFlankingNTNumber();
			int position = this.threePposition;// > 0 ? this.threePposition 
			if(position < 0 || this.fivePposition < 0) return null;
			//		: this.fivePposition;// + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
			return mafParser.getSeqsInMafFormat(contig, position, isPlusStrand, length);
		} 
		
		private String getSeedMafString(MafParser mafParser){
			int length = Math.min(7, seq.length() - 2 * FCLIP_Scorer.getFlankingNTNumber());// ;
			int position = this.threePposition;
			if(position < 0 || this.fivePposition < 0) return null;
					//: this.fivePposition;// + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
			return mafParser.getSeqsInMafFormat(contig, position + (isPlusStrand? 1 : -1), isPlusStrand, length);
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(contig); sb.append('\t');
			sb.append(threePposition); sb.append('\t');
			sb.append(fivePposition); sb.append('\t');
			sb.append(isPlusStrand? '+': '-'); sb.append('\t');
			sb.append(preLength); sb.append('\t');
			sb.append(isPaired() ? "T" : "F" ); sb.append('\t');
			sb.append(classification); sb.append('\t');
			sb.append(predictionScore); sb.append('\t');
			sb.append(threePreads == null? '_' : threePreads[0]); sb.append('\t');
			sb.append(threePreads == null? '_' : threePreads[1]); sb.append('\t');
			
			sb.append(fivePreads == null? '_' : fivePreads[0]); sb.append('\t');
			sb.append(fivePreads == null? '_' : fivePreads[1]); sb.append('\t');
			
		//	sb.append(readCount); sb.append('\t');
	//		sb.append(readLengthVariance); sb.append('\t');
			
			writeStringArray(sb, miRNAs, "_");
			sb.append(blatHits); sb.append('\t');
		//	sb.append(sci); sb.append('\t');
		//	sb.append(zscore); sb.append('\t');
		//	sb.append(dnds); sb.append('\t');
		//	sb.append(seedConservation); sb.append('\t');
			writeStringArray(sb, containingGeneAccessions, "_");
			writeStringArray(sb, containingGeneNames, "_");
			writeStringArray(sb, genomicRegions3p, "_");
			writeStringArray(sb, genomicRegions5p, "_");
			writeStringArray(sb, distToNextExon3p, " ");
			writeStringArray(sb, distToNextExon5p, " ");
			//System.out.println(sb);
			//toMappedRegionString(sb, threePcoordinate); sb.append('\t');
			//toMappedRegionString(sb, fivePcoordinate); sb.append('\t');
			sb.append(threePScore); sb.append('\t');
			sb.append(fivePScore); sb.append('\t');
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNumInStem); sb.append('\t');
			sb.append(leftPaired); sb.append('\t');
			sb.append(rightPaired); sb.append('\t');
			sb.append(overHang); sb.append('\t');
			sb.append(seqEntropy); sb.append('\t');
			sb.append(structureEntropy); sb.append('\t');
			sb.append(hetero); sb.append('\t');
			sb.append(getCleavedSequence(seq, FCLIP_Scorer.getFlankingNTNumber()));
			
			if(!misc.isEmpty()){
				sb.append('\t');
				sb.append(misc);
			}
			
			return sb.toString();
		}
		
		
		public String toSimpleString(){
			StringBuilder sb = new StringBuilder();
			sb.append(contig); sb.append('\t');
			sb.append(threePposition); sb.append('\t');
			sb.append(fivePposition); sb.append('\t');
			sb.append(isPlusStrand? '+': '-'); sb.append('\t');
			sb.append(preLength); sb.append('\t');
			sb.append(threePreads == null? '_' : threePreads[0]); sb.append('\t');
			sb.append(threePreads == null? '_' : threePreads[1]); sb.append('\t');
			
			sb.append(fivePreads == null? '_' : fivePreads[0]); sb.append('\t');
			sb.append(fivePreads == null? '_' : fivePreads[1]); sb.append('\t');
			
			writeStringArray(sb, containingGeneAccessions, "_");
			writeStringArray(sb, containingGeneNames, "_");
			writeStringArray(sb, genomicRegions3p, "_");
			writeStringArray(sb, genomicRegions5p, "_");
			writeStringArray(sb, distToNextExon3p, " ");
			writeStringArray(sb, distToNextExon5p, " ");
			
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNumInStem); sb.append('\t');
			
			sb.append(overHang);
			return sb.toString();
		}
		
		
		static ScoredPosition getScoredPositionForScoredPair(String s){ 
			ScoredPosition p = new ScoredPosition();
			String[] token = s.split("\t");
			int i = 0;
			p.contig = token[i++];
			p.threePposition = Integer.parseInt(token[i++]);
			p.fivePposition = Integer.parseInt(token[i++]);
			p.isPlusStrand = token[i++].equals("+");
			p.preLength =  Integer.parseInt(token[i++]);
			p.threePreads = new int[2];
			p.threePreads[0] = Integer.parseInt(token[i++]);
			p.threePreads[1] = Integer.parseInt(token[i++]);
			p.fivePreads = new int[2];
			p.fivePreads[0] = Integer.parseInt(token[i++]);
			p.fivePreads[1] = Integer.parseInt(token[i++]);
			
			if(!token[i].startsWith("_")){
				String[] accessions = token[i].split(";");
				p.containingGeneAccessions = new ArrayList<String>();
				for(String st : accessions){
					p.containingGeneAccessions.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] names = token[i].split(";");
				p.containingGeneNames = new ArrayList<String>();
				for(String st : names){
					p.containingGeneNames.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] regions = token[i].split(";");
				p.genomicRegions3p = new ArrayList<String>();
				for(String st : regions){
					p.genomicRegions3p.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] regions = token[i].split(";");
				p.genomicRegions5p = new ArrayList<String>();
				for(String st : regions){
					p.genomicRegions5p.add(st);
				}
			}
			i++;
			if(!token[i].startsWith(" ")){
				String[] dists = token[i].split(";");
				p.distToNextExon3p = new ArrayList<Integer>();
				for(String d : dists){
					if(d.equals(" ")) p.distToNextExon3p.add(null);
					else p.distToNextExon3p.add(Integer.parseInt(d));
				}
			}
			i++;
			if(!token[i].startsWith(" ")){
				String[] dists = token[i].split(";");
				p.distToNextExon5p = new ArrayList<Integer>();
				for(String d : dists){
					if(d.equals(" ")) p.distToNextExon5p.add(null);
					else p.distToNextExon5p.add(Integer.parseInt(d));
				}
			}
			i++;
			p.depth = Integer.parseInt(token[i++]);
			p.energy = Double.parseDouble(token[i++]);
			p.hairpinNumInStem = Integer.parseInt(token[i++]);
		
			p.overHang = Integer.parseInt(token[i++]);
			return p;
		}
		
		
		public static String getCleavedSequence(String sequence, int flankingNTNumber){
			StringBuilder sb = new StringBuilder();
			sb.append(sequence.subSequence(0, flankingNTNumber));
			sb.append(" ");
			sb.append(sequence.subSequence(flankingNTNumber, sequence.length() - flankingNTNumber));
			sb.append(" ");
			sb.append(sequence.subSequence(sequence.length() - flankingNTNumber, sequence.length()));
			return sb.toString();
		}
		
		public String toMappedRegionString(StringBuilder sb, ArrayList<Integer> coordinate){
			//StringBuilder sb = new StringBuilder();
			if(coordinate!=null){
				ArrayList<ArrayList<Integer>> mappedRegionStartsEnds = Bed12Parser.getCoordinateStartsEnds(coordinate, isPlusStrand);
				for(int j=0;j<mappedRegionStartsEnds.size()-1;j++){
					sb.append(mappedRegionStartsEnds.get(j).get(0));
					sb.append(',');
				}
				sb.append(mappedRegionStartsEnds.get(mappedRegionStartsEnds.size()-1).get(0));
			
				sb.append('\t');
				for(int j=0;j<mappedRegionStartsEnds.size()-1;j++){
					sb.append(mappedRegionStartsEnds.get(j).get(1));
					sb.append(',');
				}
				sb.append(mappedRegionStartsEnds.get(mappedRegionStartsEnds.size()-1).get(1));
			}else sb.append('\t');
			//sb.append('\t');
			return sb.toString();
		}
		
		public ScoredPosition(String contig, boolean isPlusStrand, int fivePposition, 
				int threePposition, double fivePscore, double threePscore, AnnotationFileParser annotationParser){
			this.contig = contig;
			this.isPlusStrand = isPlusStrand;
			this.fivePposition = fivePposition;
			this.threePposition = threePposition;
			this.fivePScore = fivePscore;
			this.threePScore = threePscore;
			this.preLength = Math.abs(fivePposition - threePposition) + 1;
			setGeneInfo(annotationParser);
					
		}
		
		public void setGeneInfo(AnnotationFileParser annotationParser){
			setGeneInfo(annotationParser, 0);
		}
		
		public void setGeneInfo(AnnotationFileParser annotationParser, int offset){
			if(!isPlusStrand) offset = - offset;
			if(annotationParser != null){
				HashMap<AnnotatedGene, Integer> genes5p = new HashMap<AnnotatedGene, Integer>();
				HashMap<AnnotatedGene, Integer> genes3p = new HashMap<AnnotatedGene, Integer>();
				if(fivePposition >= 0){
					ArrayList<AnnotatedGene> ag = annotationParser.getContainingGenes(contig, fivePposition + offset);
					if(ag != null){
						for(AnnotatedGene a : ag){
							genes5p.put(a, fivePposition + offset);
						}
					}
				}
				if(threePposition>= 0){
					ArrayList<AnnotatedGene> ag = annotationParser.getContainingGenes(contig, threePposition - offset);
					if(ag != null){
						for(AnnotatedGene a : ag){
							genes3p.put(a, threePposition - offset);
						}
					}
				}
				this.containingGeneAccessions = new ArrayList<String>();
				this.containingGeneNames = new ArrayList<String>();
				this.genomicRegions3p = new ArrayList<String>();
				this.genomicRegions5p = new ArrayList<String>();
				this.distToNextExon3p = new ArrayList<Integer>();
				this.distToNextExon5p = new ArrayList<Integer>();
				HashSet<AnnotatedGene> genes = new HashSet<AnnotationFileParser.AnnotatedGene>();
				genes.addAll(genes5p.keySet());
				genes.addAll(genes3p.keySet());
				for(AnnotatedGene gene : genes){
					this.containingGeneAccessions.add(gene.getAccession());
					this.containingGeneNames.add(gene.getGeneName());
					this.genomicRegions3p.add(genes3p.containsKey(gene) ? annotationParser.getGenomicRegionNameAndFrameShift(contig, isPlusStrand, genes3p.get(gene), gene, true).get(0) : "_");
					this.genomicRegions5p.add(genes5p.containsKey(gene) ? annotationParser.getGenomicRegionNameAndFrameShift(contig, isPlusStrand, genes5p.get(gene), gene, true).get(0) : "_");
					this.distToNextExon3p.add(genes3p.containsKey(gene) ? gene.getDistToNextExon(threePposition - offset) : null);
					this.distToNextExon5p.add(genes5p.containsKey(gene) ? gene.getDistToNextExon(fivePposition + offset ) : null);
					//this.genomicRegions.add(annotationParser.getGenomicRegionNameAndFrameShift(contig, isPlusStrand, fivePScore >=0? fivePposition : threePposition , gene, true).get(0));
				}
			}	
		}
		
		
		public void addGenes(ScoredPosition other){
			if(this.containingGeneAccessions == null){
				this.containingGeneAccessions = new ArrayList<String>();
				this.containingGeneNames = new ArrayList<String>();
			}
			for(String a : other.containingGeneAccessions){
				if(this.containingGeneAccessions.contains(a)) continue;
				this.containingGeneAccessions.add(a);
			}
			for(String a : other.containingGeneNames){
				if(this.containingGeneNames.contains(a)) continue;
				this.containingGeneNames.add(a);
			}
			//setGeneInfo(annotationParser);
		}
		
		public ScoredPosition setClassification(String e, double predictionScore){
			this.classification = e; 
			this.predictionScore = predictionScore; 
			return this;
		}
		
		public ScoredPosition set(RNAcofoldLauncher fold){
			depth = fold.getDepth();
			energy = fold.getEnergyPerNT();
			hairpinNumInStem = fold.getNumberOfHairPinsInStem();
			leftPaired = fold.getLeftPaired();
			rightPaired = fold.getRightPaired();
			overHang = fold.getOverHang();
			seqEntropy = fold.getSeqEntropy();
			structureEntropy = fold.getStructureEntropy();
			isFeasibleFold = fold.isFeasibleFold();
			return this;
		}
		
		public boolean isFeasibleFold() { return isFeasibleFold;}
		
		
		
		public String getSeq() {
			return seq;
		}
		
	//	public void setReadCount(int r){
	//		this.readCount = r;
	//	}
		
		public void setBlatHits(int n){
			blatHits = n;
		}
		
		public void set3pReads(int t, int f){
			threePreads = new int[2];
			threePreads[0] = t;
			threePreads[1] = f;
		}
		
		public void set5pReads(int t, int f){
			fivePreads = new int[2];
			fivePreads[0] = t;
			fivePreads[1] = f;
		}
		
		public void setSeq(String seq) {
			//this.complexity = getComplexity(seq);
			this.seq = seq;
		}

		
		public void setMiRNAs(ArrayList<MiRNA> miRNAs) {
			this.miRNAs = new ArrayList<String>();
			for(MiRNA miRNA : miRNAs){
				this.miRNAs.add(miRNA.getName());
			}
		}

		
		public void setHetero(double h){
			hetero = h;
		}
		
		public ScoredPosition setFivePposition(int i, AnnotationFileParser annotationParser) {
			fivePposition = i;
			this.preLength = Math.abs(fivePposition - threePposition) + 1;
			setGeneInfo(annotationParser);
			return this;
		}

		public ScoredPosition setThreePposition(int i, AnnotationFileParser annotationParser) {
			threePposition = i;
			this.preLength = Math.abs(fivePposition - threePposition) + 1;
			setGeneInfo(annotationParser);
			return this;
		}


		public double getFivePScore() {
			return fivePScore;
		}

		public double getThreePScore() {
			return threePScore;
		}
		
		public int[] getThreePReads(){
			return threePreads;
		}

		public int[] getFivePReads(){
			return fivePreads;
		}
		
		public int getLeftPaired() {
			return leftPaired;
		}
		
		public int getRightPaired(){
			return rightPaired;
		}
			
		
		public boolean overlapped(ScoredPosition other){ // simple overlapping not considering splices
			if(!this.contig.equals(other.contig)) return false;
			int s = this.getThreePPosition();
			int l = this.getFivePPosition();
			if(s>l){
				int t = s;
				s = l;
				l = t;
			}
			
			int so = other.getThreePPosition();
			int lo = other.getFivePPosition();
			if(so>lo){
				int t = so;
				so = lo;
				lo = t;
			}
			
			return (s>=so && s<=lo) || (l>=so && l<=lo);
			
		}
		
//		public double getSeedConservation(){
//			return seedConservation;
//		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((contig == null) ? 0 : contig.hashCode());
			result = prime * result + fivePposition;
			result = prime * result + (isPlusStrand ? 1231 : 1237);
			result = prime * result + threePposition;
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof ScoredPosition))
				return false;
			ScoredPosition other = (ScoredPosition) obj;
			if (contig == null) {
				if (other.contig != null)
					return false;
			} else if (!contig.equals(other.contig))
				return false;
			if (isPlusStrand != other.isPlusStrand)
				return false;
			if (fivePposition != other.fivePposition)
				return false;
			if (threePposition != other.threePposition)
				return false;
			return true;
		}

		public String getContig() {
			return contig;
		}

		public int compareTo(ScoredPosition o) {
			int i = new Integer(this.threePposition).compareTo(o.threePposition);
			if(i == 0) i = new Integer(this.fivePposition).compareTo(o.fivePposition);
			return i;
		}
		
		
		
	}

	private HashMap<String, ArrayList<ScoredPosition>> positionMap;
	private String header;
	
	public ScoredPositionOutputParser(String inFile){
		read(inFile);
	}
	
	public void setLongestGeneIndex(AnnotationFileParser aparser){
		for(ScoredPosition sp : getPositions()){
			sp.setLongestGeneIndex(aparser);
		}		
	}

	private void read(String outFile){
		positionMap = new HashMap<String, ArrayList<ScoredPosition>>();
		try {
			BufferedLineReader in  = new BufferedLineReader((outFile));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("Contig")){
					header = s;
					String[] token = header.split("\t");
					for(int i=0;i<token.length;i++){
						if(token[i].equals("5p rmsk")) ScoredPosition.rmsk3pHeaderColumnNumber = i;
						if(token[i].equals("3p rmsk")) ScoredPosition.rmsk5pHeaderColumnNumber = i;
					}					
					//rmsk3pHeaderColumnNumber
					
					
					continue;
				}
				ScoredPosition position = new ScoredPosition(s);				
				if(!positionMap.containsKey(position.getContig()))
					positionMap.put(position.getContig(), new ArrayList<ScoredPosition>());
				positionMap.get(position.getContig()).add(position);
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public String getHeader() { return header; }
	
	public Set<String> getContigs(){
		return positionMap.keySet();
	}
	
	public ArrayList<ScoredPosition> getPositions(String contig){
		if(positionMap.containsKey(contig)) return positionMap.get(contig);
		return new ArrayList<ScoredPosition>();
	}
	
	public ArrayList<ScoredPosition> getPositions(){
		ArrayList<ScoredPosition> positions = new ArrayList<ScoredPosition>();
		for(String key : positionMap.keySet())
			positions.addAll(positionMap.get(key));
		return positions;
	}
	
	static public double getComplexity(String seq){
		try {
		     // Encode a String into bytes
		     byte[] input = seq.getBytes("UTF-8");

		     // Compress the bytes
		     byte[] output = new byte[200];
		     Deflater compresser = new Deflater();
		     compresser.setInput(input);
		     compresser.finish();
		     int compressedDataLength = compresser.deflate(output);
		     compresser.end();
		    // System.out.println(compressedDataLength + " " + seq.length());
		     return ((double)compressedDataLength)/seq.length();
		 } catch(java.io.UnsupportedEncodingException ex) {
		     // handle
		 }
		return 0;
		 
	}
	
	static void generateArffFromCsv(String csv, String arff){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		PrintStream out;
		try {
			out = new PrintStream(arff);
			out.println(Classifier.getArffHeader());
			for(ScoredPosition sp : parser.getPositions()){
				if(sp.isPaired()) out.println(Classifier.toArffString(sp) + (sp.miRNAs != null && !sp.miRNAs.isEmpty() ? "M" : ""));
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	static public void generateBedFromCsv(String csv, String bed5p, String bed3p, boolean invertStrand){
		//if(new File(bed5p).exists() && new File(bed3p).exists()) return;
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		PrintStream out5p, out3p;
		try {
			out5p = new PrintStream(bed5p+".tmp");
			out3p = new PrintStream(bed3p+".tmp");
			
			HashSet<String> fstrs = new HashSet<String>();
			HashSet<String> tstrs = new HashSet<String>();
			
			for(ScoredPosition sp : parser.getPositions()){
				int off = sp.isPlusStrand? 0 : 1;
				boolean strand = invertStrand? !sp.isPlusStrand : sp.isPlusStrand;
				char strandChar = strand? '+' : '-';
				
				String region3p = sp.getGenomicRegions3p() == null? null : sp.getGenomicRegions3p().get(0);
				String region5p = sp.getGenomicRegions5p() == null? null : sp.getGenomicRegions5p().get(0);
				int r3 = region3p == null? 0 : (region3p.endsWith("Intron")? 1 : 2);
				int r5 = region5p == null? 0 : (region5p.endsWith("Intron")? 1 : 2);
				
				
				String fk =  sp.getContig() + sp.getThreePPosition() + strandChar;
				if(!fstrs.contains(fk)){
					String fstring = sp.getContig() + "\t" + Math.max(0, (sp.getThreePPosition()-1+off)) + "\t" + 
							Math.max(2,(sp.getThreePPosition()+1+off)) + "\t" + sp.getThreePPosition() +"\t1\t" + strandChar + "\t" + r3;
					out5p.println(fstring);
					fstrs.add(fk);
				}
				off = sp.isPlusStrand? 1 : 0;
				
				String tk =  sp.getContig() + sp.getFivePPosition() + strandChar;
				if(!tstrs.contains(tk)){
					String tstring = sp.getContig() + "\t" + Math.max(0, (sp.getFivePPosition()-1+off)) + "\t" +
							Math.max(2, (sp.getFivePPosition()+1+off)) + "\t" + sp.getFivePPosition() + "\t1\t" + strandChar + "\t" + r5;
					out3p.println(tstring);
					tstrs.add(tk);
				}
				//if(sp.isPaired()) out5p.println(sp.toArffString() + (sp.miRNAs != null && !sp.miRNAs.isEmpty() ? "M" : ""));
			}
			
			out5p.close();
			out3p.close();
			
			ShellLauncher.run("sort -k1,1 -k2,2n " + bed5p+".tmp >" + bed5p);
			ShellLauncher.run("sort -k1,1 -k2,2n " + bed3p+".tmp >" + bed3p);
			new File(bed5p+".tmp").delete();
			new File(bed3p+".tmp").delete();
			
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	static public void generateRandomizedBedFromCsv(String csv, String bed, AnnotationFileParser annotationParser, ZeroBasedFastaParser fastaParser){
		//if(new File(bed5p).exists() && new File(bed3p).exists()) return;
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		PrintStream out;
		try {
			out = new PrintStream(bed+".tmp");
			HashSet<String> accessions = new HashSet<String>();
			for(ScoredPosition sp : parser.getPositions()){
				if(sp.getContainingGeneAccessions() != null)
					accessions.addAll(sp.getContainingGeneAccessions());
			}			
			
			
			for(ScoredPosition sp : parser.getPositions()){
			//	if(sp.getGenomicRegions3p() == null || !sp.getGenomicRegions3p().get(0).endsWith("Intron")) continue;
			//	if(sp.getGenomicRegions5p() == null || !sp.getGenomicRegions5p().get(0).endsWith("Intron")) continue;
			
				boolean isInterGenic = sp.getContainingGeneAccessions() == null || sp.getContainingGeneAccessions().isEmpty(); 
			//	if(isInterGenic) continue;
				
				String region3p = sp.getGenomicRegions3p() == null? null : sp.getGenomicRegions3p().get(0);
				String region5p = sp.getGenomicRegions5p() == null? null : sp.getGenomicRegions5p().get(0);
				int r3 = region3p == null? 0 : (region3p.endsWith("Intron")? 1 : 2);
				int r5 = region5p == null? 0 : (region5p.endsWith("Intron")? 1 : 2);
				char strandChar = sp.isPlusStrand()? '+' : '-';
				
				for(int n=0;n<2;n++){
					int rd = 0;
					int regionID = 0;
					while(true){
						if(isInterGenic){
							rd = new Random().nextInt(fastaParser.getLength(sp.getContig()));
							ArrayList<AnnotatedGene> genes = annotationParser.getContainingGenes(sp.getContig(), rd);
							if(genes == null){
								regionID = 0;
								break;
							}
						}else{
							ArrayList<AnnotatedGene> genes = annotationParser.getAnnotatedGenes(sp.getContig());
							int rdi = new Random().nextInt(genes.size()); // select random gene
							AnnotatedGene gene = genes.get(rdi);
							if(accessions.contains(gene.getAccession())) continue;
							
							rd = gene.getTxStart() + new Random().nextInt(gene.getTxEnd() - gene.getTxStart());							
							String region = annotationParser.getGenomicRegionNameAndFrameShift(sp.getContig(), sp.isPlusStrand(), rd, gene, true).get(0);
							int r = region == null? 0 : (region.endsWith("Intron")? 1 : 2);
							if(r == r3 || r == r5) break;
							regionID = r;
						}
					}
					
					int off = new Random().nextBoolean()? 0 : 1;
					String fstring = sp.getContig() + "\t" + Math.max(0, (-1+off+rd)) + "\t" + 
							Math.max(2,(1+off+rd)) + "\t" + (rd) +"\t1\t" + strandChar + "\t" + regionID;
					out.println(fstring);				
				}
			}
			
			out.close();
			
			ShellLauncher.run("sort -V " + bed+".tmp >" + bed);
			new File(bed+".tmp").delete();
			
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	static public void generateFastaForMotif(String csv, ZeroBasedFastaParser fastaParser, int flankingNTNumber, String fasta5p, String fasta3p, 
			String mfile, String classification, int option, int aluOption){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		PrintStream out5p, out3p;
		PrintStream mout;
	//	double threshold = 0;// .15;
		try {
			out5p = new PrintStream(fasta5p);
			out3p = new PrintStream(fasta3p);
			mout = new PrintStream(mfile);
			int i = 0;
			int k = 0;
			int motiflenth  = 2;
			HashSet<String> fseqs = new HashSet<String>();
			HashSet<String> tseqs = new HashSet<String>();
			String[][] motifset = {{"A", "C", "G", "T", },
				{"AA", "AC", "AG", "AT","CA", "CC", "CG", "CT","GA", "GC", "GG", "GT","TA", "TC", "TG", "TT",},
				{"A..A", "A..C", "A..G", "A..T","C..A", "C..C", "C..G", "C..T","G..A", "G..C", "G..G", "G..T","T..A", "T..C", "T..G", "T..T",},
				{"AAA", "AAC", "AAG", "AAT","ACA", "ACC", "ACG", "ACT","AGA", "AGC", "AGG", "AGT","ATA", "ATC", "ATG", "ATT","CAA", "CAC", "CAG", "CAT","CCA", "CCC", "CCG", "CCT","CGA", "CGC", "CGG", "CGT","CTA", "CTC", "CTG", "CTT","GAA", "GAC", "GAG", "GAT","GCA", "GCC", "GCG", "GCT","GGA", "GGC", "GGG", "GGT","GTA", "GTC", "GTG", "GTT","TAA", "TAC", "TAG", "TAT","TCA", "TCC", "TCG", "TCT","TGA", "TGC", "TGG", "TGT","TTA", "TTC", "TTG", "TTT",},
				{"A.A", "A.C", "A.G", "A.T","C.A", "C.C", "C.G", "C.T","G.A", "G.C", "G.G", "G.T","T.A", "T.C", "T.G", "T.T",},
				
			};
			int soff = 0;
			int eoff = 30; // inclusive
			
			String[] motifs = null;
		//	soff = -15;
		//	eoff = 0;
			boolean is5pCleavage = true;
			String titleprefix = classification.equals("miRNA")? "miRNA " : "MT ";
			String titlesuffix = null;// " loop region ALU removed";
			//titlesuffix = " stem region only ALU";
			
			if(option == 0){ // 5p stem
				motifs = motifset[1];
				soff = -15;
				eoff = 0;
				motiflenth = 2;
				titlesuffix = " 5'' stem region";
			}else if(option == 1){ // 5p loop
				motifs = motifset[3];
				soff = 18;
				eoff = 25;
				motiflenth = 3;
				titlesuffix = " loop region";
			}else if(option == 2){// 3p stem
				is5pCleavage = false;
				soff = 0;
				eoff = flankingNTNumber - 5;
				motifs = motifset[2];
				motiflenth = 2;
				titlesuffix = " 3'' stem region";
			}else if(option == 3){
				is5pCleavage = false;
				soff = -10;
				eoff = 0;//flankingNTNumber - 5;
				motifs = motifset[4];
				motiflenth = 2;
				titlesuffix = " 3'' stem region";
			}
			if(aluOption == 0);// titlesuffix += 
			else if(aluOption == 1){// remove ALU
				titlesuffix += " ALU removed";
			}else{ // only ALU
				titlesuffix += " only ALU retained";
			}
			
			double[][] hm = new double[motifs.length][eoff - soff + 1];			
			
			int sum = 0;
						
			for(ScoredPosition sp : parser.getPositions()){
				boolean ismiRNA = false;
				if(sp.miRNAs != null && !sp.miRNAs.isEmpty()){
					ismiRNA = true;
				}
				if(ismiRNA){
					if(!classification.equals("miRNA")) continue;
					
				}else if(!sp.isPaired() || !sp.getClassification().equals(classification)) continue;
				if(aluOption==1 && (sp.isALU3p() && sp.isALU5p())) continue; // no alu
				if(aluOption==2 && !(sp.isALU3p() && sp.isALU5p())) continue; // only alu
					
				
				
				//else
				//System.out.println(sp);
				//	continue;
				if(sp.getOverHang() != 2) continue;
				String seq = sp.isPlusStrand()? fastaParser.getSequence(sp.getContig(), sp.getThreePPosition() - flankingNTNumber, sp.getFivePPosition() + flankingNTNumber + 1)
					    : ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(sp.getContig(), sp.getFivePPosition()-flankingNTNumber, sp.getThreePPosition() + flankingNTNumber + 1), true);
				
				String cseq = ScoredPosition.getCleavedSequence(seq, flankingNTNumber);
				//if(sp.getClassification().equals("M") && sp.getFivePScore() > threshold && sp.getThreePScore() > threshold)//  
				{
					//String seq = ScoredPosition.getCleavedSequence(sp.getSeq(), flankingNTNumber);
					String fseq = cseq.substring(0, cseq.indexOf(' '));
					//System.out.println(fseq);
					if(!fseqs.contains(fseq)){
						fseqs.add(fseq);
						out5p.println(">"+ i);
						out5p.println(fseq);
						i++;
					}
				}
				
				//if(sp.getClassification().equals("M") && sp.getFivePScore() > threshold && sp.getThreePScore() > threshold)//  && sp.getFivePScore() > 0.45
				{
					//String seq = ScoredPosition.getCleavedSequence(sp.getSeq(), flankingNTNumber);
					String tseq = cseq.substring(cseq.indexOf(' ') + 1, Math.min(cseq.length(), cseq.indexOf(' ') + 28));//seq.substring(seq.lastIndexOf(' ') + 1);
					//System.out.println(k+ " " + seq);
					//System.out.println("                          " + tseq);
					
					if(!tseqs.contains(tseq)){
						tseqs.add(tseq);
						out3p.println(">"+ k);
						out3p.println(tseq);
						k++;
					}
				}
				
				{
					/*int soff = -15;
			int eoff = 0; // inclusive
			boolean is5pCleavage = true;
			*/
					//String seq = sp.isPlusStrand()? fastaParser.getSequence(sp.getContig(), sp.getThreePPosition() - flankingNTNumber, sp.getFivePPosition() + flankingNTNumber + 1)
					//		    : ZeroBasedFastaParser.getComplementarySequence(fastaParser.getSequence(sp.getContig(), sp.getFivePPosition()-flankingNTNumber, sp.getThreePPosition() + flankingNTNumber + 1), true);
						
					String tseq = is5pCleavage? seq.toUpperCase().substring(flankingNTNumber+soff, flankingNTNumber + eoff+1+5) :
						seq.toUpperCase().substring(seq.length() - flankingNTNumber + soff , seq.length() - flankingNTNumber + eoff + 1+3);
					System.out.println(k+ " " + seq.toUpperCase());
					System.out.println("                          " + tseq);
					for(int j=0; j<motifs.length;j++){
						String motif = motifs[j];
						Pattern p = Pattern.compile(motif);
						Matcher m = p.matcher(tseq);
						while(m.find()){
							int ms = m.start();
							if(ms < hm[j].length)
								hm[j][m.start()]++;
						}						
					}					
				}
				sum++;
			}
			
			for(int j=0;j<motifs.length;j++){
				mout.println(motifs[j].replace('.', 'N') + " = [");
				for(double l : hm[j]) mout.print(l + " ");
				mout.println("];");
			}
			
			mout.print("t=[");
			for(String motif : motifs) mout.print(motif.replace('.', 'N') + "; ");
			mout.println("];");
			mout.print("motif={");
			for(String motif : motifs) mout.print("'" + motif.replace('T', 'U').replace('.', 'N') + "', ");
			mout.println("};");
			mout.print("xtick=[" );
			mout.println(is5pCleavage ? soff + ":" +eoff+"];" : "-" +soff + ":-1:-" +eoff+"];");
			mout.print("titlestr = '");
			mout.print(titleprefix);
			mout.print("motif frequency : " + titlesuffix + " (n = " + sum + ")';");
			mout.print("totaln="+sum+";");
			mout.print("motiflenth = " + motiflenth);
			mout.println(";xstr='Position from " + (is5pCleavage? "5''" : "3''") + " cleavage';");
				
			mout.close();
			
			out5p.close();
			out3p.close();
			mout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	static public void reClassify(String csv, String outFile, String arff){
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		Classifier classifier = new Classifier(arff);
		try {
			PrintStream out = new PrintStream(outFile);
			out.println(ScoredPosition.getHeader());
			for(ScoredPosition sp : parser.getPositions()){
				RNAfoldLauncher fold = new RNAfoldLauncher(sp.seq, FCLIP_Scorer.getFlankingNTNumber());
				if(!fold.isFeasibleFold()) continue;
				classifier.setClassification(sp);
				
				out.println(sp);
			}
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	static public void intersectWithBam(String csvOut, String csv, String bamfile, boolean contained){
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamfile));
				
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);		
		
		try {	
			PrintStream out = new PrintStream(csvOut);
			out.println(parser.getHeader()+"\t" + new File(bamfile).getName() + (contained? "_containedReads":"_overlappingReads"));
			for(ScoredPosition p :parser.getPositions()){
				int p1 = p.getFivePPosition();
				int p2 = p.getThreePPosition();
				if(p1 > p2){
					int t = p1;
					p1 = p2;
					p2 = t;
				}
				if(p1<0){
					p1 = p2 - FCLIP_Scorer.getMinReadDiff();
				}
				
				SAMRecordIterator iterator = reader.query(p.getContig(), p1+1, p2+1, contained);
				int n = 0;
				try{
				while(iterator.hasNext()){
					SAMRecord r = iterator.next();
					n++;
				}
				}finally{
					iterator.close();
				}
				out.println(p+"\t"+n);
				
			}
			out.close();	
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static public void main(String[] args) throws IOException{
		ScoredPositionOutputParser parser  = new ScoredPositionOutputParser("/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.mt.csv");
		int counter = 0;
		for(ScoredPosition sp : parser.getPositions()){
			boolean count = false;
			if(sp.getGenomicRegions3p() != null){
				for(String rg : sp.getGenomicRegions3p()){
					if(rg.endsWith("ORF") || rg.endsWith("UTR")){
						count = true;
						break;
					}
				}
			}
			if(sp.getGenomicRegions5p() != null){				
				for(String rg : sp.getGenomicRegions5p()){
					if(rg.endsWith("ORF") || rg.endsWith("UTR")){
						count = true;
						break;
					}
				}
			}
			if(count)counter++;
		}
		System.out.println(counter);
	}
	/*
	static public void main(String[] args) throws IOException{
		int maxOffset = 11;
		String csv = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.csv";
		String outCsv = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.regionWithOffset"+maxOffset+".txt";
		//ZeroBasedFastaParser fparser = new ZeroBasedFastaParser("/media/kyowon/Data1/fCLIP/genomes/hg38.fa");
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg38.refFlat.txt");
		PrintStream out = new PrintStream(outCsv);
		HashSet<ScoredPosition> sps = new HashSet<ScoredPositionOutputParser.ScoredPosition>();
		out.println(parser.getHeader());
		for(ScoredPosition sp : parser.getPositions()){
			if(!sp.getClassification().toUpperCase().equals("M")) continue;
			for(int o=0;o<=maxOffset;o++){
				sp.setGeneInfo(annotationParser, o);
				if(sp.getGenomicRegions3p() == null || sp.getGenomicRegions5p() == null) continue;
				boolean toPut = false;
				for(int i=0;i<sp.getGenomicRegions3p().size();i++){
					if(sp.getGenomicRegions3p().get(i).endsWith("Intron") && !sp.getGenomicRegions5p().get(i).endsWith("Intron")) toPut = true;
					if(sp.getGenomicRegions5p().get(i).endsWith("Intron") && !sp.getGenomicRegions3p().get(i).endsWith("Intron")) toPut = true;
					if(toPut) break;					
				}
				if(toPut) {
					System.out.println(sp.getGenomicRegions3p() + " " + sp.getGenomicRegions5p());
					if(!sps.contains(sp)) sps.add(sp);
					break;
				}
			}
		}
		for(ScoredPosition sp : sps){
			out.println(sp);
			
		}
		out.close();		
		*/
		/*String classification = "M";
	
		
		
		String fasta5p = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.motif5p." + classification + ".txt";
		String fasta3p = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.completeNoverlapped.motifloop." + classification + ".txt";
		int region = 0;
		ZeroBasedFastaParser fparser = new ZeroBasedFastaParser("/media/kyowon/Data1/fCLIP/genomes/hg38.fa");
		int flankingNT = 25;
		int aluOption = 2;
		for(;region<3;region++){
			String mout = "/media/kyowon/Data1/fCLIP/motif_" + classification + region + (aluOption>0? "_" + aluOption : "") +".m";
			ScoredPositionOutputParser.generateFastaForMotif(csv, fparser, flankingNT, fasta5p, fasta3p, mout, classification, region, aluOption);
		}*/


}
