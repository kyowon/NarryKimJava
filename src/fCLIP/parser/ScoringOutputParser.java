package fCLIP.parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.Deflater;

import launcher.PhyloPLauncher;
import launcher.RNAfoldLauncher;
import launcher.RNAzLauncher;
import net.sf.samtools.util.BufferedLineReader;
import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.Bed12Parser;
import parser.MafParser;
import parser.ZeroBasedFastaParser;
import util.DnDsCalculator;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import fCLIP.Classifier;
import fCLIP.MirGff3FileParser.MiRNA;
import fCLIP.Scorer;

public class ScoringOutputParser {
	public static class ScoredPosition implements Comparable<ScoredPosition>{
		private String contig;		
		private boolean isPlusStrand;
		private int fivePposition = -1;
		private int threePposition = -1;
		//private ArrayList<Integer> fivePcoordinate;
		//private ArrayList<Integer> threePcoordinate;
		
		private double sci = 0;
		private double zscore = 1;
		private double seedConservation = 1;
		private double dnds = -1;
		
		private double fivePScore = -1;
		private double threePScore = -1;
		private int depth;
		private double energy;
		private int hairpinNum;
		private int leftPaired;
		private int rightPaired;
		private int overHang;
		private String seq;
		private ArrayList<String> miRNAs = null;
		private ArrayList<String> containingGeneAccessions;
		private ArrayList<String> containingGeneNames;
		private ArrayList<String> genomicRegions5p;
		private ArrayList<String> genomicRegions3p;
		//private double complexity = 0;
		private String classification;
		private double predictionScore;
		private int[] threePreads;
		private int[] fivePreads;
		private int blatHits = 0;
		private double seqEntropy = 0;
		private double structureEntropy = 0;
		private double hetero = 0;
		
		//private int readCount = 0;
		//private double readLengthVariance = -1;
		
		public static String getHeader(){
			return "Contig\tStrand\t5p\t3p\tPaired\tClass\tPredictionScore\t5p5preads\t5p3preads\t3p5preads\t3p3preads\tMatching miRNAs\tBlatHits\tSCI\tz-score\tDnDs\tSeedConservation\tAccession\tGeneName\tGenomicRegion5p\tGenomicRegion3p\t5pSocre\t3pScore\tDepth\tEnergy\tHairpinNum\tLeftPaired\tRightPaired\tOverhang\tSeqEntropy\tStructureEntropy\tHetero\tSeq";
		}
		
		public ArrayList<String> getContainingGeneAccessions() {
			return containingGeneAccessions;
		}

		public ArrayList<String> getContainingGeneNames() {
			return containingGeneNames;
		}
		
		public double getSeqEntropy(){
			return seqEntropy;
		}
		
		public double getStructureEntropy(){
			return structureEntropy;
		}
		
		public double getZScore(){
			return zscore;
		}
		
		public int getFivePposition() {
			return fivePposition;
		}

		public int getThreePposition() {
			return threePposition;
		}

		public double getHetero(){
			return hetero;
		}
				
		public double getSCI(){
			return sci;
		}
		
		public int getOverHang(){
			return overHang;
		}
		
		public ArrayList<String> getGenomicRegions5p(){
			return genomicRegions5p;
		}
		
		public ArrayList<String> getGenomicRegions3p(){
			return genomicRegions3p;
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
		
		public int getHairpinNumber(){
			return hairpinNum;
		}
		
		public ArrayList<String> getMiRNAs(){
			return miRNAs;
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
		//public int getReadCount(){
		//	return readCount;
		//}
		
		public ScoredPosition(String s){
			String[] token = s.split("\t");
			int i = 0;
			this.contig = token[i++];
			this.isPlusStrand = token[i++].equals("+");
			this.threePposition = Integer.parseInt(token[i++]);
			this.fivePposition = Integer.parseInt(token[i++]);
		
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
				String[] miRNAsString = token[i].split(",");
				this.miRNAs = new ArrayList<String>();
				for(String st : miRNAsString){
					miRNAs.add(st);
				}
			}
			i++;
			this.blatHits = Integer.parseInt(token[i++]);
			this.sci = Double.parseDouble(token[i++]);
			this.zscore = Double.parseDouble(token[i++]);
			this.dnds = Double.parseDouble(token[i++]);
			this.seedConservation  = Double.parseDouble(token[i++]);
			
			if(!token[i].startsWith("_")){
				String[] accessions = token[i].split(",");
				this.containingGeneAccessions = new ArrayList<String>();
				for(String st : accessions){
					containingGeneAccessions.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] names = token[i].split(",");
				this.containingGeneNames = new ArrayList<String>();
				for(String st : names){
					containingGeneNames.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] regions = token[i].split(",");
				this.genomicRegions3p = new ArrayList<String>();
				for(String st : regions){
					genomicRegions3p.add(st);
				}
			}
			i++;
			if(!token[i].startsWith("_")){
				String[] regions = token[i].split(",");
				this.genomicRegions5p = new ArrayList<String>();
				for(String st : regions){
					genomicRegions5p.add(st);
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
			this.hairpinNum = Integer.parseInt(token[i++]);
			this.leftPaired = Integer.parseInt(token[i++]);
			this.rightPaired = Integer.parseInt(token[i++]);
			this.overHang = Integer.parseInt(token[i++]);
			this.seqEntropy = Double.parseDouble(token[i++]);
			this.structureEntropy = Double.parseDouble(token[i++]);
			this.hetero = Double.parseDouble(token[i++]);
			this.seq = token[i].replaceAll(" ", "");
			
		}
		
		public static void writeStringArray(StringBuilder sb, ArrayList<String> array){
			if(array != null && !array.isEmpty()){
				for(int i=0; i<array.size();i++){
					sb.append(array.get(i));sb.append(i == array.size()-1 ? '\t' : ',');
				}			
			}else sb.append("_\t");
		}
		
		public RNAzLauncher getRNAzLauncher(MafParser mafParser){
			String mafSeq = getMafString(mafParser);
			return new RNAzLauncher(mafSeq);
		}
		
		public void setRNAzScoresNSeedConservation(MafParser mafParser){
		//	String mafSeq = getMafString(mafParser);
			RNAzLauncher rna = getRNAzLauncher(mafParser);
			this.sci = rna.getSCI();
			this.zscore = rna.getZScore();
			this.seedConservation = new PhyloPLauncher(getSeedMafString(mafParser)).getPvalConservation();
		}
		
	//	public void setSeedConvervation(MafParser mafParser){
	//		String mafString = mafParser.getSeqsInMafFormat(contig, position.getCoordinate(), position.getPosition(), isPlusStrand, positionQuantityChangeLength);
	//		this.phyloP = new PhyloPLauncher(mafString).getPvalConservation();
	//	}
		
		
		public void setDnDs(MafParser mafParser){
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
		
		private String getMafString(MafParser mafParser){
			int length = seq.length() - 2 * Scorer.flankingNTNumber;
			int position = this.threePposition;// > 0 ? this.threePposition 
			if(position < 0 || this.fivePposition < 0) return null;
			//		: this.fivePposition;// + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
			return mafParser.getSeqsInMafFormat(contig, position, isPlusStrand, length);
		} 
		
		private String getSeedMafString(MafParser mafParser){
			int length = Math.min(7, seq.length() - 2 * Scorer.flankingNTNumber);// ;
			int position = this.threePposition;
			if(position < 0 || this.fivePposition < 0) return null;
					//: this.fivePposition;// + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
			return mafParser.getSeqsInMafFormat(contig, position + (isPlusStrand? 1 : -1), isPlusStrand, length);
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(contig); sb.append('\t');
			sb.append(isPlusStrand? '+': '-'); sb.append('\t');
			sb.append(threePposition); sb.append('\t');
			sb.append(fivePposition); sb.append('\t');
			sb.append(isPaired() ? "T" : "F" ); sb.append('\t');
			sb.append(classification); sb.append('\t');
			sb.append(predictionScore); sb.append('\t');
			sb.append(threePreads == null? '_' : threePreads[0]); sb.append('\t');
			sb.append(threePreads == null? '_' : threePreads[1]); sb.append('\t');
			
			sb.append(fivePreads == null? '_' : fivePreads[0]); sb.append('\t');
			sb.append(fivePreads == null? '_' : fivePreads[1]); sb.append('\t');
			
		//	sb.append(readCount); sb.append('\t');
	//		sb.append(readLengthVariance); sb.append('\t');
			
			writeStringArray(sb, miRNAs);
			sb.append(blatHits); sb.append('\t');
			sb.append(sci); sb.append('\t');
			sb.append(zscore); sb.append('\t');
			sb.append(dnds); sb.append('\t');
			sb.append(seedConservation); sb.append('\t');
			writeStringArray(sb, containingGeneAccessions);
			writeStringArray(sb, containingGeneNames);
			writeStringArray(sb, genomicRegions3p);
			writeStringArray(sb, genomicRegions5p);
			
			//toMappedRegionString(sb, threePcoordinate); sb.append('\t');
			//toMappedRegionString(sb, fivePcoordinate); sb.append('\t');
			sb.append(threePScore); sb.append('\t');
			sb.append(fivePScore); sb.append('\t');
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNum); sb.append('\t');
			sb.append(leftPaired); sb.append('\t');
			sb.append(rightPaired); sb.append('\t');
			sb.append(overHang); sb.append('\t');
			sb.append(seqEntropy); sb.append('\t');
			sb.append(structureEntropy); sb.append('\t');
			sb.append(hetero); sb.append('\t');
			sb.append(getCleavedSequence(seq, Scorer.flankingNTNumber));
			return sb.toString();
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
		
		public static String getArffHeader(){
			return "@relation w\n\n@attribute Energy numeric\n@attribute HairpinNumber numeric\n@attribute SCI numeric"
					+ "\n@attribute ZScore numeric"
					+ "\n@attribute structureEntropy numeric"
					+ "\n@attribute Hetero numeric"
					+ "\n@attribute Class {M,U}\n\n@data";
		}
		
		public String toTrainingArffString(){
			StringBuilder sb = new StringBuilder();
			sb.append(energy); sb.append(',');
			sb.append(hairpinNum); sb.append(',');
			sb.append(sci<0? "?" : sci); sb.append(',');
			sb.append(zscore==10000? "?" : zscore); sb.append(',');
		//	sb.append(is3pScored()? threePScore : "?"); sb.append(',');
		//	sb.append(is5pScored()? fivePScore : "?"); sb.append(',');
		//	sb.append(seqEntropy); sb.append(',');
			sb.append(structureEntropy); sb.append(',');
			sb.append(hetero); sb.append(',');
			if(this.miRNAs == null ||  this.miRNAs.isEmpty()){
				if(this.overHang<0 || this.overHang>5) sb.append("U");
				else sb.append("?");
			}else{
				sb.append('M');
			}
			return sb.toString();
		}
		
		public String toArffString(){
			StringBuilder sb = new StringBuilder();
			sb.append(energy); sb.append(',');
			sb.append(hairpinNum); sb.append(',');
			sb.append(sci<0? "?" : sci); sb.append(',');
			sb.append(zscore==10000? "?" : zscore); sb.append(',');
		//	sb.append(is3pScored()? threePScore : "?"); sb.append(',');
		//	sb.append(is5pScored()? fivePScore : "?"); sb.append(',');
		//	sb.append(seqEntropy); sb.append(',');
			sb.append(structureEntropy); sb.append(',');
			sb.append(hetero); sb.append(',');
			String c = (classification == null? "?" : classification.toUpperCase());
			//if(c.equals("M")){
			//	if(this.miRNAs != null && !this.miRNAs.isEmpty()) c = "M"; 
			//	else c = "m";	//}
			//if(this.isPaired()) c +="P";
			sb.append(c);
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
			setGeneInfo(annotationParser);
					
		}
		
		private void setGeneInfo(AnnotationFileParser annotationParser){
			if(annotationParser != null){
				HashMap<AnnotatedGene, Integer> genes5p = new HashMap<AnnotatedGene, Integer>();
				HashMap<AnnotatedGene, Integer> genes3p = new HashMap<AnnotatedGene, Integer>();
				if(fivePposition >= 0){
					ArrayList<AnnotatedGene> ag = annotationParser.getContainingGenes(contig, isPlusStrand, fivePposition);
					if(ag != null){
						for(AnnotatedGene a : ag){
							genes5p.put(a, fivePposition);
						}
					}
				}
				if(threePposition>= 0){
					ArrayList<AnnotatedGene> ag = annotationParser.getContainingGenes(contig, isPlusStrand, threePposition);
					if(ag != null){
						for(AnnotatedGene a : ag){
							genes3p.put(a, threePposition);
						}
					}
				}
				this.containingGeneAccessions = new ArrayList<String>();
				this.containingGeneNames = new ArrayList<String>();
				this.genomicRegions3p = new ArrayList<String>();
				this.genomicRegions5p = new ArrayList<String>();
				HashSet<AnnotatedGene> genes = new HashSet<AnnotationFileParser.AnnotatedGene>();
				genes.addAll(genes5p.keySet());
				genes.addAll(genes3p.keySet());
				for(AnnotatedGene gene : genes){
					this.containingGeneAccessions.add(gene.getAccession());
					this.containingGeneNames.add(gene.getGeneName());
					this.genomicRegions3p.add(genes3p.containsKey(gene) ? annotationParser.getGenomicRegionNameAndFrameShift(contig, isPlusStrand, genes3p.get(gene), gene, true).get(0) : "_");
					this.genomicRegions5p.add(genes5p.containsKey(gene) ? annotationParser.getGenomicRegionNameAndFrameShift(contig, isPlusStrand, genes5p.get(gene), gene, true).get(0) : "_");
					
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
		}
		
		public ScoredPosition setClassification(String e, double predictionScore){
			this.classification = e; 
			this.predictionScore = predictionScore; 
			return this;
		}
		
		public ScoredPosition set(RNAfoldLauncher fold){
			depth = fold.getDepth();
			energy = fold.getEnergyPerNT();
			hairpinNum = fold.getNumberOfHairPins();
			leftPaired = fold.getLeftPaired();
			rightPaired = fold.getRightPaired();
			overHang = fold.getOverHang();
			seqEntropy = fold.getSeqEntropy();
			structureEntropy = fold.getStructureEntropy();
			return this;
		}
		
		
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
			setGeneInfo(annotationParser);
			return this;
		}

		public ScoredPosition setThreePposition(int i, AnnotationFileParser annotationParser) {
			threePposition = i;
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
			
		public double getSeedConservation(){
			return seedConservation;
		}
		
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
	
	public ScoringOutputParser(String inFile){
		read(inFile);
	}
	

	private void read(String outFile){
		positionMap = new HashMap<String, ArrayList<ScoredPosition>>();
		try {
			BufferedLineReader in  = new BufferedLineReader(new FileInputStream(outFile));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("Contig")){
					header = s;
					continue;
				}
				ScoredPosition position = new ScoredPosition(s);				
				if(!positionMap.containsKey(position.getContig()))
					positionMap.put(position.getContig(), new ArrayList<ScoredPosition>());
				positionMap.get(position.getContig()).add(position);
			}
			in.close();
		} catch (FileNotFoundException e) {
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
		ScoringOutputParser parser = new ScoringOutputParser(csv);
		PrintStream out;
		try {
			out = new PrintStream(arff);
			out.println(ScoredPosition.getArffHeader());
			for(ScoredPosition sp : parser.getPositions()){
				if(sp.isPaired()) out.println(sp.toArffString() + (sp.miRNAs != null && !sp.miRNAs.isEmpty() ? "M" : ""));
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	static public void generateBedFromCsv(String csv, String bed5p, String bed3p, boolean invertStrand){
		ScoringOutputParser parser = new ScoringOutputParser(csv);
		PrintStream out5p, out3p;
		try {
			out5p = new PrintStream(bed5p);
			out3p = new PrintStream(bed3p);
			for(ScoredPosition sp : parser.getPositions()){
				int off = sp.isPlusStrand? 0 : 1;
				boolean strand = invertStrand? !sp.isPlusStrand : sp.isPlusStrand;
				char strandChar = strand? '+' : '-';
				 
				out5p.println(sp.getContig() + "\t" + Math.max(0, (sp.getThreePposition()-1+off)) + "\t" + Math.max(2,(sp.getThreePposition()+1+off)) + "\t" + sp.getThreePposition() +"\t1\t" + strandChar);
				off = sp.isPlusStrand? 1 : 0;
				out3p.println(sp.getContig() + "\t" + Math.max(0, (sp.getFivePposition()-1+off)) + "\t" + Math.max(2, (sp.getFivePposition()+1+off)) + "\t" + sp.getFivePposition() + "\t1\t" + strandChar);
				
				//if(sp.isPaired()) out5p.println(sp.toArffString() + (sp.miRNAs != null && !sp.miRNAs.isEmpty() ? "M" : ""));
			}
			
			out5p.close();
			out3p.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	
	static public void main(String[] args){
		ScoringOutputParser.generateBedFromCsv("/media/kyowon/Data1/Dropbox/h19x2.sorted.out.csv", "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.5p.bed",
				"/media/kyowon/Data1/Dropbox/h19x2.sorted.out.3p.bed", false);
	}
}
