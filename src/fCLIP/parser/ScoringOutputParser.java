package fCLIP.parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.Deflater;

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
import fCLIP.MirGff3FileParser.MiRNA;
import fCLIP.RNAfoldLauncher;
import fCLIP.RNAzLauncher;
import fCLIP.Scorer;

public class ScoringOutputParser {
	public static class ScoredPosition implements Comparable<ScoredPosition>{
		private String contig;		
		private boolean isPlusStrand;
		private int fivePposition = -1;
		private int threePposition = -1;
		private ArrayList<Integer> fivePcoordinate;
		private ArrayList<Integer> threePcoordinate;
		
		private double sci = -1;
		private double zscore = 10000;
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
		private ArrayList<String> genomicRegions;
		private double complexity = 0;
		private String classification;
		private double predictionScore;
		
		public static String getHeader(){
			return "Contig\tStrand\t3p\t5p\tPaired\tClass\tPredictionScore\tMatching miRNAs\tSCI\tz-score\tDnDs\tAccession\tGeneName\tGenomicRegion\t3pStarts\t3pEnds\t5pStarts\t5pEnds\t3pSocre\t5pScore\tDepth\tEnergy\tHairpinNum\tLeftPaired\tRightPaired\tOverhang\tComplexity\tSeq";
		}
		
		public ArrayList<String> getContainingGeneAccessions() {
			return containingGeneAccessions;
		}

		public ArrayList<String> getContainingGeneNames() {
			return containingGeneNames;
		}
		
		public double getZScore(){
			return zscore;
		}
		
		public double getSCI(){
			return sci;
		}
		
		public ArrayList<String> getGenomicRegions(){
			return genomicRegions;
		}
		
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

		public double getSequencComplexity(){
			return complexity;
		}
		
		public ArrayList<String> getMiRNAs(){
			return miRNAs;
		}
		
		public boolean isPlusStrand(){ return isPlusStrand;}
		
		public boolean is3pScored(){
			return threePcoordinate != null;
		}
	
		public boolean is5pScored(){
			return fivePcoordinate != null;
		}
		
		public boolean isPaired(){
			return threePcoordinate != null && fivePcoordinate !=null;
		}
		
		public String getClassification(){
			return classification;
		}
		
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
			if(!token[i].startsWith("_")){
				String[] miRNAsString = token[i].split(",");
				this.miRNAs = new ArrayList<String>();
				for(String st : miRNAsString){
					miRNAs.add(st);
				}
			}
			i++;
			
			this.sci = Double.parseDouble(token[i++]);
			this.zscore = Double.parseDouble(token[i++]);
			this.dnds = Double.parseDouble(token[i++]);
			
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
				this.genomicRegions = new ArrayList<String>();
				for(String st : regions){
					genomicRegions.add(st);
				}
			}
			i++;
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
			this.threePScore = Double.parseDouble(token[i++]);
			this.fivePScore = Double.parseDouble(token[i++]);
			this.depth = Integer.parseInt(token[i++]);
			this.energy = Double.parseDouble(token[i++]);
			this.hairpinNum = Integer.parseInt(token[i++]);
			this.leftPaired = Integer.parseInt(token[i++]);
			this.rightPaired = Integer.parseInt(token[i++]);
			this.overHang = Integer.parseInt(token[i++]);
			this.complexity = Double.parseDouble(token[i++]);
			this.seq = token[i];
			
		}
		
		public static void writeStringArray(StringBuilder sb, ArrayList<String> array){
			if(array != null && !array.isEmpty()){
				for(int i=0; i<array.size();i++){
					sb.append(array.get(i));sb.append(i == array.size()-1 ? '\t' : ',');
				}			
			}else sb.append("_\t");
		}
		
		public RNAzLauncher getRNAzLauncher(MafParser mafParser){
			String maf = getMafString(mafParser);
			//System.out.println(maf);
			//maf = "a score=1\ns s 0 1 - 1 cctagtagttgggattaca---ggtgcctgccaccatgcttgaccaatttcttgtatttttaatggcgatggggtttcacc";
			return new RNAzLauncher(maf);
		}
		
		public void setRNAzScores(MafParser mafParser){
			RNAzLauncher rna = getRNAzLauncher(mafParser);
			this.sci = rna.getSCI();
			this.zscore = rna.getZScore();
		}
		
		
		public void setDnDs(MafParser mafParser){
			int length = seq.length();
			int position = this.threePposition > 0 ? this.threePposition + (isPlusStrand? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) : this.fivePposition + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
		
			String[] seqs = mafParser.getSeqs(contig, position, isPlusStrand, length);
			if(seqs.length < 2) return;
			//
			//for(String seq : seqs) System.out.println(seq);
			
			this.dnds = DnDsCalculator.calculate(seqs, false);
			
		}
		
		private String getMafString(MafParser mafParser){
			int length = seq.length();
			int position = this.threePposition > 0 ? this.threePposition + (isPlusStrand? -Scorer.flankingNTNumber + 1 : Scorer.flankingNTNumber - 1) : this.fivePposition + (isPlusStrand? Scorer.flankingNTNumber - length : length - Scorer.flankingNTNumber);			
		//	System.out.println(contig + " " + isPlusStrand + " " + position);
			return mafParser.getSeqsInMafFormat(contig, position, isPlusStrand, length);
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
			writeStringArray(sb, miRNAs);
			sb.append(sci); sb.append('\t');
			sb.append(zscore); sb.append('\t');
			sb.append(dnds); sb.append('\t');
			writeStringArray(sb, containingGeneAccessions);
			writeStringArray(sb, containingGeneNames);
			writeStringArray(sb, genomicRegions);
			
			toMappedRegionString(sb, threePcoordinate); sb.append('\t');
			toMappedRegionString(sb, fivePcoordinate); sb.append('\t');
			sb.append(threePScore); sb.append('\t');
			sb.append(fivePScore); sb.append('\t');
			sb.append(depth); sb.append('\t');
			sb.append(energy); sb.append('\t');
			sb.append(hairpinNum); sb.append('\t');
			sb.append(leftPaired); sb.append('\t');
			sb.append(rightPaired); sb.append('\t');
			sb.append(overHang); sb.append('\t');
			sb.append(complexity); sb.append('\t');
			sb.append(seq);
			
			return sb.toString();
		}
		
		public static String getArffHeader(){
			return "@relation w\n\n@attribute Energy numeric\n@attribute HairpinNumber numeric\n@attribute SCI numeric"
					+ "\n@attribute ZScore numeric\n@attribute Overhang"
					+ " numeric\n@attribute Class {M,U}\n\n@data";
		}
		
		public String toArffString(){
			StringBuilder sb = new StringBuilder();
			//sb.append(isPaired()? "1" : "0");sb.append(',');
		//	sb.append(depth);sb.append(',');
			sb.append(energy); sb.append(',');
			sb.append(hairpinNum); sb.append(',');
			sb.append(sci<0? "?" : sci); sb.append(',');
			sb.append(zscore==10000? "?" : zscore); sb.append(',');
			//	sb.append(leftPaired); sb.append(',');
		//	sb.append(rightPaired); sb.append(',');
			sb.append(overHang); sb.append(',');
			sb.append(this.miRNAs == null || this.miRNAs.isEmpty()? 'U' : 'M');
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
		
		public ScoredPosition(String contig, boolean isPlusStrand, int fivePposition, ArrayList<Integer> fivePcoordinate, 
				int threePposition, ArrayList<Integer> threePcoordinate, AnnotationFileParser annotationParser){
			this.contig = contig;
			this.isPlusStrand = isPlusStrand;
			this.fivePposition = fivePposition;
			this.threePposition = threePposition;
			this.fivePcoordinate = fivePcoordinate;
			this.threePcoordinate = threePcoordinate;
			if(annotationParser != null){
			HashSet<AnnotatedGene> genes = new HashSet<AnnotatedGene>();
				if(fivePcoordinate != null){
					ArrayList<AnnotatedGene> ag = annotationParser.getContainingGenes(contig, isPlusStrand, fivePposition);
					if(ag != null) genes.addAll(ag);
				}
				if(threePcoordinate != null){
					ArrayList<AnnotatedGene> ag = annotationParser.getContainingGenes(contig, isPlusStrand, threePposition);
					if(ag != null) genes.addAll(ag);
				}
				this.containingGeneAccessions = new ArrayList<String>();
				this.containingGeneNames = new ArrayList<String>();
				this.genomicRegions = new ArrayList<String>();
				for(AnnotatedGene gene : genes){
					this.containingGeneAccessions.add(gene.getAccession());
					this.containingGeneNames.add(gene.getGeneName());
					this.genomicRegions.add(annotationParser.getGenomicRegionNameAndFrameShift(contig, isPlusStrand, fivePcoordinate != null?fivePposition : threePposition , gene, true).get(0));
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
		
		public ScoredPosition setFivePScore(double s){ this.fivePScore = s; return this;}
		public ScoredPosition setThreePScore(double s){ this.threePScore = s; return this;}
		public ScoredPosition setClassification(String e, double predictionScore){this.classification = e; this.predictionScore = predictionScore; return this;}
		public ScoredPosition set(RNAfoldLauncher fold){
			depth = fold.getDepth();
			energy = fold.getEnergy();
			hairpinNum = fold.getNumberOfHairPins();
			leftPaired = fold.getLeftPaired();
			rightPaired = fold.getRightPaired();
			overHang = fold.getOverHang();
			return this;
		}
		
		
		public String getSeq() {
			return seq;
		}

		public void setSeq(String seq) {
			this.complexity = getComplexity(seq);
			this.seq = seq;
		}

		
		public void setMiRNAs(ArrayList<MiRNA> miRNAs) {
			this.miRNAs = new ArrayList<String>();
			for(MiRNA miRNA : miRNAs){
				this.miRNAs.add(miRNA.getName());
			}
		}

		public int getFivePposition() {
			return fivePposition;
		}

		public int getThreePposition() {
			return threePposition;
		}

		public ScoredPosition setFivePposition(int i) {
			fivePposition = i;
			return this;
		}

		public ScoredPosition setThreePposition(int i) {
			threePposition = i;
			return this;
		}
		
		
		public ArrayList<Integer> getFivePcoordinate() {
			return fivePcoordinate;
		}

		public ArrayList<Integer> getThreePcoordinate() {
			return threePcoordinate;
		}

		public double getFivePScore() {
			return fivePScore;
		}

		public double getThreePScore() {
			return threePScore;
		}

		public int getLeftPaired() {
			return leftPaired;
		}
		
		public int getRightPaired(){
			return rightPaired;
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
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((contig == null) ? 0 : contig.hashCode());
			result = prime
					* result
					+ ((fivePcoordinate == null) ? 0 : fivePcoordinate
							.hashCode());
			result = prime * result + fivePposition;
			result = prime * result + (isPlusStrand ? 1231 : 1237);
			result = prime
					* result
					+ ((threePcoordinate == null) ? 0 : threePcoordinate
							.hashCode());
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
			if (fivePcoordinate == null) {
				if (other.fivePcoordinate != null)
					return false;
			} else if (!fivePcoordinate.equals(other.fivePcoordinate))
				return false;
			if (threePcoordinate == null) {
				if (other.threePcoordinate != null)
					return false;
			} else if (!threePcoordinate.equals(other.threePcoordinate))
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
	
	public ScoringOutputParser(String inFile){
		read(inFile);
	}
	

	private void read(String outFile){
		positionMap = new HashMap<String, ArrayList<ScoredPosition>>();
		try {
			BufferedLineReader in  = new BufferedLineReader(new FileInputStream(outFile));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("Contig")) continue;
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
	
	static public void main(String[] args){
		MafParser mafParser = new MafParser("/media/kyowon/Data1/RPF_Project/genomes/hg19/maf/");
		mafParser.readIndexFile();
		//GGTCTCAATGCCTTCATAGCCCTGTACAATGCTGCTTGATCCATATGCAACAAGGCAGCACTGTAAAGAAGCCGAGGGCAGTAAGA
		//chr5	-	852928	852863	F	M	0.8518518519	_	_	_	_	852908	852968			0.111582863	-1	27	-0.4164835165	1	8	8	0	0.5494505495	CCTCCGTGGGGTCTGTCCCCGCTGCTTGCTGTCCTCCGTGGGGTCTCTTCCCCCGGCTCGCTGTCCTCCATGGGGTCTGTCCCCCGGTTTG

		
		//TCTTACTGCCCTCGGCTTCTTTACAGTGCTGCCTTGTTGCATATGGATCAAGCAGCATTGTACAGGGCTATGAAGGCATTGAGACC
		ScoredPosition test = new ScoredPosition("chr5	-	64874930	64874877	F	U	0.5722689076	_	0	-0.96	1.6340435743	NM_015342,NM_001278927,NM_001278926,NM_001278929	PPWD1,PPWD1,PPWD1,PPWD1	NM_ORF,NM_ORF,NM_ORF,NM_ORF	64874590	64874690			0.0882978445	-1	33	-0.4769230769	1	13	11	2	0.358974359	TATATATACACATATATATGTATATATATATACACATATATATGTATATATATATACACATATATGTGTGTATATATA");
		RNAzLauncher rna = test.getRNAzLauncher(mafParser);
		test.setDnDs(mafParser);
		System.out.println(rna.getSCI() + " " + rna.getZScore());
		System.out.println(test.dnds);
		
	}
}
