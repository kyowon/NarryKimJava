package rpf;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BedCovFileParser;
import parser.ZeroBasedFastaParser;

public class Quantifier {

	private BedCovFileParser bedCovPlusFileParser;
	private BedCovFileParser bedCovMinusFileParser;
	//private AnnotationFileParser annotationFileParser = null;
	private ZeroBasedFastaParser fastaFileParser;
	//private double orfQuantity;
	//private double uorfQuantity;
	
	
	public Quantifier(String bedCovPlusFile, String bedCovMinusFile, String annotationFile, String fastaFile){
		bedCovPlusFileParser = new BedCovFileParser(bedCovPlusFile, annotationFile);
		bedCovMinusFileParser = new BedCovFileParser(bedCovMinusFile, annotationFile);	
		fastaFileParser = new ZeroBasedFastaParser(fastaFile);
	}
	
	public double getCDSQuantity(AnnotatedGene gene){
		BedCovFileParser bedCovFileParser = gene.isPlusStrand()? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCDSCoverage(gene);
	}
	
	public double getPositionQuantity(String contig, int position, boolean isPlusStrand){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position, fastaFileParser);
	}
	
	public double getPositionQuantatyChangeRatio(String contig, int position, boolean isPlusStrand, int length){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		int prevPosition = isPlusStrand? position - length : position + length;
		double qb = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, prevPosition, length) + 1;
		double qa = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, position, length) + 1; 
		return qa/qb;
	}
	
	public BedCovFileParser getBedCovPlusFileParser() {
		return bedCovPlusFileParser;
	}


	public BedCovFileParser getBedCovMinusFileParser() {
		return bedCovMinusFileParser;
	}
	
	public static void main(String[] args) {
		Quantifier test = new Quantifier("/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/RPF6_Thy_RPF_1-uncollapsed.plus.cov"
				,"/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/RPF6_Thy_RPF_1-uncollapsed.minus.cov"
				, "/media/kyowon/Data1/RPF_Project/genomes/refFlatHuman.txt",
				"/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		
		AnnotatedGene gene = new AnnotatedGene("SRXN1	NM_080725	chr20	-	627267	634014	629357	633829	2	627267,633619,	629561,634014,");
		System.out.println(test.getPositionQuantity("chr20", 306568, true));
		
		System.out.println(test.getPositionQuantity("chr20", 633829 - 1 , false));
		
		System.out.println(test.getCDSQuantity(gene));
		
	}

}
