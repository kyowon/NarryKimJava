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
	
	
	public Quantifier(String bedCovPlusFile, String bedCovMinusFile, AnnotationFileParser annotationFileParser, ZeroBasedFastaParser fastaFileParser){
		bedCovPlusFileParser = new BedCovFileParser(bedCovPlusFile, annotationFileParser);
		bedCovMinusFileParser = new BedCovFileParser(bedCovMinusFile, annotationFileParser);	
		this.fastaFileParser = fastaFileParser;
	}
	
	public double getCDSRPKM(AnnotatedGene gene){
		BedCovFileParser bedCovFileParser = gene.isPlusStrand()? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCDSCoverage(gene, true) * 1e9 / bedCovFileParser.getTotalReadCount();
	}
	
	public double getPositionRPKM(String contig, int position, boolean isPlusStrand, int maxLength){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position, maxLength, fastaFileParser, true) * 1e9 / bedCovFileParser.getTotalReadCount();
	}
	///([length of transcript]/1000)/([total reads]/10^6)
	public double getPositionQuantatyChangeRatio(String contig, int position, boolean isPlusStrand, int length){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		int prevPosition = isPlusStrand? position - length : position + length;
		double qb = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, prevPosition, length, false) + 1;
		double qa = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, position, length, false) + 1; 
		return qa/qb;
	}
	
	public double getCDSQuantity(AnnotatedGene gene){
		BedCovFileParser bedCovFileParser = gene.isPlusStrand()? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCDSCoverage(gene, false);
	}
	
	public double getPositionQuantity(String contig, int position, boolean isPlusStrand, int maxLength){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position, maxLength, fastaFileParser, false);
	}
	

	
	public BedCovFileParser getBedCovPlusFileParser() {
		return bedCovPlusFileParser;
	}


	public BedCovFileParser getBedCovMinusFileParser() {
		return bedCovMinusFileParser;
	}
	
	/*public static void main(String[] args) {
		Quantifier test = new Quantifier("/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/RPF6_Thy_RPF_1-uncollapsed.plus.cov"
				,"/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/RPF6_Thy_RPF_1-uncollapsed.minus.cov"
				, "/media/kyowon/Data1/RPF_Project/genomes/refFlatHuman.txt",
				"/media/kyowon/Data1/RPF_Project/genomes/hg19.fa");
		
		AnnotatedGene gene = new AnnotatedGene("SRXN1	NM_080725	chr20	-	627267	634014	629357	633829	2	627267,633619,	629561,634014,");
		System.out.println(test.getPositionQuantity("chr20", 306568, true));
		
		System.out.println(test.getPositionQuantity("chr20", 633829 - 1 , false));
		
		System.out.println(test.getCDSQuantity(gene));
		
	}*/

}
