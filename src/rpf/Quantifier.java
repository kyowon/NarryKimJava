package rpf;

import java.util.ArrayList;

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
		return (bedCovFileParser.getTotalCDSCoverage(gene, true) * 1e9 + 1) / (bedCovFileParser.getTotalReadCount()+1);
	}
	//10^9*C/NL
	public double getPositionRPKM(String contig, int position, boolean isPlusStrand, int maxLength){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		int offset = 30;
		return (bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position + (isPlusStrand? -offset : offset), offset, maxLength+offset, fastaFileParser, true) * 1e9 + 1) / (bedCovFileParser.getTotalReadCount()+1);
	}
	
	public ArrayList<Double> getNextStopCodonQuantityChangeRatioNStopPosition(String contig, int position, boolean isPlusStrand, int length, int maxLength){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
	//	position = isPlusStrand? position - offset : position + offset;
		ArrayList<Double>covs = bedCovFileParser.getCoverageBeforeNAfternextStopCodonNStopCodonPosition(contig, isPlusStrand, position, length, 0, maxLength, fastaFileParser, false);
		ArrayList<Double> ret = null;
		
		if(covs != null){
			ret = new ArrayList<Double>();
			if(covs.get(0)+covs.get(1)<10) ret.add(-1.0);
			else ret.add((covs.get(0)+1)/(covs.get(1)+1));
			ret.add(covs.get(2));
		}			
		return ret;
	}
	
	//  /([length of transcript]/1000)/([total reads]/10^6)
	public double getPositionQuantityChangeRatio(String contig, int position, boolean isPlusStrand, int length, int offset){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		//int position = isPlusStrand? position - offset : position + offset;
		int prevPosition = isPlusStrand? position - offset - length : position + offset + length;
		double qb = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, prevPosition, length, false);
		double qa = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, position, length, false); 
		
		if(qa + qb < 10) return -1;
		
		return (qa+1)/(qb+1);
	}
	
	public double getCDSQuantity(AnnotatedGene gene){
		BedCovFileParser bedCovFileParser = gene.isPlusStrand()? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCDSCoverage(gene, false);
	}
	
	/*public double getPositionQuantity(String contig, int position, boolean isPlusStrand, int maxLength){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position, maxLength, fastaFileParser, false);
	}*/
	
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
