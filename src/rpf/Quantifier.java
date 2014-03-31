package rpf;

import java.util.ArrayList;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BedCovFileParser;
import parser.ScoringOutputParser.ScoredPosition;
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
	
	public boolean isAbundant(ScoredPosition position, double RPKM, int maxLength, int offset){
		return RPKM < getPositionRPKM(position.getContig(), position.getPosition(), position.isPlusStrand(), maxLength, offset);
	}
	
	
	public double getCDSRPKM(AnnotatedGene gene){
		BedCovFileParser bedCovFileParser = gene.isPlusStrand()? bedCovPlusFileParser : bedCovMinusFileParser;
		return (bedCovFileParser.getTotalCDSCoverage(gene, true) * 1e9) / (bedCovFileParser.getTotalReadCount()+1);
	}
	//10^9*C/NL
	public double getPositionRPKM(String contig, int position, boolean isPlusStrand, int maxLength, int offset){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return (bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position + (isPlusStrand? -offset : offset), offset, maxLength+offset, fastaFileParser, true) * 1e9) / (bedCovFileParser.getTotalReadCount()+1);
	}
	
	public double getPositionCount(String contig, int position, boolean isPlusStrand, int maxLength, int offset){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position + (isPlusStrand? -offset : offset), offset, maxLength+offset, fastaFileParser, false);
	}
	
	public ArrayList<Double> getNextStopCodonQuantityChangeRatioNStopPosition(String contig, int position, boolean isPlusStrand, int length, int maxLength){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
	//	position = isPlusStrand? position - offset : position + offset;
		ArrayList<Double>covs = bedCovFileParser.getCoverageBeforeNAfternextStopCodonNStopCodonPosition(contig, isPlusStrand, position, length, 0, maxLength, fastaFileParser, false);
		ArrayList<Double> ret = null;
		
		if(covs != null){
			ret = new ArrayList<Double>();
			//if(covs.get(0)+covs.get(1)<2) ret.add(-1.0);
			//else 
				ret.add((covs.get(0)+2)/(covs.get(1)+2));
			ret.add(covs.get(2));
		}			
		return ret;
	}
	
	//  /([length of transcript]/1000)/([total reads]/10^6)
	public double getPositionQuantityChangeRatio(String contig, int position, boolean isPlusStrand, int length, int offset){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		//int position = isPlusStrand? position - offset : position + offset;
		int prevPosition = isPlusStrand? position - offset - length : position + offset + length;
		int afterPosition = isPlusStrand? position - offset : position + offset;
		double qb = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, prevPosition, length, false);
		double qa = bedCovFileParser.getTotalCoverage(contig, isPlusStrand, afterPosition, length, false); 
		
		//if(qa + qb < 10) return -1;
		
		return (qa+2)/(qb+2);
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
	
	public static void main(String[] args) {
		Quantifier test = new Quantifier("/media/kyowon/Data1/RPF_Project/samples/sample3/coverages/Harr_C-uncollapsed.plus.cov"
				,"/media/kyowon/Data1/RPF_Project/samples/sample3/coverages/Harr_C-uncollapsed.minus.cov"
				, new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.refFlat.txt"),
				new ZeroBasedFastaParser("/media/kyowon/Data1/RPF_Project/genomes/mm9.fa"));
		
		System.out.println(test.getPositionCount("chr7", 6324070, true, 50, 13));
		System.out.println(test.getPositionRPKM("chr7", 6324070, true, 50, 13));
		
		//AnnotatedGene gene = new AnnotatedGene("SRXN1	NM_080725	chr20	-	627267	634014	629357	633829	2	627267,633619,	629561,634014,");
		
		
	}

}
