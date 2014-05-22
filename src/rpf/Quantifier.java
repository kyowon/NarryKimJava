package rpf;

import java.util.ArrayList;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.Bed12Parser;
import parser.ZeroBasedFastaParser;
import rpf.parser.ScoringOutputParser.ScoredPosition;

public class Quantifier {

	//private BedCovFileParser bedCovPlusFileParser;
	//private BedCovFileParser bedCovMinusFileParser;
	//private AnnotationFileParser annotationFileParser = null;
	private Bed12Parser bedParser;
	private ZeroBasedFastaParser fastaFileParser;
	//private double orfQuantity;
	//private double uorfQuantity;
	
	
	public Quantifier(Bed12Parser bedParser, AnnotationFileParser annotationFileParser, ZeroBasedFastaParser fastaFileParser){
		this.bedParser = bedParser; 
		this.fastaFileParser = fastaFileParser;
	}
	
	public boolean isAbundant(ScoredPosition position, double RPKM){
		if(!position.getContig().equals(bedParser.getContig())) return false;
		return RPKM < getPositionRPKM(position.isPlusStrand(), position.getCoordinate());
	}
	
	
	public double getCDSRPKM(AnnotatedGene gene){
		double sum = 0;
		for(double c : bedParser.getCDSCoverages(gene)) sum += c;
		return (sum * 1e9) / (bedParser.getTotalReadCount()+1);
	}
	//10^9*C/NL
	public double getPositionRPKM(boolean isPlusStrand, ArrayList<Integer> coordinate){ // coordinate should take care of offset ..
		double sum = 0;
		for(double c : bedParser.get5pCoverages(isPlusStrand, coordinate)) sum += c;
		return (sum * 1e9) / (bedParser.getTotalReadCount()+1);
	}
	
	/*public double getPositionCount(String contig, int position, boolean isPlusStrand, int maxLength, int offset){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position + (isPlusStrand? -offset : offset), offset, maxLength+offset, fastaFileParser, false);
	}*/
	
	public double getReleaseScore(AnnotatedGene gene, int stopPosition, int length){
		ArrayList<Integer> coordinate = gene.getLiftOverPositions(stopPosition, length, length, false);
		double[] cov = bedParser.get5pCoverages(gene.isPlusStrand(), coordinate);
		double l = 0;
		double r = 0;
		for(int i=0;i<cov.length;i++){
			if(i < length) l+=cov[i];
			else r+=cov[i];
		}
		
		return (l+1)/(r+1);
	}
	
	
	public double getPositionQuantityChangeRatio(ScoredPosition position, int length, int offset){
		
		ArrayList<Integer> coordinate = position.getCoordinate();
		double[] cov = bedParser.get5pCoverages(position.isPlusStrand(), coordinate);
		double qa =0, qb = 0;
		int si = coordinate.indexOf(position.getPosition());
		for(int i = si - offset - length; i<si - offset; i++){
			if(i>=0 && i<cov.length)
				qb += cov[i];
		}
		
		for(int i = si + offset; i<si + offset + length; i++){
			if(i>=0 && i<cov.length)
				qa += cov[i];
		}
		//if(qa + qb < 10) return -1;
		
		return (qa+2)/(qb+2);
	}
	
	// only for NM_ORF T.. 
	public double getCDSRPKMChangeRatio(AnnotatedGene gene, int position){
		if(gene == null) return -1;
			
		ArrayList<Integer> coordinatea= gene.getLiftOverPositions(position, 0, Integer.MAX_VALUE, true);
		double qa = 0;
		for(double c : bedParser.get5pCoverages(gene.isPlusStrand(), coordinatea)) qa += c;
		
		int start = gene.isPlusStrand() ? gene.getCdsStart() : gene.getCdsEnd()-1;
		int rw = Math.abs(position - start);
		ArrayList<Integer> coordinateb= gene.getLiftOverPositions(start, 0, rw , false);
		double qb = 0;
		for(double c : bedParser.get5pCoverages(gene.isPlusStrand(), coordinateb)) qb += c;
		
		return qa/qb;
	}
	



}
