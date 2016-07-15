package rpf;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import parser.AnnotationFileParser.AnnotatedGene;
import rpf.parser.ScoringOutputParser.ScoredPosition;

public class Quantifier {

	//private BedCovFileParser bedCovPlusFileParser;
	//private BedCovFileParser bedCovMinusFileParser;
	//private AnnotationFileParser annotationFileParser = null;
	//private Bed12Parser bedParser;
	//private ZeroBasedFastaParser fastaFileParser;
	//private double orfQuantity;
	//private double uorfQuantity;
	private SamReader reader;
	private int totalReadCount = 0;
	
	public Quantifier(SamReader reader){
		this.reader = reader;
		SAMRecordIterator iterator = reader.iterator();
		while(iterator.hasNext()){
			iterator.next();
			totalReadCount++;
		}
		iterator.close();
	}
	
	public boolean isAbundant(ScoredPosition position, double RPKM){
		return RPKM < getPositionRPKM(position.getContig(), position.isPlusStrand(), position.getCoordinate());
	}
	
	
	public double getCDSRPKM(AnnotatedGene gene){
		int[] s = gene.getExonStarts(); // 0 based
		int[] e = gene.getExonEnds(); // exclusive
		ArrayList<Integer> coordinate = new ArrayList<Integer>();
		for(int i=0;i<s.length;i++){
			for(int p = s[i];p<e[i];p++){
				if(p<gene.getCdsStart() || p>=gene.getCdsEnd()) continue;
				coordinate.add(p);
			}
		}
		if(!gene.isPlusStrand()) Collections.sort(coordinate, Collections.reverseOrder());
		return getPositionRPKM(gene.getContig(), gene.isPlusStrand(), coordinate);		
	}
	//10^9*C/NL
	public double getPositionRPKM(String contig, boolean isPlusStrand, ArrayList<Integer> coordinate){ // coordinate should take care of offset ..
		double sum = 0;
		if(coordinate.isEmpty()) return 0;
		int sp = coordinate.get(0)+1; // 1 based inclusive
		int ep = coordinate.get(coordinate.size()-1)+1; // 1 based inclusive
		ArrayList<Integer> cc = new ArrayList<Integer>(coordinate);
		if(!isPlusStrand){
			int tp = sp;
			sp = ep;
			ep = tp;
			Collections.sort(cc, Collections.reverseOrder());
		}
		
		SAMRecordIterator iterator = reader.query(contig, sp, ep, false);
		
		while(iterator.hasNext()){
			SAMRecord r = iterator.next();
			
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
			ArrayList<Integer> ap = new ArrayList<Integer>();
			for(int i=0;i<blocks.size();i++){
				AlignmentBlock block = blocks.get(i);
				for(int b = block.getReferenceStart()-1; b < block.getReferenceStart()+block.getLength()-1;b++){					
					if(b < sp-1 || b > ep-1) continue;
					ap.add(b);
				}
			}
			if(ap.isEmpty()) continue;
			int i = Collections.binarySearch(cc, ap.get(0));
			if(i<0) continue;
			//System.out.println(sp + " " + ep + " " + ap);
			boolean matched = true;
			for(int j=0;i<cc.size();j++,i++){
				if(j >= ap.size()){
					matched = false;
					break;
				}
				if(ap.get(j).equals(cc.get(i))) continue;
				matched = false;
				break;
			}
			//System.out.println(matched + " " + r.getReadName());
			if(matched) sum++;
		}		
		iterator.close();	
		//System.out.println(sum);
		return ((sum * 1e9 + 1) / (cc.size()+1)) / (totalReadCount);
	}
	
	/*public double getPositionCount(String contig, int position, boolean isPlusStrand, int maxLength, int offset){
		BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;
		return bedCovFileParser.getTotalCoverageTillnextStopCodon(contig, isPlusStrand, position + (isPlusStrand? -offset : offset), offset, maxLength+offset, fastaFileParser, false);
	}*/
	
	public double getReleaseScore(AnnotatedGene gene, int stopPosition, int length){
		ArrayList<Integer> rcoordinate = gene.getLiftOverPositions(stopPosition, 0, length, false);		
		ArrayList<Integer> lcoordinate = gene.getLiftOverPositions(stopPosition + (gene.isPlusStrand()? -1 : 1), length, 0, false);		
		
		double l = getPositionRPKM(gene.getContig(), gene.isPlusStrand(), rcoordinate);
		double r = getPositionRPKM(gene.getContig(), gene.isPlusStrand(), lcoordinate);		
		
		return (l+.1)/(r+.1);
	}
	
	
	public double getPositionQuantityChangeRatio(ScoredPosition position, int length, int offset){		
		ArrayList<Integer> coordinate = position.getCoordinate();
		ArrayList<Integer> bcoordinate = new ArrayList<Integer>();
		ArrayList<Integer> acoordinate = new ArrayList<Integer>();
		
		int si = coordinate.indexOf(position.getPosition());
		//if(si<0) System.out.println(position);
		for(int i = si - offset - length; i<si - offset; i++){
			if(i>=0 && i<coordinate.size())
				bcoordinate.add(coordinate.get(i));
		}
		
		for(int i = si + offset; i<si + offset + length; i++){
			if(i>=0 && i<coordinate.size())
				acoordinate.add(coordinate.get(i));
		}
		double qa = getPositionRPKM(position.getContig(), position.isPlusStrand(), acoordinate);
		double qb = getPositionRPKM(position.getContig(), position.isPlusStrand(), bcoordinate);	
		
		return (qa+.1)/(qb+.1);
	}
	
	// only for NM_ORF T.. 
	public double getCDSRPKMChangeRatio(AnnotatedGene gene, int position){
		if(gene == null) return -1;
			
		ArrayList<Integer> coordinatea = gene.getLiftOverPositions(position, 0, Integer.MAX_VALUE, true);
		double qa =  getPositionRPKM(gene.getContig(), gene.isPlusStrand(), coordinatea);				
				
		int start = gene.isPlusStrand() ? gene.getCdsStart() : gene.getCdsEnd()-1;
		int rw = Math.abs(position - start);
		ArrayList<Integer> coordinateb= gene.getLiftOverPositions(start, 0, rw , false);
		double qb =  getPositionRPKM(gene.getContig(), gene.isPlusStrand(), coordinateb);				
		
		return qa/qb;
	}
	

	public static void main(String[] args){
		SamReader reader = SamReaderFactory.makeDefault().open(new File("/media/kyowon/Data1/RPF_Project/samples/sample5/mapped/bam/mSeq-C_1-uncollapsed.sorted.bam"));
		Quantifier test = new Quantifier(reader);
		//chr1	60139462	60140609	0_5437344-1	1	-	60139462	60140609	255,0,0	2	31,11	0,1136
		//chr1	60590217	60623422	0_3289552-2	255	-	60590217	60623422	255,0,0	2	33,14	0,33191

		ArrayList<Integer> coordinate = new ArrayList<Integer>();
		//coordinate.add(60590248);
		coordinate.add(60590249);
		coordinate.add(60623408);
		coordinate.add(60623409);
		coordinate.add(60623410);
		coordinate.add(60623411);
		coordinate.add(60623412);
		coordinate.add(60623413);
		coordinate.add(60623414);
		coordinate.add(60623415);
		coordinate.add(60623416);
		coordinate.add(60623417);
		coordinate.add(60623418);
		coordinate.add(60623419);
		coordinate.add(60623420);
		coordinate.add(60623421);
		//coordinate.add(60623422);
		test.getPositionRPKM("chr1", true, coordinate);
	}

}
