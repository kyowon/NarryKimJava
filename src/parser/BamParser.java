package parser;

import fCLIP.analysis.NewCheckCoverage;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


public class BamParser {


	
	public static void main(String[] args) {
		{
			final SamReader reader = SamReaderFactory.makeDefault().open(new File("/media/kyowon/Data1/fCLIP/samples/sample4/alignments/Set-6-siControl-1.pe.sorted.bam"));
			int p1 = 64037695-1;
			int p2 = p1+1;
				

//49012909	49012968



		//	int[] a  = NewCheckCoverage.getPassingThroughReadPairCount("chr3", 49012909, 49012968, reader, 200);
			
		//	System.out.println(a[0]);
			//System.out.println(a[1]);
			
			/*
			SAMRecordIterator iterator = reader.query("chr20", p1-300, p2+300, false);
			
			
			ArrayList<SAMRecord> srs = new ArrayList<SAMRecord>();
			//chr21 82170491 82170991
			while(iterator.hasNext()){				
				SAMRecord rec = iterator.next();
				if(!rec.getReadPairedFlag() || !rec.getFirstOfPairFlag()) continue;
								
				List<AlignmentBlock> blocks = rec.getAlignmentBlocks();
								
				
				System.out.println(" " + rec + " " + rec.getAlignmentStart() + " " + rec.getAlignmentEnd() + " " + rec.getReferenceName());
			//	for(AlignmentBlock block : blocks){
			//		System.out.println(block.getReferenceStart() + " " + (block.getReferenceStart() + block.getLength() - 1));
			//	}
				
				//	
			//	for(SAMTagAndValue t : rec.getAttributes()){
			////		System.out.println(t.tag);
			//	}
			
				srs.add(rec);
			}
			iterator.close();
			
			for(SAMRecord rec : srs){	
				System.out.println(reader.queryMate(rec).getReadName() + " " + rec.getPairedReadName());
				
			//	System.out.println(NewCheckCoverage.isPairPassingThrough(reader, rec, p1, 300));
			}
			System.out.println(srs.size());
			*/
			
		}
	}

}
