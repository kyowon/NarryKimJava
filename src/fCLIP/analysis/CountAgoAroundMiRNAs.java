package fCLIP.analysis;

import fCLIP.FCLIP_Scorer;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import parser.BufferedLineReader;
import parser.MirGff3FileParser;
import parser.MirGff3FileParser.MiRNA;

public class CountAgoAroundMiRNAs {

	
	public static void main(String[] args) throws IOException {
		String miRNAsorted = "/media/kyowon/Data1/fCLIP/pri-miRNAs_sorted.txt"; // retreive ... paired or unpaired
		String miRNAgff3 = "/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3";
		String bamfile = "/media/kyowon/Data1/fCLIP/Data/FLAG_AGO_small_RNA.bam";
		String csvOut = "/media/kyowon/Data1/fCLIP/Figure1/agoCountMiRNA.csv";
		MirGff3FileParser mirParser = new MirGff3FileParser(miRNAgff3);
		PrintStream out = new PrintStream(csvOut);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamfile));
		
		BufferedLineReader oin = new BufferedLineReader(miRNAsorted);
		HashSet<String> miRNAs = new HashSet<String>();
		String s;
		while((s=oin.readLine())!=null){
			String[] token = s.split("\t");
			if(!token[1].equals(".")) continue;
			miRNAs.add(token[0].trim());
		}
		oin.close();
		Iterator<MiRNA> miIterator = mirParser.getMiRNAIterator();
		
		while(miIterator.hasNext()){
			MiRNA mi = miIterator.next();
			
			if(!miRNAs.contains(mi.getName())) continue;
			Integer sp = mi.isPlusStrand()? mi.getPri5p() : mi.getPri3p();
			Integer lp = mi.isPlusStrand()? mi.getPri3p() : mi.getPri5p();
			
			if(sp == null) sp = lp - 160;
			if(lp == null) lp = sp + 160;
			
			if(sp >= lp) System.exit(1);
			
			SAMRecordIterator iterator = reader.query(mi.getContig(), sp+1, lp+1, true);
			int n = 0;
			try{
			while(iterator.hasNext()){
				SAMRecord r = iterator.next();
				n++;
			}
			}finally{
				iterator.close();
			}
			System.out.println(mi.getName() + "\t" + n);
			out.println(mi.getName() + "\t" + n);
		}
		out.close();
	}

}
