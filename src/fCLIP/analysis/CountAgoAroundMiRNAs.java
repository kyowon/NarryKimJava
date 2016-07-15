package fCLIP.analysis;

import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import parser.BufferedLineReader;
import parser.MirGff3FileParser;
import parser.MirGff3FileParser.MiRNA;

public class CountAgoAroundMiRNAs {

	
	public static void main(String[] args) throws IOException {
		String miRNAsorted = "/media/kyowon/Data1/fCLIP/pri-miRNAs_sorted.txt"; // retreive ... paired or unpaired
		String miRNAgff3 = "/media/kyowon/Data1/fCLIP/genomes/hsa_hg38.gff3";
		String agobamfile = "/media/kyowon/Data1/fCLIP/Data/FLAG_AGO_small_RNA.bam";
		String dicerbamfile = "/media/kyowon/Data1/fCLIP/Data/Dicer_PAR-CLIP.bam";
		String csvOut = "/media/kyowon/Data1/fCLIP/Figure4/agodicerCountMiRNA.csv";
		//String csv = "/media/kyowon/Data1/fCLIP/miRNA_ReadEnd_unpairedIncluded.csv";
		String csv = "/media/kyowon/Data1/fCLIP/samples/sample8/results/lists/cis/bamMerged.complete.csv";
		ScoredPositionOutputParser parser = new ScoredPositionOutputParser(csv);		
		
		MirGff3FileParser mirParser = new MirGff3FileParser(miRNAgff3);
		PrintStream out = new PrintStream(csvOut);
		SamReader agoReader = SamReaderFactory.makeDefault().open(new File(agobamfile));
		SamReader dicerReader = SamReaderFactory.makeDefault().open(new File(dicerbamfile));
		
		BufferedLineReader oin = new BufferedLineReader(miRNAsorted);
		HashMap<String, Integer> miRNAs = new HashMap<String, Integer>();
		String s;
		while((s=oin.readLine())!=null){
			String[] token = s.split("\t");
			int kind = 0;
			if(token[1].equals("."));
			else if(token[1].startsWith("Capped")) kind = 2;
			else if(token[1].startsWith("miRTron")) kind = 1;
			else kind = 3;
			miRNAs.put(token[0].trim(), kind);
		}
		oin.close();
		HashSet<String> miRNAsDetected = new HashSet<String>();
		
		for(ScoredPosition sp : parser.getPositions()){
			if(sp.hasMatchingMiRNA()){
				miRNAsDetected.addAll(sp.getMiRNAs());
			}
		}
		
		
		// 3139 4449 6723
		
		Iterator<MiRNA> miIterator = mirParser.getMiRNAIterator();
		out.println("miRNA\tKind(0:miRNA,1:mirtron,2:capped,3:group6)\tDicer read count\tAgo read count\tfCLIP detected");
		while(miIterator.hasNext()){
			MiRNA mi = miIterator.next();
			
			if(!miRNAs.containsKey(mi.getName())) continue;
			Integer sp = mi.isPlusStrand()? mi.get5p() : mi.get3p();
			Integer lp = mi.isPlusStrand()? mi.get3p() : mi.get5p();
			
			if(sp == null) sp = lp - 300;
			if(lp == null) lp = sp + 300;
			
			if(sp >= lp) System.exit(1);
			
			SAMRecordIterator iterator1 = agoReader.query(mi.getContig(), sp+1, lp+1, true);
			int an = 0;
			try{
			while(iterator1.hasNext()){
				iterator1.next();
				an++;
			}
			}finally{
				iterator1.close();
			}
			
			SAMRecordIterator iterator2 = dicerReader.query(mi.getContig(), sp+1, lp+1, true);
			int dn = 0;
			try{
			while(iterator2.hasNext()){
				iterator2.next();
				dn++;
			}
			}finally{
				iterator2.close();
			}
			int kind = miRNAs.get(mi.getName());
			System.out.println(mi.getName() + "\t" + kind + "\t" + dn + "\t"  + an + "\t" + (miRNAsDetected.contains(mi.getName()) ? "1" : "0"));
			out.println(mi.getName() + "\t" + kind + "\t" +  dn + "\t" + an + "\t" + (miRNAsDetected.contains(mi.getName()) ? "1" : "0"));
		}
		out.close();
	}

}
