package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BufferedLineReader;
import util.FoldChangerCalculator;
import fCLIP.Scorer;
import fCLIP.parser.ScoringOutputParser;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;

public class CheckFoldChange {
 
	
	//TODO run bedtool!!
	
	
	public static void main(String[] args) {
		String key = "siDrosha";
		//key = "siD-K-D";
		//key = "siDicer";
		String csv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.csv";
		String pairCsv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.pair.multihit.csv";
		
		String bed1 = "/media/kyowon/Data1/fCLIP/RNAseq/" + key + "_R1_Aligned_Sorted.bed";
		String bed2 = "/media/kyowon/Data1/fCLIP/RNAseq/siControl_R1_Aligned_Sorted.bed";
		String annotationFileName = "/media/kyowon/Data1/fCLIP/genomes/hg19.refFlat.txt";
		String outCsv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out." + key + ".csv";
		String outPairCsv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out. " + key + ".pair.csv";
		String backGround = "/media/kyowon/Data1/Dropbox/back." + key + ".txt";
		
		AnnotationFileParser annotationParser = new AnnotationFileParser(annotationFileName);
		ScoringOutputParser parser = new ScoringOutputParser(csv);
		FoldChangerCalculator fcc = new FoldChangerCalculator(bed1, bed2, annotationParser, true);
		HashMap<String, Double> readCCnt3pMap = new HashMap<String, Double>();
		HashMap<String, Double> readTCnt3pMap = new HashMap<String, Double>();
		
		HashMap<String, Double> readCCnt5pMap = new HashMap<String, Double>();
		HashMap<String, Double> readTCnt5pMap = new HashMap<String, Double>();
		
		int flank = 150;
		
		try {
			/*PrintStream bout = new PrintStream(backGround);
			//bout.println("b=[");
			HashMap<AnnotatedGene, double[]> map = fcc.getReadCounts();
			for(AnnotatedGene gene : map.keySet()){
				double[] rc = map.get(gene);
				bout.println(gene.getAccession() + "," + rc[0] + "," + rc[1]);
			}
			//bout.println("];");
			bout.close();
			
			System.exit(0);
			*/
			PrintStream out = new PrintStream(outCsv);
		//	out.println(fcc.getMedian());
			out.println(ScoredPosition.getHeader() + "\tControl\tTarget");
			for(ScoredPosition sp : parser.getPositions()){
				if(!sp.getContig().equals("chr1")) continue;
				//ArrayList<String> accessions = sp.getContainingGeneAccessions(); 
				out.print(sp);
				//if(accessions != null && !accessions.isEmpty()){
					Double fc3p = .0;//Double.NaN;
					Double ft3p = .0;//Double.NaN;
					Double fc5p = .0;//Double.NaN;
					Double ft5p = .0;//Double.NaN;					
					String k = sp.getContig()+"\t" + sp.isPlusStrand() + "\t" + sp.getThreePposition() + "\t" + sp.getFivePposition();
					
					int sign = sp.isPlusStrand()? 1 : -1;
					if(sp.getThreePposition() >= 0){
						fc3p = fcc.getControlReadCount(sp.getContig(), sp.getThreePposition() - sign * flank, sp.getThreePposition() - sign * 5, sp.isPlusStrand());
						ft3p = fcc.getTargetReadCount(sp.getContig(), sp.getThreePposition() - sign * flank, sp.getThreePposition() - sign * 5, sp.isPlusStrand());
						readCCnt3pMap.put(k, fc3p);
						readTCnt3pMap.put(k, ft3p);
					}
					
					if(sp.getFivePposition() >= 0){
						fc5p = fcc.getControlReadCount(sp.getContig(), sp.getFivePposition() + sign * 5, sp.getFivePposition() + sign * flank, sp.isPlusStrand());
						ft5p = fcc.getTargetReadCount(sp.getContig(), sp.getFivePposition() + sign * 5, sp.getFivePposition() + sign * flank, sp.isPlusStrand());
						readCCnt5pMap.put(k, fc5p);
						readTCnt5pMap.put(k, ft5p);
					}				
					
					/*if(sp.getThreePposition() >= 0 && sp.getFivePposition() >= 0){
					//	fc = fcc.getControlReadCount(annotationParser.getGeneByAccession(accessions.get(0)));
					//	ft = fcc.getTargetReadCount(annotationParser.getGeneByAccession(accessions.get(0)));
						fc = fcc.getControlReadCount(sp.getContig(), sp.getThreePposition() - sign * flank, sp.getFivePposition() + sign * flank, sp.isPlusStrand());
						fc -= fcc.getControlReadCount(sp.getContig(), sp.getThreePposition(), sp.getFivePposition(), sp.isPlusStrand());
						
						ft = fcc.getTargetReadCount(sp.getContig(), sp.getThreePposition() - sign * flank, sp.getFivePposition() + sign * flank, sp.isPlusStrand());
						ft -= fcc.getTargetReadCount(sp.getContig(), sp.getThreePposition(), sp.getFivePposition(), sp.isPlusStrand());
					}*/
					//	fc = fcc.getFoldChange(sp.getContig(), sp.getThreePposition() - sign * Scorer.flankingNTNumber, sp.getFivePposition() + sign * Scorer.flankingNTNumber, sp.isPlusStrand());
					out.println("\t" + (fc3p + fc5p) + "\t" + (ft3p + ft5p));
					
					
				//}else{
				//	out.println("\t%_\t%_");
				//}
			}
			out.close();
			
			BufferedLineReader in = new BufferedLineReader(pairCsv);
			PrintStream pairOut = new PrintStream(outPairCsv);
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("Contig1")){
					pairOut.println(s + "\tControl\tTarget");
					continue;
				}
				String[] token = s.split("\t");
			//	String accession1 = token[8]; 
			//	String accession2 = token[20];
				String contig1 = token[0];
				String contig2 = token[12];
				boolean isPlusStrand1 = token[1].equals("+");
				boolean isPlusStrand2 = token[13].equals("+");
				
				String k1 = contig1 + "\t" + isPlusStrand1 + "\t" + token[2] + "\t" + token[3];
				String k2 = contig2 + "\t" + isPlusStrand2 + "\t" + token[14] + "\t" + token[15];
				
				
				
				pairOut.print(s);
				double fc = readCCnt3pMap.get(k1) + readCCnt5pMap.get(k2);
				double ft = readTCnt3pMap.get(k1) + readTCnt5pMap.get(k2);
				
				pairOut.print("\t" + fc + "\t" + ft);
				
				
				/*if(!accession1.startsWith("_")){
					Double fc = Double.NaN;
					Double ft = Double.NaN;
					
					fc = fcc.getControlReadCount(contig, sp.getThreePposition() - sign * flank, sp.getFivePposition() + sign * flank, sp.isPlusStrand());
					fc -= fcc.getControlReadCount(contig, sp.getThreePposition(), sp.getFivePposition(), sp.isPlusStrand());
					
					ft = fcc.getTargetReadCount(sp.getContig(), sp.getThreePposition() - sign * flank, sp.getFivePposition() + sign * flank, sp.isPlusStrand());
					ft -= fcc.getTargetReadCount(sp.getContig(), sp.getThreePposition(), sp.getFivePposition(), sp.isPlusStrand());
					
					
					//fc = fcc.getControlReadCount(annotationParser.getGeneByAccession(accession1.split(",")[0]));	
					//ft = fcc.getTargetReadCount(annotationParser.getGeneByAccession(accession1.split(",")[0]));	
					pairOut.print("\t" + fc + "\t" + ft);
				}else
					pairOut.print("\t%_\t%_");
				
				
				if(!accession2.startsWith("_")){
					Double fc = Double.NaN;
					Double ft = Double.NaN;
					fc = fcc.getControlReadCount(annotationParser.getGeneByAccession(accession2.split(",")[0]));
					ft = fcc.getTargetReadCount(annotationParser.getGeneByAccession(accession2.split(",")[0]));	
					pairOut.print("\t" + fc + "\t" + ft);
				}else
					pairOut.print("\t%_\t%_");
				*/
				pairOut.println();
			}
			
			
			pairOut.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
