package rpf.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import parser.AnnotationFileParser;
import parser.BedCovFileParser;
import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;
import parser.ZeroBasedFastaParser;
import rpf.MakeProteinFastaAndAnnotationFile;
import rpf.Scorer;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;


public class ReadLengthDistribution
    {
    private File bamFile=null;
    private SAMFileReader inputSam=null;
    private ScoringOutputParser scoringOutputParser = null;
    private ZeroBasedFastaParser fastaParser = null;
    private BedCovFileParser bedCovPlusFileParser, bedCovMinusFileParser;
    private HashMap<String, ArrayList<ScoredPosition>> positionMap = null;
    
    public ReadLengthDistribution(String bamFileName, String scoringOutputFileName, String fasta, String bedCovPlusFile, String bedCovMinusFile, AnnotationFileParser annotationFile)
    {
        this.bamFile=new File(bamFileName);
        this.scoringOutputParser = new ScoringOutputParser(scoringOutputFileName);
        this.inputSam=new SAMFileReader(this.bamFile);
        this.fastaParser = new ZeroBasedFastaParser(fasta);
        this.inputSam.setValidationStringency(ValidationStringency.SILENT);
        bedCovPlusFileParser = new BedCovFileParser(bedCovPlusFile, annotationFile);
		bedCovMinusFileParser = new BedCovFileParser(bedCovMinusFile, annotationFile);
        getPositionMap();
    }

    public void close()
    {
        if(inputSam!=null)
            {
            inputSam.close();
            }
        this.inputSam=null;
        this.bamFile=null;
    }
  
    private void getPositionMap(){
    	/*positionMap = new HashMap<String, ArrayList<ScoredPosition>>();
    	String key = "_etc";
    	for(ScoredPosition position : scoringOutputParser.getPositions()){
	    	boolean isPlusStrand = position.isPlusStrand();
			int pos = position.getPosition(); 
						
			if(position.getGeneName() == null) key = "IG";
			else{
				if(position.isAnnotated()) key = "ANNO";
				else if(position.getGBGeneName().startsWith("LINC")) key = "LINC";
				else{
					if(pos >= position.getCdsStart() && pos < position.getCdsEnd()) key = "dORF";
					else{
						if(isPlusStrand && pos < position.getCdsStart() - 100) key = "uORF";
						else if(!isPlusStrand && pos > position.getCdsEnd() + 100) key = "uORF";
					}
				}				
			}
			if(!positionMap.containsKey(key)) positionMap.put(key, new ArrayList<ScoredPosition>());
			positionMap.get(key).add(position);
    	}*/
    }
    
    private void normalizeDistribution(double[][] dist){
    	for(int i=0;i<dist.length;i++){
    		double sum = 0;
    		for(int j=0;j<dist[i].length;j++){
    			sum += dist[i][j];
    		}
    		if(sum!=0){
    			for(int j=0;j<dist[i].length;j++){
    				dist[i][j] /= sum;
    			}
    		}
    	}
    }
    
    private double[][] getLengthDistribution(double[] signal, ScoredPosition position, int leftWindowSize, int rightWindowSize, boolean is5primeSide){    	
    	double[][] dist = new double[rightWindowSize + leftWindowSize][33];
    	String contig = position.getContig();    	
    	boolean isPlusStrand = position.isPlusStrand();
    	int start, end;
    	int pos = -1;
    	
    	if(is5primeSide){
    		pos = position.getPosition();  
    	//	if(isPlusStrand)System.out.println(fastaParser.getSequence(contig, pos, pos+3));
    	//	else System.out.println(fastaParser.getSequence(contig, pos-2, pos+1));
    	}    	
    	else{
    		HashSet<String> stopCodons = new HashSet<String>();
    		stopCodons.add("TAG");
    		stopCodons.add("TAA");
    		stopCodons.add("TGA");
    		
    		if(isPlusStrand){
    			//if(position.isAnnotated()){
    			//	pos = position.getCdsEnd() - 3;    			
    				
    			//}else
    			{
    				//System.out.println(fastaParser.getSequence(contig, position.getPosition(), position.getPosition()+3));
	    			for(int i=0;i<1500;i+=3){
	    				int t = position.getPosition() + i;
	    				if(stopCodons.contains(fastaParser.getSequence(contig, t, t+3))){
	    					pos = t;
	    					break;
	    				}
	    			}
    			}
    			//System.out.println(fastaParser.getSequence(contig, pos, pos+3));
    		}else{
    			//if(position.isAnnotated()){
    			//	pos = position.getCdsStart() + 2;    
    				//System.out.println(MakeProteinFastaAndAnnotationFile.getComplementaryCodon(fastaParser.getSequence(contig, pos-1, pos+2)));
    			//}else
    			{
    				//System.out.println(MakeProteinFastaAndAnnotationFile.getComplementaryCodon(fastaParser.getSequence(contig, position.getPosition()-2, position.getPosition()+1)));
	    			for(int i=0;i<1500;i+=3){
	    				int t = position.getPosition() - i;
	    				//System.out.println(MakeProteinFastaAndAnnotationFile.getComplementaryCodon(fastaParser.getSequence(contig, t-2, t+1)));
	    				if(stopCodons.contains(MakeProteinFastaAndAnnotationFile.getComplementaryCodon(fastaParser.getSequence(contig, t-2, t+1)))){
	    					pos = t;
	    					break;
	    				}
	    			}
	    		}   
    			
    		}
    	}
    	
    	if(pos < 0) return dist;
    	BedCovFileParser bedCovFileParser = isPlusStrand? bedCovPlusFileParser : bedCovMinusFileParser;		
    	double[] s = bedCovFileParser.getCoverages(contig, pos, leftWindowSize, rightWindowSize, isPlusStrand);
    	Scorer.normalize(s);
    	for(int i=0;i<signal.length;i++){
    		signal[i] = s[i];
    	}
    	
    	if(isPlusStrand){
			start = pos - leftWindowSize;
			end = pos + rightWindowSize;
			//System.out.println("+ " +  fastaParser.getSequence(contig, position.getPosition(), position.getPosition()+3));
		}else{
			end = pos + leftWindowSize + 1;
			start = pos - rightWindowSize + 1;
			//System.out.println("- " + MakeProteinFastaAndAnnotationFile.getComplementaryCodon(fastaParser.getSequence(contig, position.getPosition()-2, position.getPosition()+1)));
		}	       	
    	
    	SAMRecordIterator itr = inputSam.queryOverlapping(contig, start, end);
    	HashMap<Integer, ArrayList<Integer>> lengths = new HashMap<Integer, ArrayList<Integer>>();
 
    	while(itr.hasNext()){
    		SAMRecord sam = itr.next();
    		int offset = isPlusStrand? sam.getAlignmentStart() - 1 - pos :  pos - sam.getAlignmentEnd() + 1; //
    		offset += leftWindowSize;
    		if(offset < 0 || offset >= dist.length) continue;
    		if(!lengths.containsKey(offset)) lengths.put(offset, new ArrayList<Integer>());
    		lengths.get(offset).add(sam.getReadLength());
    		//if(sam.getAlignmentEnd() - sam.getAlignmentStart() + 1 > 40) System.out.println(sam.getSAMString());
    	}
    	itr.close();
    	for(int offset : lengths.keySet()){
    		for(int length : lengths.get(offset)){
    			//System.out.println(offset + " " + length);
    			dist[offset][length]++;
    		}
    	}    	
    	normalizeDistribution(dist);
    	return dist;
    }
    
    private double[][] getAverageDistribution(ArrayList<double[][]> dists){
    	double[][] avgDist = new double[dists.get(0).length][dists.get(0)[0].length];
    	for(double[][] dist : dists){
    		for(int i=0;i<avgDist.length;i++){
    			for(int j=0;j<avgDist[i].length;j++){
    				avgDist[i][j] += dist[i][j];
    			}
    		}
    	}
    	normalizeDistribution(avgDist);
    	return avgDist;
    }
    
    public double[][] getAverageDistribution(double[] avgSignal, String key, double scoreThreshold, int leftWindowSize, int rightWindowSize, boolean is5primeSide){
    	ArrayList<double[][]> dists = new ArrayList<double[][]>();
    	for(ScoredPosition position : positionMap.get(key)){
    		if(position.getScore() < scoreThreshold) continue;
    		//if(!position.getCodon().equals("CTG")) continue;
    		double[] signal = new double[avgSignal.length];
    		dists.add(getLengthDistribution(signal, position, leftWindowSize, rightWindowSize, is5primeSide));
    		for(int i=0;i<avgSignal.length;i++){
    			avgSignal[i] += signal[i];
    		}
    	}
    	Scorer.normalize(avgSignal);
    	return getAverageDistribution(dists);
    }

    
    /*
    public static void main(String[] args){
    	String key = "ANNO";
    	String sampleKey = "RPF6_Thy_RPF_1-uncollapsed";
    	//sampleKey="NS_Harr_10mHsum-uncollapsed";
    	String fasta = "/media/kyowon/Data1/RPF_Project/genomes/hg19.fa";
    	String bam = "/media/kyowon/Data1/RPF_Project/samples/sample1/mapped/bam/"+sampleKey+".sorted.bam";
    	String scored = "/media/kyowon/Data1/RPF_Project/samples/sample1/coverages5/Thy_Harr_10mHsum-uncollapsed.plus.cov.score.tsv";
    	
    	String covFilePlus = "/media/kyowon/Data1/RPF_Project/samples/sample1/coverages5/"+sampleKey+".plus.cov";
		String covFileMinus = "/media/kyowon/Data1/RPF_Project/samples/sample1/coverages5/"+sampleKey+".minus.cov";
		
		int leftWindowSize = 50;
		int rightWindowSize = 50;
		
    	double scoreThreshold = 2.0;
    	ReadLengthDistribution test = new ReadLengthDistribution(bam, scored, fasta, covFilePlus, covFileMinus);
    	double[] avgSignal = new double[leftWindowSize + rightWindowSize];
    	double[][] dist = test.getAverageDistribution(avgSignal, key, scoreThreshold, leftWindowSize, rightWindowSize, false);
    	
    	for(int i=0;i<dist.length;i++){
    		//System.out.print(i + ":\t");
    		for(int j=0;j<dist[i].length;j++){
    			System.out.print(dist[i][j] + " ");
    		}
    		System.out.println();
    	}
    	
    	System.out.println();
    	
    	for(int i=0;i<dist.length;i++){
    		System.out.print(i + "\t");
    		double sum = 0;
    		for(int j=0;j<dist[i].length;j++){
    			sum += +dist[i][j] * j;
    		}
    		System.out.println(sum);
    	}    	
    	
    	System.out.println();
    	
    	for(int i=0;i<avgSignal.length;i++){
    		System.out.println(avgSignal[i]);    		
    	}    
    	
    	test.close();
    }*/
}