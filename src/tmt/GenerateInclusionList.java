package tmt;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import sequences.ProteinFastaSequence;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Composition;
import msutil.Enzyme;

public class GenerateInclusionList {
	private static ArrayList<Annotation> getAllAnnotations(String protein, int minl, int maxl, Enzyme enzyme, int enzymeOption, int ntt, AminoAcidSet aaSet){
		ArrayList<Annotation> annotations = new ArrayList<Annotation>();
		for(int s=0;s<protein.length()-minl;s++){
			char prevAA = s>0? protein.charAt(s-1): '-';
			for(int e=s+minl;e<=Math.min(protein.length()-1, s+maxl);e++){
				char nextAA = e<protein.length()-1 ? protein.charAt(e+1) : '-';
				
				String anno = prevAA + "." + protein.substring(s, e) + "." + nextAA;
				//System.out.println(anno);
				Annotation annotation = new Annotation(anno, aaSet);
				if(enzyme != null){
					if(annotation.getPeptide() == null || annotation.getPeptide().getNumMissedCleavageSites(enzyme) > ntt) continue;
					if(enzymeOption <2){
						int n = enzymeOption;
						if(annotation.getPeptide().hasCleavageSite(enzyme)) n++;
						if(enzyme.isCTerm()){
							if(prevAA == '-') n++; 
							else if(enzyme.isCleavable(prevAA)) n++;
						}else{
							if(nextAA == '-') n++; 
							else if(enzyme.isCleavable(nextAA)) n++;
						}
						
						if(n<2) continue;
					}
				}
				annotations.add(annotation);
			}
		}
		return annotations;
	}
	
	private static HashMap<Integer,ArrayList<Float>> getExclusionMzList(String fasta, int minl, int maxl, int minc, int maxc,
			Enzyme enzyme, int enzymeOption, int ntt, float maxMass, AminoAcidSet aaSet){
		ProteinFastaSequence fs = new ProteinFastaSequence(fasta);
		long start = 0;
		HashMap<Integer,HashSet<Float>> mzs = new HashMap<Integer,HashSet<Float>>();
		for(int c=minc;c<=maxc;c++){
			mzs.put(c, new HashSet<Float>());
		}
		while(start<fs.getSize()){ // fs.getSize()	
			long str = fs.getStartPosition(start);
			String protein = fs.getMatchingEntry(str);		
			start += protein.length()+1;
			int mass = fs.getIntegerMass(str, str + protein.length());	
			
			if(mass > maxMass) continue;
			
			for(Annotation annotation : getAllAnnotations(protein, minl, maxl, enzyme, enzymeOption, ntt, aaSet)){
				for(int c=minc;c<=maxc;c++){
					mzs.get(c).add((float)(annotation.getPeptide().getParentMass()/c + Composition.PROTON));
				}
			}
		}	
	
		HashMap<Integer,ArrayList<Float>> mzlist = new HashMap<Integer,ArrayList<Float>>();
		for(int c=minc;c<=maxc;c++){
			ArrayList<Float> smzlist = new ArrayList<Float>(mzs.get(c));
			Collections.sort(smzlist);
			mzlist.put(c, smzlist);
		}
		
		//System.out.println("Number of exclusion mz values : " + mzlist.size());
		return mzlist;
	}
	
	private static HashSet<String> getFilteredProteinAnnotations(String fasta, float maxMass){
		ProteinFastaSequence fs = new ProteinFastaSequence(fasta);
		long start = 0;
		HashSet<String> annotations = new HashSet<String>();
		while(start<fs.getSize()){ // fs.getSize()	
			long str = fs.getStartPosition(start);
			String protein = fs.getMatchingEntry(str);		
			start += protein.length()+1;
			int mass = fs.getIntegerMass(str, str + protein.length());	
			
			if(mass > maxMass) continue;
			
			annotations.add(fs.getAnnotation(start));
			
		}	
		return annotations;
	}
	
	private static HashMap<String, HashMap<Integer, Float>> getInclusionMzList(String fasta, String exclusionFasta,
			int minl, int maxl, int minc, int maxc, Enzyme enzyme, int enzymeOption, int ntt, float maxMass, AminoAcidSet aaSet, int numFlankingAA){
		HashMap<String, HashMap<Integer, Float>> mzlist = new HashMap<String, HashMap<Integer, Float>>();
		
		ProteinFastaSequence fs = new ProteinFastaSequence(fasta);
		SuffixArray saExclusion = new SuffixArray(new SuffixArraySequence(exclusionFasta));
		SuffixArray saInclusion = new SuffixArray(new SuffixArraySequence(fasta));
		
		HashSet<String> annotations = getFilteredProteinAnnotations(exclusionFasta, maxMass);
		long start = 0;
		

		while(start<fs.getSize()){ // fs.getSize()	
			long str = fs.getStartPosition(start);
			String protein = fs.getMatchingEntry(str);		
			start += protein.length()+1;
				
			for(Annotation annotation : getAllAnnotations(protein, minl, maxl, enzyme, enzymeOption, ntt, aaSet)){
				ArrayList<String> matchingAnnotations = saExclusion.getAllMatchingAnnotations(annotation.getPeptide().toString());
				boolean toAdd = true;
				for(String ma : matchingAnnotations){
					if(annotations.contains(ma)){
						toAdd = false;
						break;
					}
					
				}
				
				if(!toAdd) continue;
				//annotations
				String mas = "";
				ArrayList<String> matchingAnnotationsInclusion = saInclusion.getAllMatchingAnnotations(annotation.getPeptide().toString());
				for(String ma : matchingAnnotationsInclusion){
					mas += ma + "\t";
				}
				HashMap<Integer, Float> mzPerCharge = new HashMap<Integer, Float>();							
				for(int c=minc;c<=maxc;c++){
					//ArrayList<Float> el = exclusionList.get(c);
					float mz = (float)(annotation.getPeptide().getParentMass()/c + Composition.PROTON);
					
					mzPerCharge.put(c, mz);
				}
				if(!mzPerCharge.isEmpty()){
					for(String annotationWithTwoFlankingAA : saInclusion.getAllMatchedStrings(annotation.getPeptide().toString(), numFlankingAA)){
						//V_D.MLAAEAGLLPVLCALCTCLR._T          annotationWithTwoFlankingAA
						//annotationWithTwoFlankingAA = "_T.MLAAEAGLLPVLCALCTCLR._R";
						int index1 = annotationWithTwoFlankingAA.indexOf('_');
						
						if(index1 >=1 && index1 < numFlankingAA){
							annotationWithTwoFlankingAA = annotationWithTwoFlankingAA.substring(index1);
						}
						
						int index2 = annotationWithTwoFlankingAA.lastIndexOf('_');
						
						if(index2 > annotationWithTwoFlankingAA.lastIndexOf('.')){
							annotationWithTwoFlankingAA = annotationWithTwoFlankingAA.substring(0, 1 + index2);		
						}
						//System.out.println(annotationWithTwoFlankingAA);
						//System.exit(0);
						mzlist.put(mas+annotationWithTwoFlankingAA, mzPerCharge);
						
					}
					
					
				}
			}
		}		
		
		System.out.println("Number of inclusion peptide number : " + mzlist.size());
		return mzlist;
	}
	
	public static void main(String[] args) throws IOException {
		//new Peptide("GHIFKHVTVASCPPVQIEELIEK")
		String infasta = "/media/kyowon/Data1/Dropbox/inclusion/out_0.3_smORF.fasta";
		String exfasta = "/media/kyowon/Data1/Dropbox/inclusion/MOUSE.fasta";
		String outtxt = "/media/kyowon/Data1/RPF_Project/samples/sample5/results/out_0.3_smORF.inclusion.csv"; 
		//fasta="/media/kyowon/Data1/MassSpec/HUMAN.fasta";
		Enzyme enzyme = Enzyme.TRYPSIN;
		//enzyme.isCleaved(p) getNumCleavedTermini
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
	//	Tolerance tol = new Tolerance(1f, true);
		int minc = 2;
		int maxc = 7;
		int minl = 5;
		int maxl = 40;
		int ntt = 0;
		float maxMass = 10000;
		int enzymeOption = 0; // enzymeOption 0 : full, 1 : semi, 2 : all    
		if(args.length < 5){
			System.out.println("java -jar -Xmx2g GenerateInclusionList [includsion fasta] [exclusion fasta] [output csv] [min charge] [max charge]"
					+ "[min length] [max length] [ntt] [max protein mass] [enzyme option : 0 - full try 1 - semi 2 - no enzyme]");
			System.exit(0);
		}else{
			infasta = args[0];
			exfasta = args[1];
			outtxt = args[2];
			//tol = Tolerance.parseToleranceStr(args[3]);
			minc = Integer.parseInt(args[3]);
			maxc = Integer.parseInt(args[4]);
			minl = Integer.parseInt(args[5]);
			maxl = Integer.parseInt(args[6]);
			ntt = Integer.parseInt(args[7]);
			maxMass = Float.parseFloat(args[8]);
			enzymeOption = Integer.parseInt(args[9]);
		}	
		
		
	//	HashMap<Integer, ArrayList<Float>> emzs = getExclusionMzList(exfasta, minl, maxl, minc, maxc, enzyme, 1, 3, maxMass, aaSet);
		HashMap<String, HashMap<Integer, Float>> inclusion = 
				getInclusionMzList(infasta, exfasta, minl, maxl, minc, maxc, enzyme, enzymeOption, ntt, maxMass, aaSet, 2);
		
		
		
		PrintStream out = new PrintStream(outtxt);
		out.print("#\tPeptide");
		for(int c=minc;c<=maxc;c++) out.print("\t"+c);
		out.print("\tH\tM\tC\tProtein ID");
		out.println();
		int i = 0;
		for(String key : inclusion.keySet()){
			// pep = annos + pep
			String[] token = key.split("\t");
			out.print(++i); out.print("\t"+token[token.length-1]);
			HashMap<Integer, Float> sm = inclusion.get(key);
			for(int c=minc;c<=maxc;c++){
				Float mz = sm.get(c);
				if(mz == null) out.print("\t ");
				else out.print("\t"+mz);
			}
			
			String anno = token[token.length-1];//AR.SAGPSR.SC
			String pep = anno.substring(anno.indexOf('.')+1, anno.lastIndexOf('.'));
			out.print("\t" + pep.contains("H")+"\t" + pep.contains("M")+"\t" + pep.contains("C"));
			out.print("\t");
			for(int k=0;k<token.length-1;k++){
				out.print(token[k]+";");
			}
			out.println();
		}
		
		out.close();
	}

}
