package rpf;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;
public class MakeProteinFasta {

	public static void main(String[] args) throws IOException {
		HashSet<String> genes = new HashSet<String>();
		String s;
		PrintStream out = new PrintStream("/media/kyowon/Data1/RPF_Project/proteomes/HUMAN_Thy.fasta");
		BufferedLineReader inScore = new BufferedLineReader(("/media/kyowon/Data1/RPF_Project/samples/sample1/coverages/Thy_Harr_10mHsum-uncollapsed.plus.cov.score.tsv"));
		
		while((s=inScore.readLine())!=null){
			String[] token = s.split("\t");
			if(!token[5].equals("ATG") || !token[6].equals("NM_ORF")) continue;
			if(!token[7].equals("0") || token[8].equals("T")) continue;
			if(Double.parseDouble(token[3]) < 2.1) continue;
			if(token[9].equals("_")) continue;
			genes.add(token[9]);			
		}
		
		inScore.close();
		
		BufferedLineReader inFasta = new BufferedLineReader(("/media/kyowon/Data1/RPF_Project/proteomes/HUMAN.fasta"));
		
		boolean toWrite = true;
		while((s=inFasta.readLine())!=null){
			if(s.startsWith(">")){
				String[] token = s.split(" ");
				String gene = null;
				for(String t : token){
					if(t.startsWith("GN=")){
						gene = t.substring(3);
					}
				}				
				if(gene == null || !genes.contains(gene)){
					toWrite = false;
					continue;
				}else toWrite = true;
				
			}
			if(toWrite) out.println(s);
		}
		
		inFasta.close();
		out.close();

	}

}
