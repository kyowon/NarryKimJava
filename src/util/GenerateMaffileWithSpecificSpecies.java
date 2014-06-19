package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class GenerateMaffileWithSpecificSpecies {

	public static void main(String[] args) throws IOException {
		String originalMafFolder = "/media/kyowon/Data1/RPF_Project/genomes/hg19/maf/";
		String outMafFolder = "/media/kyowon/Data1/fCLIP/genomes/maf/";
		HashSet<String> species = new HashSet<String>();
		species.add("hg19");
		species.add("mm9");
		
		for(File originalMaf : new File(originalMafFolder).listFiles()){
			if(!originalMaf.getName().endsWith(".maf")) continue;
			String outMaf = outMafFolder + originalMaf.getName();
			if(new File(outMaf).exists()) continue;
			BufferedLineReader in = new BufferedLineReader(originalMaf.getAbsolutePath());
			PrintStream out = new PrintStream(outMaf);
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("s ") || s.startsWith("q ") || s.startsWith("i ")){
					//System.out.println(s);
					String t = s.substring(2);
					for(String spec : species){
						if(t.startsWith(spec)){
							out.println(s);
							break;
						}
					}					 
				}else out.println(s);
			}
			in.close();
			out.close();
		}
		
		
		
		
	}

}
