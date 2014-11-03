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
		String outMafFolder = "/media/kyowon/Data1/fCLIP/genomes/mafSelected/";
		HashSet<String> species = new HashSet<String>();
	
		species.add("hg19");
		species.add("panTro2");// chimp
		species.add("gorGor1"); // gorilla
		species.add("ponAbe2"); // orang
		species.add("rheMac2"); // rhesus
		species.add("papHam1"); // baboon
		species.add("calJac1"); // marmoset
		species.add("tarSyr1"); // Tarsier
		species.add("micMur1"); // mouse lemur
		species.add("otoGar1");//bush baby
		
		/*
		species.add("mm10");
		species.add("rn5");
		species.add("dipOrd1");
		species.add("hetGla2");
		species.add("cavPor3");
		species.add("speTri2");
		species.add("oryCun2");
		species.add("ochPri2");
		species.add("hg19");
		species.add("panTro4");
		species.add("gorGor3");
		species.add("ponAbe2");
		species.add("nomLeu2");
		species.add("rheMac3");
		species.add("papHam1");
		species.add("calJac3");
		species.add("saiBol1");
		species.add("tarSyr1");
		species.add("micMur1");
		species.add("otoGar3");
		species.add("tupBel1");
		*/
		
		for(File originalMaf : new File(originalMafFolder).listFiles()){
			if(!originalMaf.getName().endsWith(".maf")) continue;
			String outMaf = outMafFolder + originalMaf.getName();
			//if(new File(outMaf).exists()) continue;
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
