package launcher;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Random;

public class PhyloPLauncher {
	
	private double pvalConservation = 1;
	private double pvalAcceleration = 1;
	static private String modFile = null;
	static public void setModFile(String mod){
		modFile = mod;
	}
	
	public PhyloPLauncher(String maf){
		run(maf, modFile);
	}
	
	public double getPvalConservation(){
		return pvalConservation;
	}
	
	public double getPvalAcceleration(){
		return pvalAcceleration;
	}
	
	private void run(String maf, String modFile){
		//System.out.println(modFile);
		if(maf == null || modFile == null) return;
		//System.out.println(modFile + " " + maf);
		try {
			String mafFile = System.getProperty("user.dir") + System.getProperty("file.separator") + "Tmp" + new Random().nextLong() + ".maf";
			PrintStream mafStream = new PrintStream(mafFile);
			mafStream.println(maf);
			mafStream.close();
			String[] cmd = {
					"/bin/sh",
					"-c",
					"phyloP " + modFile + " " + mafFile
					};
			ProcessBuilder pr = new ProcessBuilder(cmd);
			Process p = pr.start();
			try {
				p.waitFor();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			StringBuilder builder = new StringBuilder();
			String line = null;
			while ( (line = br.readLine()) != null) {
			   builder.append(line);
			   builder.append(System.getProperty("line.separator"));
			}
			String[] re = builder.toString().split("\n");
			
			//p-value of conservation: 3.957918e-06
			//p-value of acceleration: 9.999986e-01
			//System.out.println("phyloP " + modFile + " " + mafFile);
			for(String r : re){
			//	System.out.println(r);
				if(r.startsWith("p-value of conservation")){
					pvalConservation = Double.parseDouble(r.split(":")[1]);
				}else if(r.startsWith("p-value of acceleration")){
					pvalAcceleration = Double.parseDouble(r.split(":")[1]);
				}
			}
			
			new File(mafFile).delete();
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/*public static void main(String[] args){
		String mod = "/media/kyowon/Data1/RPF_Project/genomes/maf/phylo.mod";
		MafParser mp = new MafParser("/media/kyowon/Data1/RPF_Project/genomes/maf/");
		PhyloPLauncher.modFile = mod;
		mp.readIndexFile();
		
		String mafString = "##maf version=1 scoring=autoMZ.v1\na score=1\ns mm9.chr5 3692901 -249 + 1 "
				+ "CTGGAAGAGAATATTCTGGCCAGGGAAGAGTTGTCTCAAGCTGGTGACAGT\ns "
				+ "oryCun1.chr5 3692901 -249 - 1 TCGGGAGGAAGCCTTCCGGCTGAGGAAGCGGCGCCTCCAGCCGGTGACAGT\ns "
				+ "rn4.chr5 3692901 -249 - 1 CTGGAAGAGAACATTCTGGCTGGGGAAGAGTTGTCTCAAGCCGACGACAGT\ns "
				+ "otoGar1.chr5 3692901 -249 - 1 TTGGAAGACAATATTCTGGCTGGGGAAGCAGTGTCTCAAGCTGGTGACAGT\ns "
				+ "cavPor2.chr5 3692901 -249 - 1 TTGGAAGAAAATACTCTGGCTGGAGAGGTGTCATTGCAAGGTAGTGACAGT\ns "
			//	+ "tupBel1.chr5 3692901 -249 - 1 TTGGAAGAAAATACTGTGGTCGGGGAAGCAGCATCTCAAGTTGGTGACAGT\ns "
				+ "tupBel1.chr5 3692901 -249 - 1 TGGGTTTTTTATACTTTTTTCTTTTAAGCAGCATCTCAAGAAAATGACAGT\ns "
				+ "calJac1.chr5 3692901 -249 - 1 TTGGAAGAAAATATTCTAACTGAGGAAGCAGCATCTCAAGCTGGTGACAGT";
		//System.out.println(mp.getSeqsInMafFormat("chrX", 150158242, true, 69*3));
		
		System.out.println(mafString + "\n" + new PhyloPLauncher(mafString).getPvalConservation());
	}*/
}
