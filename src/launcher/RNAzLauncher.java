package launcher;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;

public class RNAzLauncher {
	
	private double sci = 0;
	private double zscore = 1;
	
	public RNAzLauncher(String maf){
		run(maf);
	}
	
	public double getSCI(){ return sci; }
	public double getZScore() { return zscore; }
	
	private void run(String maf){
		if(maf == null) return;
		try {
			String mafFile = System.getProperty("user.dir") + System.getProperty("file.separator") + "Tmp.maf";
			PrintStream mafStream = new PrintStream(mafFile);
			mafStream.println(maf);
			mafStream.close();
			String[] cmd = {
					"/bin/sh",
					"-c",
					"RNAz " + mafFile
					};
			ProcessBuilder pr = new ProcessBuilder(cmd);		 
			Process p = pr.start();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			StringBuilder builder = new StringBuilder();
			String line = null;
			while ( (line = br.readLine()) != null) {
			   builder.append(line);
			   builder.append(System.getProperty("line.separator"));
			}
			String[] re = builder.toString().split("\n");
			
			for(String r : re){
				//System.out.println(r);
				if(r.startsWith(" Structure conservation index")){
					sci = Double.parseDouble(r.split(":")[1]);
				}else if(r.startsWith(" Mean z-score")){
					zscore = Double.parseDouble(r.split(":")[1]);
				}
			}
			new File(mafFile).delete();
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
}
