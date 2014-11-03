package launcher;

import java.io.IOException;

public class BlatLauncher {
	static public void run(String databaseFasta, String queryFasta, String outputPsl, double percentIdentity){
		try {
			String[] cmd = {
					"/bin/sh",
					"-c",
					"blat " + databaseFasta + " " + queryFasta + " " + outputPsl + " -minIdentity=" + Double.toString(percentIdentity),
					};
			ProcessBuilder pr = new ProcessBuilder(cmd);		 
			Process p = pr.start();
			try {
				p.waitFor();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
}
