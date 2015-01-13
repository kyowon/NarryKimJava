package launcher;

import java.io.IOException;

public class ShellLauncher {
	static public void run(String c){
		try {
			//System.out.println(c);
			String[] cmd = {
					"/bin/sh",
					"-c",
					c,
					};
			ProcessBuilder pr = new ProcessBuilder(cmd);		 
			Process p = pr.start();
			try {
				p.waitFor();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
}
