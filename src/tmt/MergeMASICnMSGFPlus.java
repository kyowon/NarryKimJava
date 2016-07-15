package tmt;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import parser.BufferedLineReader;


public class MergeMASICnMSGFPlus {

	public static void main(String[] args) {
		String masic = args[0];
		String msgf = args[1];
		String out = args[2];
		
		BufferedLineReader inMasic;
		BufferedLineReader inMsgf;
		try {
			inMasic = new BufferedLineReader((masic));
			String s;
			HashMap<Integer, String> masicStringMap = new HashMap<Integer, String>();
			while((s=inMasic.readLine())!=null){
				int sn = -1;
				String[] token = s.split("\t");				
				if(!s.startsWith("Dataset")) sn = Integer.parseInt(token[1]);
				
				StringBuffer sb = new StringBuffer();
				for(int i=4;i<token.length;i++){
					sb.append(token[i]);
					if(i<token.length-1)sb.append('\t');
				}
				masicStringMap.put(sn, sb.toString());
			}
			
			inMasic.close();
			
			inMsgf = new BufferedLineReader((msgf));
			PrintStream op = new PrintStream(out);
			while((s=inMsgf.readLine())!=null){
				int sn = -1;
				String[] token = s.split("\t");				
				if(!s.startsWith("#SpecFile")) sn = Integer.parseInt(token[2]);
				op.println(s + "\t" + masicStringMap.get(sn));
			}
			op.close();
			inMsgf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
