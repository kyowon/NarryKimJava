package rpf.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class tmp2 {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		BufferedLineReader in = new BufferedLineReader("/home/kyowon/Desktop/pair");
		PrintStream out = new PrintStream("/home/kyowon/Desktop/pair.fa");
		String s;
		HashSet<String> seqs = new HashSet<String>();
		int i = 0;
		while((s=in.readLine()) != null){
			if(s.isEmpty()) continue;
			int j = s.indexOf(' ');
			seqs.add(s.substring(0, j<0 ? 20 : j));
			
		}
		
		in.close();
		
		for(String seq :seqs){
			out.println(">seq" + i++);
			out.println(seq);
		}

		out.close();
	}

}
