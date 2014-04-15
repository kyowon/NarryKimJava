package rpf.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class tmp2 {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		BufferedLineReader in = new BufferedLineReader("/home/kyowon/Desktop/target.csv");
		BufferedLineReader in2 = new BufferedLineReader("/media/kyowon/Data1/RPF_Project/samples/sample4/results/out1_0.3.csv");
		PrintStream out = new PrintStream("/home/kyowon/Desktop/DrYou.txt");
		String s;
		HashSet<String> targetAcc = new HashSet<String>();
		HashSet<String> diffAcc = new HashSet<String>();
		while((s=in.readLine())!=null){
			targetAcc.add(s);
		}
		//15
		while((s=in2.readLine())!=null){
			String acc = s.split("\t")[5];
			if(s.split("\t")[10].equals("true")) continue;
			if(targetAcc.contains(acc)){
				System.out.println(s);
				out.println(s);
				diffAcc.add(acc);
			}
		}
		System.out.println(diffAcc.size());
		out.close();
		in2.close();
		in.close();
	}

}
