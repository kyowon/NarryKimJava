package util;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import parser.BufferedLineReader;

public class Arff2Mfile {

	public static void generate(String arff, String m) {
		try {
			BufferedLineReader in = new BufferedLineReader(arff);
			PrintStream out = new PrintStream(m);
			String s;
			out.print("att={");
			ArrayList<String> classes = null;
			while((s=in.readLine())!=null){
				
				if(s.startsWith("@attribute")){
					out.print("'"+ s.split(" ")[1]+"'; ");
				}else if(s.startsWith("@data")){
					out.println("};");
					out.println("value = [");
					classes = new ArrayList<String>();
				}else if(classes != null){
					String[] vals = s.split(",");
					if(vals[vals.length-1].equals("?")) continue;
					for(int i=0;i<vals.length-1;i++){
						String val = vals[i];
						if(val.equals("?")) val = "NaN";
						out.print(val + ",");
					}
					out.println();
					classes.add(vals[vals.length-1]);
				}
			}
			out.println("];");
			out.print("group={");
			for(String c : classes){
				out.println("'"+c+"'");
			}
			out.println("};");
			out.close();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}

}
