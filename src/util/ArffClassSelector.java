package util;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import parser.BufferedLineReader;

public class ArffClassSelector {
	
	
	public static void main(String[] args) {
		String arff = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out.arff";
		String arffOut = arff + ".arff";
		ArrayList<String> classNames = new ArrayList<String>();
		classNames.add("MP");
		classNames.add("mP");
		classNames.add("UP");
		
		
		try {
			BufferedLineReader in = new BufferedLineReader(arff);
			PrintStream out = new PrintStream(arffOut);
			String s;
			boolean filter = false;
			//@attribute Class
			while((s=in.readLine())!=null){
				if(filter){
					//System.out.println("hh" + s);
					for(String c : classNames){
						if(s.endsWith(c)){
							out.println(s);
							break;
						}
					}
				}else if(s.startsWith("@attribute Class")){
					out.print("@attribute Class {");
					for(String c : classNames){
						out.print(c + ",");
					}
					out.println("}");
				}else if(s.startsWith("@data")){
					filter = true;
					out.println(s);
				}else{
					out.println(s);
				}
			}
			
			in.close();
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}
