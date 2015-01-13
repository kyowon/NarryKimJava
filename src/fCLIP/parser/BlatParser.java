package fCLIP.parser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;

public class BlatParser {
	public class BlatResult{// TODO fill later
		private String line;
		private String qname;
		public BlatResult(String s){
			line = s;
			String[] token = s.split("\t");
			qname = token[9];
		}
		public String getQname() { return qname; }
	}
	private HashMap<String, ArrayList<BlatResult>> resultMap;
	private String file;
	
	private void run(){
		try {
			BufferedLineReader in = new BufferedLineReader(file);
			String s;
			boolean start = false;
			resultMap = new HashMap<String, ArrayList<BlatResult>>();
			while((s=in.readLine())!=null){
				if(s.startsWith("-----------")){
					start = true;
					continue;
				}
				if(!start) continue;
				BlatResult r = new BlatResult(s);
				String qname = r.getQname();
				if(!resultMap.containsKey(qname)) resultMap.put(qname, new ArrayList<BlatParser.BlatResult>());
				resultMap.get(qname).add(r);
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public BlatParser(String file){
		this.file=  file;
		run();
	}
	public ArrayList<BlatResult> getResults(String qname){
		return resultMap.get(qname);
	}
}
