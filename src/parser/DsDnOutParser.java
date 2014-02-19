package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import net.sf.samtools.util.BufferedLineReader;

public class DsDnOutParser {
	private static HashMap<String, HashMap<Boolean, HashMap<Integer, ArrayList<double[]>>>> map; // integer -> mapped to smap
	private static HashMap<String, HashMap<Boolean, ArrayList<Integer>>> smap;
	private static HashMap<String, HashMap<Boolean, ArrayList<Integer>>> emap;// smap << emap
	
	private static void read(String filename){
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(filename));
			String s;
			String contig = "";
			boolean isPlusStrand = false;;
			int sp = 0;
			int ep = 0;
			ArrayList<double[]> val = null;
			map = new HashMap<String, HashMap<Boolean, HashMap<Integer, ArrayList<double[]>>>>();
			smap = new HashMap<String, HashMap<Boolean, ArrayList<Integer>>>();
			emap = new HashMap<String, HashMap<Boolean, ArrayList<Integer>>>();
			
			while((s=in.readLine())!=null){
				String[] token = s.split("\t");
				if(s.startsWith("CONTIG")){
					contig = token[1];
					continue;
				}else if(s.startsWith("PLUS") || s.startsWith("MINUS")){
					if(val != null){
						if(!map.containsKey(contig)){
							map.put(contig, new HashMap<Boolean, HashMap<Integer, ArrayList<double[]>>>());
							smap.put(contig, new HashMap<Boolean, ArrayList<Integer>>());
							emap.put(contig, new HashMap<Boolean, ArrayList<Integer>>());
						}
						if(!map.get(contig).containsKey(isPlusStrand)){
							map.get(contig).put(isPlusStrand, new HashMap<Integer, ArrayList<double[]>>());
							smap.get(contig).put(isPlusStrand, new ArrayList<Integer>());
							emap.get(contig).put(isPlusStrand, new ArrayList<Integer>());
						}
						map.get(contig).get(isPlusStrand).put(sp, val);
						smap.get(contig).get(isPlusStrand).add(sp);
						emap.get(contig).get(isPlusStrand).add(ep);
					}					
					isPlusStrand = s.startsWith("PLUS");
					sp = Integer.parseInt(token[1]);
					ep = Integer.parseInt(token[2]);
					val = new ArrayList<double[]>();
					continue;
				}
				val.add(new double[]{Double.parseDouble(token[0]), Double.parseDouble(token[1])});
			}			
			in.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static double getDsDn(String contig, boolean isPlusStrand, ArrayList<Integer> positions){
		if(!smap.containsKey(contig)) return 1;
		ArrayList<Integer> sps = smap.get(contig).get(isPlusStrand);
		ArrayList<Integer> eps = emap.get(contig).get(isPlusStrand);
		HashMap<Integer, ArrayList<double[]>> ss = map.get(contig).get(isPlusStrand);
		int sp = -1;
		int ep = -1;
		ArrayList<double[]> sss = null; 
		double nNSSubs=0, nNSSites=0;
		for(int i=0;i<positions.size();i+=3){
			int cp = positions.get(i);
			if(!(cp>=sp && cp <= ep)){
				int index = Collections.binarySearch((isPlusStrand? sps : eps), cp);
				if(index < 0) index = - index - 2;
				sp = sps.get(index);
				ep = eps.get(index);
				sss = ss.get(sp);
			}				
			if(cp>=sp && cp <= ep && sss!=null){
				double[] nsss = sss.get(isPlusStrand? cp-sp : ep-cp);
				nNSSubs += nsss[0];
				nNSSites += nsss[1];
				//System.out.println(nsss[0] + " " + nsss[1]);
			}else{
				nNSSubs += 0.9404296875; 
				nNSSites += 2.28125;
				// poor score
			}				
		}//(nNSSubs * nSSites) / (nNSSites * nSSubs);
		return (nNSSubs * (positions.size() - nNSSites)) / (nNSSites * (positions.size()/3 - nNSSubs));
	}
	
	static public void main(String[] args){
		DsDnOutParser.read("/media/kyowon/Data1/RPF_Project/genomes/mm9/maf/dsdnOutTest.txt");
		ArrayList<Integer> positions = new ArrayList<Integer>();
		positions.add(3015970);
		positions.add(3015970+1);
		positions.add(3015970+2);
		positions.add(3015970+3);
		System.out.println(DsDnOutParser.getDsDn("chr18", true, positions));
	}
	
	
}
