package fCLIP.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import uninovo.parser.BufferedLineReader;

public class CheckGeneFoldChange {

	public static void main(String[] args) {
		String key = "siD-K-D";
		String csv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out." + key + ".csv";
		String txt = "/media/kyowon/Data1/Dropbox/back." + key + ".txt";
		
		try {
			BufferedLineReader inTxt = new BufferedLineReader(txt);
			HashMap<String, double[]> cntMap = new HashMap<String, double[]>();
			String s;
			
			while((s=inTxt.readLine())!=null){
				String[] token = s.split(",");
				String acc = token[0];
				double[] cnt = new double[2];
				cnt[0] = Double.parseDouble(token[1]);
				cnt[1] = Double.parseDouble(token[2]);
				cntMap.put(acc, cnt);
			}
			inTxt.close();
			BufferedLineReader inCsv = new BufferedLineReader(csv);
			PrintStream out = new PrintStream(csv + ".gene.csv");
			while((s=inCsv.readLine())!=null){
				if(s.startsWith("Contig")){
					out.println(s + "\tGcontrol\tGtarget" );
				}else if(s.startsWith("chr")){
					out.print(s);
					if(s.equals("_")){
						out.print("\t%\t%");
					}else{
						String[] token = s.split("\t");
						int index = -1;
						String[] accs = token[17].split(",");
						String[] gns = token[18].split(",");
						String[] reg1 = token[19].split(",");
						String[] reg2 = token[20].split(",");
						
						for(int i=0; i<gns.length;i++){
							String gn = gns[i];
							if(!gn.startsWith("MIR") && ( (reg1.length>i &&  !reg1[i].contains("Intron")) && (reg2.length>i && !reg2[i].contains("Intron")))){
								index = i;
								break;
							}
						}
						if(index >=0){
							
							
							double[] t = cntMap.get(accs[index]);
							if(t != null)
								out.print("\t"+t[0]+"\t"+t[1]);
							else out.print("\t%\t%");
						}else{
							out.print("\t%\t%");
						}
					}
					out.println();
				}
			}
			out.close();
			inCsv.close();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}

}
