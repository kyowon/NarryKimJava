package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import uninovo.parser.BufferedLineReader;

public class CheckGeneFoldChange {

	public static void main(String[] args) { // TODO
		String key = "DKD";
		String csv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out." + key + ".csv";
		String txt = "/media/kyowon/Data1/Dropbox/back.si" + key + ".txt";
		String mfile = "/media/kyowon/Data1/Dropbox/gene_" + key + ".m";;
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
			PrintStream outm = new PrintStream(mfile);
			HashMap<String, StringBuilder> mStringMap = new HashMap<String, StringBuilder>();
			
			while((s=inCsv.readLine())!=null){
				if(s.startsWith("Contig")){
					out.println(s + "\tGcontrol\tGtarget\tGIntronContol\tGIntronTarget\tGsemiIntronControl\tGsemiIntronTarget" );
				}else if(s.startsWith("chr")){
					out.print(s);
					String mKeyPrefix = key + "Gene";
					
					String[] token = s.split("\t");
				//	int index = -1;
					if(!token[11].equals("_")) mKeyPrefix += "MiRNA";
					else{
						mKeyPrefix += token[5] + token[4];
					}
					
					String[] accs = token[17].split(",");
					String[] gns = token[18].split(",");
					String[] reg1 = token[19].split(",");
					String[] reg2 = token[20].split(",");
					String c1c = ""; // orfs
					String c2c = ""; // introns
					String c3c = ""; // semi-introns
					String c1t = ""; // orfs
					String c2t = ""; // introns
					String c3t = ""; // semi-introns
					
					for(int i=0; i<gns.length;i++){
						String gn = gns[i];
						if(gn.startsWith("MIR")){
							// write nothing
						}else if(reg1.length>i && !reg1[i].contains("Intron") && reg2.length>i && !reg2[i].contains("Intron")){ 
							double[] t = cntMap.get(accs[i]);
							if(t != null){
								c1c += t[0] + ";";
								c1t += t[1] + ";";
								
								String mKey = mKeyPrefix + "ORFs";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t[0] + " " + t[1] + ";");
								
							}
							//else out.print("\t%\t%");
						}else if((reg1.length<=i || reg1[i].contains("Intron")) && (reg2.length<=i || reg2[i].contains("Intron"))){
							double[] t = cntMap.get(accs[i]);
							if(t != null){
								c2c += t[0] + ";";
								c2t += t[1] + ";";
								
								String mKey = mKeyPrefix + "Introns";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t[0] + " " + t[1] + ";");
							}
						}else{
							double[] t = cntMap.get(accs[i]);
							if(t != null){
								c3c += t[0] + ";";
								c3t += t[1] + ";";
								
								String mKey = mKeyPrefix + "SemiIntrons";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t[0] + " " + t[1] + ";");
							}
						}										
						
					}
					out.print("\t"+c1c + "\t" + c1t + "\t"+c2c + "\t" + c2t + "\t"+c3c + "\t" + c3t);
				
					out.println();
				}
			}
			out.close();
			
			for(String mkey : mStringMap.keySet()){
				outm.println(mkey + "= [");
				outm.println(mStringMap.get(mkey));
				outm.println("];");
				outm.println(mkey + "FC = " + mkey + "(find(" + mkey + "(:,1)>=30 & " + mkey + "(:,2)>=30),:);" + mkey + "FC = log("+mkey+"FC(:,2)./"+ mkey + "FC(:,1))/log(2);" );
			}
			
			
			outm.close();
			inCsv.close();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}

}
