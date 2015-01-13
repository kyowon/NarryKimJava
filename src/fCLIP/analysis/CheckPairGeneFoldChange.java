package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import uninovo.parser.BufferedLineReader;

public class CheckPairGeneFoldChange {

	public static void main(String[] args) {
		String key = "Drosha";
		key = "Dicer";
		key = "DKD";
		String csv = "/media/kyowon/Data1/Dropbox/h19x2.sorted.out." + key + ".U.pair.csv";
		String txt = "/media/kyowon/Data1/Dropbox/back.si" + key + ".txt";
		String mfile = "/media/kyowon/Data1/Dropbox/gene_pair" + key + ".m";;
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
					out.println(s + "\tGcontrol1\tGtarget1\tGIntronContol1\tGIntronTarget1\tGcontrol2\tGtarget2\tGIntronContol2\tGIntronTarget2" );
				}else if(s.startsWith("chr")){
					out.print(s);
					String mKeyPrefix = key + "PairGene";
					
					String[] token = s.split("\t");
				//	int index = -1;
					
					mKeyPrefix += token[32]; // class TODO
					
					String[] accs1 = token[8].split(","); 
					String[] accs2 = token[14].split(",");
					
					//String[] gns = token[18].split(",");
					String[] reg1 = token[24].split(",");
					String[] reg2 = token[26].split(",");
					
					String c1c1 = ""; // orfs
					String c2c1 = ""; // introns
				//	String c3c1 = ""; // semi-introns
					String c1t1 = ""; // orfs
					String c2t1 = ""; // introns
				//	String c3t1 = ""; // semi-introns
					
					String c1c2 = ""; // orfs
					String c2c2 = ""; // introns
				//	String c3c2 = ""; // semi-introns
					String c1t2 = ""; // orfs
					String c2t2 = ""; // introns
				//	String c3t2 = ""; // semi-introns
					for(int i=0; i<accs1.length;i++){
						if(reg1.length>i && !reg1[i].contains("Intron")){// && reg2.length>j && !reg2[j].contains("Intron")){ 
							double[] t1 = cntMap.get(accs1[i]);
							if(t1 != null){
								c1c1 += t1[0] + ";";
								c1t1 += t1[1] + ";";
								
								String mKey = mKeyPrefix + "ORFs";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t1[0] + " " + t1[1] + ";");
							}	
						}else{ // intron
							double[] t1 = cntMap.get(accs1[i]);
							if(t1 != null){
								c2c1 += t1[0] + ";";
								c2t1 += t1[1] + ";";
								
								String mKey = mKeyPrefix + "Introns";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t1[0] + " " + t1[1] + ";");
							}	
						}
					}
					
					for(int i=0; i<accs2.length;i++){
						if(reg2.length>i && !reg2[i].contains("Intron")){// && reg2.length>j && !reg2[j].contains("Intron")){ 
							double[] t1 = cntMap.get(accs2[i]);
							if(t1 != null){
								c1c2 += t1[0] + ";";
								c1t2 += t1[1] + ";";
								
								String mKey = mKeyPrefix + "ORFs";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t1[0] + " " + t1[1] + ";");
							}	
						}else{ // intron
							double[] t1 = cntMap.get(accs2[i]);
							if(t1 != null){
								c2c2 += t1[0] + ";";
								c2t2 += t1[1] + ";";
								
								String mKey = mKeyPrefix + "Introns";
								if(!mStringMap.containsKey(mKey)) mStringMap.put(mKey, new StringBuilder());
								mStringMap.get(mKey).append(t1[0] + " " + t1[1] + ";");
							}	
						}
					}
					
					out.print("\t"+c1c1 + "\t" + c1t1 + "\t"+c2c1 + "\t" + c2t1 + "\t"+c1c2 + "\t" + c1t2 + "\t"+c2c2 + "\t" + c2t2);
				
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
