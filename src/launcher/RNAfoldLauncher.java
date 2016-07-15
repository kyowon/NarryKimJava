package launcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class RNAfoldLauncher {
	
	//private String seq;
	private double energy = 0;
	private double energyPerNT = 0;
	
	private int flankingNTnumber = 15;
	private int depth = 0;
	private int numberOfHairPins = 0;
	private int overHang = -1000;
	private int leftPaired = 0;
	private int rightPaired = 0;
	private double seqEntropy = 0;
	private double structureEntropy = 0;
	public final static int[] overhanglimit = {0,5};
	public final static int[] pairedNumberlimit = {8,6};
	private String structure;
	private boolean useCentroid = false;
	
	
	public double getSeqEntropy() {
		return seqEntropy;
	}
	
	public double getStructureEntropy(){
		return structureEntropy;
	}
	
	public double getEnergy() {
		return energy;
	}

	public double getEnergyPerNT(){
		return energyPerNT;
	}

	public int getDepth() {
		return depth;
	}

	public int getLeftPaired(){
		return leftPaired;
	}
	
	public int getRightPaired(){
		return rightPaired;
	}
	
	public int getOverHang() {
		return overHang;
	}
	
	public boolean isFeasibleFold(){
		if(overHang < overhanglimit[0] || overHang > overhanglimit[1]) return false;
	//	if(leftPaired < pairedNumberlimit[0] || rightPaired < pairedNumberlimit[1]) return false;
		return true;
	}
			

	public int getNumberOfHairPins(){
		return numberOfHairPins;
	}
	public RNAfoldLauncher(String seq, int flankingNTnumber){
		this(seq, flankingNTnumber, false);
	}
	
	public RNAfoldLauncher(String seq, int flankingNTnumber, boolean useCentroid){
		this.flankingNTnumber = flankingNTnumber;
		this.useCentroid = useCentroid;
		if(seq != null && !seq.isEmpty())
			run(seq);
	}
	
	public String getStructureString(){
		return structure;
	}
	
	private void run(String seq){//seq = pre + flanking.
		try {
			String pseq = seq;//.replace(' ', '\0');
			
			String[] cmd = {
					"/bin/sh",
					"-c",
					"echo "+ pseq + " | RNAfold --noPS -p"
					};
			ProcessBuilder pr = new ProcessBuilder(cmd);		 
			Process p = pr.start();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			StringBuilder builder = new StringBuilder();
			String line = null;
			while ( (line = br.readLine()) != null) {
			//	System.out.println(line);
			   builder.append(line);
			   builder.append(System.getProperty("line.separator"));
			}
			String[] re = builder.toString().split("\n");
		//	String[] result = re[1].split(" ");
			
		//	for(String kk : re) 
		//		System.out.println(builder.toString());
			
			int maxDepth = 0;
		//	int maxCenterDepth = 0;
			String st = re[1].substring(0, re[1].indexOf(' ')); 
			if(useCentroid) st = re[3].substring(0, re[3].indexOf(' '));
			
			int d = 0;
			int l = 0;
			int r = 0;
			int nh = 0;
			boolean left = true;
			for(int i=0;i<st.length();i++){
				char c = st.charAt(i);
				if(c == '('){
					d ++;
					left = true;
				}
				else if(c == ')'){
					d --;
					if(left) nh++;
					left = false;
				}//else if(c == '.'){
				//	if(i >= flankingNTnumber) r++;
				//	if(i< st.length() - flankingNTnumber) l++;
				//}
				maxDepth = maxDepth > d ? maxDepth : d;
			}
			
			for(int i=0;i<flankingNTnumber;i++){
			//	System.out.println("l");
				if(st.charAt(i) == '(')//break;
					l++;
			}
			
			for(int i=st.length()-1;i>st.length() - 1 - flankingNTnumber;i--){
			//	System.out.println("r");
				if(st.charAt(i) == ')') //break;
					r++;
			}
			structure = st;
			//System.out.println(st);
		//	if(nh == 1){	
			
			
			/*
			
				this.overHang = 0;
				
				String fst = st.substring(0, st.length()/2);
				//System.out.println(fst);
				String tst = st.substring(st.length()/2+1);
				
				//fst tst
				int fivePaired = 0;
				for(int i=fst.length()-1; i>=flankingNTnumber; i--){
					if(fst.charAt(i) == '(') fivePaired ++;
					else if(fst.charAt(i) == ')') fivePaired--;
				}
				
				int threePaired = 0;
				for(int i=0; i<tst.length() - flankingNTnumber; i++){
					if(tst.charAt(i) == ')') threePaired ++;
					else if(tst.charAt(i) == '(') threePaired--;
				}
				
				//int diff = fivePaired - threePaired;
				
			//	int fp = 0;
				for(int i=flankingNTnumber;i<fst.length();i++){
					if(fst.charAt(i)=='.') this.overHang--;
					else{
						if(threePaired >= fivePaired) break;
						if(fst.charAt(i) == '('){
							fivePaired --;
						}else if(fst.charAt(i) == ')'){
							fivePaired++;
						}					
						this.overHang--;
						
					}
				}
				
			//	int tp = 0;
				for(int i=tst.length() - flankingNTnumber-1;i>=0;i--){
					if(tst.charAt(i)=='.') this.overHang++;
					else{
						if(threePaired <= fivePaired) break;	
						if(tst.charAt(i) == ')'){
							threePaired --;
						}else if(tst.charAt(i) == '('){
							threePaired++;
						}
						this.overHang++;
										
					}
				}
				*/
			this.overHang = 0;
			
			String fst = st.substring(0, st.length()/2);
			//fst = fst.replace("&", "");
			//System.out.println(st);
			String tst = st.substring(st.length()/2);
			//System.out.println(fst);
			//System.out.println(tst);
			//fst tst
			int fivePaired = 0;
			for(int i=fst.length()-1; i>=flankingNTnumber; i--){
				if(fst.charAt(i) == '(') fivePaired ++;
				else if(fst.charAt(i) == ')') fivePaired--;
			}
			//.((((..(((((. .(((((.(((..(..((((((((.(((.((.(.
			//...).))))).)))))))).)))).))))).  .))))))))).....
			int threePaired = 0;
			for(int i=0; i<tst.length() - flankingNTnumber; i++){
				if(tst.charAt(i) == ')') threePaired ++;
				else if(tst.charAt(i) == '(') threePaired--;
			}
			
			//int diff = fivePaired - threePaired;
		//	int fp = 0;
			//System.out.println(threePaired + " " + fivePaired);
			for(int i=flankingNTnumber;i<fst.length();i++){
				if(fst.charAt(i)=='.') this.overHang--;
				else{
					if(threePaired >= fivePaired) break;
					if(fst.charAt(i) == '('){
						fivePaired --;
					}else if(fst.charAt(i) == ')'){
						fivePaired++;
					}					
					this.overHang--;					
				}
			}
		//	System.out.println(this.overHang);
		//	int tp = 0;
			for(int i=tst.length() - flankingNTnumber-1;i>=0;i--){
				if(tst.charAt(i)=='.') this.overHang++;
				else{
					if(threePaired <= fivePaired) break;	
					if(tst.charAt(i) == ')'){
						threePaired --;
					}else if(tst.charAt(i) == '('){
						threePaired++;
					}
					this.overHang++;									
				}
			}
			//System.out.println(this.overHang);
				/*int[] threePPairedIndices = new int[500];
				int[] fivePPairedIndices = new int[500];
				
				int prev3p = 0;
				for(int i=0; i<flankingNTnumber;i++){
					if(st.charAt(i) == '('){
						threePPairedIndices[prev3p++] = i;
					}
				}
				
				int prev5p = 0;
				for(int i=st.length() - 1; i>= st.length() - flankingNTnumber ;i--){
					if(st.charAt(i) == ')'){
						fivePPairedIndices[prev5p++] = i;
					}
				}
				//System.out.println(seq);
				int prev = Math.min(prev3p, prev5p);
				for(int i=prev > 0 ? threePPairedIndices[prev-1]+1 : 0; i<flankingNTnumber;i++){
					this.overHang ++;
				}
				for(int i=prev > 0 ? fivePPairedIndices[prev-1]-1 : st.length() - 1; i>= st.length() - flankingNTnumber ;i--){
					this.overHang --;
				}
			
				*/
				//for(int i=first-1;  i>= st.length() - flankingNTnumber ;i--){
				//	this.overHang--;
			//	}
				/*int td = 0;
				boolean startFolding = false;
				boolean endFolding = false;
				for(int i=flankingNTnumber;i<st.length() - flankingNTnumber;i++){
					//System.out.print(st.charAt(i));
					if(endFolding) this.overHang++;
					if(st.charAt(i) == '('){
						startFolding = true;
						td++;
					}
					if(st.charAt(i) == ')'){
						td--;
					}
					if(startFolding && td == 0){
						endFolding = true;
					}
					if(!startFolding) this.overHang--;
				}
				
				if(!endFolding){
					startFolding = false;
					this.overHang = 0;
					td = 0;
					for(int i=st.length() - flankingNTnumber-1;i>=flankingNTnumber;i--){
						if(endFolding) this.overHang--;
						if(st.charAt(i) == ')'){
							startFolding = true;
							td++;
						}
						if(st.charAt(i) == '('){
							td--;
						}
						if(startFolding && td == 0){
							endFolding = true;
						}
						if(!startFolding) this.overHang++;
						//System.out.println(this.overHang);
					}
				}*/
				
				
				
		//	}
				// this.overHang -= td;
			/*for(int i=st.length()/2; i<st.length();i++){
				//System.out.println(i + " " + result2.charAt(i));
				if(st.charAt(i) == '(') break;
				if(st.charAt(st.length() - i - 1) == ')') break;
				if(st.charAt(i) == ')' || st.charAt(st.length() - i - 1) == '(')		
					maxCenterDepth++;
			}*/
			//......................
			this.depth = maxDepth;
			this.numberOfHairPins = nh;
		//	this.centerDepth = maxCenterDepth;
			this.energy = Double.parseDouble(re[1].substring(re[1].indexOf(' ')+2, re[1].length()-1).trim());
			this.energyPerNT = energy/pseq.length();
			this.seqEntropy = getEntropy(pseq.substring(flankingNTnumber, pseq.length() - flankingNTnumber), 2);
			this.structureEntropy = getEntropy(st.substring(flankingNTnumber, st.length() - flankingNTnumber), 4);
			
		//	System.out.println(seq.substring(flankingNTnumber, seq.length() - flankingNTnumber));
			//this.overHang = l-r;
			leftPaired = l;
			rightPaired = r;
		//	int i = result.lastIndexOf(')');
		//	if(ei >= 0 && i >= 0)
		//		ret.add(Double.parseDouble(result.substring(ei + 1 , i)));
		//	else ret.add(0.0);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	//GUGAUCUCGACAGGAAACGUGCCACAGAGCAGUAGUGCGCAGGCAAGACUUUUCAGUGACGCCUUGUGGAACGCAGUUCAUGAUGUCCUAGCAGCUCUCACUAAGGGAACUGUACAUUCU
	//(((...((.(((((....(((.(((.(((.(((..(((....)))..))).))).))).)))))))).)).))).......(((((....((((.((((.....)))).)))))))))..
	 
	
	
	public String getBalancedStructure(char cleavageChar){
		String ustr = getUnWindedStructureStr(structure);
		StringBuilder sbl = new StringBuilder();
		StringBuilder sbr = new StringBuilder();
		
		int lindex = 0;
		int rindex = ustr.length();
		
		while(ustr.charAt(lindex) != '('){
			sbl.append(ustr.charAt(lindex++));
		}
		
		while(ustr.charAt(rindex) != ')'){
			sbr.insert(0, ustr.charAt(rindex--));
		}
		
		boolean[] cleavageSpecified = new boolean[2];
		while(lindex < rindex){
			char l = ustr.charAt(lindex);
			char r = ustr.charAt(rindex);
			if(l == r){
				sbl.append(l);
				sbr.insert(0,r);
				lindex++;
				rindex--;
			}else if(l == '('){
				rindex--;
			}else if(r == ')'){
				lindex++;
			}
			
			if(!cleavageSpecified[0] && lindex == flankingNTnumber){
				sbl.append(cleavageChar);
				cleavageSpecified[0] = true;
			}
			if(!cleavageSpecified[1] && rindex == ustr.length() - flankingNTnumber){
				sbr.insert(0,r);
				cleavageSpecified[1] = true;
			}
			
		}
		if(lindex == rindex) sbl.append(ustr.charAt(lindex));
		
		return sbl.toString();
	}
	
	
	
	

	public static String getUnWindedStructureStr(String str){
		//"(((...((.(((((....(((.(((.(((.(((..(((....)))..))).))).))).)))))))).)).))).......(((((....((((.((((.....)))).))))))))).."
		StringBuffer unWindedStr = new StringBuffer();
		boolean startScan = false;
		int c = 0;
		for(int i=str.length()/2-1;i>=0;i--){
			char toput = str.charAt(i);
			if(!startScan && toput == '('){
				startScan = true;
			}
			if(startScan){
				if(toput == ')') // local hairpin
				{
					toput = '.';
					c++;
				}
				if(c>0 && toput == '('){
					toput = '.';
					c--;
				}
			}
			unWindedStr.insert(0, toput);
		}
		
		c = 0;
		startScan = false;
		
		for(int i=str.length()/2;i<str.length();i++){
			char toput = str.charAt(i);
			if(!startScan && toput == ')'){
				startScan = true;
			}
			if(startScan){
				if(toput == '(') // local hairpin
				{
					toput = '.';
					c++;
				}
				if(c>0 && toput == ')'){
					toput = '.';
					c--;
				}
			}
			unWindedStr.append(toput);
		}
		
		
		return unWindedStr.toString();
	}
	

	
	public int[] getBulgeMatrix5p(int l, int r, boolean unwind){ 
		String str = structure;
		if(unwind) str = getUnWindedStructureStr(str);
		//if(!structure.equals(str)) System.out.println(structure+"\n"+str);
		int[] unpaired5p = new int[l+r+1];
		
		for(int i=0;i<unpaired5p.length-1;i++){
			unpaired5p[i] = 1;
		}
		
		int index = l-1; // index for paired5p
		boolean unpairedExist = false;
		boolean countStart = false;
		for(int i=0; i<l;i++){//5p flnaking region
			if(index < 0) break;			
			if(str.charAt(l-i-1) == '('){
				if(countStart){
					unpaired5p[index] = unpairedExist? 1 : 0;
					
				}
				index--;
				unpairedExist = false;
				countStart = true;
			}else{
				unpairedExist = true;
			}			
		}
		
		index = l-1; // right side
		unpairedExist = false;
		
		for(int i=0;i<str.length()/2-l-1;i++){
		//	System.out.println((l+i) + " " + str.charAt(l+i));	
			if(index >= unpaired5p.length-1) break;
			if(str.charAt(l+i) == '('){
				unpaired5p[index++] = unpairedExist? 1 : 0;
				unpairedExist = false;
			
			}else{
				unpairedExist = true;
			}	
		}
		
		return unpaired5p;
		
	}
	
	public int[] getBulgeMatrix3p(int l, int r, boolean unwind){ 
		int[] unpaired3p = new int[l+r+1];
		String str = structure;
		if(unwind) str = getUnWindedStructureStr(str);
		//System.out.println(l+ " " + r);
		
		for(int i=0;i<unpaired3p.length-1;i++){
			unpaired3p[i] = 1;
		}
		
		int paircntr = 0;
		for(int i=0; i<l ; i++){
			if(str.charAt(i) == '(') paircntr ++;
			else if(str.charAt(i) == ')') paircntr --;
		}
		int pairoffset = 0;
		int indexInitValue = 0;
		for(; paircntr>0 && indexInitValue<str.length(); indexInitValue++){
			char c = str.charAt(str.length()-indexInitValue-1);
			if(c == '('){
				paircntr++;
				if(indexInitValue >= l) pairoffset++;
			}
			else if(c == ')'){
				paircntr--;
				if(indexInitValue >= l) pairoffset--;
			}
		}
			
		int index = l - 1; // index for paired5p
		//System.out.println(pairoffset);
		
		boolean unpairedExist = false;
		boolean countStart = false;
		for(int i=0; i<indexInitValue;i++){//5p flnaking region
			if(index < 0) break;
			if(str.charAt(str.length() - indexInitValue + i) == ')'){
				if(countStart){
					unpaired3p[index] = unpairedExist? 1 : 0;
					
				}
				index--;
				unpairedExist = false;
				countStart = true;
			}else{
				unpairedExist = true;
			}			
		}
		
		index = l - 1; // right side
		unpairedExist = false;
		
		for(int i=0;i<str.length()/2-indexInitValue-1;i++){
		//	System.out.println((l+i) + " " + str.charAt(l+i));		
			if(index >= unpaired3p.length-1) break;
			if(str.charAt(str.length() - indexInitValue - 1 - i) == ')'){
				unpaired3p[index++] = unpairedExist? 1 : 0;
				unpairedExist = false;
			}else{
				unpairedExist = true;
			}	
		}	
		unpaired3p[unpaired3p.length-1] = pairoffset;
		return unpaired3p;
		
	}
	
	private double getEntropy(String seq, int stepSize){
		HashMap<String, Double> freqMap = new HashMap<String, Double>();
		for(int i=0;i<seq.length() - stepSize + 1;i+=stepSize){
			String s = seq.substring(i, i+stepSize);
			if(!freqMap.containsKey(s)) freqMap.put(s, .0);
			freqMap.put(s, freqMap.get(s) + 1);
		}
		double e = 0;
		for(String s : freqMap.keySet()){
			double freq = freqMap.get(s) / (seq.length() / stepSize);
			e -= freq * Math.log(freq) / Math.log(2);
		}
		return e;
	}
	
	
	static public void main(String[] args){
		String seq = "TGTGTGCACATGTGC CCAGGGCCCGGGACAGCGCCACGGAAGAGGACGCACCCGGCTGTGTGC ACATGTGCCCAGGGC"; // f
		seq = "AAAGGACCCTTCCAG AGGGCCCCCCCTCAATCCTGTTGTGCCTAATTCAGAGGGTTGGGTGGAGGCTCTCC TGAAGGGCTCTGAAG";
		seq = "TGAGAGCTGCAGGGT GGGTGGGCCCCTGGCTGTGCTGGGCCCTTCTCGATCATCTGAGGACCTGGCCGGCCCCCTCCC TTCCTCAGTCTCTTC";
		seq = seq.replaceAll(" ", "");
		String seqw = "";
		for(char s : seq.toCharArray()){
			seqw = s + seqw;
		}//ACGGAGGGGTAGAGC AGCCCTCGGCGGCCCGGGGGGCGGGCGGCGGTGCCCGTCCCGGGGCTGCGCGAGGCACAGGC GCCGCGCCCAGCGGA

//		System.out.println(seq);  
//		System.out.println(seqw);
		RNAfoldLauncher test = new RNAfoldLauncher(seq , 15);
		
		
		//ACCTATCCATTTATCATCCATCTACGAGGTGTCTGGGATGTAATGGATGG
		//ACCUAUCCAUUUAUCAUCCAUCUACAGGUGGACCAAAAUGUAGGAGAUCU
		//TGCTCAGGCATAACCCCTCG GAACGGCTGCCCCTGGCCCAGGTCTCAGCCCACCCTTGGGTCCGGGCCAACTCTCGG AGGGTGCTGCCTCCCTCTGC
		System.out.println(test.energyPerNT);
		System.out.println(test.numberOfHairPins);
		System.out.println(test.leftPaired);
		System.out.println(test.rightPaired);
		System.out.println(test.seqEntropy);
		System.out.println(test.structureEntropy);
		System.out.println(test.numberOfHairPins);
		//System.out.println(test.centerDepth);
		System.out.println(test.overHang);
		System.out.println(seq);
		System.out.println(test.getStructureString());
		for(int b : test.getBulgeMatrix5p(25, 30, true))
			System.out.print(b + "");
		System.out.println();
		for(int b : test.getBulgeMatrix3p(25, 30, true))
			System.out.print(b + "");
		System.out.println();
	//	1 1 1 1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 0 
	//	17 15
	//	1 1 1 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 0 
		
		
	//	(.((.((((((((((.((.((((.(((...(((((........)).))).))))))).)).)))))))))).)).)..
	//	1 1 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0 1 0 0 0 0 1 1 1 1 1 1 1 0 
	//	17 15
	//	1 1 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 1 0 1 1 1 1 1 1 1 0 
	//	1 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 
	//	1 1 1 1 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 -1 
	}
	
}
