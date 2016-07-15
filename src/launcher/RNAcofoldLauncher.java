package launcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

public class RNAcofoldLauncher {
	public final int upperStemLength = 22;
	public final int lowerStemLength = 11;
	private double energy = 0;
	private double energyPerNT = 0;
	
	private int flankingNTnumber = 25;
	private int depth = 0;
	private int numberOfHairPinsInStem = 0;
	private int overHang = -1000;
	private int leftPaired = 0;
	private int rightPaired = 0;
	private int loopLength = 0;
	private int leftUpstreamPaired = 0;
	private int rightUpstreamPaired = 0;
	private String structure;
	private String pstructure;
	
	private double seqEntropy = 0;
	private double structureEntropy = 0;
	private boolean has11Burge3pUp = false;
	private boolean has11Burge3pDown = false;
	private boolean has11Burge5pUp = false;
	private boolean has11Burge5pDown = false;
	private boolean has22Burge3pUp = false;
	private boolean has22Burge3pDown = false;
	private boolean has22Burge5pUp = false;
	private boolean has22Burge5pDown = false;
	
	static private int seqLength = 100000;
	
	
	public boolean has11Burge3pUp() { return has11Burge3pUp; }
	public boolean has11Burge3pDown() { return has11Burge3pDown; }
	public boolean has11Burge5pUp() { return has11Burge5pUp; }
	public boolean has11Burge5pDown() { return has11Burge5pDown; }

	public boolean has22Burge3pUp() { return has22Burge3pUp; }
	public boolean has22Burge3pDown() { return has22Burge3pDown; }
	public boolean has22Burge5pUp() { return has22Burge5pUp; }
	public boolean has22Burge5pDown() { return has22Burge5pDown; }

	public double getSeqEntropy() {
		return seqEntropy;
	}
	
	public String getStructureString(){
		return structure;
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
	
	public int getLeftUpstreamPaired(){
		return leftUpstreamPaired;
	}
	
	public int getRightUpstreamPaired(){
		return rightUpstreamPaired;
	}
	
	public int getOverHang() {
		return overHang;
	}

	public int getNumberOfHairPinsInStem(){
		return numberOfHairPinsInStem;
	}
	
	public int getLoopLength(){
		return loopLength;
	}
	
	public RNAcofoldLauncher(String seq, int flankingNTnumber){
		this.flankingNTnumber = flankingNTnumber;
		if(seq != null && !seq.isEmpty())
			run(seq);
		//System.out.println(seq + " " + this.numberOfHairPinsInStem + " " + flankingNTnumber);
	}
	
	public boolean isFeasibleFold(){
		if(overHang < RNAfoldLauncher.overhanglimit[0] || overHang > RNAfoldLauncher.overhanglimit[1]) return false;
		//if(leftPaired < RNAfoldLauncher.pairedNumberlimit[0] || rightPaired < RNAfoldLauncher.pairedNumberlimit[1]) return false;
		return true;
	}
	
	public static void setSeqLength(int s){ seqLength = s;}
	public static int getSeqLength() {return seqLength;}
	
	private void run(String seq){//seq = pre + flanking.
		try {
			StringBuilder pseqb = new StringBuilder();//
			String fivePString = seq.replace(' ', '\0').substring(0, Math.min(seq.length()/2, flankingNTnumber + seqLength/2));
			pseqb.append(fivePString);
			
			
			pseqb.append("\\&");
			String threePString = seq.substring(Math.max(seq.length() - flankingNTnumber - seqLength/2 , seq.length()/2)); 
			
			
			
			pseqb.append(threePString);
			
			//System.out.println(pseqb);
			
			String pseq = pseqb.toString();
		//	System.out.println(pseq);
			//.replace("&", "\\&");
			String[] cmd = {
				"/bin/sh",
				"-c",
				"echo "+ pseq + " | RNAcofold --noPS -p"
			};
			ProcessBuilder pr = new ProcessBuilder(cmd);		 
			Process p = pr.start();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			StringBuilder builder = new StringBuilder();
			
			String line = null;
			while ( (line = br.readLine()) != null) {
			   builder.append(line);
			   builder.append(System.getProperty("line.separator"));
			}
			String[] re = builder.toString().split("\n");
		
			int maxDepth = 0;
			String st = re[1].substring(0, re[1].indexOf(' ')).replace("&", "");
			structure = st;
			pstructure = re[2].substring(0, re[2].indexOf(' ')).replace("&", "");
			if(pstructure.contains("|")) pstructure = structure; // TODO
			
			String fst = st.substring(0, st.length()/2);
			//fst = fst.replace("&", "");
			//System.out.println(fst);
			String tst = st.substring(st.length()/2);
			//System.out.println(tst);
			for(int i=flankingNTnumber+9;i<Math.min(flankingNTnumber+12, fst.length());i++){
				if(fst.charAt(i) == '.'){
					has11Burge5pUp = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber+20;i<Math.min(flankingNTnumber+23, fst.length());i++){
				if(fst.charAt(i) == '.'){
					has22Burge5pUp = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber-11;i>=Math.max(flankingNTnumber-13, 0);i--){
				if(fst.charAt(i) == '.'){
					has11Burge5pDown = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber-22;i>=Math.max(flankingNTnumber-24, 0);i--){
				if(fst.charAt(i) == '.'){
					has22Burge5pDown = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber+9;i<Math.min(flankingNTnumber+12, tst.length());i++){
				if(tst.charAt(i) == '.'){
					has11Burge3pDown = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber+20;i<Math.min(flankingNTnumber+23, tst.length());i++){
				if(tst.charAt(i) == '.'){
					has22Burge3pDown = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber-11;i>=Math.max(flankingNTnumber-13, 0);i--){
				if(tst.charAt(i) == '.'){
					has11Burge3pUp = true;
					break;
				}
			}
			
			for(int i=flankingNTnumber-22;i>=Math.max(flankingNTnumber-24, 0);i--){
				if(tst.charAt(i) == '.'){
					has22Burge3pUp = true;
					break;
				}
			}
			
			int d = 0;
			int l = 0;
			int r = 0;
			int lu = 0;
			int ru = 0;
			
			int tnh = 0;
			
			int ln = 0;
			boolean left = false;
			for(int i=0;i<st.length();i++){
				char c = st.charAt(i);
				if(c == '('){
					d ++;
					left = true;
					ln = 0;
				}else if(c == ')'){
					d --;
					if(left){
						loopLength = loopLength > ln ? loopLength : ln;
					}
					ln = 0;
					left = false;
				}else if(c == '.'){
					ln++;
				}
				maxDepth = maxDepth > d ? maxDepth : d;
			}
			//System.out.println(tnh + " " + hnh);
			
			
			left = false;
			for(int i=flankingNTnumber-lowerStemLength;i<flankingNTnumber+upperStemLength;i++){
				char c = st.charAt(i);
				if(c == '('){
					left = true;
				}else if(c == ')'){
					if(left){
						tnh++;
					}
					left = false;
				}
			}
			left = false;
			for(int i=st.length() - flankingNTnumber - upperStemLength;i<st.length() - flankingNTnumber + lowerStemLength;i++){
				char c = st.charAt(i);
				if(c == '('){
					left = true;
				}else if(c == ')'){
					if(left){
						tnh++;
					}
					left = false;
				}
			}
			
			
			for(int i=flankingNTnumber-1;i>=0;i--){
			//	System.out.println("l");
				if(st.charAt(i) == '(')//break;
					l++;
				//if(st.charAt(i) == '.' && st.charAt(st.length() - 1 - i) == '.') break;
			}
			
			for(int i=flankingNTnumber;i<st.length()/2;i++){
				if(st.charAt(i) == '(')//break;
					lu++;
			}
			
			for(int i=st.length() - flankingNTnumber;i<st.length();i++){ //st.length()-1
				if(st.charAt(i) == ')') //break;
					r++;
				//if(st.charAt(i) == '.' && st.charAt(st.length() - 1 - i) == '.') break;
			}
			
			for(int i=st.length() - 1 - flankingNTnumber;i>st.length()/2;i--){
				if(st.charAt(i) == ')') //break;
					ru++;
			}
			
		//	if(nh == 1){			
			this.overHang = 0;
			
			
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
		//	System.out.println(fivePaired + " " + threePaired);
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
			
			/*
			int[] threePPairedIndices = new int[500];
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
			
			this.depth = maxDepth;
			this.numberOfHairPinsInStem = tnh;
		//	this.centerDepth = maxCenterDepth;
			this.energy = Double.parseDouble(re[1].substring(re[1].indexOf(' ')+2, re[1].length()-1).trim());
			
			pseq = pseq.replace("\\&", "");
			this.energyPerNT = energy/pseq.length();
		//	System.out.println(pseq.length());
			if( pseq.length() - flankingNTnumber <=0) System.out.println(pseq);
			
			this.seqEntropy = getEntropy(pseq.substring(flankingNTnumber, pseq.length() - flankingNTnumber), 2);
			this.structureEntropy = getEntropy(st.substring(flankingNTnumber, st.length() - flankingNTnumber), 4);
			
		//	System.out.println(seq.substring(flankingNTnumber, seq.length() - flankingNTnumber));
			//this.overHang = l-r;
			leftPaired = l;
			rightPaired = r;
			
			leftUpstreamPaired = lu;
			rightUpstreamPaired = ru;
		//	int i = result.lastIndexOf(')');
		//	if(ei >= 0 && i >= 0)
		//		ret.add(Double.parseDouble(result.substring(ei + 1 , i)));
		//	else ret.add(0.0);
		} catch (IOException e) {
			e.printStackTrace();
		}
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
	
	
	
	
	
	
	public static String getUnWindedStructureStr(String str){
		//"(((...((.(((((....(((.(((.(((.(((..(((....)))..))).))).))).)))))))).)).))).......(((((..,.(({(.((((.....)))).))))})))).."
		StringBuffer unWindedStr = new StringBuffer();
		boolean startScan = false;
		int c = 0;
		for(int i=str.length()/2-1;i>=0;i--){
			char toput = str.charAt(i);
			if(!startScan && (toput == '(')){
				startScan = true;
			}
			if(startScan){
				if(toput == ')') // local hairpin
				{
					toput = '.';
					c++;
				}else if(c>0){					
					if(toput == '(') c--;
					toput = '.';
				}
			}
			unWindedStr.insert(0, toput);
		}
		
		
		startScan = false;
		c = 0;
		for(int i=str.length()/2;i<str.length();i++){
			char toput = str.charAt(i);
			if(!startScan && (toput == ')')){
				startScan = true;
			}
			if(startScan){
				if(toput == '(') // local hairpin
				{
					toput = '.';
					c++;
				}else if(c>0){					
					if(toput == ')' ) c--;
					toput = '.';
				}
			}
			unWindedStr.append(toput);
		}
		
		
		return unWindedStr.toString();
	}
	
	
	public float[] getBulgeMatrix5p(int l, int r, boolean unwind, boolean correspondTo3p){ 
		String str = structure;
		if(unwind) str = getUnWindedStructureStr(str);
		float[] unpaired5p = new float[l+r+1];
		for(int i=0;i<unpaired5p.length-1;i++){
			unpaired5p[i] = 1;
		}
		
		int paircntr = 0;
		for(int i=str.length()-1; i>=str.length()-l ; i--){
			if(str.charAt(i) == ')') paircntr ++;
			else if(str.charAt(i) == '(') paircntr --;
		}
		//int indexInitValue = 0;
		int pairOffset = 0;
		
		for(int i=0; i<l || paircntr>0; i++){
			char c = str.charAt(i);
			if(c == ')'){
				paircntr++;				
			}
			else if(c == '('){
				paircntr--;
			}
			if(i == l-1) pairOffset = paircntr;
			//if(paircntr>=0) indexInitValue++;
		}
		
		ArrayList<Float> tmp = new ArrayList<Float>();
		float pairProb  = 0;
		int offset = 0;
		for(int i=0;i<r+l;i++){
			char pc = str.charAt(i);
			if(pc == '('){
				tmp.add(pairProb);	
				pairProb = 0;
			}
			if(pc != '('){
				if(pc == '.') pairProb = Math.max(pairProb, 1);
				else if(pc == ',') pairProb = Math.max(pairProb, 2f/3f);
				else pairProb = Math.max(pairProb, 1f/3f);
			}
			if(i == l-1) offset = tmp.size() + (correspondTo3p ? pairOffset : 0);
		}
		//System.out.println(offset);
		for(int i=0;i<Math.min(tmp.size(), r+offset);i++){
			unpaired5p[i+l-offset] =  tmp.get(i);
		}
		
		if(correspondTo3p) unpaired5p[unpaired5p.length-1] = -pairOffset;		
		return unpaired5p;
		
	}
	/*public float[] getBulgeMatrix5p(int l, int r, boolean unwind, boolean correspondTo3p){ 
		String str = pstructure;
		if(unwind) str = getUnWindedStructureStr(str);
		float[] unpaired5p = new float[l+r+1];
		
		for(int i=0;i<unpaired5p.length-1;i++){
			unpaired5p[i] = 1;
		}
		
		int paircntr = 0;
		for(int i=str.length()-1; i>=str.length()-l ; i--){
			if(str.charAt(i) == ')') paircntr ++;
			else if(str.charAt(i) == '(') paircntr --;
		}
		int indexInitValue = 0;
		int pairOffset = 0;
		//System.out.println(structure);
		//System.out.println(getUnWindedStructureStr(structure);
		//System.out.println(pstructure);
		//System.out.println(str);
		for(int i=0; i<l || paircntr>0; i++){
			char c = str.charAt(i);
			if(c == ')'){
				paircntr++;				
			}
			else if(c == '('){
				paircntr--;
			}
			if(i == l-1) pairOffset = paircntr;
			if(paircntr>=0) indexInitValue++;
		}
		
		int index = l-1; // index for paired5p
		float pairProb = 0;
		boolean countStart = false;
		int ll = correspondTo3p? indexInitValue : l;
		//System.out.println(ll);
		for(int i=0; i<ll;i++){//5p flnaking region
			if(index < 0) break;			
			if(str.charAt(ll-i-1) == '('){
				if(countStart){
					unpaired5p[index] = pairProb;
					
				}
				index--;
				pairProb = 0;
				countStart = true;
			}else{
				char pc = str.charAt(ll-i-1);
				if(pc == '.') pairProb = Math.max(pairProb, 1);
				else if(pc == ',') pairProb = Math.max(pairProb, 2f/3f);
				else pairProb = Math.max(pairProb, 1f/3f);
			}			
		}
		
		index = l-1; // right side
		pairProb = 0;
		//System.out.println(str.length()/2-ll-1);
		for(int i=0;i<str.length()/2-ll-1;i++){
			//System.out.println((ll+i) + " " + str.charAt(ll+i));	
			if(index >= unpaired5p.length-1) break;
			if(str.charAt(ll+i) == '('){
				unpaired5p[index++] = pairProb;	
				pairProb = 0;
			}else{
				char pc = str.charAt(ll+i);
				if(pc == '.') pairProb = Math.max(pairProb, 1);
				else if(pc == ',') pairProb = Math.max(pairProb, 2f/3f);
				else pairProb = Math.max(pairProb, 1f/3f);
			}	
		}
		if(correspondTo3p) unpaired5p[unpaired5p.length-1] = -pairOffset;
		return unpaired5p;
		
	}*/


	
	public String getBalancedStructure(char cleavageChar, boolean unwind){
		String ustr = unwind ? getUnWindedStructureStr(structure) : structure;
		StringBuilder sbl = new StringBuilder();
		StringBuilder sbr = new StringBuilder();
		//System.out.println(ustr);
		int lindex = 0;
		int rindex = ustr.length()-1;
		boolean[] cleavageSpecified = new boolean[2];
		
		while(ustr.charAt(lindex) != '(' && lindex < ustr.length()){
			sbl.append(ustr.charAt(lindex++));
			if(!cleavageSpecified[0] && lindex == flankingNTnumber){
				sbl.append(cleavageChar);
				cleavageSpecified[0] = true;
			}
		}
		
		while(ustr.charAt(rindex) != ')' && rindex >=0){
			sbr.insert(0, ustr.charAt(rindex--));
			if(!cleavageSpecified[1] && rindex == ustr.length() - flankingNTnumber - 1){
				sbr.insert(0,cleavageChar);
				cleavageSpecified[1] = true;
			}
		}
		//System.out.println(sbl + " " + sbr);
		while(lindex < rindex){
			char l = ustr.charAt(lindex);
			char r = ustr.charAt(rindex);
			//System.out.println(lindex + " " + rindex + " " + l + " " + r ); //48 113 . (
			if(l == '(' && r == ')' || l == r){
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
			if(!cleavageSpecified[1] && rindex == ustr.length() - flankingNTnumber - 1){
				sbr.insert(0,cleavageChar);
				cleavageSpecified[1] = true;
			}
			
		}
		if(lindex == rindex) sbl.append(ustr.charAt(lindex));
		
		return sbl.toString() + sbr.toString();
	}
	
	
	
	public int[] getBulgeNTcountsAround(int l, int r){
		int[] c = new int[2];
		
		for(int i=l-1;i>=0;i--){
			if(structure.charAt(i) == '.') c[0]++;
			else break;
		}
		
		for(int i=l;i<structure.length();i++){
			if(structure.charAt(i) == '.') c[0]++;
			else break;
		}
		
		for(int i=r-1;i>=0;i--){
			if(structure.charAt(i) == '.') c[1]++;
			else break;
		}
		
		for(int i=r;i<structure.length();i++){
			if(structure.charAt(i) == '.') c[1]++;
			else break;
		}
		
		return c;
	}
	
	public float[][] getBulgeMatrices(int l, int r, boolean unwind){
		float[][] bms = new float[2][l+r];		
		String bs = getBalancedStructure(' ', unwind);
		
		for(int i=0;i<bms.length;i++){
			for(int j=0;j<bms[i].length;j++){
				bms[i][j] = 1;
			}
		}
		
		int lc = bs.indexOf(' ');
		int rc = bs.lastIndexOf(' ');
		
		int offset = 0;
		char prevC = '.';
		for(int i=0;i<bms[0].length;i++){
			int index = i - l + lc;
			if(index <0) continue;
			if(index >= bs.length()) break;
			if(bs.charAt(index) == ' '){
				offset = 1;
				continue;
			}
			if(bs.charAt(index) == '(' && prevC == '(') bms[0][i-offset] = 0;	
			prevC = bs.charAt(index);
		}
		
		prevC = '.';
		offset = -1;
		for(int i=0;i<bms[1].length;i++){
			int index = rc + l - i - 1; //bs.length() - 
			
			if(index <0) break;
			if(index >= bs.length()) continue;;	
			if(bs.charAt(index) == ' '){
				offset = 0;
				continue;
			}
			//System.out.println(index);
			if(bs.charAt(index) == ')' && prevC == ')') bms[1][(i-offset)] = 0;
			prevC = bs.charAt(index);
		}
		
		return bms;
		
	}
	
	
	public float[] getBulgeMatrix3p(int l, int r, boolean unwind, boolean correspondTo5p){ 
		float[] unpaired3p = new float[l+r+1];
		String str = structure;
		if(unwind) str = getUnWindedStructureStr(str);
		
		for(int i=0;i<unpaired3p.length-1;i++){
			unpaired3p[i] = 1;
		}
		
		int paircntr = 0;
		for(int i=0; i<l ; i++){
			if(str.charAt(i) == '(') paircntr ++;
			else if(str.charAt(i) == ')') paircntr --;
		}
		//int indexInitValue = 0;
		int pairOffset = 0;
		for(int i=0; i<l || paircntr > 0 ; i++){
			char c = str.charAt(str.length()-i-1);
			if(c == '('){
				paircntr++;				
			}
			else if(c == ')'){
				paircntr--;				
			}
			if(i == l-1) pairOffset = paircntr;
		//	if(paircntr>=0) indexInitValue++;
		}
	
		ArrayList<Float> tmp = new ArrayList<Float>();
		float pairProb  = 0;
		int offset = 0;
		for(int i=0;i<r+l;i++){
			char pc = str.charAt(str.length() - i - 1);
			if(pc == ')'){
				tmp.add(pairProb);	
				pairProb = 0;
			}
			if(pc != ')'){
				if(pc == '.') pairProb = Math.max(pairProb, 1);
				else if(pc == ',') pairProb = Math.max(pairProb, 2f/3f);
				else pairProb = Math.max(pairProb, 1f/3f);
			}
			if(i == l-1) offset = tmp.size() + (correspondTo5p ? pairOffset : 0);
		}
		//System.out.println(offset);
		for(int i=0;i<Math.min(tmp.size(), r+offset);i++){
			unpaired3p[i+l-offset] =  tmp.get(i);
		}
		
		if(correspondTo5p) unpaired3p[unpaired3p.length-1] = -pairOffset;
		return unpaired3p;
		
	}
	
	/* 1111111110000000000010100   1000000000000000000100111111110
	public float[] getBulgeMatrix3p(int l, int r, boolean unwind, boolean correspondTo5p){ 
		float[] unpaired3p = new float[l+r+1];
		String str = pstructure;
		if(unwind) str = getUnWindedStructureStr(str);
		
		for(int i=0;i<unpaired3p.length-1;i++){
			unpaired3p[i] = 1;
		}
		
		int paircntr = 0;
		for(int i=0; i<l ; i++){
			if(str.charAt(i) == '(') paircntr ++;
			else if(str.charAt(i) == ')') paircntr --;
		}
		int indexInitValue = 0;
		int pairOffset = 0;
	//	System.out.println(pstructure);
	//	System.out.println(str);
		for(int i=0; i<l || paircntr > 0 ; i++){
			char c = str.charAt(str.length()-i-1);
			if(c == '('){
				paircntr++;				
			}
			else if(c == ')'){
				paircntr--;				
			}
			if(i == l-1) pairOffset = paircntr;
			if(paircntr>=0) indexInitValue++;
		}
		//System.out.println(indexInitValue);

		
		int index = l - 1; // index for paired5p
		
		float pairProb = 0;
		boolean countStart = false;
		int ll = correspondTo5p? indexInitValue : l;
		for(int i=0; i<ll;i++){//5p flnaking region
			if(index < 0) break;
			if(str.charAt(str.length() - ll + i) == ')'){
				if(countStart){
					unpaired3p[index] = pairProb;
					
				}
				index--;
				pairProb = 0;
				countStart = true;
			}else{
				char pc = str.charAt(str.length() - ll + i);
				if(pc == '.') pairProb = Math.max(pairProb, 1);
				else if(pc == ',') pairProb = Math.max(pairProb, 2f/3f);
				else pairProb = Math.max(pairProb, 1f/3f);
			}			
		}
		
		index = l - 1; // right side
		pairProb = 0;
		
		for(int i=0;i<str.length()/2-ll-1;i++){
		//	System.out.println((l+i) + " " + str.charAt(l+i));		
			if(index >= unpaired3p.length-1) break;
			if(str.charAt(str.length() - ll - 1 - i) == ')'){
				unpaired3p[index++] = pairProb;
				pairProb = 0;
			}else{
				char pc = str.charAt(str.length() - ll - 1 - i);
				if(pc == '.') pairProb = Math.max(pairProb, 1);
				else if(pc == ',') pairProb = Math.max(pairProb, 2f/3f);
				else pairProb = Math.max(pairProb, 1f/3f);
			}	
		}	
		if(correspondTo5p) unpaired3p[unpaired3p.length-1] = -pairOffset;
		return unpaired3p;
		
	}
	
	TTGTTCAGAAAGTCTGTTGTTGTAAACATCCCCGACTGGAAGCTGTAAGACACAGCTAAGCTTTCAGTCAGATGTTTGCTGCTACCGGCTATTCACAGACAT
TGATGATGCTGCTGATGCTG GCGGTGATCCCGATGGTGTGAGCTGGAAATGGGGTGCTACGtcatcgttgtcatcgtca tcatcatcatcCGAGCagcc
TGATGATGCTGCTGATGCTGGCGGTGATCCCGATGGTGTGAGCTGGAAATGGGGTGCTACGtcatcgttgtcatcgtcatcatcatcatcCGAGCagcc
TATGTGGGCAGGGCCCTGGG GAGCTGAGGCTCTGGGGGTGGCCGGGGCTGACCCTGGGCCTCTGCTCCC CAGTGTCTGACCGCGACCGC
TATGTGGGCAGGGCCCTGGGGAGCTGAGGCTCTGGGGGTGGCCGGGGCTGACCCTGGGCCTCTGCTCCCCAGTGTCTGACCGCGACCGC
AAAAAAATGGGTTCCTAGGA AGAGGTAGTAGGTTGCATAGTTTTAGGGCAGGGATTTTGCCCACAAGGAGGTAACTATACGACCTGCTGCCTTTC TTAGGGCCTTATTATTCACC
AAAAAAATGGGTTCCTAGGAAGAGGTAGTAGGTTGCATAGTTTTAGGGCAGGGATTTTGCCCACAAGGAGGTAACTATACGACCTGCTGCCTTTCTTAGGGCCTTATTATTCACC
TGACCTCTCTAACAAGGTGC AGAGCTTAGCTGATTGGTGAACAGTGATTGGTTTCCGCTTTGTTCACAGTGGCTAAGTTCTGC ACCTGAAGAGAAGGTGAGAT
TGACCTCTCTAACAAGGTGCAGAGCTTAGCTGATTGGTGAACAGTGATTGGTTTCCGCTTTGTTCACAGTGGCTAAGTTCTGCACCTGAAGAGAAGGTGAGAT
TCGATTGGACCCGCCCTCCG GTGCCTACTGAGCTGATATCAGTTCTCATTTTACACACTGGCTCAGTTCAGCAGGAACAG GAGTCGAGCCCTTGAGCAAA
	*/
	static public void main(String[] args){
		String seq = "cgagccgccgccgcccggg ccgatgcccccggcgccgcggcgcggcggaggtgtgggcgtgggcggcggcggcacgggcg tgggcggcggcgATCGCGAC"; // f
		//seq = "AAAGGACCCTTCCAG AGGGCCCCCCCTCAATCCTGTTGTGCCTAATTCAGAGGGTTGGGTGGAGGCTCTCC TGAAGGGCTCTGAAG";
		//System.out.println(seq);
		//setSeqLength(20);
		//seq = "GCUGGAAGGUGUAGGUACCCUCAAU GGCUCAGUAGCCAGUGUAGAUCCUGUCUUUCGUAAUCAGCAGCUACAUCUGGCUACUGGGUCUC UGAUGGCAUCUUCUAGCUUCUGCUU";
		System.out.println(seq);
		seq = seq.replaceAll(" ", "");
		int flanking = 20;
		
//		System.out.println(seq);
//		System.out.println(seqw);
		RNAcofoldLauncher test = new RNAcofoldLauncher(seq , flanking);
		
		
	//	Classifier classifier = new Classifier("/media/kyowon/Data1/fCLIP/samples/sample4/results/training/merged.arff");
		
		//ScoredPosition sp = new ScoredPosition("chr12	57111441	57111384	-	58	T	U	0.5383901862	12	9	9	9	_	0	NM_001178081;NM_003153;NM_001178080;NR_033659;NM_001178078	STAT6;STAT6;STAT6;STAT6;STAT6	_;_;_;_;_	NM_5_UTR;NM_5_UTR;NM_5_UTR;NR;NM_5_UTR	 ; ; ; ; 	0;0;0;0;0	0.4782023273	0.444002939	28	-0.4114	1	12	10	2	3.7532696895	3.2516291674	-0.5744490595	cttcctcttcctcct ccagcccactttctcttctctGTGTCGTCAGAGCTCCAGGGAGGGACCTGGGTAGAAG GAGAAGCCGGAAACA");
		
	//	System.out.println(classifier.classify(classifier.toInstance(sp)));
		
		//ACCTATCCATTTATCATCCATCTACGAGGTGTCTGGGATGTAATGGATGG
		//ACCUAUCCAUUUAUCAUCCAUCUACAGGUGGACCAAAAUGUAGGAGAUCU
		//TGCTCAGGCATAACCCCTCG GAACGGCTGCCCCTGGCCCAGGTCTCAGCCCACCCTTGGGTCCGGGCCAACTCTCGG AGGGTGCTGCCTCCCTCTGC
		//System.out.println(test.energyPerNT);
		//System.out.println(test.numberOfHairPins);
		//System.out.println(test.leftPaired);
		//System.out.println(test.rightPaired);
		//System.out.println(test.seqEntropy);
		//System.out.println(test.structureEntropy);
		//System.out.println(test.numberOfHairPins);
		//System.out.println(test.centerDepth);
		System.out.println(test.energyPerNT);
		System.out.println(test.numberOfHairPinsInStem);
		System.out.println(test.leftPaired);
		System.out.println(test.rightPaired);
		System.out.println(test.seqEntropy);
		System.out.println(test.structureEntropy);
		System.out.println(test.getNumberOfHairPinsInStem());
		//System.out.println(test.centerDepth);
		System.out.println(test.overHang);
		System.out.println(seq);
		System.out.println(test.getStructureString());
		//System.out.println(test.pstructure);
		
		System.out.println(test.getBalancedStructure(' ', true));
	
//  	11111111000000000001    01001  0000000000000000001001111111110
//	    11111111110000000000    00100  1000000000000010000100111111110

	//  1111111110000000000110010000000000000000001001111111111-1
	//  11111111100000000000010010000000000000100000111111111110
		for(float b : test.getBulgeMatrices(flanking, 30, true)[0])
			System.out.print((int)b + "");
		System.out.println();
		for(float b : test.getBulgeMatrices(flanking, 30, true)[1])
			System.out.print((int)b + "");
		System.out.println();
		//System.out.println(test.pstructure);
		/*System.out.println(test.has11Burge3pDown());
		System.out.println(test.has11Burge3pUp());
		System.out.println(test.has22Burge3pDown());
		System.out.println(test.has22Burge3pUp());
		System.out.println(test.has11Burge5pDown());
		System.out.println(test.has11Burge5pUp());
		System.out.println(test.has22Burge5pDown());
		System.out.println(test.has22Burge5pUp());*/
	}
}
