package launcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class RNAfoldLauncher {
	
	//private String seq;
	private double energy = 0;
	private double energyPerNT = 0;
	
	private int flankingNTnumber = 11;
	private int depth = 0;
	private int numberOfHairPins = 0;
	private int overHang = -1000;
	private int leftPaired = 0;
	private int rightPaired = 0;
	private double seqEntropy = 0;
	private double structureEntropy = 0;
	public final static int[] overhanglimit = {0,5};
	public final static int[] pairedNumberlimit = {8,6};
	
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
		this.flankingNTnumber = flankingNTnumber;
		if(seq != null && !seq.isEmpty())
			run(seq);
	}
	
	private void run(String seq){//seq = pre + flanking.
		try {
			String pseq = seq;//.replace(' ', '\0');
			
			String[] cmd = {
					"/bin/sh",
					"-c",
					"echo "+ pseq + " | RNAfold --noPS"
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
			
			//for(String kk : result) 
			//System.out.println(builder.toString());
			
			int maxDepth = 0;
		//	int maxCenterDepth = 0;
			String st = re[1].substring(0, re[1].indexOf(' '));
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
			
		//	if(nh == 1){			
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
		String seq = "CGCCACGCCTCCGCTGGCGACGGGACATTATTACTTTTGGTACGCGCTGTGACACTTCAAACTCGTACCGTGAGTAATAATGCGCCGTCCACGGCACCGCATCGAAAAC";
		seq = seq.replaceAll(" ", "");
		String seqw = "";
		for(char s : seq.toCharArray()){
			seqw = s + seqw;
		}
//		System.out.println(seq);
//		System.out.println(seqw);
		RNAfoldLauncher test = new RNAfoldLauncher(seq , 20);
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
	}
	
}
