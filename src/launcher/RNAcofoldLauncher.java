package launcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

import fCLIP.Classifier;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class RNAcofoldLauncher {
	private double energy = 0;
	private double energyPerNT = 0;
	
	private int flankingNTnumber = 25;
	private int depth = 0;
	private int numberOfHairPins = 0;
	private int overHang = -1000;
	private int leftPaired = 0;
	private int rightPaired = 0;
	private int loopLength = 0;
	private int leftUpstreamPaired = 0;
	private int rightUpstreamPaired = 0;
	
	
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

	public int getNumberOfHairPins(){
		return numberOfHairPins;
	}
	
	public int getLoopLength(){
		return loopLength;
	}
	
	public RNAcofoldLauncher(String seq, int flankingNTnumber){
		this.flankingNTnumber = flankingNTnumber;
		if(seq != null && !seq.isEmpty())
			run(seq);
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
			
			
			
			String pseq = pseqb.toString();
		//	System.out.println(pseq);
			//.replace("&", "\\&");
			String[] cmd = {
				"/bin/sh",
				"-c",
				"echo "+ pseq + " | RNAcofold --noPS"
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
			//for(String kk : re) 
			//	System.out.println(builder.toString());
			
			int maxDepth = 0;
			String st = re[1].substring(0, re[1].indexOf(' ')).replace("&", "");;
		//	System.out.println(st);
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
			
			int nh = 0;
			int ln = 0;
			boolean left = true;
			for(int i=0;i<st.length();i++){
				char c = st.charAt(i);
				if(c == '('){
					d ++;
					left = true;
					ln = 0;
				}else if(c == ')'){
					d --;
					if(left){
						nh++;
						loopLength = loopLength > ln ? loopLength : ln;
					}
					ln = 0;
					left = false;
				}else if(c == '.'){
					ln++;
				}
				//else if(c == '.'){
				//	if(i >= flankingNTnumber) r++;
				//	if(i< st.length() - flankingNTnumber) l++;
				//}
				maxDepth = maxDepth > d ? maxDepth : d;
			}
			
			for(int i=flankingNTnumber-1;i>=0;i--){
			//	System.out.println("l");
				if(st.charAt(i) == '(')//break;
					l++;
				if(st.charAt(i) == '.' && st.charAt(st.length() - 1 - i) == '.') break;
			}
			
			for(int i=flankingNTnumber;i<st.length()/2;i++){
				if(st.charAt(i) == '(')//break;
					lu++;
			}
			
			for(int i=st.length() - flankingNTnumber;i<st.length();i++){ //st.length()-1
				if(st.charAt(i) == ')') //break;
					r++;
				if(st.charAt(i) == '.' && st.charAt(st.length() - 1 - i) == '.') break;
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
			this.numberOfHairPins = nh;
		//	this.centerDepth = maxCenterDepth;
			this.energy = Double.parseDouble(re[1].substring(re[1].indexOf(' ')+2, re[1].length()-1).trim());
			
			pseq = pseq.replace("\\&", "");
			this.energyPerNT = energy/pseq.length();
		//	System.out.println(pseq.length());
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
	
	
	static public void main(String[] args){
		String seq = "CTTCAGCAAACATTAGGAGAGTATC TTCTCTGTTTTGGCCATGTGTGTACTCACAGCCCCTCACACATGGCCGAAACAGAGAAGT TACTTTCCTAATATTtgcctccttg";
		//System.out.println(seq);
		//setSeqLength(20);
		seq = seq.replaceAll(" ", "");
		String seqw = "";
		for(char s : seq.toCharArray()){
			seqw = s + seqw;
		}
		
//		System.out.println(seq);
//		System.out.println(seqw);
		RNAcofoldLauncher test = new RNAcofoldLauncher(seq , 25);
		
		
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
		System.out.println(test.overHang);
		System.out.println(test.loopLength);
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
