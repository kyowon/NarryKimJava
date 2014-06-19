package fCLIP;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

public class RNAfoldLauncher {
	
	//private String seq;
	double energy = 0;
	int flankingNTnumber = 11;
	int depth = 0;
	int numberOfHairPins = 0;
	int overHang = -1000;
	int leftPaired = 0;
	int rightPaired = 0;
	
	public double getEnergy() {
		return energy;
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

	public int getNumberOfHairPins(){
		return numberOfHairPins;
	}
	
	public RNAfoldLauncher(String seq, int flankingNTnumber){
		this.flankingNTnumber = flankingNTnumber;
		if(seq != null && !seq.isEmpty())
			run(seq);
	}
	
	public Instance toInstance(Instances insts){
		Instance inst = new DenseInstance(6);
		inst.setDataset(insts);
	//	inst.setValue(0, getDepth());
		inst.setValue(0, getEnergy());
		inst.setValue(1, getNumberOfHairPins());
		
		inst.setMissing(2);
		inst.setMissing(3);
	//	inst.setValue(3, getLeftPaired());
	//	inst.setValue(4, getRightPaired());
		inst.setValue(4, getOverHang());
		inst.setClassMissing();
		return inst;
	}
	
	private void run(String seq){
		try {
			String[] cmd = {
					"/bin/sh",
					"-c",
					"echo "+ seq + " | RNAfold --noPS"
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
				int td = 0;
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
				}
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
			this.energy = Double.parseDouble(re[1].substring(re[1].indexOf(' ')+2, re[1].length()-1).trim()) / seq.length();
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
	 
	static public void main(String[] args){
		RNAfoldLauncher test = new RNAfoldLauncher("CCTCCGTGGGGTCTGTCCCCGCTGCTTGCTGTCCTCCGTGGGGTCTCTTCCCCCGGCTCGCTGTCCTCCATGGGGTCTGTCCCCCGGTTTG" , Scorer.flankingNTNumber);
		//ACCTATCCATTTATCATCCATCTACGAGGTGTCTGGGATGTAATGGATGG
		//ACCUAUCCAUUUAUCAUCCAUCUACAGGUGGACCAAAAUGUAGGAGAUCU
		System.out.println(test.energy);
		System.out.println(test.numberOfHairPins);
		System.out.println(test.leftPaired);
		System.out.println(test.rightPaired);
		
		//System.out.println(test.centerDepth);
		System.out.println(test.overHang);
	}
	
}
