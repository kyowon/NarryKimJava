package fCLIP;

import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;
import weka.classifiers.trees.RandomForest;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

public class Classifier {
	private RandomForest rf;
	private Instances trainData;
	public Classifier(String train){
		try {
			trainData = DataSource.read(train);
			trainData.setClassIndex(trainData.numAttributes() - 1);
			//rf.setNumExecutionSlots(Thread.activeCount());
			
			rf = new RandomForest();
			//rf.setMaxDepth(2);
			
			rf.setNumTrees(30);

			//rf.setMinNumObj(5);
			rf.buildClassifier(trainData);
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public Instances getDataset(){
		return trainData;
	}
	public double classify(Instance inst){
		try {
			double clsLabel = rf.classifyInstance(inst);
			double[] dist = rf.distributionForInstance(inst);
			//System.out.println(clsLabel);
			//String c = inst.classAttribute().value((int) clsLabel); 
			double prediction = dist[(int) clsLabel];
			return prediction * (clsLabel> 0 ? -1 : 1);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return 0;
	}
	
	public void setClassification(ScoredPosition p){
		if(!p.isFeasibleFold()) p.setClassification("U", 1);
		else{
		//if(!p.isPaired() && (p.getOverHang() < 0 || p.getOverHang() > 5)) p.setClassification("T", 1);
		//else{
			double prediction = classify(toInstance(p));
			String classificaition = prediction>=0.0? (prediction >= 0.7?"M" : "m") : "U";
			
			//if(p.getSCI()<0) 
			//	classificaition = classificaition.toLowerCase();
			p.setClassification(classificaition, Math.abs(prediction));
		}
	}
	
	public void setClassification(ScoredPair p){
		if(!p.isFeasibleFold()) p.setClassification("U", 1);
		//if(p.getOverHang() < 0 || p.getOverHang() > 4) p.setClassification("u", 1);
		else{
			double prediction = classify(toInstance(p));
			p.setClassification(prediction>=0.0? (prediction >= 0.7?"M" : "m") : "U", Math.abs(prediction));
		}
	}
	
	public Instance toInstance(ScoredPosition p){
		int i = 0;
		Instance inst = new DenseInstance(7);
		inst.setDataset(this.getDataset());
		inst.setValue(i++, p.getEnergy());
		inst.setValue(i++, p.getHairpinNumber());
		//if(p.getSCI()<0) inst.setMissing(i);
		//else inst.setValue(i, p.getSCI());
		//i++;
		//if(p.getZScore() == 10000) inst.setMissing(i);
		//else inst.setValue(i, p.getZScore());
		//i++;
		if(p.is3pScored()) inst.setValue(i, p.getThreePScore());
		else inst.setMissing(i);
		i++;
		if(p.is5pScored()) inst.setValue(i, p.getFivePScore());
		else inst.setMissing(i);
		i++;
		inst.setValue(i++, p.getStructureEntropy());
		if(p.getHetero()<0) inst.setMissing(i);
		else inst.setValue(i, p.getHetero());
		inst.setClassMissing();
		return inst;
	}
	
	public Instance toInstance(ScoredPair p){
		int i = 0;
		Instance inst = new DenseInstance(7);
		inst.setDataset(this.getDataset());
		inst.setValue(i++, p.getEnergy());
		inst.setValue(i++, p.getHairpinNumber());
		//if(p.getSCI()<0) inst.setMissing(i);
		//else inst.setValue(i, p.getSCI());
		//i++;
		//if(p.getZScore() == 10000) inst.setMissing(i);
		//else inst.setValue(i, p.getZScore());
		//i++;
		inst.setValue(i++, p.getThreePScore());
		inst.setValue(i++, p.getFivePScore());
		inst.setValue(i++, p.getStructureEntropy());
		if(p.getHetero()<0) inst.setMissing(i);
		else inst.setValue(i, p.getHetero());
		inst.setClassMissing();
		return inst;
	}
	
	static public String toArffString(ScoredPair pair){
		
		StringBuilder sb = new StringBuilder();
		sb.append(pair.getEnergy()); sb.append(',');
		sb.append(pair.getHairpinNumber()); sb.append(',');
		//sb.append(pair.getSCI()<0? "?" : pair.getSCI()); sb.append(',');
		//sb.append(pair.getZScore()==10000? "?" : pair.getZScore()); sb.append(',');
		sb.append(pair.getThreePScore()); sb.append(',');
		sb.append(pair.getFivePScore()); sb.append(',');
	//	sb.append(seqEntropy); sb.append(',');
		sb.append(pair.getStructureEntropy()); sb.append(',');
		sb.append("?"); sb.append(',');
		//sb.append(fivePhetero); sb.append(',');
		
		sb.append(pair.getClassification() == null? "?" : pair.getClassification().toUpperCase());
		return sb.toString();
	}
	
	//if(overHang < overhanglimit[0] || overHang > overhanglimit[1]) return false;
	
	static public String toTrainingArffString(ScoredPosition p){
		StringBuilder sb = new StringBuilder();
		sb.append(p.getEnergy()); sb.append(',');
		sb.append(p.getHairpinNumber()); sb.append(',');
		//sb.append(p.getSCI()<0? "?" : p.getSCI()); sb.append(',');
		//sb.append(p.getZScore()==10000? "?" : p.getZScore()); sb.append(',');
		sb.append(p.is3pScored()? p.getThreePScore() : "?"); sb.append(',');
		sb.append(p.is5pScored()? p.getFivePScore() : "?"); sb.append(',');
	//	sb.append(seqEntropy); sb.append(',');
		sb.append(p.getStructureEntropy()); sb.append(',');
		sb.append(p.getHetero()); sb.append(',');
		//boolean isOverHangOK = p.getOverHang() >= 0 && p.getOverHang() <= 5 && ;
		if(!p.hasMatchingMiRNA()){
			if(!p.isFeasibleFold()) sb.append("U");
			else sb.append("?");
		}else{
			if(!p.isFeasibleFold()) sb.append("?");
			else sb.append('M');
		}
		return sb.toString();
	}
	
	static public String toArffString(ScoredPosition p){
		StringBuilder sb = new StringBuilder();
		sb.append(p.getEnergy()); sb.append(',');
		sb.append(p.getHairpinNumber()); sb.append(',');
		//sb.append(p.getSCI()<0? "?" : p.getSCI()); sb.append(',');
		//sb.append(p.getZScore()==10000? "?" : p.getZScore()); sb.append(',');
		sb.append(p.is3pScored()? p.getThreePScore() : "?"); sb.append(',');
		sb.append(p.is5pScored()? p.getFivePScore() : "?"); sb.append(',');
	//	sb.append(seqEntropy); sb.append(',');
		sb.append(p.getStructureEntropy()); sb.append(',');
		sb.append(p.getHetero()); sb.append(',');
		String c = (p.getClassification() == null? "?" : p.getClassification().toUpperCase());
		//if(c.equals("M")){
		//	if(this.miRNAs != null && !this.miRNAs.isEmpty()) c = "M"; 
		//	else c = "m";	//}
		//if(this.isPaired()) c +="P";
		sb.append(c);
		return sb.toString();
	}
	
	public static String getArffHeader(){
		return "@relation w\n"
				+ "\n@attribute Energy numeric"
				+ "\n@attribute HairpinNumber numeric"
				//+ "\n@attribute ZScore numeric"
				+ "\n@attribute 5PScore numeric"
				+ "\n@attribute 3PScore numeric"
				+ "\n@attribute StructureEntropy numeric"
				+ "\n@attribute Hetero numeric"
				+ "\n@attribute Class {M,U}\n\n@data";
	}
	/*	public String toArffString(){
			StringBuilder sb = new StringBuilder();
			sb.append(energy); sb.append(',');
			sb.append(hairpinNum); sb.append(',');
			sb.append(sci<0? "?" : sci); sb.append(',');
			sb.append(zscore==10000? "?" : zscore); sb.append(',');
			sb.append(is3pScored()? threePScore : "?"); sb.append(',');
			sb.append(is5pScored()? fivePScore : "?"); sb.append(',');
			sb.append(classification == null? "?" : classification);
			return sb.toString();
		}
		*/
}
