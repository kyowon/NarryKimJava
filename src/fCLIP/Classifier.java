package fCLIP;

import java.util.ArrayList;

import fCLIP.parser.PairScoringOutputParser.ScoredPair;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;
import weka.classifiers.trees.J48;
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
			
			rf = new RandomForest();
			//rf.setMaxDepth(2);
			rf.setNumTrees(50);
			//dt.setMinNumObj(5);
			//rf.setNumExecutionSlots(Thread.activeCount());
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
		//if(!p.isPaired() && (p.getOverHang() < 0 || p.getOverHang() > 5)) p.setClassification("T", 1);
		//else{
			double prediction = classify(toInstance(p));
			String classificaition = prediction>=0? "M" : "U";
			if(p.getSCI()<0) classificaition = classificaition.toLowerCase();
			p.setClassification(classificaition, Math.abs(prediction));
		//}
	}
	
	public void setClassification(ScoredPair p){
		//if(p.getOverHang() < 0 || p.getOverHang() > 4) p.setClassification("u", 1);
		//else{
			double prediction = classify(toInstance(p));
			p.setClassification(prediction>=0? "M" : "U", Math.abs(prediction));
		//}
	}
	
	/*sb.append(seqEntropy); sb.append(',');
			sb.append(structureEntropy); sb.append(',');
			sb.append(threePhetero); sb.append(',');
			sb.append(fivePhetero); sb.append(',');*/
	public Instance toInstance(ScoredPosition p){
		int i = 0;
		Instance inst = new DenseInstance(7);
		inst.setDataset(this.getDataset());
		inst.setValue(i++, p.getEnergy());
		inst.setValue(i++, p.getHairpinNumber());
		if(p.getSCI()<0) inst.setMissing(i);
		else inst.setValue(i, p.getSCI());
		i++;
		if(p.getZScore() == 10000) inst.setMissing(i);
		else inst.setValue(i, p.getZScore());
		i++;
		//if(p.is3pScored()) inst.setValue(4, p.getThreePScore());
		//else inst.setMissing(4);
		//if(p.is5pScored()) inst.setValue(5, p.getFivePScore());
		//else inst.setMissing(5);
	//	inst.setValue(i++, p.getSeqEntropy());
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
		if(p.getSCI()<0) inst.setMissing(i);
		else inst.setValue(i, p.getSCI());
		i++;
		if(p.getZScore() == 10000) inst.setMissing(i);
		else inst.setValue(i, p.getZScore());
		i++;
		//inst.setValue(4, p.getThreePScore());
		//inst.setValue(5, p.getFivePScore());
	//	inst.setValue(i++, p.getSeqEntropy());
		inst.setValue(i++, p.getStructureEntropy());
		if(p.getHetero()<0) inst.setMissing(i);
		else inst.setValue(i, p.getHetero());
		inst.setClassMissing();
		return inst;
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
