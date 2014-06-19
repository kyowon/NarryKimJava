package fCLIP;

import java.util.ArrayList;

import fCLIP.parser.PairScoringOutputParser.ScoredPair;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;
import weka.classifiers.trees.J48;
import weka.classifiers.trees.RandomForest;
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
			rf.setNumTrees(100);
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 0;
	}
	
	public void setClassification(ScoredPosition p){
		double prediction = classify(p.toInstance(trainData));
		p.setClassification(prediction>=0? "M" : "U", Math.abs(prediction));
	}
	
	public void setClassification(ScoredPair p){
		double prediction = classify(p.toInstance(trainData));
		p.setClassification(prediction>=0? "M" : "U", Math.abs(prediction));
	}
}
