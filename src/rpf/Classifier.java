package rpf;

import java.util.ArrayList;
import java.util.Random;

import parser.MergedFileParser.MergedResult;
import weka.classifiers.Evaluation;
import weka.classifiers.trees.RandomForest;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

public class Classifier {
	
	static public void classify(String train, String test, ArrayList<MergedResult> mrs, int sampleIndex, boolean removeLengthAttribute, double predictionThreshold){
		try {
			Instances testData = DataSource.read(test);
			Instances trainData = DataSource.read(train);		
			Instances filteredTestData = testData;
			Instances filteredTrainData = trainData;
			
			if(removeLengthAttribute){
				String[] options = new String[2];
				options[0] = "-R";
				options[1] = "14";
				Remove rm = new Remove();
				rm.setOptions(options);
				rm.setInputFormat(testData);
				filteredTestData = Filter.useFilter(testData, rm);
				filteredTrainData = Filter.useFilter(trainData, rm);
			}
			
			filteredTestData.setClassIndex(filteredTestData.numAttributes() - 1);
			filteredTrainData.setClassIndex(filteredTrainData.numAttributes() - 1);
			
			RandomForest rf = new RandomForest();
			
			rf.setNumTrees(100);
			rf.setMaxDepth(10);
			rf.setNumFeatures(0);
			//rf.setNumExecutionSlots(Thread.activeCount());
			rf.buildClassifier(filteredTrainData);
			
			
			// create copy
		//	Instances labeled = new Instances(filteredTestData);
			// label instances
			System.out.println("Classifying " + filteredTestData.size() + " instances");
			for (int i = 0; i < filteredTestData.numInstances(); i++) {
				double clsLabel = rf.classifyInstance(filteredTestData.instance(i));
				double[] dist = rf.distributionForInstance(filteredTestData.instance(i));
				String c = filteredTestData.classAttribute().value((int) clsLabel);
				double prediction = dist[(int) clsLabel];
				if(prediction < predictionThreshold) c = c.toLowerCase();
				mrs.get(i).setPredictedClasses(c, prediction, sampleIndex);
				
				//System.out.print((i+1) + " - ");
				//System.out.print(filteredTestData.instance(i).toString(filteredTestData.classIndex()) + " - ");
				//System.out.println( + " - " + );
			}
			// save newly labeled data
			
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static public void evaluate(String train1, String train2, boolean removeLengthAttribute, double predictionThreshold){
		try {
			Instances testData = DataSource.read(train1);
			Instances trainData = DataSource.read(train2);		
			Instances filteredTestData = testData;
			Instances filteredTrainData = trainData;
			
			if(removeLengthAttribute){
				String[] options = new String[2];
				options[0] = "-R";
				options[1] = "14";
				Remove rm = new Remove();
				rm.setOptions(options);
				rm.setInputFormat(testData);
				
				filteredTestData = Filter.useFilter(testData, rm);
				filteredTrainData = Filter.useFilter(trainData, rm);
			}
			
			/*
			RemoveWithValues rmv = new RemoveWithValues();
			rmv.setAttributeIndex("11");
			rmv.setSplitPoint(6);
			//rmv.setInvertSelection(true);
			rmv.setInputFormat(filteredTestData);
			
			
		//	filteredTestData = Filter.useFilter(filteredTestData, rmv);
		//	filteredTrainData = Filter.useFilter(filteredTrainData, rmv);
			*/
			filteredTestData.setClassIndex(filteredTestData.numAttributes() - 1);
			filteredTrainData.setClassIndex(filteredTrainData.numAttributes() - 1);
			
			
			RandomForest rf = new RandomForest();
			
			rf.setNumTrees(100);
			rf.setMaxDepth(10);
			rf.setNumFeatures(0);
			//rf.setNumExecutionSlots(Thread.activeCount());
			rf.buildClassifier(filteredTrainData);
			
			System.out.println("Number of Instances before prediction score filtering\t"+filteredTestData.size());
			for (int i = 0; i < filteredTestData.numInstances(); i++) {
				//System.out.println(filteredTestData.instance(i));
				double clsLabel = rf.classifyInstance(filteredTestData.instance(i));
				double[] dist = rf.distributionForInstance(filteredTestData.instance(i));
				if(dist[(int) clsLabel] < predictionThreshold){
					filteredTestData.remove(i);
					//System.out.println(clsLabel + " " + dist[(int) clsLabel]);
				}
			}
			System.out.println("Number of Instances after prediction score filtering\t"+filteredTestData.size());			
			
			// evaluate classifier and print some statistics
			Evaluation eval = new Evaluation(filteredTrainData);
			eval.evaluateModel(rf, filteredTestData);
			
			System.out.println(eval.toSummaryString("\nResults\n\n", true));
			System.out.println(eval.toMatrixString());
			System.out.println(eval.toClassDetailsString());
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static public void evaluate(String train1, boolean removeLengthAttribute, double predictionThreshold){
		try {
			Instances testData = DataSource.read(train1);
			//Instances trainData = DataSource.read(train2);		
			Instances filteredTestData = testData;
			//Instances filteredTrainData = trainData;
			
			if(removeLengthAttribute){
				String[] options = new String[2];
				options[0] = "-R";
				options[1] = "14";
				Remove rm = new Remove();
				rm.setOptions(options);
				rm.setInputFormat(testData);
				filteredTestData = Filter.useFilter(testData, rm);
			//	filteredTrainData = Filter.useFilter(trainData, rm);
			}
			filteredTestData.setClassIndex(filteredTestData.numAttributes() - 1);
			RandomForest rf = new RandomForest();
			
			rf.setNumTrees(100);
			rf.setMaxDepth(10);
			rf.setNumFeatures(0);
			//rf.setNumExecutionSlots(Thread.activeCount());
			rf.buildClassifier(filteredTestData);
			
			System.out.println("Number of Instances before prediction score filtering\t"+filteredTestData.size());
			for (int i = 0; i < filteredTestData.numInstances(); i++) {
				double clsLabel = rf.classifyInstance(filteredTestData.instance(i));
				double[] dist = rf.distributionForInstance(filteredTestData.instance(i));
				if(dist[(int) clsLabel] < predictionThreshold) filteredTestData.remove(i);
			}
			System.out.println("Number of Instances after prediction score filtering\t"+filteredTestData.size());
			
		//	filteredTrainData.setClassIndex(filteredTrainData.numAttributes() - 1);
			
			Evaluation eval = new Evaluation(filteredTestData);
	
			eval.crossValidateModel(rf, filteredTestData, 10, new Random(1));
			System.out.println(eval.toSummaryString("\nResults\n\n", true));
			System.out.println(eval.toMatrixString());
			System.out.println(eval.toClassDetailsString());
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static public void main(String[] args){
		String train1 = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out1_0.3.csv_classification/train_0.arff";
		String train2 = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out2_0.3.csv_classification/train_0.arff";
		//String test = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out_0.3.csv_test/test_0.arff";
		//String output = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out_0.3.csv_test/test_0_classified.arff";
		boolean removeLengthAttribute = true;
		
		Classifier.evaluate(train1, train2, removeLengthAttribute, .75);
	//	Classifier.classify(train1, test, output, removeLengthAttribute);
	}
	
}
