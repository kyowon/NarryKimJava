package rpf;

import java.util.ArrayList;
import java.util.Random;

import parser.MergedFileParser.MergedResult;
import weka.classifiers.Evaluation;
import weka.classifiers.trees.RandomForest;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSink;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

public class Classifier {
	
	static public void classify(String train, String test, ArrayList<MergedResult> mrs, int sampleIndex, boolean removeLengthAttribute){
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
			rf.buildClassifier(filteredTrainData);
			
			
			// create copy
		//	Instances labeled = new Instances(filteredTestData);
			// label instances
			System.out.println("Classifying " + filteredTestData.size() + " instances");
			for (int i = 0; i < filteredTestData.numInstances(); i++) {
				double clsLabel = rf.classifyInstance(filteredTestData.instance(i));
				double[] dist = rf.distributionForInstance(filteredTestData.instance(i));
				
				mrs.get(i).setPredictedClasses(filteredTestData.classAttribute().value((int) clsLabel), dist[(int) clsLabel], sampleIndex);
				
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
	
	static public void evaluate(String train1, String train2, boolean removeLengthAttribute){
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
			
			filteredTestData.setClassIndex(filteredTestData.numAttributes() - 1);
			filteredTrainData.setClassIndex(filteredTrainData.numAttributes() - 1);
			
			RandomForest rf = new RandomForest();
			rf.setNumTrees(50);
			rf.setNumFeatures(0);
			rf.buildClassifier(filteredTrainData);

			// evaluate classifier and print some statistics
			Evaluation eval = new Evaluation(filteredTrainData);
			eval.evaluateModel(rf, filteredTestData);
			System.out.println(eval.toSummaryString("\nResults\n\n", false));
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static public void evaluate(String train1, boolean removeLengthAttribute){
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
		//	filteredTrainData.setClassIndex(filteredTrainData.numAttributes() - 1);
			
			RandomForest rf = new RandomForest();
			rf.setNumTrees(100);
			rf.setMaxDepth(10);
			rf.setNumFeatures(0);
			rf.buildClassifier(filteredTestData);

			Evaluation eval = new Evaluation(filteredTestData);
	
			eval.crossValidateModel(rf, filteredTestData, 10, new Random(1));
			System.out.println(eval.toSummaryString("\nResults\n\n", false));

			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static public void main(String[] args){
		String train1 = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out_0.3.csv_train/train_0.arff";
		String train2 = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out2_0.3.csv_train/train_0.arff";
		String test = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out_0.3.csv_test/test_0.arff";
		String output = "/media/kyowon/Data1/RPF_Project/samples/sample3/results/out_0.3.csv_test/test_0_classified.arff";
		boolean removeLengthAttribute = true;
		
		Classifier.evaluate(train1, train2, removeLengthAttribute);
	//	Classifier.classify(train1, test, output, removeLengthAttribute);
	}
	
}
