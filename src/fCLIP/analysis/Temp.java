package fCLIP.analysis;

import java.util.ArrayList;
import java.util.HashSet;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class Temp {

	public static void main(String[] args) {
		ScoredPairOutputParser parser = new ScoredPairOutputParser("/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/trans/merged.trans.complete.M.csv");
		int threshold = 600;
		HashSet<ScoredPosition> tpositions = new HashSet<ScoredPosition>();
		HashSet<ScoredPosition> fpositions = new HashSet<ScoredPosition>();
		
		for(ScoredPair pair1 : parser.getPairs()){
			tpositions.add(pair1.getPairedScoredPositions()[0]);
			fpositions.add(pair1.getPairedScoredPositions()[1]);
		}
		
		HashSet<ScoredPosition> list = new HashSet<ScoredPosition>();
		for(ScoredPosition p1 : tpositions){
			for(ScoredPosition p2 : tpositions){
				if(p1.equals(p2)) continue;
				if(!p1.getContig().equals(p2.getContig())) continue;
				if(Math.abs(p1.getFivePPosition() - p2.getFivePPosition()) < threshold)
					list.add(p1);
				if(Math.abs(p1.getThreePPosition() - p2.getThreePPosition()) < threshold)
					list.add(p1);
			}	
		}
		System.out.println("3p");
		for(ScoredPosition p : list)
			System.out.println(p);
		list.clear();
		
		for(ScoredPosition p1 : fpositions){
			for(ScoredPosition p2 : fpositions){
				if(p1.equals(p2)) continue;
				if(!p1.getContig().equals(p2.getContig())) continue;
				if(Math.abs(p1.getFivePPosition() - p2.getFivePPosition()) < threshold)
					list.add(p1);
				if(Math.abs(p1.getThreePPosition() - p2.getThreePPosition()) < threshold)
					list.add(p1);
			}	
		}
		System.out.println("5p");
		for(ScoredPosition p : list)
			System.out.println(p);
		
	}

}
