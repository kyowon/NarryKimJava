package fCLIP.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import util.Arff2Mfile;

public class GenerateMfile {
    
	public static void generateMfileForAttFigure(String mOut, String arff, String attm, String outFolder){
		Arff2Mfile.generate(arff, attm);
		try {
			PrintStream out = new PrintStream(mOut);
			String attmName = new File(attm).getName();
			out.println(attmName.substring(0, attmName.length()-2));
			out.println("for i=1:size(att,1)-1");
			out.println("for j=i+1:size(att,1)-1");
			out.println("createAttFigure(value, att, i, j, group, '" + outFolder + "');");
			out.println("end\nend");
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void generateMfileForFigure(String mOut, String figFolder, String prefix, String[] cks, String[] kks, String mFolder, int numThreshold, String[] legends){
		try {
			ArrayList<String> mfiles = new ArrayList<String>();
			for(File f : new File(mFolder).listFiles()){
				if(f.getName().endsWith(".m") && (!f.getName().contains("gene"))) mfiles.add(f.getName());
			}
			PrintStream out = new PrintStream(mOut);
			out.println("close all;");
			for(String mfile : mfiles){
				if(!mfile.startsWith(prefix) && !mfile.startsWith("background")) continue;
				out.println(mfile.substring(0, mfile.length()-2));
			}
			//for(int r=0;r<6;r++){
			for(int i=0;i<cks.length;i++){
				String cko = cks[i];
				String ck = cko.replace('-', '_');
				String kko = kks[i];
				String kk = kko.replace('-', '_');
				String figOut = figFolder + ck + kk;
				out.println("clear cs" + (i+1));
				out.println("cs" + (i+1) + " = {");
				for(String legend : legends){
					out.println(legend + "_" + ck +";");
				}
				out.println("};");
				out.println("clear ts" + (i+1));
				out.println("ts" + (i+1) + " = {");
				for(String legend : legends){
					out.println(legend + "_" + kk + ";");
				}
				out.println("};");
				
				out.println("legends = {");
				for(String legend : legends){
					out.println("'" + legend +"';");
				}
				out.println("};");
				out.println("%createFC(cs" + (i+1)+ ", ts"  + (i+1) + ", legends, " + numThreshold + ",'" + figOut + "','" + cko + " vs "+ kko + "');");
			}
			out.println("clear cs;");
			out.println("clear ts;");
			
			out.println("for i=1:size(cs1,1)");
			out.println("\tcs{i}=cs1{i}+cs2{i}+cs3{i};");
			out.println("\tts{i}=ts1{i}+ts2{i}+ts3{i};");
			out.println("end");

			out.println("for i=1:size(cs,2)");
			out.println("\tfor j=3:size(cs{i},1)");
			out.println("\t\tcs{i}(j,:)=cs{i}(j,:)/3;");
			out.println("\t\tts{i}(j,:)=ts{i}(j,:)/3;");
			out.println("\tend");
			out.println("end");
			
			out.println("cs=cs';");
			out.println("ts=ts';");
			
			//out.println("createFC(cs, ts, legends, " + numThreshold + ",'" + figFolder  + "','" + cko + " vs "+ kko + "');");
			
			out.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
	}
	
	public static void generateGeneMfileForFigure(String mOut, String figFolder, String prefix, String[] cks, String[] kks, String mFolder, int numThreshold, String[] legends){
		try {
			ArrayList<String> mfiles = new ArrayList<String>();
			for(File f : new File(mFolder).listFiles()){
				if(f.getName().endsWith("gene.m")) mfiles.add(f.getName());
			}
			PrintStream out = new PrintStream(mOut);
			out.println("close all;");
			for(String mfile : mfiles){
			//	if(!mfile.startsWith(prefix) && !mfile.startsWith("background")) continue;
				out.println(mfile.substring(0, mfile.length()-2));
			}
			//for(int r=0;r<6;r++){
			for(int i=0;i<cks.length;i++){
				String cko = cks[i];
				String ck = cko.replace('-', '_');
				String kko = kks[i];
				String kk = kko.replace('-', '_');
				String figOut = figFolder + ck + kk;
				
				out.println("clear cs" + (i+1));
				out.println("cs" + (i+1) + " = {");
				for(String legend : legends){
					out.println("Gene" + legend + "_" + ck +";");
				}
				out.println("};");
				out.println("clear ts" + (i+1));
				out.println("ts" + (i+1) + " = {");
				for(String legend : legends){
					out.println("Gene" + legend + "_" + kk + ";");
				}
				out.println("};");
				
				out.println("legends = {");
				for(String legend : legends){
					out.println("'" + legend +"';");
				}
				out.println("};");
				out.println("%createGeneFC(cs" + (i+1)+ ", ts"  + (i+1) + ", legends, " + numThreshold + ",'" + figOut + "','" + cko + " vs "+ kko + "');");
				
				out.println("clear acccs" + (i+1));
				out.println("acccs" + (i+1) + " = {");
				for(String legend : legends){
					out.println("Acc" + legend + "_" + ck +";");
				}
				out.println("};");
				out.println("clear accts" + (i+1));
				out.println("accts" + (i+1) + " = {");
				for(String legend : legends){
					out.println("Acc" + legend + "_" + kk + ";");
				}
				out.println("};");
				
				//out.println("%createGeneFC(cs" + (i+1)+ ", ts"  + (i+1) + ", legends, " + numThreshold + ",'" + figOut + "','" + cko + " vs "+ kko + "');");
			}
			out.println("clear cs;");
			out.println("clear ts;");
			
			
			out.println("for i=1:size(cs1,1)");
			out.println("\tcs{i}=cs1{i}+cs2{i}+cs3{i};");
			out.println("\tts{i}=ts1{i}+ts2{i}+ts3{i};");
			out.println("end");

			out.println("for i=1:size(cs,2)");
			out.println("\tfor j=2:size(cs{i},1)");
			out.println("\t\tcs{i}(j,:)=cs{i}(j,:)/3;");
			out.println("\t\tts{i}(j,:)=ts{i}(j,:)/3;");
			out.println("\tend");
			out.println("end");
			
		
			out.println("cs=cs';");
			out.println("ts=ts';");
			
			out.println("clear acccs;");
			out.println("clear accts;");
			

			
		
			out.println("acccs=acccs1';");
			out.println("accts=accts1';");

			//out.println("createFC(cs, ts, legends, " + numThreshold + ",'" + figFolder  + "','" + cko + " vs "+ kko + "');");
			
			out.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		if(args[args.length - 1].equals("ATT")){
			String mOut = args[0];
			String arff = args[1];
			String attm = args[2];
			String outFolder = args[3];
			generateMfileForAttFigure(mOut, arff, attm, outFolder);
		}
		if(args[args.length - 1].equals("FIG")){
			String mOut = args[0];
			String figFolder = args[1];
			String prefix = args[2];
			String[] controlKeys = args[3].split(",");
			String[] kdKeys = args[4].split(",");
			String mFolder = args[5];
			int numThreshold = Integer.parseInt(args[6]);
			String[] legends = args[7].split(",");
			generateMfileForFigure(mOut, figFolder, prefix, controlKeys, kdKeys, mFolder, numThreshold, legends);
		}
		if(args[args.length - 1].equals("GENEFIG")){
			String mOut = args[0];
			String figFolder = args[1];
			String prefix = args[2];
			String[] controlKeys = args[3].split(",");
			String[] kdKeys = args[4].split(",");
			String mFolder = args[5];
			int numThreshold = Integer.parseInt(args[6]);
			String[] legends = args[7].split(",");
			generateGeneMfileForFigure(mOut, figFolder, prefix, controlKeys, kdKeys, mFolder, numThreshold, legends);
		}
	}

}
