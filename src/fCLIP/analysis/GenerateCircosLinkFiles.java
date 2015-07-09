package fCLIP.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class GenerateCircosLinkFiles {
	static void run(String pairCsv, String link, String histo5p, String histo3p, String gene5p, String gene3p, String scatter5p, String scatter3p){
		ScoredPairOutputParser parser = new ScoredPairOutputParser(pairCsv);
		PrintStream linkOut;
		try {
			linkOut = new PrintStream(link);
			
			HashMap<ScoredPosition, String> sMap = new HashMap<ScoredPosition, String>();
			HashMap<ScoredPosition, String> eMap = new HashMap<ScoredPosition, String>();
			
			HashMap<ScoredPosition, HashMap<String, Integer>> shistMap = new HashMap<ScoredPosition, HashMap<String, Integer>>();
			HashMap<ScoredPosition, HashMap<String, Integer>> ehistMap = new HashMap<ScoredPosition, HashMap<String, Integer>>();
			ArrayList<String> linkColors = new ArrayList<String>();
			ArrayList<String> colors = new ArrayList<String>();
						
			colors.add("normal");colors.add("etc");colors.add("line");colors.add("sine");
			
			String[] oheader = ScoredPair.getHeader().split("\t");
			String[] header = parser.getHeader().split("\t");
			int index = 0;
			for(int i=oheader.length;i<header.length;i++){
				if(header[i].equals("5p rmsk")){
					index = i;
					break;
				}
			}
			index -= oheader.length;
			
			for(ScoredPair p : parser.getPairs()){
				String[] misc = p.getMisc().split("\t");				
				String linkColor, linkStartColor, linkEndColor;
				if(misc[index].startsWith("SINE")){
					linkStartColor = "sine";
					if(misc[index+1].startsWith("SINE")){
						linkEndColor = "sine";
						linkColor = "ss";
					}else if(misc[index+1].startsWith("LINE")){
						linkEndColor = "line";
						linkColor = "sl";
					}else if(misc[index+1].startsWith("_")){
						linkEndColor = "normal";
						linkColor = "sn";
					}else{
						linkEndColor = "etc";
						linkColor = "se";
					}
				}else if(misc[index].startsWith("LINE")){
					linkStartColor = "line";
					if(misc[index+1].startsWith("SINE")){
						linkEndColor = "sine";
						linkColor = "sl";
					}else if(misc[index+1].startsWith("LINE")){
						linkEndColor = "line";
						linkColor = "ll";
					}else if(misc[index+1].startsWith("_")){
						linkEndColor = "normal";
						linkColor = "ln";
					}else{
						linkEndColor = "etc";
						linkColor = "le";
					}
				}else if(misc[index].startsWith("_")){
					linkStartColor = "normal";
					if(misc[index+1].startsWith("SINE")){
						linkEndColor = "sine";
						linkColor = "sn";
					}else if(misc[index+1].startsWith("LINE")){
						linkEndColor = "line";
						linkColor = "ln";
					}else if(misc[index+1].startsWith("_")){
						linkEndColor = "normal";
						linkColor = "nn";
					}else{
						linkEndColor = "etc";
						linkColor = "ne";
					}
				}else{
					linkStartColor = "etc";
					if(misc[index+1].startsWith("SINE")){
						linkEndColor = "sine";
						linkColor = "se";
					}else if(misc[index+1].startsWith("LINE")){
						linkEndColor = "line";
						linkColor = "le";
					}else if(misc[index+1].startsWith("_")){
						linkEndColor = "normal";
						linkColor = "ne";
					}else{
						linkEndColor = "etc";
						linkColor = "ee";
					}
				}
				
				ScoredPosition[] sps = p.getPairedScoredPositions();
				if(!shistMap.containsKey(sps[1])){
					shistMap.put(sps[1], new HashMap<String, Integer>());
				}
					
				if(!ehistMap.containsKey(sps[0])){
					ehistMap.put(sps[0], new HashMap<String, Integer>());
				}
				
				if(!shistMap.get(sps[1]).containsKey(linkStartColor))
					shistMap.get(sps[1]).put(linkStartColor, 0);
				if(!ehistMap.get(sps[0]).containsKey(linkEndColor))
					ehistMap.get(sps[0]).put(linkEndColor, 0);
					
				sMap.put(sps[0], linkStartColor);
				eMap.put(sps[1], linkEndColor);
				
				shistMap.get(sps[1]).put(linkStartColor,shistMap.get(sps[1]).get(linkStartColor) + 1);
				ehistMap.get(sps[0]).put(linkEndColor,ehistMap.get(sps[0]).get(linkEndColor) + 1);
				//colors.add(linkStartColor);
				//colors.add(linkEndColor);
				linkColors.add(linkColor);
			}
			
		//	HashMap<ScoredPosition, Integer> cummulativeS = new HashMap<ScoredPosition, Integer>();
		//	HashMap<ScoredPosition, Integer> cummulativeE = new HashMap<ScoredPosition, Integer>();
			
			PrintStream outs = new PrintStream(scatter3p);
			PrintStream oute = new PrintStream(scatter5p);
			
			PrintStream outsGN = new PrintStream(gene3p);
			PrintStream outeGN = new PrintStream(gene5p);
			
			for(ScoredPosition p : sMap.keySet()){
				StringBuilder sb = new StringBuilder();
				String contig = p.getContig();
				int pos5 = p.getFivePPosition();
				int pos3 = p.getThreePPosition();
				sb.append("hs");
				sb.append(contig.substring(contig.indexOf("chr") + 3));
				sb.append(' ');
				sb.append(Math.min(pos3, pos5));
				sb.append(' ');
				sb.append(Math.max(pos3, pos5));
				sb.append(' ');
				if(p.getContainingGeneNames() != null && !p.getContainingGeneNames().isEmpty()){
					outsGN.println(sb.toString() + p.getContainingGeneNames().get(0) + " color=" + sMap.get(p));
				}
				outs.println(sb.toString() + "0 color=" + sMap.get(p));
			}
			
			for(ScoredPosition p : eMap.keySet()){
				StringBuilder sb = new StringBuilder();
				String contig = p.getContig();
				int pos5 = p.getFivePPosition();
				int pos3 = p.getThreePPosition();
				sb.append("hs");
				sb.append(contig.substring(contig.indexOf("chr") + 3));
				sb.append(' ');
				sb.append(Math.min(pos3, pos5));
				sb.append(' ');
				sb.append(Math.max(pos3, pos5));
				sb.append(' ');
				if(p.getContainingGeneNames() != null && !p.getContainingGeneNames().isEmpty()){
					outeGN.println(sb.toString() + p.getContainingGeneNames().get(0) + " color=" + eMap.get(p));
				}
				oute.println(sb.toString() + "0 color=" + eMap.get(p));
			}
			
			outs.close();
			oute.close();
			outsGN.close();
			outeGN.close();
			
			int width = 4000000;
			PrintStream outsc = new PrintStream(histo5p);
			PrintStream outec = new PrintStream(histo3p);
			for(ScoredPosition p : shistMap.keySet()){
				StringBuilder sb = new StringBuilder();
				String contig = p.getContig();
				int pos5 = p.getFivePPosition();
				int pos3 = p.getThreePPosition();
				sb.append("hs");
				sb.append(contig.substring(contig.indexOf("chr") + 3));
				sb.append(' ');
				sb.append(Math.min(pos3, pos5) - width);
				sb.append(' ');
				sb.append(Math.max(pos3, pos5) + width);
				sb.append(' ');
				
				for(String color : colors){
					Integer n = shistMap.get(p).get(color);
					if(n == null) n = 0;
					sb.append(n+",");
				}
				
				outsc.println(sb.toString());	
			}
				
			for(ScoredPosition p : ehistMap.keySet()){
				StringBuilder sb = new StringBuilder();
				String contig = p.getContig();
				int pos5 = p.getFivePPosition();
				int pos3 = p.getThreePPosition();
				sb.append("hs");
				sb.append(contig.substring(contig.indexOf("chr") + 3));
				sb.append(' ');
				sb.append(Math.min(pos3, pos5) - width);
				sb.append(' ');
				sb.append(Math.max(pos3, pos5) + width);
				sb.append(' ');
				
				for(String color : colors){
					Integer n = ehistMap.get(p).get(color);
					if(n == null) n = 0;
					sb.append(n+",");
				}
			
				outec.println(sb.toString());	
			}
			outsc.close();
			outec.close();
			
			
			for(int n=0; n< parser.getPairs().size(); n++){
				ScoredPair p = parser.getPairs().get(n);
				ScoredPosition[] sps = p.getPairedScoredPositions();	
				for(int i = 0; i<sps.length;i++){
					StringBuilder sb = new StringBuilder();
					String contig = sps[i].getContig();
					int pos5 = sps[i].getFivePPosition();
					int pos3 = sps[i].getThreePPosition();
					sb.append("hs");
					sb.append(contig.substring(contig.indexOf("chr") + 3));
					sb.append(' ');
					sb.append(Math.min(pos3, pos5));
					sb.append(' ');
					sb.append(Math.max(pos3, pos5));
					sb.append(' ');
					linkOut.print(sb.toString());
					//	if(i == 0) linkStartOut.print(sb.toString());
				//	else linkEndOut.print(sb.toString());
				}
				linkOut.println("color="+linkColors.get(n));
			//	linkStartOut.println("0 fill_color="+linkStartColor + "_a2");
			//	linkEndOut.println("0.5 fill_color="+linkEndColor + "_a2");
				
			}
			
			linkOut.close();
		//	linkStartOut.close();
		//	linkEndOut.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static void generateConf(String confout, String confin, String dir, String  outImage, String link, String histo5p, String histo3p, String gene5p, String gene3p, String scatter5p, String scatter3p){
		try {
			BufferedLineReader in = new BufferedLineReader(confin);
			String s;
			PrintStream out = new PrintStream(confout);
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String r = s;
				if(s.contains("$DIR"))
					r = s.replace("$DIR", dir);
				if(s.contains("$OUTPUT"))
					r = s.replace("$OUTPUT", outImage);
				if(s.contains("$H5"))
					r = s.replace("$H5", histo5p);
				if(s.contains("$H3"))
					r = s.replace("$H3", histo3p);
				if(s.contains("$G5"))
					r = s.replace("$G5", gene5p);
				if(s.contains("$G3"))
					r = s.replace("$G3", gene3p);
				if(s.contains("$L"))
					r = s.replace("$L", link);
				if(s.contains("$S5"))
					r = s.replace("$S5", scatter5p);
				if(s.contains("$S3"))
					r = s.replace("$S3", scatter3p);
				
				
				out.println(r);
			}			
			out.close();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		
		String confout = args[0];
		String link = args[1];
		String histo5p = args[2];
		String histo3p = args[3];
		String gene5p = args[4];
		String gene3p = args[5];
		String scatter5p = args[6];
		String scatter3p = args[7];
		String outDir = args[8];
		String outImage = args[9];
		String csv = args[10];
		String confin = args[11];
		
		run(csv, link, histo5p, histo3p, gene5p, gene3p, scatter5p, scatter3p);
		generateConf(confout, confin, outDir, outImage, link, histo5p, histo3p, gene5p, gene3p, scatter5p, scatter3p);
	}
	
}
