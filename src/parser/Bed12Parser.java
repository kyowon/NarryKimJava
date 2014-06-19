package parser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import parser.AnnotationFileParser.AnnotatedGene;

public class Bed12Parser {
	
	final int maxReadLength = 300;
	
	public static class Bed12Read{
		private static boolean read3p = false;
		//private String contig;
		private int[] fivePs;// inclusive
		private int[] threePs; // inclusive
		private short[] depth = null;
		
		//private boolean isPlusStrand;
		
		public Bed12Read(String s){
			String[] token = s.split("\t");
			//this.contig = token[0];
			int start = Integer.parseInt(token[1]);
			boolean isPlusStrand = token[5].equals("+");
			
			String[] lenString = token[10].split(",");
			int[] lens = new int[lenString.length]; 
			for(int i=0;i<lenString.length;i++){
				lens[i] = Integer.parseInt(lenString[i]);
			}
			String[] startString = token[11].split(",");
			int[] starts = new int[startString.length];
			for(int i=0;i<startString.length;i++){
				starts[i] = Integer.parseInt(startString[i]);
			}
			fivePs = new int[starts.length];
			threePs = new int[starts.length];
			if(isPlusStrand){
				for(int i=0;i<starts.length;i++){
					fivePs[i] = start + starts[i];
					threePs[i] = fivePs[i] + lens[i] - 1;
				}
			}else{
				for(int i=0;i<starts.length;i++){
					fivePs[fivePs.length - i - 1] = start + starts[i] + lens[i] - 1;
					threePs[threePs.length - i - 1] = fivePs[fivePs.length - i - 1] - lens[i] + 1;
				}
			}			
		}
		
		public boolean isSpliced(){
			return fivePs.length > 1;
		}
		
		public short get5pDepth() {return depth == null? 1 : depth[0];}
		public short get3pDepth() {return (depth == null || depth.length < 2) ? 1: depth[1];}
		
		public int get5p(){
			return fivePs[0];
		}
		
		public int get3p(){
			return threePs[threePs.length - 1];
		}

		public int[] get5ps(){
			return fivePs;
		}
		
		public int[] get3ps(){
			return threePs;
		}
		
		private void initDepth(){
			if(depth != null) return;
			depth = new short[read3p? 2 : 1];
			for(int i=0;i<depth.length;i++){
				depth[i] = 1;
			}
		}
		
		public void set5pDepth(short i){
			initDepth();
			depth[0] = i;
		}
		
		public void set3pDepth(short i){
			initDepth();
			depth[1] = i;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Arrays.hashCode(fivePs);
			result = prime * result + Arrays.hashCode(threePs);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Bed12Read other = (Bed12Read) obj;
			if (!Arrays.equals(fivePs, other.fivePs))
				return false;
			if (!Arrays.equals(threePs, other.threePs))
				return false;
			return true;
		}

		public boolean is5pcontained(ArrayList<Integer> coordinates, ArrayList<ArrayList<Integer>> splices, boolean isPlusStrand){
			if(isPlusStrand){
				if(fivePs[0] < coordinates.get(0)) return false;		// || threePs[threePs.length - 1] > coordinates.get(coordinates.size()-1)			
			}else{
				if(fivePs[0] > coordinates.get(0)) return false; // || threePs[threePs.length - 1] < coordinates.get(coordinates.size()-1)
			}
		
			for(ArrayList<Integer> splice : splices){
			//	if(isPlusStrand){
				if(fivePs[0] >= splice.get(0) && fivePs[0] <= splice.get(1)|| threePs[threePs.length - 1] >= splice.get(0) && threePs[threePs.length - 1] <= splice.get(1)) return false;					
			//	}else{
			//		if(fivePs[0] <= splice.get(1) || threePs[threePs.length - 1] >= splice.get(0)) return false;
			//	}
			}				
			if(isSpliced()){
				HashSet<ArrayList<Integer>> rs = new HashSet<ArrayList<Integer>>();
				for(int i=0;i<fivePs.length-1;i++){
					ArrayList<Integer> s = new ArrayList<Integer>();
					if(isPlusStrand){
						s.add(threePs[i] + 1);
						s.add(fivePs[i+1] - 1);
					}else{
						s.add(fivePs[i+1] + 1);
						s.add(threePs[i] - 1);
					}
					rs.add(s);
				}
				//System.out.println(splices + "\n*" + rs + " " + splices.containsAll(rs));
				return splices.containsAll(rs);
			
			}
			return true;
		}
		
		public boolean is3pcontained(ArrayList<Integer> coordinates, ArrayList<ArrayList<Integer>> splices, boolean isPlusStrand){
			if(isPlusStrand){
				if(threePs[0] < coordinates.get(0)) return false;		// || threePs[threePs.length - 1] > coordinates.get(coordinates.size()-1)			
			}else{
				if(threePs[0] > coordinates.get(0)) return false; // || threePs[threePs.length - 1] < coordinates.get(coordinates.size()-1)
			}
			
			for(ArrayList<Integer> splice : splices){
				if(fivePs[0] >= splice.get(0) && fivePs[0] <= splice.get(1)|| 
						threePs[threePs.length - 1] >= splice.get(0) && threePs[threePs.length - 1] <= splice.get(1)) return false;					
			}				
			if(isSpliced()){
				HashSet<ArrayList<Integer>> rs = new HashSet<ArrayList<Integer>>();
				for(int i=0;i<fivePs.length-1;i++){
					ArrayList<Integer> s = new ArrayList<Integer>();
					if(isPlusStrand){
						s.add(threePs[i] + 1);
						s.add(fivePs[i+1] - 1);
					}else{
						s.add(fivePs[i+1] + 1);
						s.add(threePs[i] - 1);
					}
					rs.add(s);
				}
				//System.out.println(splices + "\n*" + rs + " " + splices.containsAll(rs));
				return splices.containsAll(rs);
			}
			
			return true;
		}
		
		private boolean containing(int position, boolean isPlusStrand){
			if(isPlusStrand){
				for(int i=0;i<fivePs.length;i++){
					if(fivePs[i] <= position && threePs[i] >= position) return true;
				}
			}else{
				for(int i=0;i<fivePs.length;i++){
					if(fivePs[i] >= position && threePs[i] <= position) return true;
				}
			}
			
			return false;
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			//sb.append(contig);sb.append('\t');
			//sb.append(isPlusStrand? '+': '-');sb.append('\t');
			//if(isPlusStrand){
				for(int i=0;i<fivePs.length-1;i++){
					sb.append(fivePs[i]);sb.append(',');
				}	
				sb.append(fivePs[fivePs.length-1]);
				sb.append('\t');
				for(int i=0;i<threePs.length-1;i++){
					sb.append(threePs[i]);sb.append(',');
				}	
				sb.append(threePs[threePs.length-1]);
			/*}else{
				for(int i=fivePs.length-1;i>0;i--){
					sb.append(fivePs[i]);sb.append(',');
				}	
				sb.append(fivePs[0]);
				sb.append('\t');
				for(int i=threePs.length-1;i>0;i--){
					sb.append(threePs[i]);sb.append(',');
				}	
				sb.append(threePs[0]);
			}*/
			return sb.toString();
			
			
		}
		
	}

	private HashMap<Boolean, HashMap<Integer, ArrayList<Bed12Read>>> read5pMap;
	private HashMap<Boolean, HashMap<Integer, ArrayList<Bed12Read>>> read3pMap;
	private HashMap<Boolean, HashMap<Integer, Short>> readDepthMap;
	private HashMap<Boolean, HashMap<Integer, HashSet<Integer>>> splices3pDirection; // inclusive
	private HashMap<Boolean, HashMap<Integer, HashSet<Integer>>> splices5pDirection; // inclusive
	private int totalReadCount = 0;
	
	/**
	 * @return the totalReadCount
	 */
	public int getTotalReadCount() {
		return totalReadCount;
	}
	
	private String filename;
	private String contig;
	
	public String getFileName() {return filename; }
	public String getContig() { return contig; }
	
	public Bed12Parser(String filename, AnnotationFileParser annotationParser, String contig){
		this(filename, annotationParser, contig, false);
	}
	
	public Bed12Parser(String filename, AnnotationFileParser annotationParser, String contig, boolean read3p){
		Bed12Read.read3p = read3p;
		this.filename = filename;
		this.contig = contig;
		read5pMap = new HashMap<Boolean, HashMap<Integer,ArrayList<Bed12Read>>>(); 
		if(read3p){ 
			read3pMap = new HashMap<Boolean, HashMap<Integer,ArrayList<Bed12Read>>>();
			readDepthMap = new HashMap<Boolean, HashMap<Integer,Short>>();
		}
		splices3pDirection = new HashMap<Boolean, HashMap<Integer,HashSet<Integer>>>();
		splices5pDirection = new HashMap<Boolean, HashMap<Integer,HashSet<Integer>>>();
		//HashMap<String, HashMap<Boolean, ArrayList<Bed12Read>>> readList = new HashMap<String, HashMap<Boolean,ArrayList<Bed12Read>>>(); 
		try {
			BufferedLineReader in = new BufferedLineReader(filename);
			String s;
			while((s=in.readLine())!=null){
				totalReadCount++;
				String[] token = s.split("\t");
				if(!token[0].equals(this.contig)) continue;
				Bed12Read read = new Bed12Read(s);
				
				boolean isPlusStrand = token[5].equals("+");//l//read.isPlusStrand();
				
			//	HashMap<Boolean, ArrayList<Bed12Read>> sreadList = readList.get(contig);
				
				if(!read5pMap.containsKey(isPlusStrand)){
					read5pMap.put(isPlusStrand, new HashMap<Integer,ArrayList<Bed12Read>>());
					if(read3p){
						read3pMap.put(isPlusStrand, new HashMap<Integer,ArrayList<Bed12Read>>());
						readDepthMap.put(isPlusStrand, new HashMap<Integer, Short>());
					}
					splices3pDirection.put(isPlusStrand, new HashMap<Integer, HashSet<Integer>>());
					splices5pDirection.put(isPlusStrand, new HashMap<Integer, HashSet<Integer>>());
			//		sreadList.put(isPlusStrand, new ArrayList<Bed12Read>());
				}
				int position = read.get5p();
			//	ArrayList<Bed12Read> ssreadList = sreadList.get(isPlusStrand);
				HashMap<Integer,ArrayList<Bed12Read>> sreads = read5pMap.get(isPlusStrand);
				if(!sreads.containsKey(position)) sreads.put(position,new ArrayList<Bed12Read>());
				ArrayList<Bed12Read> ssreads = sreads.get(position);

				Bed12Read lastRead = null;
				if(ssreads.size() > 0) lastRead = ssreads.get(ssreads.size()-1);
				if(lastRead != null &&read.equals(lastRead)){
					lastRead.set5pDepth((short)(lastRead.get5pDepth()+1));
				}else
					ssreads.add(read);
				
				if(read3p){
					int position3p = read.get3p();
					//	ArrayList<Bed12Read> ssreadList = sreadList.get(isPlusStrand);
					HashMap<Integer,ArrayList<Bed12Read>> sreads3p = read3pMap.get(isPlusStrand);
					if(!sreads3p.containsKey(position3p)) sreads3p.put(position3p, new ArrayList<Bed12Read>());
					ArrayList<Bed12Read> ssreads3p = sreads3p.get(position3p);

					Bed12Read lastRead3p = null;
					if(ssreads3p.size() > 0) lastRead3p = ssreads3p.get(ssreads3p.size()-1);
					if(lastRead3p != null && read.equals(lastRead3p)){
						lastRead3p.set3pDepth((short)(lastRead3p.get3pDepth()+1));
					}else
						ssreads3p.add(read);
					
				//	if(position3p == 133680428-5) System.out.println(ssreads3p.size());
					
					HashMap<Integer, Short> sreadDepth = readDepthMap.get(isPlusStrand);
					if(isPlusStrand){
						for(int p=read.get5p();p<=read.get3p();p++){
							Short  n = sreadDepth.get(p);
							n = n==null? 0 : n;
							sreadDepth.put(p, (short)(n + 1));
						}
					}else{
						for(int p=read.get3p();p<=read.get5p();p++){
							Short  n = sreadDepth.get(p);
							n = n==null? 0 : n;
							sreadDepth.put(p, (short)(n + 1));
						}
					}
				}
				
				HashMap<Integer,HashSet<Integer>> s3 = splices3pDirection.get(isPlusStrand);
				HashMap<Integer,HashSet<Integer>> s5 = splices5pDirection.get(isPlusStrand);
				
				//if(read.getSplices() != null){
				if(read.isSpliced()){
					for(int i=0;i<read.get5ps().length-1;i++){
						int fiveP = read.get5ps()[i+1];
						int threeP = read.get3ps()[i];
						if(!s5.containsKey(fiveP)) s5.put(fiveP, new HashSet<Integer>());
						s5.get(fiveP).add(threeP);
						if(!s3.containsKey(threeP)) s3.put(threeP, new HashSet<Integer>());
						s3.get(threeP).add(fiveP);
					}
				}
			}
			in.close();
		
			if(annotationParser != null){
				Iterator<AnnotatedGene> geneIterator = annotationParser.getAnnotatedGeneIterator(contig);
				while(geneIterator.hasNext()){
					AnnotatedGene gene = geneIterator.next();
					HashMap<Integer,HashSet<Integer>> ss3 = splices3pDirection.get(gene.isPlusStrand());
					HashMap<Integer,HashSet<Integer>> ss5 = splices5pDirection.get(gene.isPlusStrand());
					if(ss3 != null){
						if(gene.isPlusStrand()){
							for(int[] intron : gene.getIntrons()){
								if(!ss3.containsKey(intron[0]-1)) ss3.put(intron[0]-1, new HashSet<Integer>());
								ss3.get(intron[0]-1).add(intron[1]+1);
								if(!ss5.containsKey(intron[1]+1)) ss5.put(intron[1]+1, new HashSet<Integer>());
								ss5.get(intron[1]+1).add(intron[0]-1);
							}
						}else{
							for(int[] intron : gene.getIntrons()){
								if(!ss5.containsKey(intron[0]-1)) ss5.put(intron[0]-1, new HashSet<Integer>());
								ss5.get(intron[0]-1).add(intron[1]+1);
								if(!ss3.containsKey(intron[1]+1)) ss3.put(intron[1]+1, new HashSet<Integer>());
								ss3.get(intron[1]+1).add(intron[0]-1);
							}
						}
					}
					
				}
			}
			for(boolean isPlusStrand : splices3pDirection.keySet()){
				HashMap<Integer,ArrayList<Bed12Read>> ssreads = read5pMap.get(isPlusStrand);
				HashMap<Integer,HashSet<Integer>> ss3 = splices3pDirection.get(isPlusStrand);
				for(int position : ss3.keySet()){
					if(isPlusStrand){
						for(int offset = 0; offset >= -maxReadLength ; offset --){
							int p = position + offset;
							ArrayList<Bed12Read> sssreads = ssreads.get(p);
							if(sssreads == null) continue;
							for(Bed12Read read : sssreads){
								if(read.containing(position + 1, isPlusStrand)){
									ss3.get(position).add(position + 1);
									//System.out.println(contig + " " + position);
									break;
								}
							}
						}	
					}else{
						for(int offset = 0; offset <= maxReadLength ; offset ++){
							int p = position + offset;
							ArrayList<Bed12Read> sssreads = ssreads.get(p);
							if(sssreads == null) continue;
							for(Bed12Read read : sssreads){
								if(read.containing(position - 1, isPlusStrand)){
									ss3.get(position).add(position - 1);
									//System.out.println(contig + " " + position);
									break;
								}
							}
						}
					}
					
				}
			}
		
			for(boolean isPlusStrand : splices5pDirection.keySet()){
				HashMap<Integer,ArrayList<Bed12Read>> ssreads = read5pMap.get(isPlusStrand);
				HashMap<Integer,HashSet<Integer>> ss5 = splices5pDirection.get(isPlusStrand);
				for(int position : ss5.keySet()){
					if(isPlusStrand){
						for(int offset = 0; offset >= -maxReadLength ; offset --){
							int p = position + offset;
							ArrayList<Bed12Read> sssreads = ssreads.get(p);
							if(sssreads == null) continue;
							for(Bed12Read read : sssreads){
								if(read.containing(position - 1, isPlusStrand)){
									ss5.get(position).add(position - 1);
									//System.out.println(contig + " " + position);
									break;
								}
							}
						}	
					}else{
						for(int offset = 0; offset <= maxReadLength ; offset ++){
							int p = position + offset;
							ArrayList<Bed12Read> sssreads = ssreads.get(p);
							if(sssreads == null) continue;
							for(Bed12Read read : sssreads){
								if(read.containing(position + 1, isPlusStrand)){
									ss5.get(position).add(position + 1);
									//System.out.println(contig + " " + position);
									break;
								}
							}
						}
					}
					
				}
			
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Iterator<Integer> getNonZero3pPositionIterator(boolean isPlusStrand){
		ArrayList<Integer> positions = new ArrayList<Integer>();
		if(read3pMap == null || !read3pMap.containsKey(isPlusStrand)) return positions.iterator();
		positions.addAll(read3pMap.get(isPlusStrand).keySet());
		Collections.sort(positions);
		return positions.iterator();
	}
	
	public Iterator<Integer> getNonZero5pPositionIterator(boolean isPlusStrand){
		ArrayList<Integer> positions = new ArrayList<Integer>();
		if(!read5pMap.containsKey(isPlusStrand)) return positions.iterator();
		positions.addAll(read5pMap.get(isPlusStrand).keySet());
		Collections.sort(positions);
		return positions.iterator();
		//HashMap<String, HashMap<Boolean, HashMap<Integer, HashMap<Bed12Read, Short>>>> reads;
	}
	
	public ArrayList<ArrayList<Integer>> getCoordinates(int position, int length, boolean isPlusStrand, boolean is3pDirection){
		boolean sign = isPlusStrand ? is3pDirection : !is3pDirection;
		//System.out.println(sign);
		HashMap<Boolean, HashMap<Integer,HashSet<Integer>>> splices = is3pDirection? splices3pDirection : splices5pDirection;
		HashMap<Integer,HashSet<Integer>> ss = splices.get(isPlusStrand);
		
		ArrayList<ArrayList<Integer>> init = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> sinit = new ArrayList<Integer>();
		if(length > 0 )
			sinit.add(position);
		init.add(sinit);
		ArrayList<ArrayList<Integer>> coordinates = updateCoordinates(init, ss, length - 1, sign);
		if(!is3pDirection){
			for(ArrayList<Integer> coordinate : coordinates){
				Collections.reverse(coordinate);
			}
		}
		return coordinates;
	}
	
	public ArrayList<ArrayList<Integer>> getCoordinates(int position, int len5pDirection, int len3pDirection, boolean isPlusStrand){
		ArrayList<ArrayList<Integer>> threePCoordinates = getCoordinates(position, len3pDirection, isPlusStrand, true); 
		ArrayList<ArrayList<Integer>> fivePCoordinates = getCoordinates(position, len5pDirection + 1, isPlusStrand, false);
				
		ArrayList<ArrayList<Integer>> coordinates = new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> threePCoordinate : threePCoordinates){
			for(ArrayList<Integer> fivePCoordinate : fivePCoordinates){
				ArrayList<Integer> coordinate = new ArrayList<Integer>(fivePCoordinate);
				coordinate.remove(coordinate.size()-1);
				coordinate.addAll(threePCoordinate);
				coordinates.add(coordinate);
			}
		}
		
		return coordinates;
	}
	
	// get splice such that int[0] < int[1]
	static public ArrayList<ArrayList<Integer>> getSplices(ArrayList<Integer> coordinate){
		ArrayList<ArrayList<Integer>> splices = new ArrayList<ArrayList<Integer>>();
		for(int i=0;i<coordinate.size()-1;i++){
			int v1 = coordinate.get(i);
			int v2 = coordinate.get(i+1);
			if(Math.abs(v1 - v2) > 1){
				ArrayList<Integer> s = new ArrayList<Integer>();
				s.add(Math.min(v1, v2) + 1);
				s.add(Math.max(v1, v2) - 1);
				splices.add(s);
			}
		}
		return splices;
	}
	
	
	// from coordiante to coordinate start/end..
	static public ArrayList<ArrayList<Integer>> getCoordinateStartsEnds(ArrayList<Integer> coordinate, boolean isPlusStrand){
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();	
		//int end = 0;
		if(isPlusStrand){
			int start = coordinate.get(0);
			for(int i=0;i<coordinate.size()-1;i++){
				int v1 = coordinate.get(i);
				int v2 = coordinate.get(i+1);
				if(Math.abs(v1 - v2) > 1){
					ArrayList<Integer> s = new ArrayList<Integer>();
					
					s.add(start);
					s.add(v1+1);
			//	end = v1;
					start = v2;
					ret.add(s);
				}
			}
			
			ArrayList<Integer> s = new ArrayList<Integer>();
			s.add(start);
			s.add(coordinate.get(coordinate.size()-1)+1);
		
			ret.add(s);
		}else{
			int start = coordinate.get(coordinate.size()-1);
			for(int i=coordinate.size()-1;i>0;i--){
				int v1 = coordinate.get(i);
				int v2 = coordinate.get(i-1);
				if(Math.abs(v1 - v2) > 1){
					ArrayList<Integer> s = new ArrayList<Integer>();
					
					s.add(start);
					s.add(v1+1);
			//	end = v1;
					start = v2;
					ret.add(s);
				}
			}
			
			ArrayList<Integer> s = new ArrayList<Integer>();
			s.add(start);
			s.add(coordinate.get(0)+1);
		
			ret.add(s);
			
		}
		return ret;
	}
	
	static public ArrayList<Integer> getCoordinate(ArrayList<ArrayList<Integer>> coordinateStartsEnds, boolean isPlusStrand){
		ArrayList<Integer> coordinate = new ArrayList<Integer>();
		for(ArrayList<Integer> se : coordinateStartsEnds){
			for(int i=se.get(0);i<se.get(1);i++){
				coordinate.add(i);
			}
		}
		if(!isPlusStrand) Collections.sort(coordinate, Collections.reverseOrder());
		return coordinate;
	}
	
	
	public double[] get3pSignalForfCLIP(boolean isPlusStrand, ArrayList<Integer> coordinate){
		double[] covs5p  = get5pCoverages(isPlusStrand, coordinate);
		double[] covs3p = get3pCoverages(isPlusStrand, coordinate);
		double[] depths = getDepths(isPlusStrand, coordinate);
		double[] signal = new double[depths.length];
		
		for(int i=0;i<signal.length;i++){
		//	if(i%2 == 0){
				signal[i] = depths[i] -  (covs5p[i] / (Math.max(depths[i] - covs5p[i] + 1, 1)) * (i>0 ? covs3p[i-1] : 0));
				signal[i] = Math.sqrt(Math.abs(signal[i])) * Math.signum(signal[i]);
		//	}else{
			//	signal[i] = depths[i/2];// +  covs5p[i/2] ;
		//	}
		}
		return signal;
		
	}
	
	
	public double[] get5pSignalForfCLIP(boolean isPlusStrand, ArrayList<Integer> coordinate){
		double[] covs5p  = get5pCoverages(isPlusStrand, coordinate);
		double[] covs3p = get3pCoverages(isPlusStrand, coordinate);
		double[] depths = getDepths(isPlusStrand, coordinate);
		double[] signal = new double[depths.length];
		
		for(int i=0;i<signal.length;i++){
		//	if(i%2 == 0){
				signal[i] = depths[i] -  (covs3p[i] / (Math.max(depths[i] - covs3p[i] + 1, 1)) * (i<covs5p.length-1 ? covs5p[i+1] : 0));
				signal[i] = Math.sqrt(Math.abs(signal[i])) * Math.signum(signal[i]);
		//	}else{
			//	signal[i] = depths[i/2];// +  covs5p[i/2] ;
		//	}
		}
		return signal;
		
	}
	
	//HashMap<String, HashMap<Boolean, HashMap<Integer, HashMap<Bed12Read, Short>>>> reads;
	public double[] get5pCoverages(boolean isPlusStrand, ArrayList<Integer> coordinate){
		double[] covs = new double[coordinate.size()];	
		HashMap<Integer, ArrayList<Bed12Read>> ssreads = read5pMap.get(isPlusStrand);
		if(ssreads == null) return covs;
		
		//HashSet<Integer> coordinateSet = new HashSet<Integer>(coordinates);
		ArrayList<ArrayList<Integer>> splices = getSplices(coordinate);
		for(int i =0; i<coordinate.size(); i++){
			int position = coordinate.get(i);
			ArrayList<Bed12Read> c = ssreads.get(position);
			if(c == null) continue;
			int j=0;
			for(Bed12Read read : c){
				if(read.is5pcontained(coordinate, splices, isPlusStrand)){					
					j+=read.get5pDepth();
				}//else System.out.println(splices + " " + read + " " + coordinates.get(0) + " " + coordinates.get(coordinates.size()-1));
			}
			covs[i] = j;
		}
		
		return covs;
	}
	
	// do not consider mapped reads to the coordinate...
	public double[] getDepths(boolean isPlusStrand, ArrayList<Integer> coordinate){
		double[] covs = new double[coordinate.size()];	
		HashMap<Integer, Short> sreadDepth = readDepthMap.get(isPlusStrand);
		if(sreadDepth == null) return covs;
		for(int i =0; i<coordinate.size(); i++){
			int position = coordinate.get(i);
			Short s = sreadDepth.get(position);
			if(s == null) s = 0;
			
			covs[i] = s;
		}
		
		return covs;
		//sreadDepth
	}
	
	public double[] get3pCoverages(boolean isPlusStrand, ArrayList<Integer> coordinate){
		if(read3pMap == null) return null;
		double[] covs = new double[coordinate.size()];	
		HashMap<Integer, ArrayList<Bed12Read>> ssreads = read3pMap.get(isPlusStrand);
		if(ssreads == null) return covs;
		
		ArrayList<ArrayList<Integer>> splices = getSplices(coordinate);
		for(int i =0; i<coordinate.size(); i++){
			int position = coordinate.get(i);
			
			
			ArrayList<Bed12Read> c = ssreads.get(position);
			
			if(c == null) continue;
			int j=0;
		//	if(position == 133680428 - 5) System.out.println(c.size());
			for(Bed12Read read : c){
				if(read.is3pcontained(coordinate, splices, isPlusStrand)){					
					j+=read.get3pDepth();
				//	System.out.println(read.get3pDepth());
				}//else System.out.println(splices + " " + read + " " + coordinates.get(0) + " " + coordinates.get(coordinates.size()-1));
			}
			covs[i] = j;
		}
		
		return covs;
	}
	
	
	public double[] getCDSCoverages(AnnotatedGene gene){
		if(!contig.equals(getContig())) return null;
		ArrayList<Integer> coordinate = gene.getLiftOverPositions(gene.isPlusStrand() ? gene.getCdsStart() : gene.getCdsEnd() - 1, 0, Integer.MAX_VALUE, true);
		return get5pCoverages(gene.isPlusStrand(), coordinate);
	}
	
	private ArrayList<ArrayList<Integer>> updateCoordinates(ArrayList<ArrayList<Integer>> coordinates, HashMap<Integer,HashSet<Integer>> ss, int len, boolean sign){
		if(len <= 0) return coordinates;		
		ArrayList<ArrayList<Integer>> newCoordinates = new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> coordinate : coordinates){
			int lp = coordinate.get(coordinate.size()-1);
			HashSet<Integer> ps = ss.get(lp);
			if(ps == null){
				ArrayList<Integer> newCoordinate = new ArrayList<Integer>(coordinate);
				int np = lp + (sign? 1 : -1);
				newCoordinate.add(np);
				newCoordinates.add(newCoordinate);
			}else{
				for(int np : ps){
					//System.out.println(lp + " " + ps);
					ArrayList<Integer> newCoordinate = new ArrayList<Integer>(coordinate);
					newCoordinate.add(np);
					newCoordinates.add(newCoordinate);
				}
			}
		}
		return updateCoordinates(newCoordinates, ss, len-1, sign);
	}
	
	static public void main(String[] args){//565858,413662
		//Bed12Read read = new Bed12Read("chr17	15142909	15165755	0_349741-3	255	-	15142909	15165755	255,0,0	2	19,10	0,22836");
		//System.out.println(read);
		
		//
		Bed12Parser bedParser = new Bed12Parser("/media/kyowon/Data1/fCLIP/samples/sample2/bed/x2.sorted.bed", new AnnotationFileParser("/media/kyowon/Data1/fCLIP/genomes/hg19.refFlat.txt"),
				"chrX", true);
			
		// chrX: 133680358-133680428 [
		ArrayList<Integer> coordinate = new ArrayList<Integer>();
		for(int i= 133680428; i>=133680358; i--) coordinate.add(i);
		
		
		
		for(double d : bedParser.get5pSignalForfCLIP(false, coordinate))
			System.out.println(d);
//		Nkx6-2	NM_183248	chr7	-	146767118	146768696	146767332	146768438	3	146767118,146767779,146768032,	146767587,146767952,146768696,
//		Zfp78	NM_001112805	chr7	+	6316014	6335315	6324070	6332201	4	6316014,6324012,6325833,6330918,	6316054,6324103,6325960,6335315,
		
		//Bed12Parser test = new Bed12Parser("/media/kyowon/Data1/fCLIP/samples/sample1/bed/Drosha2.sorted.bed", parser, "chr5", true);
	
	//	Collections.reverse(coordinate);
	}
		
	
}
