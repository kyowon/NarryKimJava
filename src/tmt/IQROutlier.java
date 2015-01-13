package tmt;

public class IQROutlier{
	/*
	private int numChannel = 6;
	private MSGFPlusParser parser;
	private HashMap<HashSet<String>, ArrayList<MSGFPlusPSM>> groupedPSMs;
	private ArrayList<Integer> ionIntensityColumnNumbers; 
	private String outfile;
	
	public IQROutlier(String file, String outfile){
		parser = new MSGFPlusParser(file);
		this.outfile = outfile;
	}
	
	
	
    public static double quartile(double[] values, double lowerPercent) {

        if (values == null || values.length == 0) {
            throw new IllegalArgumentException("The data array either is null or does not contain any data.");
        }

        // Rank order the values
        double[] v = new double[values.length];
        System.arraycopy(values, 0, v, 0, values.length);
        Arrays.sort(v);

        int n = (int) Math.round(v.length * lowerPercent / 100);
        n = Math.min(n, v.length-1);
        return v[n];

    }
    
    private double[][] getIonIntensities(ArrayList<MSGFPlusPSM> psms){
    	double[][] intensities = new double[psms.size()][numChannel];
    	for(int i=0;i<psms.size();i++){
    		MSGFPlusPSM psm = psms.get(i);
    		String[] token = psm.toString().split("\t");
    		for(int j=0;j<ionIntensityColumnNumbers.size();j++){
    			intensities[i][j] = Double.parseDouble(token[ionIntensityColumnNumbers.get(j)]);
    		}
    	}
    	return intensities;
    }
    
    private double[][] getRatios(double[][] ionIntensities){
    	double[][] ratios = new double[ionIntensities.length][];
    	for(int i=0;i<ionIntensities.length;i++){
    		ratios[i] = new double[ionIntensities[i].length-1];
    		for(int j=1;j<ionIntensities[i].length;j++){
    			ratios[i][j-1] = ionIntensities[i][j]/ionIntensities[i][0];
    		}
    	}
    	return ratios;
    }
    
    private double getUpperOutlierThreshold(double[] tratio){
    	if(tratio.length < 2) return 1000000;
    	double upper = quartile(tratio, 75);
    	double lower = quartile(tratio, 25);
    	
    	return upper + (upper - lower) * 1.5;
    }
    
    private double getLowerOutlierThreshold(double[] tratio){
    	if(tratio.length < 2) return 0;
    	double upper = quartile(tratio, 75);
    	double lower = quartile(tratio, 25);
    	
    	return lower - (upper - lower) * 1.5;
    }
    
    private double[] getOutlierThresholds(double[][] ratios){
    	double[][] tratios = MC.transpose(ratios);
    	double[] thresholds = new double[tratios.length * 2];
    	for(int i=0;i<tratios.length;i++){
    		thresholds[i] = getUpperOutlierThreshold(tratios[i]);
    		thresholds[i + tratios.length] = getLowerOutlierThreshold(tratios[i]);
    	}
    	return thresholds;
    }
    
    private int[][] getOutlierIndices(double[][] ratios, double[] outliers){
    	int[][] indices = new int[ratios.length][];
    	for(int i=0;i<ratios.length;i++){
        	indices[i] = new int[ratios[i].length];
        	for(int j=0; j<ratios[i].length; j++){
        		if(ratios[i][j] <=0 || Double.isNaN(ratios[i][j]) || Double.isInfinite(ratios[i][j]) || Double.isInfinite(-ratios[i][j])) indices[i][j]=2;
        		else if(ratios[i][j] > outliers[j]) indices[i][j] = 1;
        		else if(ratios[i][j] < outliers[ratios[i].length + j]) indices[i][j] = -1;
        	}
    	}
    	return indices;
    }
    
    private double[] getOutlierRemovedRatio(double[][] ionIntensities, int[][] outlierIndices){
    	double[] ratio = new double[ionIntensities[0].length - 1];
    	double sum = 0;
    	for(int i=0;i<ionIntensities.length; i++){
    		boolean toSum = true;
    		for(int j : outlierIndices[i]){
    			if(j!=0){
    				toSum = false;
    				break;
    			}
    		}
    		if(!toSum) continue;
    		for(int j=0;j<ionIntensities[i].length;j++){
    			if(j == 0) sum += ionIntensities[i][j];
    			else ratio[j - 1] += ionIntensities[i][j];
    		}
    	}
    	
    	for(int i=0;i<ratio.length;i++){
    		ratio[i] /= sum;
    	}
    	
    	return ratio;
    }
    
    private void updateGroupedPSMs(){
    	groupedPSMs = new HashMap<HashSet<String>, ArrayList<MSGFPlusPSM>>();
    	ionIntensityColumnNumbers = new ArrayList<Integer>();
    	String[] headerToken = parser.getHeader().split("\t");
    	for(int i=0;i<headerToken.length;i++){
    		if(headerToken[i].startsWith("Ion_") && !headerToken[i].endsWith("MZ")){
    			ionIntensityColumnNumbers.add(i);
    		}
    	}
    	numChannel = ionIntensityColumnNumbers.size();
    	
    	for(MSGFPlusPSM psm : parser.getPsms()){
    		if(psm.getMiscString() == null || psm.getMiscString().isEmpty()) continue; 
    		HashSet<String> proteins = new HashSet<String>();
    		for(String p : psm.getProteins()) proteins.add(p);
    		if(!groupedPSMs.containsKey(proteins)) groupedPSMs.put(proteins, new ArrayList<MSGFPlusParser.MSGFPlusPSM>());
    		groupedPSMs.get(proteins).add(psm);
    	}
    }
    
    public void run(){
    	PrintStream out;
		try {
			out = new PrintStream(outfile);
			updateGroupedPSMs();
	    	out.println(getHeader());
	    	for(ArrayList<MSGFPlusPSM> psms : groupedPSMs.values()){
	    		out.print(getString(psms));
	    	}    	
	    	out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    private String getHeader(){
    	String header = parser.getHeader() + "\t";
    	for(int i=0;i<numChannel-1;i++){
    		header += "R_" + (i+2) + "/R_1\t";
    	}
    	for(int i=0;i<numChannel-1;i++){
    		header += "Outlier_Index_R_" + (i+2) + "/R_1\t";
    	}
    	for(int i=0;i<numChannel-1;i++){
    		header += "UpperThreshold_R_" + (i+2) + "/R_1\t";
    	}
    	for(int i=0;i<numChannel-1;i++){
    		header += "LowerThreshold_R_" + (i+2) + "/R_1\t";
    	}
    	for(int i=0;i<numChannel-1;i++){
    		header += "Outlier_Removed_Ratio" + (i+2) + "/R_1";
    		if(i<numChannel-1) header +="\t";
    	}
    	return header;
    }
    
    private String getString(ArrayList<MSGFPlusPSM> psms){
    	double[][] ionIntensities = getIonIntensities(psms);
    	double[][] ratios = getRatios(ionIntensities);
    	double[] outlierThresholds = getOutlierThresholds(ratios);
    	int[][] outlierIndices = getOutlierIndices(ratios, outlierThresholds);
    	double[] ratio = getOutlierRemovedRatio(ionIntensities, outlierIndices);
    	StringBuilder sb = new StringBuilder();
    	for(int i=0; i<psms.size();i++){
    		MSGFPlusPSM psm = psms.get(i);
    		sb.append(psm); sb.append('\t');
    		for(int j=0; j<ratios[i].length;j++){
    			sb.append(ratios[i][j]); sb.append('\t');
    		}
    		for(int j=0; j<outlierIndices[i].length;j++){
    			sb.append(outlierIndices[i][j]); sb.append('\t');
    		}
			for(int j=0; j<outlierThresholds.length;j++){
    			if(i==0) sb.append(outlierThresholds[j]); 
    			sb.append('\t');
    		}
			for(int j=0; j<ratio.length;j++){
    			if(i==0)sb.append(ratio[j]); 
    			sb.append(j<ratio.length - 1 ? '\t' : '\n');
    		}
    	}
    	return sb.toString();
    }
    
    static public void main(String[] args){
    	if(args == null || args.length < 2){
    		System.out.println("java -jar IQROutlier.jar [input.tsv] [output.tsv]");
    		System.exit(0);
    	}
    	IQROutlier test = new IQROutlier(args[0], args[1]);
    	test.run();
    }*/
    
}
