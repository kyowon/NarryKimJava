package tmt;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.util.ArithmeticUtils;

import msutil.Composition;
import msutil.Peptide;

public class IsotopeDistributionGenerator {
	 
	 
	//private static final double[] ProbC13 = {1, 0, 0, 0};
	//private static final double[] ProbN15 = {1, 0, 0, 0};
	private static final double[] ProbC = { .9893, 0.0107, 0, 0 };
    private static final double[] ProbH = { .999885, .000115, 0, 0 };
    private static final double[] ProbN = { 0.99632, 0.00368, 0, 0 };
    private static final double[] ProbO = { 0.99757, 0.00038, 0.00205, 0 };
    private static final double[] ProbS = { 0.9493, 0.0076, 0.0429, 0.0002 };
    private static int[][][] _possibleIsotopeCombinations;

     private double[] dist ;
     private int mostIntenseIsotopomerIndex;
     private float monoIsotopeMass; 
     
     private static void ComputePossibleIsotopeCombinations(int max) // called just once. 
     {   
    	 int maxIsotopeNumberInElement = ProbC.length - 1;
    	 boolean[] node = new boolean[max + 1];
         boolean[][] edge = new boolean[max + 1][maxIsotopeNumberInElement + 1];
    	 node[0] = true;
         
         for(int n=1;n<=max;n++){
        	 for(int j=1;j<=maxIsotopeNumberInElement;j++){
        		 int index = n-j;
        		 if(index < 0) continue;
        		 if(node[index]){
        			 node[n] = true;
        			 edge[n][j] = true;
        		 }
        	 }
         }
         
         ArrayList<ArrayList<ArrayList<Integer>>> possibleIsotopeCombinations = new ArrayList<ArrayList<ArrayList<Integer>>>();
         ArrayList<ArrayList<Integer>> initcomb = new ArrayList<ArrayList<Integer>>();
         initcomb.add(new ArrayList<Integer>());
         possibleIsotopeCombinations.add(initcomb);
            //Arrays.equals(ary1,ary2)      
         for (int i = 1; i <= max ; i++)
         {
        	ArrayList<ArrayList<Integer>> comb = new ArrayList<ArrayList<Integer>>();
        	for(int j=1; j<edge[i].length;j++){
        		boolean e = edge[i][j];
        		if(!e) continue;
        		ArrayList<ArrayList<Integer>> pc =  possibleIsotopeCombinations.get(i-j);
        		for(ArrayList<Integer> path : pc){
        			ArrayList<Integer> newPath = new ArrayList<Integer>(path);
        			newPath.add(j);
        			Collections.sort(newPath);
        			if(!comb.contains(newPath)) comb.add(newPath);
        		}
        	}
        	possibleIsotopeCombinations.add(comb);
        	//System.out.println(comb);
         }
         
         _possibleIsotopeCombinations = new int[max + 1][][];
         for(int i=0;i<=max;i++){
        	 _possibleIsotopeCombinations[i] = new int[possibleIsotopeCombinations.get(i).size()][maxIsotopeNumberInElement];
        	 for(int j=0;j<_possibleIsotopeCombinations[i].length;j++){
        		 for(int k : possibleIsotopeCombinations.get(i).get(j)){
        			 _possibleIsotopeCombinations[i][j][k-1]++;
        		 }
        		 //for(int k=0;k<maxIsotopeNumberInElement;k++)
        		 //    System.out.print(_possibleIsotopeCombinations[i][j][k] + " ");
        		 //System.out.print(" ; ");
        	 }
        	 //System.out.println();
         }
       //  _possibleIsotopeCombinations = possibleIsotopeCombinations;
     }
     
     
     private static double GetIsotopeProbability(int[] number, double[] means, double[] exps)
     {
         double prob = 1.0;
         for (int i = 1; i <= Math.min(ProbC.length - 1, number.length); i++)
         {
             double mean = means[i];
             double exp = exps[i];
             if (number[i - 1] == 0) prob *= exp;
             else
                 prob *=
                     (Math.pow(mean, number[i - 1]) * exp /  ArithmeticUtils.factorial(number[i - 1]));
         }
         return prob;
     }
     
   //  public static int MaxNumIsotopes = 100;
     //public static double IsotopeRelativeIntensityThreshold = 0.0;
     
     public double[] GetIsotopomerEnvelop(int c, int h, int n, int o, int s, int tmtTagCount) 
     {
    	 //System.out.println("??" + c);
         double[] dist = new double[_possibleIsotopeCombinations.length];
         double[] means = new double[ProbC.length];
         double[] exps = new double[means.length];
         for (int i = 0; i < means.length; i++) // precalculate means and thier exps
         {
             means[i] = c * ProbC[i] + h * ProbH[i] + n * ProbN[i] + o * ProbO[i] + s * ProbS[i] + (i == 0 ? tmtTagCount : 0);
             exps[i] = Math.exp(means[i]);
         }

         // This assumes that the envelop is unimodal.
         double maxHeight = 0.0;
         int isotopeIndex = 0;
         setMostIntenseIsotopomerIndex(-1);
         for (; isotopeIndex < _possibleIsotopeCombinations.length; isotopeIndex++)
         {
             for (int[] isopeCombinations : _possibleIsotopeCombinations[isotopeIndex])
             {
                 dist[isotopeIndex] += GetIsotopeProbability(isopeCombinations, means, exps);
             }
         //    if (Double.isInfinite(dist[isotopeIndex]))
          //   {
                 
           //  }
             if (dist[isotopeIndex] > maxHeight)
             {
                 maxHeight = dist[isotopeIndex];
                 setMostIntenseIsotopomerIndex(isotopeIndex);
             }
            // else if (dist[isotopeIndex] < maxHeight * IsotopeRelativeIntensityThreshold)
            // {
            //     break;
            // }
         }

         double[] normalizedDist = new double[dist.length];
         double norm = 0;
         for (int i = 0; i < dist.length; i++) norm += dist[i] * dist[i];
         norm = Math.sqrt(norm);        
         
         
         for (int i = 0; i < normalizedDist.length; i++)
         {
             normalizedDist[i] = dist[i] / norm;
         }

         return normalizedDist;
     }
     
     
	public IsotopeDistributionGenerator(Composition composition, float monoIsotopeMass, int maxBin){
		this(composition, monoIsotopeMass, maxBin, 0);
	}
	
	public IsotopeDistributionGenerator(Composition composition, float monoIsotopeMass, int maxBin, int tmtTagCount){ // composition should exclude TMT compositions
		ComputePossibleIsotopeCombinations(maxBin);
		//H(20) C(8) 13C(4) N 15N O(2)	
		dist = GetIsotopomerEnvelop(composition.getC(), composition.getH(), 
				composition.getN(), composition.getO(), composition.getS(), tmtTagCount); // TODO make TMT class..
		
		//double tmtMass = 20 * Composition.H + 8 * Composition.C + 4 * Composition.C13 + Composition.N + Composition.N15 + 2 * Composition.O;
		//System.out.println(tmtMass);
		setMonoIsotopeMass(monoIsotopeMass);
	}
	
	
	public double[] getDist() { return dist; }
	
	
	public HashMap<Float, Double> getDistMap(int charge){
		HashMap<Float, Double> dm = new HashMap<Float, Double>();
		float mass = monoIsotopeMass;
		
		for(int i=0;i<dist.length;i++){
			double d = dist[i];
			dm.put(mass/charge, d);
			mass += Composition.ISOTOPE;
		}
		return dm;
	}
	
	
	static public void main(String[] args){
		for(double d : new IsotopeDistributionGenerator(new Peptide("LRGM+15.995NEHQC+57.021ELR").getComposition(),1560.7356f, 6, 0).getDist()){
			System.out.println(d);
		}
	}

	public int getMostIntenseIsotopomerIndex() {
		return mostIntenseIsotopomerIndex;
	}


	public void setMostIntenseIsotopomerIndex(int mostIntenseIsotopomerIndex) {
		this.mostIntenseIsotopomerIndex = mostIntenseIsotopomerIndex;
	}


	public float getMonoIsotopeMass() {
		return monoIsotopeMass;
	}

	public float getMaxMass(){
		return (float) (monoIsotopeMass + dist.length * Composition.ISOTOPE);
	}

	public void setMonoIsotopeMass(float monoIsotopeMass) {
		this.monoIsotopeMass = monoIsotopeMass;
	}
	
}
