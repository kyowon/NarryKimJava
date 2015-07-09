package tmt;


import java.io.*;

import util.Codon;
import jspp.SearchMatrix;
import jspp.SignalPeptidePredictor;



public class JSPP {
  public JSPP() {
  }


  public static void main(String[] args) {
    
     
    
 //   System.out.println(pmatrix.getAbsolutePath());
  //  System.exit(1);
    
    
    
      try {
      File pmatrix=new File("extlib/jspp/matrices(SearchMatrix-objects)/eukarya.smx");
      ObjectInputStream in = new ObjectInputStream(new FileInputStream(pmatrix));
      SearchMatrix smp = (SearchMatrix) in.readObject();
      in.close();
       
      SignalPeptidePredictor pd = new SignalPeptidePredictor(smp);
      String sequence = "ATGCTAAGTTCCCGGGCGCAGGCGGCGAGGACGGCGGCCGACAAGGCCCT";
      sequence = Codon.getAminoAcids(sequence);
      sequence = "MAVMAPRTLLLLLSGALALTQTWAGSHSMR";
      System.out.println(sequence);
      
        int cpos = pd.predictEnhancedPosition(sequence);
        String yn;
        if (pd.isSignalPeptide()) {
          yn = "Y";
        }
        else {
          yn = "N";
        }
        double score=pd.getScore();
        System.out.println(cpos + "\t" + yn + "\t"+score+"\n");
        
        
        sequence = "PEPTIDE";
        pd.predictEnhancedPosition(sequence);
        
        
        System.out.println(cpos + "\t" + yn + "\t"+pd.getScore()+"\n");
     // }
      //se.close();
      //out.close();
    }
    catch (IOException ex) {
      System.err.println("IO-Error: "+ex);
    }
    catch (ClassNotFoundException ex) {
      System.err.println("Error reading matrices (wrong file ?)");
    }


  }

}

