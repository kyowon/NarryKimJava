package rpf.deprecated;

import java.io.File;
import java.io.IOException;

import net.sf.picard.reference.IndexedFastaSequenceFile;

public class test {

	public static void main(String[] args) throws IOException {
		IndexedFastaSequenceFile f = new IndexedFastaSequenceFile(new File("/media/kyowon/Data1/RPF_Project/data/hg19.fa"));
		
		for(byte b : f.getSubsequenceAt("chr14", (long)39565322-2, 39565322).getBases())
			System.out.println((char)b);//39565322 - 3 : 39565321 written : 39565322
		
		for(byte b : f.getSubsequenceAt("chr9", 136325189, 136325191).getBases())
			System.out.println((char)b);//136325188 : 136325190 written 136325188

		
		f.close();
	}

}
