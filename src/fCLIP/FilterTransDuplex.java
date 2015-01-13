package fCLIP;

public class FilterTransDuplex {

	public static void main(String[] args) {
		PipeLine.filterTransPairs(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]));
	}

}
