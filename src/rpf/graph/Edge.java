package rpf.graph;


public class Edge {
	private Node lNode;
	private Node rNode;
	private double score;
	
	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public Node getlNode() {
		return lNode;
	}

	public Node getrNode() {
		return rNode;
	}

	public Edge(Node lNode, Node rNode){
		this.lNode = lNode;
		this.rNode = rNode;
	}

	@Override
	public int hashCode() {
		return lNode.hashCode() * rNode.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(super.equals(obj)) return true;
		if(obj instanceof Edge){
			Edge o = (Edge)obj;
			return lNode.equals(o.lNode) && rNode.equals(o.rNode);
		}
		return false;
	}

	@Override
	public String toString() {
		return lNode + " " + rNode;
	}
	
	
}
	
	