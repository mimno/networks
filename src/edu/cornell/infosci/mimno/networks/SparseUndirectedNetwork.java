package edu.cornell.infosci.mimno.networks;

import gnu.trove.*;
import cc.mallet.types.*;
import cc.mallet.util.Randoms;
import java.io.*;
import java.util.zip.*;
import java.util.*;

public class SparseUndirectedNetwork {

	HashMap<String, HashSet<String>> edges;
	ArrayList<String> nodes;
	HashMap<String,Integer> dictionary;

	// Held-out pairs combines both testing and validation pairs.
	// Inside the training loop, we don't care about the distinction,
	//  we just need to know which pairs to ignore.
	HashMap<String, HashSet<String>> heldOutPairs = null;
	int[][] heldOutPairIDs = null;
	int[][] trainingPairIDs = null;
	
	HashMap<String, HashSet<String>> validationPairs = null;
	HashMap<String, HashSet<String>> testingPairs = null;

	int numNodes;
	int numEdges;

	public SparseUndirectedNetwork (String filename) throws IOException {
		numEdges = 0;
		
		edges = new HashMap<String, HashSet<String>>();
		nodes = new ArrayList<String>();
		dictionary = new HashMap<String,Integer>();

		numEdges = fill(filename, edges);

		numNodes = nodes.size();

		trainingPairIDs = new int[numNodes][];
		for (int leftNode = 0; leftNode < numNodes; leftNode++) {
			HashSet<String> trainingRightNodes = edges.get(nodes.get(leftNode));
			trainingPairIDs[leftNode] = new int[ trainingRightNodes.size() ];
			int i = 0;
			for (String rightNode: trainingRightNodes) {
				trainingPairIDs[leftNode][i] = dictionary.get(rightNode);
				i++;
			}
		}
	}

	public int fill(String filename, HashMap<String, HashSet<String>> pairs) throws IOException {
		BufferedReader in = null;
		int numPairs = 0;

		if (filename.endsWith(".gz")) {
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))));
		}
		else {
			in = new BufferedReader(new FileReader(filename));
		}

		String line = null;
		while ((line = in.readLine()) != null) {
			if (! line.startsWith("#") && ! line.equals("")) {
				String[] fields = line.split("\\s+");
				assert(fields.length >= 2) : line;

				String node1 = fields[0];
				String node2 = fields[1];

				lookup(node1);
				lookup(node2);

				if (! pairs.containsKey(node1)) {				
					pairs.put(node1, new HashSet<String>());
				}
				pairs.get(node1).add(node2);

				if (! pairs.containsKey(node2)) {
					pairs.put(node2, new HashSet<String>());
				}
				pairs.get(node2).add(node1);

				numPairs++;
				if (numPairs % 100000 == 0) { System.out.println(numPairs); }
			}
		}
		in.close();

		return numPairs;
	}

	public void lookup(String nodeName) {
		if (! dictionary.containsKey(nodeName)) {
			dictionary.put(nodeName, nodes.size());
			nodes.add(nodeName);
		}
	}

	public int getDegree(int node) {
		return getDegree(nodes.get(node));
	}

	public int getDegree(String nodeName) {
		return edges.get(nodeName).size();
	}

	public int getNumEdges() {
		return numEdges;
	}

	public void loadTesting(String filename) throws IOException {
		testingPairs = new HashMap<String, HashSet<String>>();
		fill(filename, testingPairs);		

		if (heldOutPairs == null) {
			heldOutPairs = new HashMap<String, HashSet<String>>();
		}
		heldOutPairs.putAll(testingPairs);
		cacheHeldOutPairIDs();
	}

	public void loadValidation(String filename) throws IOException {
		validationPairs = new HashMap<String, HashSet<String>>();
		fill(filename, validationPairs);		

		if (heldOutPairs == null) {
			heldOutPairs = new HashMap<String, HashSet<String>>();
		}
		heldOutPairs.putAll(validationPairs);
		cacheHeldOutPairIDs();
	}

	public void cacheHeldOutPairIDs() {
		heldOutPairIDs = new int[numNodes][];

		for (int leftNode = 0; leftNode < numNodes; leftNode++) {
			HashSet<String> heldOutRightNodes = getHeldOutNodes(nodes.get(leftNode));
			if (heldOutRightNodes != null) {
				// Cache the list of heldout nodes
				heldOutPairIDs[leftNode] = new int[ heldOutRightNodes.size() ];
				int i = 0;
				for (String rightNode: heldOutRightNodes) {
					heldOutPairIDs[leftNode][i] = dictionary.get(rightNode);
					i++;
				}
				
				// Also recache the training nodes
				HashSet<String> trainingRightNodes = (HashSet<String>) edges.get(nodes.get(leftNode)).clone();
				trainingRightNodes.removeAll(heldOutRightNodes);
				trainingPairIDs[leftNode] = new int[ trainingRightNodes.size() ];
				i = 0;
				for (String rightNode: trainingRightNodes) {
					trainingPairIDs[leftNode][i] = dictionary.get(rightNode);
					i++;
				}
			}
		}
	}

	public HashSet<String> getHeldOutNodes(String node) {
		if (heldOutPairs == null) { return null; }
		if (heldOutPairs.containsKey(node)) { return heldOutPairs.get(node); }
		return null;
	}

	public int[] getHeldOutPairIDs(String node) {
		if (heldOutPairs == null || heldOutPairIDs == null) { return null; }
		return heldOutPairIDs[ dictionary.get(node) ];
	}

	public int[] getHeldOutPairIDs(int nodeID) {
		if (heldOutPairs == null || heldOutPairIDs == null) { return null; }
		return heldOutPairIDs[nodeID];
	}

	public int[] getTrainingPairIDs(String node) {
		return trainingPairIDs[ dictionary.get(node) ];
	}

	public int[] getTrainingPairIDs(int nodeID) {
		return trainingPairIDs[nodeID];
	}

	public double getAveragePrecision(double[][] nodeCommunityWeights, HashMap<String, HashSet<String>> evaluationSets, int limit) {
		int numCommunities = nodeCommunityWeights[0].length;
		ArrayList<IDSorter> weightedPairs = new ArrayList<IDSorter>();
		double averagePrecision = 0.0;

		if (evaluationSets == null) { return averagePrecision; }
		for (String leftNodeName: evaluationSets.keySet()) {
			int leftNode = dictionary.get(leftNodeName);
			double[] leftWeights = nodeCommunityWeights[leftNode];
			for (String rightNodeName: evaluationSets.get(leftNodeName)) {
				int rightNode = dictionary.get(rightNodeName);
				
				// the list is symmetric, so only process each link once
				if (rightNode < leftNode) { continue; }

				double[] rightWeights = nodeCommunityWeights[rightNode];
				double innerProduct = 0.0;
				for (int community = 0; community < numCommunities; community++) {
					innerProduct += leftWeights[community] * rightWeights[community];
				}

				int edgeCount = 0;
				if (edges.get(leftNodeName).contains(rightNodeName)) {
					edgeCount = 1;
				}
				
				// Note that the log probability of a zero is the *negative* inner product
				weightedPairs.add(new IDSorter(edgeCount, innerProduct));
			}
		}

		Collections.sort(weightedPairs);

		double numEdgesSoFar = 0.0; // double to avoid casting to floating point
		if (limit == 0) { limit = weightedPairs.size(); }
		for (int rank = 0; rank < limit; rank++) {
			numEdgesSoFar += weightedPairs.get(rank).getID();
			averagePrecision += numEdgesSoFar / (rank + 1);
		}

		averagePrecision /= limit;

		return averagePrecision;
	}

	public double getHeldOutLogProbability(double[][] nodeCommunityWeights, HashMap<String, HashSet<String>> evaluationSets) {
		return getHeldOutLogProbability(nodeCommunityWeights, evaluationSets, null);

	}

	public double getHeldOutLogProbability(double[][] nodeCommunityWeights, HashMap<String, HashSet<String>> evaluationSets, PrintWriter out) {
		int numCommunities = nodeCommunityWeights[0].length;

		double heldOutLogProbability = 0.0;
		if (evaluationSets == null) { return heldOutLogProbability; }
		for (String leftNodeName: evaluationSets.keySet()) {
			int leftNode = dictionary.get(leftNodeName);
			double[] leftWeights = nodeCommunityWeights[leftNode];
			for (String rightNodeName: evaluationSets.get(leftNodeName)) {
				int rightNode = dictionary.get(rightNodeName);
				
				// the list is symmetric, so only process each link once
				if (rightNode < leftNode) { continue; }

				double[] rightWeights = nodeCommunityWeights[rightNode];
				double innerProduct = 0.0;
				for (int community = 0; community < numCommunities; community++) {
					innerProduct += leftWeights[community] * rightWeights[community];
				}
				heldOutLogProbability -= innerProduct;

				// Is there a link between these nodes?				
				boolean isEdge = edges.get(leftNodeName).contains(rightNodeName);

				if (isEdge) {
					heldOutLogProbability += Math.log(innerProduct);
				}

				if (out != null) {
					out.format("%s\t%s\t%f\t%f\n", leftNodeName, rightNodeName, -innerProduct, isEdge ? Math.log(innerProduct) : 0.0);
				}
			}
		}
		
		return heldOutLogProbability;
	}
	
	public Object[] getNodes() {
		return nodes.toArray();
	}

	public int getID(String s) {
		return dictionary.get(s);
	}

	public static void main(String[] args) throws Exception {
		SparseUndirectedNetwork network = new SparseUndirectedNetwork(args[0]);

		for (String node2 : network.edges.get("1")) {
			System.out.println(node2);
		}

		Object[] array = network.getNodes();
		for (int i = 0; i < 30; i++) {
			System.out.println(" " + array[i]);
		}
	}

}