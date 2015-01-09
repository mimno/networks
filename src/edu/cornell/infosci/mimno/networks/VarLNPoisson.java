package edu.cornell.infosci.mimno.networks;

import edu.cornell.infosci.mimno.UnicodeBarplot;
import java.util.*;
import cc.mallet.types.*;
import cc.mallet.util.*;
import java.io.*;

public class VarLNPoisson {

	static cc.mallet.util.CommandOption.String inputFile = new cc.mallet.util.CommandOption.String(VarLNPoisson.class, "input", "FILENAME", true, null,
		 "The filename of a network file.", null);
	
	static cc.mallet.util.CommandOption.String testingFile = new cc.mallet.util.CommandOption.String(VarLNPoisson.class, "testing-file", "FILENAME", true, null,
		 "The filename of a set of held-out pairs.", null);
	
	static cc.mallet.util.CommandOption.String validationFile = new cc.mallet.util.CommandOption.String(VarLNPoisson.class, "validation-file", "FILENAME", true, null,
		 "The filename of a set of pairs used to diagnose convergence.", null);
	
	static cc.mallet.util.CommandOption.String testingProbsFile = new cc.mallet.util.CommandOption.String(VarLNPoisson.class, "testing-results-file", "FILENAME", true, "heldout.txt",
		 "The filename for the probabilities of a set of held-out pairs.", null);
	
	static cc.mallet.util.CommandOption.String modelFile = new cc.mallet.util.CommandOption.String(VarLNPoisson.class, "model-file", "FILENAME", true, null,
		 "The filename to write node-community weights (tab delimited).", null);
	
	static cc.mallet.util.CommandOption.String communityNodesFile = new cc.mallet.util.CommandOption.String(VarLNPoisson.class, "community-nodes-file", "FILENAME", true, null,
		 "The filename to write top nodes for each community.", null);
	
	static cc.mallet.util.CommandOption.Integer numCommunitiesOption = new cc.mallet.util.CommandOption.Integer(VarLNPoisson.class, "num-communities", "INTEGER", true, 10,
		 "The number of communities to find.", null);

	static cc.mallet.util.CommandOption.Integer numIterationsOption = new cc.mallet.util.CommandOption.Integer(VarLNPoisson.class, "num-iterations", "INTEGER", true, 200,
		 "The number of gradient steps.", null);

	static cc.mallet.util.CommandOption.Double meanPriorOption = new cc.mallet.util.CommandOption.Double(VarLNPoisson.class, "mean-prior", "NUMBER", true, 0.0001,
		 "The exponentiated mean (must be positive).", null);
	
	static cc.mallet.util.CommandOption.Double precisionPriorOption = new cc.mallet.util.CommandOption.Double(VarLNPoisson.class, "precision-prior", "NUMBER", true, 0.05,
		 "The precision (inverse variance) of the distribution.", null);

	static cc.mallet.util.CommandOption.Boolean variationalOption = new cc.mallet.util.CommandOption.Boolean(VarLNPoisson.class, "variational", "TRUE/FALSE", true, true,
		 "Whether to use a Laplace variational update.", null);

	static cc.mallet.util.CommandOption.Double learningRateOption = new cc.mallet.util.CommandOption.Double(VarLNPoisson.class, "learning-rate", "NUMBER", true, 0.9,
		 "A scalar factor for the learning rate.", null);

	int numCommunities;
	int numNodes;
	Object[] nodes;
	SparseUndirectedNetwork network;
	
	double[] nodeMeans;

	double[][] nodeCommunityWeights;
	double[][] nodeCommunityGradients;

	double[] communitySums;

	boolean heldOutPairsExist = false;

	Randoms random;
	
	public VarLNPoisson(String networkFilename, int numCommunities) throws Exception {
		network = new SparseUndirectedNetwork(networkFilename);
		nodes = network.getNodes();
		
		this.numNodes = nodes.length;
		this.numCommunities = numCommunities;

		nodeMeans = new double[numNodes];
		Arrays.fill(nodeMeans, 0.0);

		nodeCommunityWeights = new double[numNodes][numCommunities];
		nodeCommunityGradients = new double[numNodes][numCommunities];

		communitySums = new double[numCommunities];

		random = new Randoms();

		for (int node = 0; node < numNodes; node++) {
			for (int community = 0; community < numCommunities; community++) {
				nodeCommunityWeights[node][community] = 0.01 * random.nextUniform();
			}
		}
	}

	public void loadTesting(String filename) throws IOException {
		heldOutPairsExist = true;
		network.loadTesting(filename);
	}

	public void loadValidation(String filename) throws IOException {
		heldOutPairsExist = true;
		network.loadValidation(filename);
	}

	public double iterate(double learningRate, double precisionPrior) {

		double logLikelihood = 0.0;
		double entropy = 0.0;

		Arrays.fill(communitySums, 0.0);

		double[] parameterSums = new double[numCommunities];
		for (int node = 0; node < numNodes; node++) {
			Arrays.fill(nodeCommunityGradients[node], 0.0);
			for (int community = 0; community < numCommunities; community++) {
				parameterSums[community] += nodeCommunityWeights[node][community];
			}
		}

		double[] distribution = new double[numCommunities];

		int links = 0;

		HashSet<String> heldOutNodes = null;

		for (int leftNode = 0; leftNode < numNodes; leftNode++) {
			for (int rightNode: network.getTrainingPairIDs(leftNode)) {
				links++;
				
				if (leftNode > rightNode) { continue; }
				
				double sum = 0.0;
				for (int community = 0; community < numCommunities; community++) {
					distribution[community] = nodeCommunityWeights[leftNode][community] * 
						nodeCommunityWeights[rightNode][community];
					sum += distribution[community];
				}

				logLikelihood += Math.log(sum) - sum;

				double normalizer = 1.0 / sum;
				for (int community = 0; community < numCommunities; community++) {
					double prob = distribution[community] * normalizer;
					//if (prob > 0.0) { entropy -= prob * Math.log(prob); }
					nodeCommunityGradients[leftNode][community] += prob;
					nodeCommunityGradients[rightNode][community] += prob;
					communitySums[community] += 2.0 * prob;
					assert(! Double.isNaN(communitySums[community])) : "NaN in link part"; 
				}
			}

			for (int community = 0; community < numCommunities; community++) {
				//nodeCommunityWeights[leftNode][community] = Math.max(nodeCommunityWeights[leftNode][community], 0.000001);
				double value = precisionPrior * (nodeMeans[leftNode] - Math.log(nodeCommunityWeights[leftNode][community]));
				nodeCommunityGradients[leftNode][community] += value;
				communitySums[community] += value;
				assert(! Double.isNaN(communitySums[community])) : "NaN in prior part: " + nodeMeans[leftNode] + " - Math.log(" + nodeCommunityWeights[leftNode][community] + "), " + nodeMeans[leftNode]; 
			}
		}

		double[] communityNormalizers = new double[numCommunities];
		for (int community = 0; community < numCommunities; community++) {
			communityNormalizers[community] = 1.0 / Math.sqrt(communitySums[community]);
		}

		for (int node = 0; node < numNodes; node++) {
			for (int community = 0; community < numCommunities; community++) {
				double before = nodeCommunityWeights[node][community];
				nodeCommunityWeights[node][community] *= (1.0 - learningRate);
				nodeCommunityWeights[node][community] += learningRate * nodeCommunityGradients[node][community] * communityNormalizers[community];
				assert(! Double.isNaN(nodeCommunityWeights[node][community])) : "NaN: " + before + " * " +  nodeCommunityGradients[node][community] + " / sqrt( " + communitySums[community] + " )";
				assert(! Double.isInfinite(nodeCommunityWeights[node][community])) : "infinite: " + before + " * " +  nodeCommunityGradients[node][community] + " * " + communityNormalizers[community];
			}
		}

		if (variationalOption.value) {
			double[] communitySums = new double[numCommunities];
			for (int node = 0; node < numNodes; node++) {
				for (int community = 0; community < numCommunities; community++) {
					communitySums[community] += nodeCommunityWeights[node][community];
				}
			}
			
			for (int node = 0; node < numNodes; node++) {
				//if (node == 0) { System.out.format("%f %d ", learningRate, network.edges.get(nodes[0]).size()); }
				for (int community = 0; community < numCommunities; community++) {
					double value = 
						Math.exp(0.05 * learningRate * 0.5 / (precisionPrior + nodeCommunityWeights[node][community] * 
															  (communitySums[community] -
															   nodeCommunityWeights[node][community])));
					nodeCommunityWeights[node][community] *= value;
					
					if (value > 150.0) { System.out.format("%.4f\t%f\t%f\n", value, nodeCommunityWeights[node][community], communitySums[community]); }
				}
			}
		}
			
		return logLikelihood / links; 
	}

	public void setNodeMeanPriors() {
		for (int node = 0; node < numNodes; node++) {
			double sum = 0.0; //meanPriorOption.value;
			for (int community = 0; community < numCommunities; community++) {
				sum += Math.log(nodeCommunityWeights[node][community] + 0.000001);
			}
			nodeMeans[node] = sum / (2.0 * numCommunities);
		}
	}	

	public static void main (String[] args) throws Exception {
		CommandOption.setSummary (VarLNPoisson.class, "Batch MLE inference for a Poisson community model");
		CommandOption.process (VarLNPoisson.class, args);
		
		VarLNPoisson poisson = new VarLNPoisson(inputFile.value, numCommunitiesOption.value);

		if (testingFile.wasInvoked()) {
			poisson.loadTesting(testingFile.value);
		}

		if (validationFile.wasInvoked()) {
			poisson.loadValidation(validationFile.value);
		}

		CommunityReports reports = new CommunityReports(poisson.network, poisson.numCommunities, poisson.nodeCommunityWeights);

		double previousScore = Double.NEGATIVE_INFINITY;

		double previousValidation = Double.NEGATIVE_INFINITY;
		int numIterationsWithTinyChange = 0;
		int numIterationsWithNegativeChange = 0;
	
		double precisionPrior = 0.001;

		for (int iteration = 0; iteration < numIterationsOption.value; iteration++) {
			//double score = poisson.iterate(learningRateOption.value * Math.pow(1.0 / (10 + iteration), 0.6));

			long startTime = System.currentTimeMillis();
			double score = poisson.iterate(learningRateOption.value, precisionPrior);
			long iterationTime = System.currentTimeMillis() - startTime;

			poisson.setNodeMeanPriors();

			if (validationFile.wasInvoked()) {
				double currentValidation =
					poisson.network.getHeldOutLogProbability(poisson.nodeCommunityWeights, poisson.network.validationPairs, null);

				double averagePrecision =
					poisson.network.getAveragePrecision(poisson.nodeCommunityWeights, poisson.network.validationPairs, 0);
				
				System.out.format("%f\t%f\t%f\t%f\n", currentValidation, score, averagePrecision, Math.abs(currentValidation - previousValidation));

				if (Math.abs(currentValidation - previousValidation) < 0.001) {
					numIterationsWithTinyChange++;
				}
				if (currentValidation < previousValidation) {
					numIterationsWithNegativeChange++;
				}
				else {
					numIterationsWithNegativeChange = 0;
				}

				previousValidation = currentValidation;
			}

			if (iteration % 10 == 0 || numIterationsWithTinyChange >= 5 || numIterationsWithNegativeChange >= 3) {
				reports.setNodeCommunityWeights(poisson.nodeCommunityWeights);
				double heldout = 0.0;
				if (testingProbsFile.wasInvoked()) {
					PrintWriter out = new PrintWriter(testingProbsFile.value);
					heldout = reports.getHeldOutLogProbability(out);
					out.close();
				}

				IDSorter[][] sortedNodes = reports.getSortedNodes(30);
				//poisson.printCommunityNodes(sortedNodes, 10);
				System.out.format("%d\t%f\t%f\t%d\t", iteration, heldout, reports.totalDensity(sortedNodes, 10), iterationTime);
				reports.printCommunitySizes();
			}

			//System.out.format("%d\t%.2f\n", iteration, score);
			//if (score < previousScore) { break; }
			previousScore = score;
			
			if (numIterationsWithTinyChange >= 5 || numIterationsWithNegativeChange >= 3) {
				break;
			}
		}
		//poisson.printNodeCommunities();
		IDSorter[][] sortedNodes = reports.getSortedNodes(30);
		//reports.printCommunityNodes(sortedNodes, 10);

		if (communityNodesFile.wasInvoked()) {
			PrintWriter out = new PrintWriter(communityNodesFile.value);
			reports.communityNodeJSON(sortedNodes, 30, out);
			out.close();
		}

		if (modelFile.wasInvoked()) {
			PrintWriter out = new PrintWriter(modelFile.value);
			reports.writeNodeCommunities(out);
			out.close();
		}
	}
	
}