package edu.cornell.infosci.mimno.networks;

import edu.cornell.infosci.mimno.UnicodeBarplot;
import gnu.trove.*;
import cc.mallet.types.*;
import cc.mallet.util.Randoms;
import java.io.*;
import java.util.zip.*;
import java.util.*;

public class CommunityReports {

	SparseUndirectedNetwork network;
	int numNodes;
	int numCommunities;
	public double[][] nodeCommunityWeights;
	double[] communitySums;
	int numUpdates = 0;

	public CommunityReports(SparseUndirectedNetwork network, int numCommunities, double[][] initialWeights) {
		this.network = network;
		this.numNodes = network.numNodes;
		this.numCommunities = numCommunities;

		setNodeCommunityWeights(initialWeights);
	}

	/** Create a new copy of the node community weights */
	public void setNodeCommunityWeights(double[][] weights) {
		this.nodeCommunityWeights = new double[numNodes][numCommunities];
		this.communitySums = new double[numCommunities];
		for (int node = 0; node < numNodes; node++) {
			System.arraycopy(weights[node], 0, nodeCommunityWeights[node], 0, numCommunities);
			for (int community = 0; community < numCommunities; community++) {
				communitySums[community] += nodeCommunityWeights[node][community];
			}
		}
		numUpdates = 1;
	}

	public void averageNodeCommunityWeights(double[][] weights) {
		numUpdates++;
		double oldFactor = (numUpdates - 1.0) / numUpdates;
		double newFactor = 1.0 / numUpdates;
		
		Arrays.fill(communitySums, 0.0);

		for (int node = 0; node < numNodes; node++) {
			for (int community = 0; community < numCommunities; community++) {
				nodeCommunityWeights[node][community] *= oldFactor;
				nodeCommunityWeights[node][community] += newFactor * weights[node][community];
				communitySums[community] += nodeCommunityWeights[node][community];
			}
		}
	}

	/** A lower-memory version of getSortedNodes() */
	public IDSorter[][] getSortedNodes(int limit) {

		if (limit > numNodes) { limit = numNodes; }

		IDSorter[][] sortedNodes = new IDSorter[numCommunities][limit];

		double[] communityMinimums = new double[numCommunities];
		Arrays.fill(communityMinimums, Double.POSITIVE_INFINITY);
		for (int community = 0; community < numCommunities; community++) {
			for (int i = 0; i < limit; i++) {
				sortedNodes[community][i] = new IDSorter(i, nodeCommunityWeights[i][community]);
				if (nodeCommunityWeights[i][community] < communityMinimums[community]) {
					communityMinimums[community] = nodeCommunityWeights[i][community];
				}
			}
			Arrays.sort(sortedNodes[community]);
		}

		// We've already recorded the first [limit] nodes, so proceed from there
		for (int node = limit; node < numNodes; node++) {
			for (int community = 0; community < numCommunities; community++) {
				// Is this node worth even considering?
				if (nodeCommunityWeights[node][community] > communityMinimums[community]) {
					// The last element is guaranteed to be the smallest, so replace it
					sortedNodes[community][limit-1] = new IDSorter(node, nodeCommunityWeights[node][community]);
					// The new last node may belong elsewhere in the sorted list...
					Arrays.sort(sortedNodes[community]);
					// ... but now we can guarantee that the last is the smallest.
					communityMinimums[community] = sortedNodes[community][limit - 1].getWeight();
				}
			}
		}

		return sortedNodes;
	}

    public double getDensity(IDSorter[] sortedNodes, int n) {
		int links = 0;
		int pairs = 0;

		if (n > sortedNodes.length) { n = sortedNodes.length; }
		
		int i = 0;
        Object[] topNodes = new Object[n];
        for (IDSorter sorter : sortedNodes) {
            topNodes[i] = network.nodes.get(sorter.getID());
            i++;
			if (i == n) { break; }
        }

        for (int row = 0; row < n-1; row++) {
            for (int col = row+1; col < n; col++) {
				if (network.edges.get(topNodes[row]).contains(topNodes[col])) {
                    links++;
                }
                pairs++;
            }
		}

        return (double) links / pairs;
    }

	public double totalDensity(IDSorter[][] communitySortedNodes, int numNodesToDisplay) {
        double sumDensity = 0.0;

        for (int community = 0; community < numCommunities; community++) {
            double density = getDensity(communitySortedNodes[community], numNodesToDisplay);
            sumDensity += density;
        }

        return sumDensity / numCommunities;
    }

	public void printCommunitySizes() {
		System.out.println(UnicodeBarplot.getBars(communitySums));
	}

	public void printCommunityNodes(IDSorter[][] communitySortedNodes, int numNodesToDisplay) {
		for (int community = 0; community < numCommunities; community++) {
			IDSorter[] sortedNodes = communitySortedNodes[community];
                        
			double density = getDensity(sortedNodes, numNodesToDisplay);
			System.out.format("%d\t%f\t%f\n", community, communitySums[community], density);
			//System.out.format("%d\t%f\t%f\n", community, parameterSum, density);

			int position = 0;
			for (IDSorter sorter : sortedNodes) {
				//System.out.format(" %s:%.3f", nodes[sorter.getID()], sorter.getWeight());
				System.out.format(" %s", network.nodes.get(sorter.getID()));

				position++;
				if (position == numNodesToDisplay) {
					break;
				}
			}

			System.out.println();
		}
	}

	public void printCommunityNodes(IDSorter[][] communitySortedNodes, int numNodesToDisplay, double[] shapes, double[] rates) {
		for (int community = 0; community < numCommunities; community++) {
			IDSorter[] sortedNodes = communitySortedNodes[community];
                        
			double density = getDensity(sortedNodes, numNodesToDisplay);
			System.out.format("%d\t%f\t%f\t%f/%f\n", community, communitySums[community], density, shapes[community], rates[community]);
			//System.out.format("%d\t%f\t%f\n", community, parameterSum, density);

			int position = 0;
			for (IDSorter sorter : sortedNodes) {
				//System.out.format(" %s:%.3f", nodes[sorter.getID()], sorter.getWeight());
				System.out.format(" %s", network.nodes.get(sorter.getID()));

				position++;
				if (position == numNodesToDisplay) {
					break;
				}
			}

			System.out.println();
		}
	}

	public void writeNodeCommunities(PrintWriter out) {
		for (int node = 0; node < numNodes; node++) {
			Formatter formatter = new Formatter();
			formatter.format("%s\t%d", network.nodes.get(node), node);
			for (int community = 0; community < numCommunities; community++) {
				formatter.format("\t%f", nodeCommunityWeights[node][community]);
				//formatter.format("\t%f / %f", nodeCommunityShapes[node][community], nodeCommunityRates[node][community]);
			}
			out.println(formatter);
		}
	}

	public void communityNodeJSON(IDSorter[][] communitySortedNodes, int numNodesToDisplay, PrintWriter out) {
		out.print("[");

		for (int community = 0; community < numCommunities; community++) {
			IDSorter[] sortedNodes = communitySortedNodes[community];
			
			Formatter buffer = new Formatter();
			
			double density = getDensity(sortedNodes, numNodesToDisplay);
			buffer.format("{\"id\":%d, \"sum\": %f, \"density\": %f, \"nodes\":[", community, communitySums[community], density);

			int position = 0;
			for (IDSorter sorter : sortedNodes) {
				//System.out.format(" %s:%.3f", nodes[sorter.getID()], sorter.getWeight());
				buffer.format("{ \"name\": \"%s\", \"weight\": %f }", network.nodes.get(sorter.getID()), sorter.getWeight());

				position++;
				if (position == numNodesToDisplay) {
					break;
				}
				else { buffer.format(","); }
			}

			buffer.format("]}");
			if (community < numCommunities - 1) { out.print(buffer + ","); }
			else { out.print(buffer); }
		}
		out.println("]");
	}

	public double getHeldOutLogProbability() {
		return network.getHeldOutLogProbability(nodeCommunityWeights, network.testingPairs, null);
	}

	public double getHeldOutLogProbability(PrintWriter out) {
		return network.getHeldOutLogProbability(nodeCommunityWeights, network.testingPairs, out);
	}
}