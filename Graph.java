package com.graph;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import org.jgrapht.*;
import org.jgrapht.graph.*;

class ValueComparator implements Comparator<DefaultWeightedEdge> {
    private HashMap<DefaultWeightedEdge, Double> map;
    public ValueComparator(HashMap<DefaultWeightedEdge, Double> map) {
        this.map = map;
    } 
    public int compare(DefaultWeightedEdge a, DefaultWeightedEdge b) {
        if (map.get(a) >= map.get(b)) {
            return 1;
        } else {
            return -1;
        }
    }
}
class QueueComparator implements Comparator<Integer> {
    private double[] dist;
    public QueueComparator(double[] dist) {
        this.dist = dist;
    } 
	public int compare (Integer v1, Integer v2) {
		double dist1 = dist[v1];
		double dist2 = dist[v2];
		if (dist1 >= dist2) 
			return 1;
		else
			return -1;
	}
}

public class Graph {
	public static void main (String[] args) throws Exception {
		String input = args[0];
		fileIO(input);		
	}
	//Generate a undirected weighted graph from input file
	public static SimpleWeightedGraph<String, DefaultWeightedEdge> generateGraph(String fileName) throws Exception {
		SimpleWeightedGraph<String, DefaultWeightedEdge> graph = 
				new SimpleWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		File file = new File(fileName);
		Scanner scanner = new Scanner(file);
		String first = scanner.nextLine();
		int n = Integer.valueOf(first);
		for (int i=0;i<n;i++){
			graph.addVertex(Integer.toString(i));
		}
		while(scanner.hasNext()){
			String line = scanner.nextLine();
			String set[] = line.split("\\s+");
			graph.setEdgeWeight(graph.addEdge(set[0], set[1]), Double.parseDouble(set[2]));
//			System.out.println("edge is "+graph.getEdge(set[0], set[1])+", weight is "+graph.getEdgeWeight(graph.getEdge(set[0],set[1])));
		}
		scanner.close();
		return graph;
	}
	//Breadth first traverse by Queue
	public static ArrayList<String> breadthFT (SimpleWeightedGraph<String, DefaultWeightedEdge> graph){
		ArrayList<String> result = new ArrayList<String>();
		ArrayList<String> white = new ArrayList<String>(graph.vertexSet());
		ArrayList<String> gray = new ArrayList<String>();
		ArrayList<String> black = new ArrayList<String>();
		ArrayList<String> vertexSet = new ArrayList<String>(graph.vertexSet());
		for(String vertex : vertexSet){			
			if(white.contains(vertex)){
				result.add(vertex);
				white.remove(vertex);
				gray.add(vertex);
				Queue<String> bfs = new LinkedList<String>();
				bfs.add(vertex);
				while(!bfs.isEmpty()){
					String head = bfs.poll();
					ArrayList<String> neighbors = new ArrayList<String>(Graphs.neighborListOf(graph, head));
					for(String neighbor : neighbors){
						if(white.contains(neighbor)){
							white.remove(neighbor);
							gray.add(neighbor);
							bfs.add(neighbor);
							result.add(neighbor);
						}
					}
					gray.remove(head);
					black.add(head);
				}
			}
		}
		return result;
	}
	//Dijkstra's Shortest Paths Algorithm
	public static void dijkstraSP (SimpleWeightedGraph<String, DefaultWeightedEdge> graph, PrintWriter output){
		Set<String> vertexSet = graph.vertexSet();
		for(String start:vertexSet){
			double [] dist = new double[vertexSet.size()];
			Integer [] path = new Integer [vertexSet.size()];
			PriorityQueue<Integer> vertexQueue = new PriorityQueue<Integer>(vertexSet.size(), new QueueComparator(dist));
			for(String newVertex:vertexSet){
				if(!newVertex.equals(start)){
					int i = Integer.valueOf(newVertex);
					dist[i] = Double.MAX_VALUE;					
					vertexQueue.add(i);
				}
			}
			int s = Integer.valueOf(start);
			dist[s] = 0;
			vertexQueue.add(s);
			Set<Integer> Vt = new HashSet<Integer>();
			for(int j=0;j<vertexSet.size();j++){
				int u = vertexQueue.poll();
				Vt.add(u);
				List<String> neighbors = Graphs.neighborListOf(graph, Integer.toString(u));
				for(String neighbor : neighbors){
					int v = Integer.valueOf(neighbor);
					if(!Vt.contains(v)){
						if((dist[u]+graph.getEdgeWeight(graph.getEdge(Integer.toString(u), neighbor)))<dist[v]){
							dist[v] = dist[u]+graph.getEdgeWeight(graph.getEdge(Integer.toString(u), neighbor));
							path[v] = u;
							vertexQueue.remove(v);
							vertexQueue.add(v);
						}
					}
				}
			}
			ArrayList<Integer> pp = new ArrayList<Integer>();
			
			for(int i = 0;i<path.length;i++){
				pp.add(path[i]);
			}
//			System.out.println(pp);
			//Print result into file
			for(int i =s+1;i<vertexSet.size();i++){
				System.out.print(s+" -> "+i+" = ");
				output.print(s+" -> "+i+" = ");
				if(path[i] != null){
					ArrayList<Integer> SP = new ArrayList<Integer>();
					int check = i;
					SP.add(i);
					
					while(path[check] !=s) {
						SP.add(path[check]);
						check = path[check];
					}
					SP.add(s);
					Collections.reverse(SP);
					double weight;
					for(int index=0;index<SP.size()-2;index++){
						
						weight = graph.getEdgeWeight(graph.getEdge(Integer.toString(SP.get(index)),Integer.toString(SP.get(index+1))));
					
						System.out.print("("+SP.get(index)+", "+SP.get(index+1)+", "+new DecimalFormat("##.0").format(weight)+")"+" -> ");
						output.print("("+SP.get(index)+", "+SP.get(index+1)+", "+new DecimalFormat("##.0").format(weight)+")"+" -> ");
					}
					weight = graph.getEdgeWeight(graph.getEdge(Integer.toString(SP.get(SP.size()-2)), Integer.toString(SP.get(SP.size()-1))));
					System.out.println("("+SP.get(SP.size()-2)+", "+SP.get(SP.size()-1)+", "+new DecimalFormat("##.0").format(weight)+")");
					System.out.println("	path weight = " + new DecimalFormat("##.00").format(dist[i]));
					output.println("("+SP.get(SP.size()-2)+", "+SP.get(SP.size()-1)+", "+new DecimalFormat("##.0").format(weight)+")");
					output.println("	path weight = " + new DecimalFormat("##.00").format(dist[i]));
				}else{
					System.out.println("Disconnected graph (no path)");
					output.println("Disconnected graph (no path)");
				}
			}
			
		}
	}
	//Depth first traverse by Stack
	public static ArrayList<String> depthFT (SimpleWeightedGraph<String, DefaultWeightedEdge> graph){
		ArrayList<String> result = new ArrayList<String>();
		ArrayList<String> white = new ArrayList<String>(graph.vertexSet());
		ArrayList<String> gray = new ArrayList<String>();
		ArrayList<String> black = new ArrayList<String>();
		ArrayList<String> vertexSet = new ArrayList<String>(graph.vertexSet());
		for(String vertex : vertexSet){			
			if(white.contains(vertex)){
				result.add(vertex);
				Stack<String> dfs = new Stack<String>();
				dfs.push(vertex);
				white.remove(vertex);
				gray.add(vertex);
				while(!dfs.isEmpty()){
					String visited = dfs.peek();
					ArrayList<String> neighbors = new ArrayList<String>(Graphs.neighborListOf(graph, visited));
					String next = null;
					for(String neighbor : neighbors){
						if(white.contains(neighbor)){
							next = neighbor;
							break;
						}
					}
					if(next == null){
						dfs.pop();
						gray.remove(visited);
						black.add(visited);
					} else {
						dfs.push(next);
						white.remove(next);
						gray.add(next);
						result.add(next);
					}
				}
			}
		}
		return result;
	}
	//Minimal spanning tree by Kruskal's algorithm
	public static UndirectedGraph<String, DefaultEdge> minST (SimpleWeightedGraph<String, DefaultWeightedEdge> graph){
		UndirectedGraph<String, DefaultEdge> result = 
				new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);
		for(String vertex:graph.vertexSet()){
			result.addVertex(vertex);
		}
		Set<DefaultWeightedEdge> edgeset = graph.edgeSet();
		final HashMap<DefaultWeightedEdge, Double> map = new HashMap<DefaultWeightedEdge, Double>();
		for(DefaultWeightedEdge edge:edgeset){
			double weight = graph.getEdgeWeight(edge);
			map.put(edge,weight);
		}
		TreeMap<DefaultWeightedEdge, Double> sortedMap = new TreeMap<DefaultWeightedEdge, Double>(new ValueComparator(map));
		sortedMap.putAll(map);
		for(DefaultWeightedEdge edge:sortedMap.keySet()){
			result.addEdge(graph.getEdgeSource(edge),graph.getEdgeTarget(edge));
			if(acyclicTest(result) == false){
				result.removeEdge(graph.getEdgeSource(edge),graph.getEdgeTarget(edge));
			}
			if(result.edgeSet().size()==result.vertexSet().size()-1){
				break;
			}
		}
		return result;
	}
	//Acyclic Test by DFS (Back edge detection)
	public static boolean acyclicTest (UndirectedGraph<String, DefaultEdge> graph){
		ArrayList<String> white = new ArrayList<String>(graph.vertexSet());
		ArrayList<String> vertexSet = new ArrayList<String>(graph.vertexSet());
		Set<DefaultEdge> edgeset = graph.edgeSet();
		Set<DefaultEdge> spanedge = new HashSet<DefaultEdge>();
		for(String vertex : vertexSet){			
			if(white.contains(vertex)){
				Stack<String> dfs = new Stack<String>();
				dfs.push(vertex);
				white.remove(vertex);
				while(!dfs.isEmpty()){
					String visited = dfs.peek();
					ArrayList<String> neighbors = new ArrayList<String>(Graphs.neighborListOf(graph, visited));
					String next = null;
					for(String neighbor : neighbors){
						if(white.contains(neighbor)){
							next = neighbor;
							break;
						}
					}					
					if(next == null){
						dfs.pop();
					} else {
						spanedge.add(graph.getEdge(visited, next));
						dfs.push(next);
						white.remove(next);
					}
				}
			}
		}
		for(DefaultEdge edge : edgeset){
			if (!spanedge.contains(edge)){
				return false;
			}
		}
		return true;
	}
	public static void fileIO(String input) throws Exception {
		String outputname = input.replace("input", "output").replace("_a", "");
		PrintWriter output = new PrintWriter(outputname);
		SimpleWeightedGraph<String, DefaultWeightedEdge> graph = generateGraph(input);	
		//DFS
		ArrayList<String> dfs = depthFT(graph);
		System.out.println("Depth First Search Traversal: "+dfs+"\n");
		output.println("Depth First Search Traversal: "+dfs+"\n");
		//BFS
		ArrayList<String> bfs = breadthFT(graph);
		System.out.println("Breadth First Search Traversal: "+bfs+"\n");		
		output.println("Breadth First Search Traversal: "+bfs+"\n");
		//MST
		UndirectedGraph<String, DefaultEdge> minSTree = minST(graph);
		ArrayList<String> mstTraversal = new ArrayList<String>();
		for(DefaultEdge edge : minSTree.edgeSet()){
			if(!mstTraversal.contains(minSTree.getEdgeSource(edge))){
				mstTraversal.add(minSTree.getEdgeSource(edge));
			}
			if(!mstTraversal.contains(minSTree.getEdgeTarget(edge))){
				mstTraversal.add(minSTree.getEdgeTarget(edge));
			}
		}
		System.out.println("Minimum Spanning Tree:\n"+"V = "+mstTraversal);
		output.println("Minimum Spanning Tree:\n"+"V = "+mstTraversal);		
		System.out.println("E = " + minSTree.edgeSet());
		output.println("E = " + minSTree.edgeSet());		
		System.out.print("Weight = ");
		output.print("Weigh t= ");
		double totalWeight = 0;
		List<Double> weights = new ArrayList<Double>();
		for(DefaultEdge edge : minSTree.edgeSet()){
			double weight = graph.getEdgeWeight(graph.getEdge(minSTree.getEdgeSource(edge),minSTree.getEdgeTarget(edge)));
			totalWeight+=weight;
			weights.add(weight);
		}
		System.out.println(weights);
		output.println(weights);
		System.out.println("Total weight is "+totalWeight+"\n");
		output.println("Total weight is "+totalWeight+"\n");
		//Shortest Paths
		System.out.println("Shortest Paths:");
		output.println("Shortest Paths:");
		dijkstraSP(graph,output);		
		output.close();
	}
}