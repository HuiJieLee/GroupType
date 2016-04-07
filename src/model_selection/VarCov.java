/**
 * 
 */
package model_selection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * This class calculates the variance-covariance matrix of mu_kj (estimated by MCEM)  
 * It produces the input files for Multidivtime
 * 

 * @author Hui-Jie Lee
 *
 */
public class VarCov {
	/** Number of iterations / trees in a file*/
	private int C;
	/** Number of branches */
	private int branchNum;
	/** The number of types of changes. */
	private int numTypeChanges;
	/** The number of types of states. */
	private int numTypeStates;
	/** Store number of changes */
	private int[][][] numberOfChanges;
	/** Store proportion of states*/
	private double[][][] propStates;
	/** Store mu_kj */
	private double[][] mu;
	/** Store c_k */
	private double[] c_k;
	/** Tree structure for output */
	private Tree tree;
	/** Store corresponding state b(k) for each type k*/
	private Hashtable<Integer, Integer> map;
	/** Store variance-covariance of mu_kj */
	private double[][][] cov;
	/** Store number of subset */
	private int N_G;


	/**
	 * Constructor for ungrouped case
	 * @param s
	 * @param mu
	 * @param c_k
	 * @param w_c
	 * @param tree
	 */
	public VarCov(SufficientStatistics s, double[][] mu, double[] c_k, Tree tree){
		this.C = s.getC();
		this.branchNum = s.getNumberOfBranches();
		this.numTypeChanges = s.getNumberOfTypes();
		this.numTypeStates = s.getNumberOfStates();
		this.numberOfChanges = s.getNumberOfChanges();
		this.propStates = s.getPropStates();
		this.mu = mu;
		this.c_k = c_k;
		this.tree = tree;
		this.N_G = s.getNumberOfTypes();
		this.cov = new double[numTypeChanges][branchNum][branchNum];
		this.map = s.getStartingState();
		setCov(calculateInformation());
	}
	
	/**
	 * Constructor for grouped case
	 * @param group 
	 * @param s
	 * @param mu
	 * @param c_k
	 * @param w_c
	 * @param tree
	 */
	public VarCov(ArrayList<ArrayList<Integer>> group, SufficientStatistics s, double[][] mu, double[] c_k, Tree tree){
		this.C = s.getC();
		this.branchNum = s.getNumberOfBranches();
		this.numTypeChanges = s.getNumberOfTypes();
		this.numTypeStates = s.getNumberOfStates();
		this.numberOfChanges = s.getNumberOfChanges();
		this.propStates = s.getPropStates();
		this.mu = mu;
		this.c_k = c_k;
		this.tree = tree;
		this.N_G = group.size();
		this.cov = new double[group.size()][branchNum][branchNum];
		this.map = s.getStartingState();
		setCov(calculateInformation(group));
	}
	
	/**
	 * This method computes the information matrix (ungrouped)
	 * term 1: -1/C * \sum   d^2 log P(mu|X, M^{(c)} 
	 * 						------------------------
	 * 							   d mu^2.
	 * 
	 * term 2: -1/ C *  \sum  ( d log P(mu|X, M^{(c)} )  ( d log P(mu|X, M^{(c)} )T
	 * 					     (---------------------  )  (---------------------  )
	 * 						 ( 	d mu  			   )  ( 	d mu  			  )
	 * 
	 * term 3: ( 1/C * \sum   d log P(mu|X, M^{(c)}  ) ( 1/C * \sum  d log P(mu|X, M^{(c)})T
	 * 		   (			 ----------------------- ) (	      ----------------------- )
	 * 		   (				   d mu			     ) (			     d mu			   )
	 * 
	 * @return info
	 */
	private double[][][] calculateInformation() {
		double[][][] info = new double[numTypeChanges][branchNum][branchNum];
		
		
		// calculate term 1, only on diagonal
		for (int k = 0; k < numTypeChanges; k++) {
			for (int j = 0; j < branchNum; j++) {
				double sum = 0;
				for (int c = 0; c < C; c++) {//numberOfChanges = new int[branchNum][N_G][C]; mu = new double[branchNum][N_G];
					sum += numberOfChanges[j][k][c] / (mu[j][k] * mu[j][k]);					
				}
				info[k][j][j] = sum/C;
			}
		}
		
		//calculate term 2
		for (int k = 0; k < numTypeChanges; k++) {
			int bk = map.get(k);
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < branchNum; l++) {
					double sum = 0;
					for (int c = 0; c < C; c++) {	//propStates = new double[branchNum][numState][C];
						sum += (-1) * (numberOfChanges[j][k][c]/mu[j][k] - propStates[j][bk][c]) * (numberOfChanges[l][k][c]/mu[l][k] - propStates[l][bk][c]);
					}
					info[k][j][l] += sum/C;
				}
			}
		}
		
		//calculate term 3
		double[][] term3 = new double[numTypeChanges][branchNum];
		for (int k = 0; k < numTypeChanges; k++) {
			int bk = map.get(k);
			for (int j = 0; j < branchNum; j++) {
				double sum = 0;
				for (int c = 0; c < C; c++) {
					sum +=  numberOfChanges[j][k][c]/mu[j][k] - propStates[j][bk][c];
				}
				term3[k][j] = sum/C;
			}
		}
		for (int k = 0; k < numTypeChanges; k++) {
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < branchNum; l++) {
					info[k][j][l] += term3[k][j] * term3[k][l];
				}
			}
		}
		
		/*print info
		for (int k = 0; k < numTypeChanges; k++) {
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < branchNum; l++) {
					System.out.println(info[k][j][l]);
				}
			}
		}*/
		
		return info;
	}
	
	/**
	 * This method computes the information matrix (grouped)
	 * term 1: -1/C * \sum   d^2 log P(mu|X, M^{(c)} 
	 * 						------------------------
	 * 							  d mu^2.
	 * 
	 * term 2: -1/C *  \sum  ( d log P(mu|X, M^{(c)} )  ( d log P(mu|X, M^{(c)} )T
	 * 					    (---------------------  )  (---------------------  )
	 * 						( 	d mu  			   )  ( 	d mu  			  )
	 * 
	 * term 3: ( 1/C * \sum   d log P(mu|X, M^{(c)}  ) ( 1/C * \sum    d log P(mu|X, M^{(c)})T
	 * 		   (			  ---------------------- ) (		    ----------------------- )
	 * 		   (				  d mu			     ) (				d mu			   )
	 * 
	 * @param group
	 * @return info
	 */
	private double[][][] calculateInformation(ArrayList<ArrayList<Integer>> group) {
		double[][][] info = new double[group.size()][branchNum][branchNum];
		
		
		// calculate sufficient statistics sum N_kj over all k in group g and sum c_k * phi_b(k)j over all k in group g
		int[][][] Ngjc = new int[group.size()][branchNum][C];
		double[][][] cphi = new double[group.size()][branchNum][C];

		for (int g = 0; g < group.size(); g++) {
			ArrayList<Integer> currentList = group.get(g);
			for (int j = 0; j < branchNum; j++) {
				for (int c = 0; c < C; c++) {//numberOfChanges = new int[branchNum][N_G][C]; mu = new double[branchNum][N_G];
					for (int k = 0; k < currentList.size(); k++) {	
						int index = currentList.get(k);
						int bk = map.get(index);
						Ngjc[g][j][c] += numberOfChanges[j][index][c];		//c_k = new double[s.getNumberOfTypes()]; propStates = new double[branchNum][numState][C]
						cphi[g][j][c] += c_k[index] * propStates[j][bk][c];
					}			
				}
			}
		}
				
		// calculate term 1, diagonal only
		for (int g = 0; g < group.size(); g++) {
			for (int j = 0; j < branchNum; j++) {
				double sum = 0;
				for (int c = 0; c < C; c++) {//numberOfChanges = new int[branchNum][N_G][C]; mu = new double[branchNum][N_G];
					sum +=  Ngjc[g][j][c]/(mu[j][g] * mu[j][g]);				
				}
				info[g][j][j] = sum/C;
			}
		}
		
		//calculate term 2
		for (int g = 0; g < group.size(); g++) {
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < branchNum; l++) {
					double sum = 0;
					for (int c = 0; c < C; c++) {	//propStates = new double[branchNum][numState][C];
						sum += (-1) * (Ngjc[g][j][c]/mu[j][g] - cphi[g][j][c]) * (Ngjc[g][l][c]/mu[l][g] - cphi[g][l][c]);
					}
					info[g][j][l] += sum/C;
				}
			}
		}
		
		//calculate term 3
		double[][] term3 = new double[group.size()][branchNum];
		for (int g = 0; g < group.size(); g++) {
			for (int j = 0; j < branchNum; j++) {
				double sum = 0;
				for (int c = 0; c < C; c++) {
					sum += (Ngjc[g][j][c]/mu[j][g] - cphi[g][j][c]);
				}
				term3[g][j] = sum/C;
			}
		}
		for (int g = 0; g < group.size(); g++) {
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < branchNum; l++) {
					info[g][j][l] += term3[g][j] * term3[g][l];
				}
			}
		}
		
		/*print info
		for (int g = 0; g < group.size(); g++) {
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < branchNum; l++) {
					System.out.println(info[g][j][l]);
				}
			}
		}*/
		
		return info;
	}
	
	
	/**
	 * This function calls calculateInformation() and inverse the information matrix.
	 */
	private void setCov(double[][][] info) {
		for (int k = 0; k < info.length; k++) {
			RealMatrix a = new Array2DRowRealMatrix(info[k]);
			a = MatrixUtils.inverse(a);	//inverse matrix
			cov[k] = a.getData();
		}		
	}
	
	
	
	
	

	
	
	/** 
	 * Print covariance matrix for type g substitution
	 * @param g index of the type of substitution
	 * @return covariance matrix
	 */
	public String printCovarianceMatrix(int g) {
		DecimalFormat formatter = new DecimalFormat("#.#################");
		String s = "variance-covariance matrix follows:\n";
		double[][] covMatrix = cov[g];
		for (int i = 0; i < branchNum; i ++) {
			for (int j = 0; j < branchNum; j ++) {
				s += formatter.format(covMatrix[i][j]) + " ";
			}
			s += "\n";
		}
		
		return s;
	}
	
	
    /** 
     * Create output files: o.estb.type, o.estb.group and substitutionLength.txt
     * o.estb.type and o.estb.group can be used as the inputs to Multidivtime.
     *
     * o.estb.type for single nucleotide substitutions are for 12 different substitution types
     * o.estb.group for single nucleotide substitutions are for 6 different substitution types (combine strand symmetric substitution types)
     *
     * o.estb.type for triple site nucleotide substitutions are for 4 different substitution types (non-CpG/CpG transversion/transition)
     * o.estb.group for triplet site nucleotide substitutions are for 9 different substitution types (defined in manuscript)
     *
     * substitutionLength.txt file stores the estimated substitution lengths for each strand symmetric substitution type.
     * This makes it easier to report results.
     **/
	public void printOutput() {
		for (int g = 0; g < N_G; g++) {
			File output;
			output = new File("o.estb.type"+g);
			//produce mu_g: double[branchNum], mu = new double[branchNum][N_G];
			double[] mu_g = new double[branchNum];
			for (int j = 0; j < branchNum; j++) {
				mu_g[j] = mu[j][g];
			}
			String newick = getNewickTree(mu_g);
			PrintStream print = null;
			
		  	  try{
		  		print = new PrintStream(output);	 
		  		//print tree in newick format with branch lengths averaged over iterations
			  	  print.print(newick);
			  	  print.println();
			  	  print.print(printNodeInfo());
			  	  print.print(printCovarianceMatrix(g));
	  	  	  } catch(FileNotFoundException e){
	  	  	  	  System.out.println("Problem creating file!");	  
	  	  	  } finally {
	  	        if (print != null) print.close();
	  	    }
		  	  
		}
				
		File output;
		output = new File("substitutionLength.txt");
		PrintStream print = null;
		try {
			print = new PrintStream(output);
			// print substitution lengths for branch 0 to the last branch to file
			for (int g = 0; g < N_G; g++) {
				print.println("Substitution length type "+ g);
				for (int j = 0; j < branchNum; j++) {
					print.println(mu[j][g]);
				}
			}
		} catch (FileNotFoundException e) {
			System.out.println("Problem creating substitution length file!");
		} finally {
  	        if (print != null) print.close();
  	    }
			
	}
	
    /** 
     * This function print some output to screen. 
     * It prints the number of substitutions per type per branch, the proportion of time of each type per branch,
     * and the sum of mu for each type per branch to screen.
     **/
	public void printToScreen(){
		for (int k = 0; k < numTypeChanges; k++) {
			for (int j = 0; j < branchNum; j++) {
				int count = 0;
				for (int c = 0; c < C; c++) {
					//numberOfChanges = new int[branchNum][N_G][C];
					count += numberOfChanges[j][k][c];
				}
				System.out.println("Type "+k +" changes on branch "+j+": "+count);
			}
		}
		
		
		for (int k = 0; k < numTypeStates; k++) {
			for (int j = 0; j < branchNum; j++) {
				double prop = 0.0;
				for (int c = 0; c < C; c++) {
					//propStates = double[branchNum][numState][C];
					prop += propStates[j][k][c];
				}
				System.out.println("Proportion of time in state "+k+" on branch "+j+": "+prop);
			}
		}
		//mu = new double[branchNum][N_G];
		
		System.out.println("=============================================");
		for (int g = 0; g < N_G; g++) {
			for (int j  = 0; j < branchNum; j++) {
				System.out.println("Type "+g+" branch "+j+" mu "+ mu[j][g]);
			}
		}
		
		// c_k = new double[s.getNumberOfTypes()];
		for (int g = 0; g < numTypeChanges; g++) {
			System.out.println("C_k for "+g+" is "+ c_k[g]);
		}
		
	}

	/** 
	 * Store newick tree with mu for a specific type g
	 * @param mu_g, double[branchNum], mu of a specific type g
	 */
	public String getNewickTree(double[] mu_g) {	
		String newick = inOrderNewick(tree.root, mu_g);		
		
		return newick;
	}
	
	
	/** 
	 * Recursive print newick format
	 * @param root cursor to traverse the tree
	 * @return newick tree
	 */
	public String inOrderNewick(TreeNode root, double[] mu_g) {
		if (!root.isLeaf()) {
	          String output = "";
	          output += "(";
	          output += inOrderNewick(root.firstChild(), mu_g);
	          output += ",";
	          output += inOrderNewick(root.lastChild(), mu_g);
	          output += ")";
	          if (root.isRoot()) {
	        	  output += ";";
	          } else {
	        	output += ":" + mu_g[root.getNodeNum()];  
	          }
	          return output;
	      } else {
	    	  
	          return root.getName() + ":" + mu_g[root.getNodeNum()];
	      }
	}
	
	/** 
	 * Print node information. Used for Multidivtime input.
	 * 1. list of names and tip numbers
	 * 2. list of child1, child2, ..., parent
	 * @return String containing the info
	 */
	public String printNodeInfo() {
		String info = "list of names and tip numbers follows:\n";
		for (int i = 0; i < tree.getNumLeaves(); i++) {
			info += tree.getNodeByNodeNum(i).getName();
			info += "  " + i + "\n";
		}
		
		info += "list of child1, child2, ..., parent follows:\n";
		for (int i = tree.getNumLeaves(); i < tree.getTotalNodeCount(); i++) {
			info += tree.getNodeByNodeNum(i).getChild(0).getNodeNum() + " ";
			info += tree.getNodeByNodeNum(i).getChild(1).getNodeNum() + " ";
			info += i + "\n";
		}
		return info;
	}
	

}
