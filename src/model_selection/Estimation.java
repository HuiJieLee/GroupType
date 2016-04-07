/**
 * 
 */
package model_selection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;


/**
 * @author Hui-Jie Lee
 *
 */
public class Estimation {
	
	/** Number of iterations / trees in a file*/
	private int C;
	/** Number of subsets */
	private int N_G;
    /** Number of branches */
	private int branchNum;
	/** Number of states */
	private int numStates;
	/** Store number of changes */
	private int[][][] numberOfChanges;
	/** Store proportion of states*/
	private double[][][] propStates;
	/** Store mu_kj */
	private double[][][] mu;
	/** Store c_k */
	private double[][] c_k;
	/** Store corresponding state b(k) for each type k*/
	private Hashtable<Integer, Integer> map;
	/** Store corresponding subset number for each type k */
	private Hashtable<Integer, Integer> subset;

	/**
	 * Constructor for 9 context-dependent substitution types (ungrouped)
	 * @param s
	 * @param w_GTR
	 */
	public Estimation(SufficientStatistics s) {
		this.C = s.getC();
		this.N_G = s.getNumberOfTypes();
		this.branchNum = s.getNumberOfBranches();
		this.numStates = s.getNumberOfStates();
		this.numberOfChanges = s.getNumberOfChanges();
		this.propStates = s.getPropStates();
		this.mu = new double[branchNum][N_G][C];
		this.c_k = new double[N_G][C];
		for (double[] row: c_k) {
			Arrays.fill(row, 1.0); //c_k = 1 for ungrouped substitution types
		}
		this.map = s.getStartingState();
		createSubsetTable();
		initialize();
		estimatesPerMap();
		/*
		System.out.println("Print mu:");
		//print mu
		for(int j = 0; j < branchNum; j++) {
			for (int k = 0; k < N_G; k++) {
				System.out.println("Mu " + k + " on branch " + j + " = " + mu[j][k]);
			}
		}
		
		System.out.println("=============================================");
		*/
	}
	
	/**
	 * Constructor for grouped context-dependent substitution types
	 * @param group store grouping info, each ArrayList contains types in one subset
	 * @param s
	 * @param w_GTR
	 */
	public Estimation(ArrayList<ArrayList<Integer>> group, SufficientStatistics s){
		this.C = s.getC();
		this.N_G = group.size(); //number of subsets
		this.branchNum = s.getNumberOfBranches();
		this.numStates = s.getNumberOfStates();
		this.numberOfChanges = s.getNumberOfChanges();
		this.propStates = s.getPropStates();
		this.mu = new double[branchNum][N_G][C];
		for (double[][] square : mu) {
	        for (double[] line : square) {
	            Arrays.fill(line, 0.2); //initialize with 0.2
	        }
	    }
		
		this.c_k = new double[s.getNumberOfTypes()][C]; // dim = 9 x 1
		for (double[] row: c_k) {
			Arrays.fill(row, 0.2); //c_k = 1 for ungrouped substitution types
		}
		for(int i = 0 ; i < group.size() ; i++) {
			ArrayList<Integer> currentList = group.get(i);
		    //now reorder the list if necessary
			reoderList(currentList, s.getNumberOfChanges());
			int index = currentList.get(0);
			Arrays.fill(c_k[index], 1.0);	//set c_k for the first type in each subset = 1.		
		}		
		this.map = s.getStartingState();
		createSubsetTable(group);
		System.out.println("Group!");
		initialize(group);

		/*
		System.out.println("Print mu_g:");
		//print mu
		for(int j = 0; j < branchNum; j++) {
			for (int k = 0; k < N_G; k++) {
				System.out.println("Mu " + k + " on branch " + j + " = " + mu[j][k]);
			}
		}
		
		System.out.println("Print c_k:");
		for (int k = 0; k < s.getNumberOfTypes(); k++) {
			System.out.println("c_" + k +  " = " + c_k[k]);
		}
		
		System.out.println("=============================================");
		*/
	}
	
	/**
	 * This method checks the currentList and reorder the list if needed.
	 * Case 1: Group has only one type or the first type has count > 0. =>no need to reorder.
	 * Case 2: Group has more than one type and every type in this group has count 0. => no need to reorder.
	 * Case 3: Group has more than one type and there is at least one type (not the first one) has count > 0.
	 * 			=> swap first element with the non-zero count element.
	 * @param currentList
	 * @param numberOfChanges
	 * @return
	 */
	private ArrayList<Integer> reoderList(ArrayList<Integer> currentList, int[][][] count) {
		int type = currentList.get(0);
		if (currentList.size() > 1 && isZeroCount(count, type)) {
			// if all types have zero count, no need to swap
			for (int i = 1; i < currentList.size(); i++) {
				int temp = currentList.get(i);
				if (!isZeroCount(count, temp)) {
					Collections.swap(currentList, 0, i);
					return currentList;
				}
			}
		}
		return currentList;
	}
	
	/**
	 * Check whether numberOfChanges summing over all branches and all iterations = 0.
	 * If so, then this type cannot be the first type in a group.
	 * If every type in one group has zero count, then keep the original order.
	 * In this case, mu_g(k)j = 0 and c_k = 1.
	 * @param count new int[branchNum][N_G][C];
	 * @param type
	 * @return
	 */
	private boolean isZeroCount(int[][][] count, int type) {
		int sum = 0;
		for (int j = 0; j < branchNum; j++) {
			for (int c = 0; c < C; c++) {
				sum += count[j][type][c];
			}
		}
		
		if (sum == 0) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * For ungrouped case, every type is in Subset 0.
	 */
	private void createSubsetTable() {
		subset = new Hashtable<Integer, Integer>();
		for (int k = 0; k < N_G; k++) {
			subset.put(k, 0);
		}
	}
	
	/**
	 * For grouped case
	 * key = type
	 * value = subset number
	 * @param group
	 */
	private void createSubsetTable(ArrayList<ArrayList<Integer>> group) {
		subset = new Hashtable<Integer, Integer>();
		for (int i = 0; i < group.size(); i++) {
			ArrayList<Integer> currentList = group.get(i);
			for (int j = 0; j < currentList.size(); j++) {
				subset.put(currentList.get(j), i);
			}
		}
	}
	
	/**
	 * Step 1: initialize mu_{kj}^{(0)} for ungrouped
	 * 	this.mu = new double[branchNum][N_G][C];
	 *	this.c_k = new double[N_G][C];
	 */
	private void initialize() {
		for (int j = 0; j < branchNum; j++) {
			for (int k = 0; k < N_G; k++) {
				for (int c = 0; c < C; c++) {
					int temp1 = numberOfChanges[j][k][c];
					double temp2 = propStates[j][map.get(k)][c];
					if(temp1 != 0 && temp2 != 0) {
						mu[j][k][c] = numberOfChanges[j][k][c]/propStates[j][map.get(k)][c];
					} else {
						mu[j][k][c] = 0;
					}
				}
			}
		}
		
	}
	
	private void estimatesPerMap() {		
		File output;
		PrintStream print = null;
		for (int k = 0; k < N_G; k++) {
			for (int j = 0; j < branchNum; j++) {
				output = new File("Mu"+k+j+".txt");
				try{
					print = new PrintStream(output);
					for (int c = 0; c < C; c++) {
						print.print(mu[j][k][c]);
						print.print(" ");
					}
				} catch (FileNotFoundException e) {
					System.out.println("Problem creating file!");
				} finally {
					if (print != null) print.close();
				}
			}
			
		}
		
	}
	
	
	/**
	 * Step 1: initialize mu_{g(k)j}^{(0)} for grouped
	 * @param group
	 */
	private void initialize(ArrayList<ArrayList<Integer>> group) {		
		Object[] obj = successiveSubstitution(c_k, mu, group, numberOfChanges, propStates);
		c_k = (double[][]) obj[0];
		mu = (double[][][]) obj[1];
		
	}
	
	/**
	 * Find max value in a 2D array
	 * @param array
	 * @return max
	 */
	public static double findMax(double[][] array) {
		double max = array[0][0];
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++){
				if (array[i][j] > max) {
					max = array[i][j];
				}
			}
		}
		return max;
	}
	
	/**
	 * Find max value in an array
	 * @param array
	 * @return max
	 */
	public static double findMax(double[] array) {
		double max = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
			}
		}
		return max;
	}
	
	/**
	 * Successive substitutions solving for c_k and mu for the first time 
	 * 		this.mu = new double[branchNum][N_G][C];
	 *	this.c_k = new double[N_G][C];
	 * @param c_k 
	 * @param mu
	 * @param group
	 * @param n_kj
	 * @param phi_bkj
	 * @return
	 */
	public Object[] successiveSubstitution(double[][] c_k_in, double[][][] mu_in, ArrayList<ArrayList<Integer>> group, int[][][] n_kj, double[][][] phi_bkj) {
		double[][] mu_diff = new double[branchNum][N_G];
		double[][] mu_temp = new double[branchNum][N_G];
		
		for (int c = 0; c < C; c++) {
			System.out.println(c);
			for (int j = 0; j < branchNum; j++) {
				for (int k = 0; k < N_G; k++) {
					mu_diff[j][k] = mu_in[j][k][c];
					mu_temp[j][k] = mu_in[j][k][c];
				}
			}
			
			double[] c_diff = new double[c_k.length];	//initialize with all 0s. c_k = new double[N_G][C];
			//store mu_g(k)j^{(r)}
			double[] c_temp = new double[c_k.length];
			for (int k = 0; k < c_k.length; k++) {
				c_temp[k] = c_k_in[k][c];
			}
			
			int i = 0;
			while (findMax(c_diff) > 0.001 || findMax(mu_diff) > 0.001) {
				System.out.println("iteration" + i);
				i++;
				// update mu_g(k)j first
				for (int g = 0; g < group.size(); g++) {
					ArrayList<Integer> currentList = group.get(g);
					for (int j = 0; j < branchNum; j++) {
						int numerator = 0;
						double denominator = 0;		
						//now iterate on the current list
						for (int k = 0; k < currentList.size(); k++) {
							int index = currentList.get(k);
							numerator += n_kj[j][index][c];
							denominator += phi_bkj[j][map.get(index)][c] * c_k_in[index][c];
						}
						if (numerator != 0 && denominator != 0) {
							mu_in[j][g][c] = numerator / denominator;
						} else {
							mu_in[j][g][c] = 0;
						}
						mu_diff[j][g] = Math.abs(mu_temp[j][g] - mu_in[j][g][c]);
						mu_temp[j][g] = mu_in[j][g][c];
					}			    
				}
				
				//update c_k 
				for (int g = 0; g < group.size(); g++) {
					ArrayList<Integer> currentList = group.get(g);
					for (int l = 1; l < currentList.size(); l++) {	// don't need to update first c_k in each subset since it is 1.
						int k = currentList.get(l); //get type index
						int numerator = 0;
						double denominator = 0;
						for (int j = 0; j < branchNum; j ++) {
							numerator += n_kj[j][k][c];
							denominator += phi_bkj[j][map.get(k)][c] * mu_in[j][g][c];
						}
						if (numerator != 0 && denominator != 0) {
							c_k_in[k][c] = numerator/denominator;		
						} 
						c_diff[k] = Math.abs(c_temp[k]-c_k_in[k][c]);
						c_temp[k] = c_k_in[k][c];
						
					}
				}
							
			}
		}
		
		
		
		
		return new Object[]{c_k_in, mu_in};
	}
	
	
	/**
	 * Calculate log P(M^{(c)}, X|mu^{(r)}) for each mapping at iteration r
	 * @param mu double[branchNum][N_G];
	 * @return log weight (numerator)
	 */
	private double[] calculateLogW(double[][] mu) {
		double[] w_c = new double[C]; 
		for (int c = 0; c < C; c++) {
			double sum = 0;
			for (int j = 0; j < branchNum; j++) { 
				for (int k = 0; k < 9; k++) { //numberOfChanges [branchNum][numTypeChanges][C]
					if(mu[j][k] == 0 || numberOfChanges[j][k][c] == 0 || propStates[j][map.get(k)][c] == 0) {
						sum += 0;
					} else {
						sum += numberOfChanges[j][k][c] * Math.log(mu[j][k]); //N_kj * log (mu_kj)
						sum -= mu[j][k] * propStates[j][map.get(k)][c]; // - mu_kj * phi_{b(k)j}
					}
				}
			}
			w_c[c] = sum;
		}
		
		return w_c;
	}
	
	/**
	 * Recover mu_kj given mu_g(k)j and c_k. mu_kj = c_k * mu_g(k)j.
	 * @param mu_gkj [branchNum][N_G] where N_G = # of subsets
	 * @param c_k [9]
	 * @return mu_kj
	 */
	private double[][] convertMu(double[][] mu_gkj, double[] c_k) {
		double[][] mu_k = new double[branchNum][9];
		for (int j = 0; j < branchNum; j++) {
			for (int k = 0; k < 9; k++) {
				mu_k[j][k] = c_k[k] * mu_gkj[j][subset.get(k)];
			}
		}
		
		return mu_k;
	}
	
	/**
	 * Return mu_kj or mu_g(k)j 
	 * @return
	 */
	public double[][][] getMu() {
		return mu;
	}
	
	/**
	 * Return c_k
	 * @return
	 */
	public double[][] getCk() {
		return c_k;
	}
	
	
}
