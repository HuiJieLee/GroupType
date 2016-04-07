package model_selection;

import java.util.*;

/**
 * This class combined number of substitution types and proportion times in states based on whether it is context-dependent.
 * @author Hui-Jie
 *
 */
public class SufficientStatistics {
    /** Store proportion of states*/
	private double[][][] propStates;
	/** Store number of changes */
	private int[][][] numberOfChanges;
	/** Number of types, GTR = 12, CpG = 9 */
	private int N_G;
	/** Number of states, GTR = 4, CpG = 3 */
	private int numState;
	/** Number of branches */
	private int branchNum;
	/** Number of iterations */
	private int C;
	/** Store corresponding state b(k) for each type k*/
	private Hashtable<Integer, Integer> map;
	
	
	/**
	 * Constructor for GTR model or 9 context-dependent substitution types
	 * @param isGTR
	 * @param branchNum
	 * @param C
	 * @param propStates	
	 * @param numberOfChanges
	 */
	public SufficientStatistics(boolean isGTR, int branchNum, int C, double[][][] propStates, int[][][] numberOfChanges) {
		this.branchNum = branchNum;
		this.C = C;
		if(isGTR) {
			this.numState = 4;
			this.N_G = 12;
			this.propStates = new double[branchNum][numState][C]; //4 possible states, A, T, C, G
			this.numberOfChanges = new int[branchNum][N_G][C];
			setGTR(propStates, numberOfChanges);
			createHashTable(true);
		} else {
			this.numState = 3;
			this.N_G = 9;
			this.propStates = new double[branchNum][numState][C]; //3 possible states, non-CpG G+C site, non-CpG A+T site, CpG C+G site 
			this.numberOfChanges = new int[branchNum][N_G][C];
			setCpG(propStates, numberOfChanges);
			createHashTable(false);
			/*print sufficient stats
			System.out.println("Prop states: ");
			for (int j = 0; j < branchNum; j++) {
				System.out.println("j = " + j);
				for (int l = 0; l < numState; l++) {
					System.out.println("l = " + l);
					for (int c = 0; c < C; c++) {
						System.out.println("Phi_jb(k)c = " + propStates[j][l][c]);
					}
				}
			}
			
			System.out.println("Number of changes: ");
			for (int j = 0; j < branchNum; j++) {
				System.out.println("j = " + j);
				for (int k = 0; k < N_G; k++) {
					System.out.println("k = " + k);
					for (int c = 0; c < C; c++) {
						System.out.println("c = " + c);
						System.out.println("N_jkc = " + numberOfChanges[j][k][c]);
					}
				}
			}*/
			
		}
		
	}
	
	/**
	 * Key = type k
	 * Value = state b(k)
	 * @param isGTR
	 */
	private void createHashTable(boolean isGTR) {
		map = new Hashtable<Integer, Integer>();
		if (isGTR) {
			map.put(0, 0);
			map.put(1, 0);
			map.put(2, 0);
			map.put(3, 1);
			map.put(4, 1);
			map.put(5, 1);
			map.put(6, 2);
			map.put(7, 2);
			map.put(8, 2);
			map.put(9, 3);
			map.put(10, 3);
			map.put(11, 3);			
		} else {
			map.put(0, 0);
			map.put(1, 0);
			map.put(2, 1);
			map.put(3, 1);
			map.put(4, 0);
			map.put(5, 1);
			map.put(6, 2);
			map.put(7, 2);
			map.put(8, 2);
		}
	}
	
	/**
	 * Set propStates and numberOfChanges for GTR model
	 * @param propStates
	 * @param numberOfChanges
	 */
	private void setGTR(double[][][] propStates, int[][][] numberOfChanges) {
		for(int j = 0; j < branchNum; j++) {
			for (int c = 0; c < C; c++) {
				//set propStates
				this.propStates[j][0][c] = propStates[j][0][c] + propStates[j][4][c]; //G
				this.propStates[j][1][c] = propStates[j][1][c] + propStates[j][5][c]; //C
				this.propStates[j][2][c] = propStates[j][2][c];	//T
				this.propStates[j][3][c] = propStates[j][3][c];	//A
				
				//set numberOfChanges
				this.numberOfChanges[j][0][c] = numberOfChanges[j][0][c] + numberOfChanges[j][12][c]; //G->C
				this.numberOfChanges[j][1][c] = numberOfChanges[j][2][c] + numberOfChanges[j][14][c]; //G->T
				this.numberOfChanges[j][2][c] = numberOfChanges[j][8][c] + numberOfChanges[j][16][c]; //G->A
				this.numberOfChanges[j][3][c] = numberOfChanges[j][1][c] + numberOfChanges[j][13][c]; //C->G
				this.numberOfChanges[j][4][c] = numberOfChanges[j][9][c] + numberOfChanges[j][17][c]; //C->T
				this.numberOfChanges[j][5][c] = numberOfChanges[j][3][c] + numberOfChanges[j][15][c]; //C->A
				this.numberOfChanges[j][6][c] = numberOfChanges[j][6][c];  //T->G
				this.numberOfChanges[j][7][c] = numberOfChanges[j][11][c]; //T->C
				this.numberOfChanges[j][8][c] = numberOfChanges[j][4][c];  //T->A
				this.numberOfChanges[j][9][c] = numberOfChanges[j][10][c]; //A->G
				this.numberOfChanges[j][10][c] = numberOfChanges[j][7][c]; //A->C
				this.numberOfChanges[j][11][c] = numberOfChanges[j][5][c]; //A->T		
			}
		}
	}
	
	/**
	 * Set propStates and numberOfChanges for 9 context-dependent substitution types (without grouping)
	 * @param propStates
	 * @param numberOfChanges
	 */
	private void setCpG(double[][][] propStates, int[][][] numberOfChanges) {
		for (int j = 0; j < branchNum; j++) {
			for (int c = 0; c < C; c++) {
				//set propStates
				this.propStates[j][0][c] = propStates[j][0][c] + propStates[j][1][c]; //non-CpG C+G
				this.propStates[j][1][c] = propStates[j][2][c] + propStates[j][3][c]; //non-CpG A+T
				this.propStates[j][2][c] = propStates[j][4][c] + propStates[j][4][c]; //CpG
				//set numberOfChanges
				this.numberOfChanges[j][0][c] = numberOfChanges[j][0][c] + numberOfChanges[j][1][c]; //non-CpG G->C, C->G
				this.numberOfChanges[j][1][c] = numberOfChanges[j][2][c] + numberOfChanges[j][3][c]; //non-CpG G->T, C->A
				this.numberOfChanges[j][2][c] = numberOfChanges[j][4][c] + numberOfChanges[j][5][c]; //non-CpG T->A, A->T
				this.numberOfChanges[j][3][c] = numberOfChanges[j][6][c] + numberOfChanges[j][7][c]; //non-CpG T->G, A->C
				this.numberOfChanges[j][4][c] = numberOfChanges[j][8][c] + numberOfChanges[j][9][c]; //non-CpG G->A, C->T
				this.numberOfChanges[j][5][c] = numberOfChanges[j][10][c] + numberOfChanges[j][11][c]; //non-CpG A->G, T->C
				this.numberOfChanges[j][6][c] = numberOfChanges[j][12][c] + numberOfChanges[j][13][c]; //CpG G->C, C->G
				this.numberOfChanges[j][7][c] = numberOfChanges[j][14][c] + numberOfChanges[j][15][c]; //CpG G->T, C->A
				this.numberOfChanges[j][8][c] = numberOfChanges[j][16][c] + numberOfChanges[j][17][c]; //CpG G->A, C->T
			}
		}
	}
	
	/**
	 * Get propStates
	 * @return propStates
	 */
	public double[][][] getPropStates() {
		return propStates;
	}
	
	/**
	 * Get numberOfChanges
	 * @return numberOfChanges
	 */
	public int[][][] getNumberOfChanges() {
		return numberOfChanges;
	}
	
	/**
	 * Get number of types
	 * @return
	 */
	public int getNumberOfTypes() {
		return N_G;
	}
	
	/**
	 * Get number of states
	 * @return
	 */
	public int getNumberOfStates() {
		return numState;
	}

	/**
	 * Return C
	 * @return
	 */
	public int getC() {
		return C;
	}
	
	/**
	 * Return number of branches
	 * @return
	 */
	public int getNumberOfBranches() {
		return branchNum;
	}
	
	/**
	 * Return map: key = type, value = starting state
	 * @return map
	 */
	public Hashtable<Integer, Integer> getStartingState() {
		return map;
	}
}
