package model_selection;

import java.util.ArrayList;
import java.util.Hashtable;

/**
 * This class store the sufficient statistics for GTR model and the parameter values mu_k.
 * It also calculates the log P(M^{(c)}, X|mu_GTR) for each mapping.
 * @author Hui-Jie Lee
 *
 */
public class GTR {
	/** Store mu_k */
	private double[] mu_k;
	/** Store rate away from state*/
	private double[] R_k;
	/** Store sufficient statistics n_k */
	private int[][][] numberOfChanges;
	/** Store sufficient statistics phi_b(k) */
	private double[][][] propStates;
	/** Store sufficient statistics number of states at root */
	private int[][] rootState;
	/** Store log P(M^{(c)}, X|mu_GTR) */
	private double[] logL;
	/** Store branch length */
	private double[] branchLength;	
	/** Number of mappings */
	private int C;
	/** Store corresponding state b(k) for each type k*/
	private Hashtable<Integer, Integer> map;
	
	/**
	 * Constructor
	 * @param pi pi_G, pi_C, pi_T, pi_A
	 * @param R R_GC/R_CG, R_GT/R_TG, R_GA/R_AG, R_CT/R_TC, R_CA/R_AC, R_TA/R_AT
	 * @param propStates double[branchNum][numState][C] 4 possible states, G, C, T, A
	 * @param numberOfChanges int[branchNum][N_G][C] N_G = 12, GC, GT, GA, CG, CT, CA, TG, TC, TA, AG, AC, AT
	 * @param rootStateCount state count at root int[4][C] G,C,T,A
	 * @param branch length double[branchNum]
	 * @param C number of mappings
	 * @param branchNum number of branches
	 */
	public GTR(double[] pi, double[] R, SufficientStatistics suff, int[][] rootStateCount, double[] branchLength) {
		this.C = suff.getC();
		this.mu_k = new double[12];
		this.numberOfChanges = suff.getNumberOfChanges();
		this.propStates = suff.getPropStates();
		this.setRootState(rootStateCount);
		this.branchLength = branchLength;
		this.logL = new double[C];
		this.map = suff.getStartingState();
		this.R_k = new double[4];
		
		setMu_k(pi, R);
		setR_k();
		calculateLogL(pi, rootStateCount);
		
	}
	
	
	/**
	 * Set mu_k
	 * @param pi
	 * @param R
	 */
	private void setMu_k(double[] pi, double[] R) {
		mu_k[0] = R[0]*pi[1]; //G->C = mu_1 * pi_C
		mu_k[1] = R[1]*pi[2]; //G->T = mu_2 * pi_T
		mu_k[2] = R[2]*pi[3]; //G->A = mu_3 * pi_A
		mu_k[3] = R[0]*pi[0]; //C->G = mu_1 * pi_G
		mu_k[4] = R[3]*pi[2]; //C->T = mu_4 * pi_T
		mu_k[5] = R[4]*pi[3]; //C->A = mu_5 * pi_A
		mu_k[6] = R[1]*pi[0]; //T->G = mu_2 * pi_G
		mu_k[7] = R[3]*pi[1]; //T->C = mu_4 * pi_C
		mu_k[8] = R[5]*pi[3]; //T->A = mu_6 * pi_A
		mu_k[9] = R[1]*pi[0]; //A->G = mu_2 * pi_G
		mu_k[10] = R[4]*pi[1];//A->C = mu_5 * pi_C
		mu_k[11] = R[5]*pi[2];//A->T = mu_6 * pi_T
	}
	
	private void setR_k(){
		for (int i = 0; i < 4; i++) {
			R_k[i] = mu_k[i*3] + mu_k[i*3+1] + mu_k[i*3+2];
		}		
	}
	
	
	/**
	 * This method calculates log P(M^{(c)}, X|mu_GTR) for each mapping
	 * @param pi
	 * @param rootStateCount
	 */
	private void calculateLogL(double[] pi, int[][] rootStateCount) {
		
		for (int c = 0; c < C; c++) {
			double sum = 0;
			//root state
			for (int l = 0; l < 4; l++) {
				sum += rootStateCount[l][c] * Math.log(pi[l]);
			}
			//branch 
			for (int j = 0; j < branchLength.length; j++) {
				for (int k = 0; k < 12; k++) { //numberOfChanges [branchNum][numTypeChanges][C]
					if(numberOfChanges[j][k][c] == 0 || propStates[j][map.get(k)][c] == 0) {
						sum += 0;
					} else {
						sum += numberOfChanges[j][k][c] * Math.log(mu_k[k]*branchLength[j]); //N_kj * log (mu_kj)
						sum -= (mu_k[k]*branchLength[j]) * propStates[j][map.get(k)][c]; // - mu_kj * phi_{b(k)j}
						
						//sum += -propStates[j][map.get(k)][c]*R_k[map.get(k)]*branchLength[j];
						
					}
				}
			}
			logL[c] = sum;
			System.out.println("GTR weight logL[" + c + "] = " + sum);
			
		}
		/*
		for (int c = 0; c < C; c++) {
			double sum = 0;
			//root state 
			for (int l = 0; l < 4; l++) {
				sum += rootStateCount[l][c] * Math.log(pi[l]);
			}
					
			for (int k = 0; k < 12; k++) {
				sum += n_k[k][c] * Math.log(mu_k[k]);  // sum N_k * log(mu_k) over all k
			}
			//state G: -3 * (mu_1 + mu_2 + mu_3) * phi_G
			sum += -3.0 * (mu_k[0] + mu_k[1] + mu_k[2]) * phi_bk[0][c];
			//state C: -3 * (mu_4 + mu_5 + mu_6) * phi_C
			sum += -3.0 * (mu_k[3] + mu_k[4] + mu_k[5]) * phi_bk[1][c];
			//state T: -3 * (mu_7 + mu_8 + mu_9) * phi_T
			sum += -3.0 * (mu_k[6] + mu_k[7] + mu_k[8]) * phi_bk[2][c];
			//state A: -3 * (mu_10 + mu_11 + mu_12) * phi_A
			sum += -3.0 * (mu_k[9] + mu_k[10] + mu_k[11]) * phi_bk[3][c];
			
			logL[c] = sum;
			System.out.println("GTR weight logL[" + c + "] = " + sum);
		}
		*/
	}
	
	/**
	 * Return log P(M^{(c)}, X|mu_GTR) for each mapping
	 * @return logL
	 */
	public double[] getLogL() {
		return logL;
	}

	public int[][] getRootState() {
		return rootState;
	}

	public void setRootState(int[][] rootState) {
		this.rootState = rootState;
	}
}
