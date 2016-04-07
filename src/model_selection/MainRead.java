package model_selection;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.ArrayList;


public class MainRead {
	
	/** Number of iterations / trees in a file*/
	private int C;
	/** Number of branches */
	private int branchNum;
	/** Store proportion of states*/
	private double[][][] propStates;
	/** Store number of changes */
	private int[][][] numberOfChanges;

	
	
	public MainRead(int C, int branchNum) {
		this.C = C;
		this.branchNum = branchNum;
	}
	
	private double[][][] getProp(String State, int branchNum, int C) {
		double[][][] prop = new double[branchNum][6][C];
		
		InputStream inStream = this.getClass().getResourceAsStream(new File("../" + State).getPath().toString());
    	BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
    	String line;
		try {
			for (int j = 0; j < branchNum; j++) {
				for (int k = 0; k < 6; k++) {
					line = r.readLine();
					String[] stringArray = line.split(" ");
					for (int c = 0; c < C; c++) {
						prop[j][k][c] = Double.parseDouble(stringArray[c]);
					}
				}
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
            try {
                inStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
		
		return prop;
	}
	
	private int[][][] getCount(String numChange, int branchNum, int C) {
		InputStream inStream = this.getClass().getResourceAsStream(new File("../" + numChange).getPath().toString());
    	BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
    	String line;
    	int[][][] count = new int[branchNum][18][C];
    	    	
		try {
			for (int j = 0; j < branchNum; j++) {
				for (int k = 0; k < 18; k++) {
					line = r.readLine();
					String[] stringArray = line.split(" ");
					for (int c = 0; c < C; c++) {
						count[j][k][c] = Integer.parseInt(stringArray[c]);
					}
				}
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
            try {
                inStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
		
		return count;
	}
	

	
	private Tree getTree(String outgroupFile, String newick) {
		StringReader str = new StringReader(newick);
		TreeParser tp = new TreeParser(str, outgroupFile);
		Tree tree = tp.tokenize(); 
		
		return tree; 
	}
	
	/**
	 * This function parse the file that stores substitution types grouping info.
	 * @param filename
	 * @return an arraylist of arraylist storing grouping info. each list contains a list of elements in that subset.
	 */
	public ArrayList<ArrayList<Integer>> grouping(String filename){
	    ArrayList<ArrayList<Integer>> outer = new ArrayList<ArrayList<Integer>>();
	    ArrayList<Integer> inner = new ArrayList<Integer>();   
	    int N_G;
		    
	    InputStream inStream = this.getClass().getResourceAsStream(new File("../" + filename).getPath().toString());
    	BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
    	String line;
		try {
			line = r.readLine();
			//first line contains the number of subsets
             N_G = Integer.parseInt(line);             
             for (int i = 0; i < N_G; i++) {
             	line = r.readLine(); // each line contains numbers in that subset
             	String[] elements = line.split(" ");
             	for (int j = 0; j < elements.length; j++) {
             		int element = Integer.parseInt(elements[j]);
             		inner.add(element-1);	//index = element - 1, i.e. type = 0, ..., 8 instead of 1,..., 9
             	}
             	outer.add(inner);
             	inner = new ArrayList<Integer>();
             }
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
            try {
                inStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
	
	    return outer;
	}

	/**
	 * Read sufficient statistics from files
	 * @param args[0] NumChange file name
	 * @param args[1] PropState file name
	 * @param args[2] number of branches
	 * @param args[3] number of iterations
	 * @param args[4] outgroupFileName
	 * @param args[5] grouping files
	 */
	public static void main(String[] args) {
		String numChange = args[0];
		String State = args[1];
		int branchNum = Integer.parseInt(args[2]);
		int C = Integer.parseInt(args[3]);
		String outgroup = args[4];
		
		if (args.length == 5) {
			
			MainRead main = new MainRead(C, branchNum);
	    	double[][][] prop = new double[branchNum][6][C];
	    	int[][][] count = new int[branchNum][18][C];
			prop = main.getProp(State, branchNum, C);
			count = main.getCount(numChange, branchNum, C);
			
			String newick = "((((human_G:0.017228:G,orangutan_G:0.018799:G)_G:0.002182:G,gibbon_G:0.021454:G)_G:0.011141:G,((rhesus_G:0.008145:G,baboon_G:0.008297:G)_G:0.00478:G,greenMonkey_G:0.013339:G)_G:0.027207:G)_G:0.017203:G,(marmoset_G:0.034144:G,squirrelMonkey_G:0.033447:G)_G:0.04928:G,bushbaby_G:0.244601:G)_G;";
			Tree tree = main.getTree(outgroup, newick); 
            
			SufficientStatistics contextdept = new SufficientStatistics(false, branchNum, C, prop, count);		
			int[][][] num = contextdept.getNumberOfChanges();

			Estimation mc = new Estimation(contextdept);
			// call VarCov to estimate var-cov matrix of mu_kj, output for Multidivtime
			//VarCov vc = new VarCov(contextdept, mc.getMu(), mc.getCk(), tree);
			MCVar vc = new MCVar(contextdept, mc.getMu(), mc.getCk(), tree);
			vc.printOutput();
			vc.printToScreen();
			
			
		} else if (args.length == 6) {
			String groupFile = args[5];
			
			MainRead main = new MainRead(C, branchNum);
	    	double[][][] prop = new double[branchNum][6][C];
	    	int[][][] count = new int[branchNum][18][C];
			prop = main.getProp(State, branchNum, C);
			count = main.getCount(numChange, branchNum, C);
	    	
			String newick = "((((human_G:0.017228:G,orangutan_G:0.018799:G)_G:0.002182:G,gibbon_G:0.021454:G)_G:0.011141:G,((rhesus_G:0.008145:G,baboon_G:0.008297:G)_G:0.00478:G,greenMonkey_G:0.013339:G)_G:0.027207:G)_G:0.017203:G,(marmoset_G:0.034144:G,squirrelMonkey_G:0.033447:G)_G:0.04928:G,bushbaby_G:0.244601:G)_G;";
			Tree tree = main.getTree(outgroup, newick); 
			
			ArrayList<ArrayList<Integer>> group = main.grouping(groupFile);
            
			SufficientStatistics contextdept = new SufficientStatistics(false, branchNum, C, prop, count);			
			Estimation mc = new Estimation(group, contextdept);
			//VarCov vc = new VarCov(group, contextdept, mc.getMu(), mc.getCk(), tree);
			MCVar vc = new MCVar(group, contextdept, mc.getMu(), mc.getCk(), tree);
			vc.printOutput();
			vc.printToScreen();		
			
		} else {
			System.out.println("Argument error.");
		}

	}

}
